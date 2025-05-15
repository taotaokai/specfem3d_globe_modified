!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


  subroutine locate_point(x_target,y_target,z_target,lat_target,lon_target,ispec_selected,xi,eta,gamma, &
                          x,y,z,distmin_not_squared,POINT_CAN_BE_BURIED)

! locates target point inside best mesh element

  use constants_solver, only: &
    NGLLX,NGLLY,NGLLZ,MIDX,MIDY,MIDZ,HUGEVAL, &
    USE_DISTANCE_CRITERION

  use shared_parameters, only: R_PLANET_KM

  use specfem_par, only: &
    nspec => NSPEC_CRUST_MANTLE

  use specfem_par_crustmantle, only: &
    ibool => ibool_crust_mantle, &
    xstore => xstore_crust_mantle, &
    ystore => ystore_crust_mantle, &
    zstore => zstore_crust_mantle

  ! for point search
  use specfem_par, only: &
    typical_size_squared, &
    lat_min,lat_max,lon_min,lon_max,xyz_midpoints

  use kdtree_search, only: kdtree_find_nearest_neighbor

!debug: sync
!  use specfem_par, only: NPROCTOT_VAL,myrank

  implicit none

  double precision,intent(in) :: x_target,y_target,z_target
  double precision,intent(in) :: lat_target,lon_target

  integer,intent(out) :: ispec_selected
  double precision,intent(out) :: xi,eta,gamma
  double precision,intent(out) :: x,y,z
  double precision,intent(out) :: distmin_not_squared

  logical,intent(in) :: POINT_CAN_BE_BURIED

  ! local parameters
  integer :: ix_initial_guess,iy_initial_guess,iz_initial_guess
  integer :: ispec,i,j,k,iglob

  double precision :: lat,lon
  double precision :: dx,dy,dz
  double precision :: distmin_squared,dist_squared

  logical :: target_located

  ! nodes search
  double precision,dimension(3) :: xyz_target
  double precision :: dist_min
  integer :: inode_min

  ! brute-force search for closest element
  ! if set to .false. (default), a kd-tree search for initial guess element is used
  logical,parameter :: DO_BRUTE_FORCE_SEARCH = .false.

  ! looks for closer estimates in neighbor elements if needed
  logical,parameter :: DO_ADJACENT_SEARCH = .true.

  ! initialize location
  ! if we have not located a target element, the receiver is not in this slice
  ! therefore use first element only for fictitious iterative search
  ispec_selected = 1
  ix_initial_guess = MIDX
  iy_initial_guess = MIDY
  iz_initial_guess = MIDZ

  ! puts initial xi/eta/gamma outside of element
  xi = 2.d0
  eta = 2.d0
  gamma = 2.d0

  ! set distance to huge initial value
  distmin_squared = HUGEVAL

  ! limits latitude to [-90.0,90.0]
  lat = lat_target
  if (lat < -90.d0) lat = -90.d0
  if (lat > 90.d0) lat = 90.d0

  ! limits longitude to [0.0,360.0]
  lon = lon_target
  if (lon < 0.d0) lon = lon + 360.d0
  if (lon > 360.d0) lon = lon - 360.d0

  ! checks if receiver in this slice
  if (lat >= lat_min .and. lat <= lat_max .and. &
      lon >= lon_min .and. lon <= lon_max) then
    target_located = .true.
  else
    target_located = .false.
  endif

  ! debug
  !print *,'target located:',target_located,'lat',sngl(lat),sngl(lat_min),sngl(lat_max),'lon',sngl(lon),sngl(lon_min),sngl(lon_max)

  ! checks if anything to do
  if (.not. target_located) then
    ! point not found in this slice
    !
    ! returns initial guess point
    iglob = ibool(ix_initial_guess,iy_initial_guess,iz_initial_guess,ispec_selected)
    x = dble(xstore(iglob))
    y = dble(xstore(iglob))
    z = dble(xstore(iglob))

    ! compute final distance between asked and found (converted to km)
    dx = x_target-x
    dy = y_target-y
    dz = z_target-z
    distmin_not_squared = dsqrt(dx*dx + dy*dy + dz*dz) * R_PLANET_KM

    ! all done
    return
  endif

  ! point in this slice
  if (DO_BRUTE_FORCE_SEARCH) then
    ! loops over elements to find nearest location
    ! searches closest GLL point
    if (USE_DISTANCE_CRITERION) then
      ! loops over all elements
      do ispec = 1,nspec
        ! exclude elements that are too far from target
        dx = x_target - xyz_midpoints(1,ispec)
        dy = y_target - xyz_midpoints(2,ispec)
        dz = z_target - xyz_midpoints(3,ispec)
        dist_squared = dx*dx + dy*dy + dz*dz

        !  we compare squared distances instead of distances themselves to significantly speed up calculations
        if (dist_squared > typical_size_squared) cycle

        ! loop only on points inside the element
        ! exclude edges to ensure this point is not shared with other elements
        do k = 2,NGLLZ-1
          do j = 2,NGLLY-1
            do i = 2,NGLLX-1
              iglob = ibool(i,j,k,ispec)

              dx = x_target - dble(xstore(iglob))
              dy = y_target - dble(ystore(iglob))
              dz = z_target - dble(zstore(iglob))
              dist_squared = dx*dx + dy*dy + dz*dz

              !  keep this point if it is closer to the receiver
              !  we compare squared distances instead of distances themselves to significantly speed up calculations
              if (dist_squared < distmin_squared) then
                distmin_squared = dist_squared
                ispec_selected = ispec
                ix_initial_guess = i
                iy_initial_guess = j
                iz_initial_guess = k
              endif

            enddo
          enddo
        enddo

      ! end of loop on all the spectral elements in current slice
      enddo
    else
      ! searches through all elements
      do ispec = 1,nspec
        ! loop only on points inside the element
        ! exclude edges to ensure this point is not shared with other elements
        do k = 2,NGLLZ-1
          do j = 2,NGLLY-1
            do i = 2,NGLLX-1
              iglob = ibool(i,j,k,ispec)

              dx = x_target - dble(xstore(iglob))
              dy = y_target - dble(ystore(iglob))
              dz = z_target - dble(zstore(iglob))
              dist_squared = dx*dx + dy*dy + dz*dz

              !  keep this point if it is closer to the receiver
              !  we compare squared distances instead of distances themselves to significantly speed up calculations
              if (dist_squared < distmin_squared) then
                distmin_squared = dist_squared
                ispec_selected = ispec
                ix_initial_guess = i
                iy_initial_guess = j
                iz_initial_guess = k
                target_located = .true.
              endif
            enddo
          enddo
        enddo
      ! end of loop on all the spectral elements in current slice
      enddo
    endif ! USE_DISTANCE_CRITERION

    ! to get a good starting point for the iterative refinement,
    ! updates final guess for closest GLL point in the best selected element

    ! note: assumes that before we used looping range from [2, NGLL - 1]
    ! top/bottom surfaces
    do k = 1,NGLLZ,NGLLZ-1
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec_selected)

          dx = x_target - dble(xstore(iglob))
          dy = y_target - dble(ystore(iglob))
          dz = z_target - dble(zstore(iglob))
          dist_squared = dx*dx + dy*dy + dz*dz

          ! take this point if it is closer to the receiver
          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          if (dist_squared < distmin_squared) then
            distmin_squared = dist_squared
            ix_initial_guess = i
            iy_initial_guess = j
            iz_initial_guess = k
          endif
        enddo
      enddo
    enddo
    ! front/back surfaces
    do k = 1,NGLLZ
      do j = 1,NGLLY,NGLLY-1
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec_selected)

          dx = x_target - dble(xstore(iglob))
          dy = y_target - dble(ystore(iglob))
          dz = z_target - dble(zstore(iglob))
          dist_squared = dx*dx + dy*dy + dz*dz

          ! take this point if it is closer to the receiver
          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          if (dist_squared < distmin_squared) then
            distmin_squared = dist_squared
            ix_initial_guess = i
            iy_initial_guess = j
            iz_initial_guess = k
          endif
        enddo
      enddo
    enddo
    ! left/right surfaces
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX,NGLLX-1
          iglob = ibool(i,j,k,ispec_selected)

          dx = x_target - dble(xstore(iglob))
          dy = y_target - dble(ystore(iglob))
          dz = z_target - dble(zstore(iglob))
          dist_squared = dx*dx + dy*dy + dz*dz

          ! take this point if it is closer to the receiver
          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          if (dist_squared < distmin_squared) then
            distmin_squared = dist_squared
            ix_initial_guess = i
            iy_initial_guess = j
            iz_initial_guess = k
          endif
        enddo
      enddo
    enddo

  else
    ! kdtree search
    xyz_target(1) = x_target
    xyz_target(2) = y_target
    xyz_target(3) = z_target

    ! finds closest point (inside GLL points) in this chunk
    call kdtree_find_nearest_neighbor(xyz_target,inode_min,dist_min)

    ! debug
    !print *,'kd-tree found location :',inode_min

    ! sets element
    ispec_selected = inode_min

    ! loops over GLL points in this element to get (i,j,k) for initial guess
    do k = 2,NGLLZ-1
      do j = 2,NGLLY-1
        do i = 2,NGLLX-1
          iglob = ibool(i,j,k,ispec_selected)

          dx = x_target - dble(xstore(iglob))
          dy = y_target - dble(ystore(iglob))
          dz = z_target - dble(zstore(iglob))
          dist_squared = dx*dx + dy*dy + dz*dz

          ! take this point if it is closer to the receiver
          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          if (dist_squared < distmin_squared) then
            distmin_squared = dist_squared
            ix_initial_guess = i
            iy_initial_guess = j
            iz_initial_guess = k
          endif
        enddo
      enddo
    enddo
  endif ! DO_BRUTE_FORCE_SEARCH

  !debug
  !print *,'initial guess ',ix_initial_guess,iy_initial_guess,iz_initial_guess,ispec_selected

  ! ****************************************
  ! find the best (xi,eta,gamma)
  ! ****************************************

!deprecated, but might be useful if one wishes to get an exact GLL point location:
!
!      ! for point sources, the location will be exactly at a GLL point
!      ! otherwise this tries to find best location
!
!      if (USE_FORCE_POINT_SOURCE) then
!        ! store xi,eta,gamma and x,y,z of point found
!        ! note: they have range [1.0d0,NGLLX/Y/Z], used for point sources
!        !          see e.g. in compute_add_sources.f90
!        xi_subset(isource_in_this_subset) = dble(ix_initial_guess)
!        eta_subset(isource_in_this_subset) = dble(iy_initial_guess)
!        gamma_subset(isource_in_this_subset) = dble(iz_initial_guess)
!
!        iglob = ibool(ix_initial_guess,iy_initial_guess, &
!            iz_initial_guess,ispec_selected)
!        xyz_found_subset(1,isource_in_this_subset) = xstore(iglob)
!        xyz_found_subset(2,isource_in_this_subset) = ystore(iglob)
!        xyz_found_subset(3,isource_in_this_subset) = zstore(iglob)
!
!        ! compute final distance between asked and found (converted to km)
!        final_distance_subset(isource_in_this_subset) = &
!          dsqrt((x_target-xyz_found_subset(1,isource_in_this_subset))**2 + &
!                (y_target-xyz_found_subset(2,isource_in_this_subset))**2 + &
!                (z_target-xyz_found_subset(3,isource_in_this_subset))**2)*R_PLANET/1000.d0
!      endif

  ! gets xi/eta/gamma and corresponding x/y/z coordinates
  call find_local_coordinates(x_target,y_target,z_target,xi,eta,gamma,x,y,z, &
                              ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess, &
                              POINT_CAN_BE_BURIED)

  ! loops over neighbors and try to find better location
  if (DO_ADJACENT_SEARCH) then
    ! checks if position lies on an element boundary
    if (abs(xi) > 1.099d0 .or. abs(eta) > 1.099d0 .or. abs(gamma) > 1.099d0) then
      ! searches for better position in neighboring elements
      call find_best_neighbor(x_target,y_target,z_target,xi,eta,gamma,x,y,z,ispec_selected,distmin_squared,POINT_CAN_BE_BURIED)
    endif
  endif ! DO_ADJACENT_SEARCH

  ! compute final distance between asked and found (converted to km)
  dx = x_target-x
  dy = y_target-y
  dz = z_target-z
  distmin_not_squared = dsqrt(dx*dx + dy*dy + dz*dz) * R_PLANET_KM

  end subroutine locate_point

!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_local_coordinates(x_target,y_target,z_target,xi,eta,gamma,x,y,z, &
                                    ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess, &
                                    POINT_CAN_BE_BURIED)

  use constants_solver, only: &
    NGNOD,HUGEVAL,NUM_ITER

  use specfem_par_crustmantle, only: &
    ibool => ibool_crust_mantle, &
    xstore => xstore_crust_mantle,ystore => ystore_crust_mantle,zstore => zstore_crust_mantle

  ! for point search
  use specfem_par, only: &
    anchor_iax,anchor_iay,anchor_iaz, &
    xigll,yigll,zigll

  implicit none

  double precision,intent(in) :: x_target,y_target,z_target
  double precision,intent(out) :: xi,eta,gamma
  double precision,intent(out) :: x,y,z

  integer,intent(in) :: ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess

  logical,intent(in) :: POINT_CAN_BE_BURIED

  ! local parameters
  integer :: ia,iter_loop
  integer :: iglob

  ! coordinates of the control points of the surface element
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  double precision :: dx,dy,dz
  double precision :: d_sq,d_min_sq
  double precision :: dxi,deta,dgamma
  double precision :: xix,xiy,xiz
  double precision :: etax,etay,etaz
  double precision :: gammax,gammay,gammaz

  ! define coordinates of the control points of the element
  do ia = 1,NGNOD
    iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec_selected)
    xelm(ia) = dble(xstore(iglob))
    yelm(ia) = dble(ystore(iglob))
    zelm(ia) = dble(zstore(iglob))
  enddo

  ! use initial guess in xi and eta
  xi = xigll(ix_initial_guess)
  eta = yigll(iy_initial_guess)
  gamma = zigll(iz_initial_guess)

  ! impose receiver exactly at the surface
  if (.not. POINT_CAN_BE_BURIED) gamma = 1.d0

  d_min_sq = HUGEVAL

  ! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

    ! recompute Jacobian for the new point
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                            xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

    ! compute distance to target location
    dx = - (x - x_target)
    dy = - (y - y_target)
    dz = - (z - z_target)

    !debug
    !print *,'  iter ',iter_loop,'dx',sngl(dx),sngl(dx_min),'dy',sngl(dy),sngl(dy_min),'dz',sngl(dz),sngl(dz_min),d_min_sq

    ! distance squared
    d_sq = dx*dx + dy*dy + dz*dz

    ! compute increments
    if (d_sq < d_min_sq) then
      d_min_sq = d_sq

      dxi = xix*dx + xiy*dy + xiz*dz
      deta = etax*dx + etay*dy + etaz*dz
      dgamma = gammax*dx + gammay*dy + gammaz*dz
    else
      ! new position is worse than old one, no change necessary
      ! stop, no further improvements
      !dxi = 0.d0
      !deta = 0.d0
      !dgamma = 0.d0
      exit
    endif

    ! decreases step length if step is large
    if ((dxi*dxi + deta*deta + dgamma*dgamma) > 1.0d0) then
      dxi = dxi * 0.33333333333d0
      deta = deta * 0.33333333333d0
      dgamma = dgamma * 0.33333333333d0
    endif
    ! alternative: impose limit on increments (seems to result in slightly less accurate locations)
    !if (abs(dxi) > 0.3d0 ) dxi = sign(1.0d0,dxi)*0.3d0
    !if (abs(deta) > 0.3d0 ) deta = sign(1.0d0,deta)*0.3d0
    !if (abs(dgamma) > 0.3d0 ) dgamma = sign(1.0d0,dgamma)*0.3d0

    !debug
    !print *,'  dxi/..',(dxi**2 + deta**2 + dgamma**2),dxi,deta,dgamma

    ! update values
    xi = xi + dxi
    eta = eta + deta
    gamma = gamma + dgamma

    ! impose that we stay in that element
    ! (useful if user gives a receiver outside the mesh for instance)
    ! we can go slightly outside the [1,1] segment since with finite elements
    ! the polynomial solution is defined everywhere
    ! can be useful for convergence of iterative scheme with distorted elements
    !if (xi > 1.10d0) xi = 1.10d0
    !if (xi <-1.10d0) xi = -1.10d0
    if (abs(xi) > 1.10d0) xi = sign(1.10d0,xi)
    !if (eta > 1.10d0) eta = 1.10d0
    !if (eta <-1.10d0) eta = -1.10d0
    if (abs(eta) > 1.10d0) eta = sign(1.10d0,eta)
    !if (gamma > 1.10d0) gamma = 1.10d0
    !if (gamma <-1.10d0) gamma = -1.10d0
    if (abs(gamma) > 1.10d0) gamma = sign(1.10d0,gamma)

  ! end of non linear iterations
  enddo

  ! impose receiver exactly at the surface
  if (.not. POINT_CAN_BE_BURIED) gamma = 1.d0

  ! compute final coordinates of point found
  call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  end subroutine find_local_coordinates

!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_best_neighbor(x_target,y_target,z_target, &
                                xi_found,eta_found,gamma_found,x_found,y_found,z_found, &
                                ispec_selected,final_distmin_squared,POINT_CAN_BE_BURIED)

  use constants_solver, only: &
    NGLLX,NGLLY,NGLLZ,MIDX,MIDY,MIDZ

  use shared_parameters, only: R_PLANET_KM

  use specfem_par, only: &
    nspec => NSPEC_CRUST_MANTLE

  use specfem_par_crustmantle, only: &
    ibool => ibool_crust_mantle, &
    xstore => xstore_crust_mantle,ystore => ystore_crust_mantle,zstore => zstore_crust_mantle

  ! for point search
  use specfem_par, only: &
    xadj,adjncy

  use kdtree_search, only: kdtree_find_nearest_neighbor

  implicit none

  double precision,intent(in) :: x_target,y_target,z_target
  double precision,intent(inout) :: xi_found,eta_found,gamma_found
  double precision,intent(inout) :: x_found,y_found,z_found

  integer,intent(inout) :: ispec_selected
  double precision,intent(inout) :: final_distmin_squared
  logical,intent(in) :: POINT_CAN_BE_BURIED

  ! local parameters
  integer :: ix_initial_guess,iy_initial_guess,iz_initial_guess
  integer :: ispec,i,j,k,iglob

  double precision :: dx,dy,dz
  double precision :: dist_squared
  double precision :: distmin_squared_guess

  ! nodes search
  double precision :: xi_n,eta_n,gamma_n,x_n,y_n,z_n ! neighbor position result

  ! neighbor elements
  integer,parameter :: MAX_NEIGHBORS = 50   ! maximum number of neighbors (around 37 should be sufficient for crust/mantle)
  integer :: index_neighbors(MAX_NEIGHBORS*MAX_NEIGHBORS) ! including neighbors of neighbors
  integer :: num_neighbors

  integer :: ii,jj,i_n,ientry,ispec_ref
  logical :: do_neighbor
  logical :: is_better_location

  ! for close points
  double precision, parameter :: TOL_POINT_DISTANCE = 1.d-10
  ! for close points relative to each other
  double precision, parameter :: TOL_RELATIVE_POINT_DISTANCE = 1.d-15

  ! verbose output
  logical,parameter :: DEBUG = .false.

  ! best distance to target .. so far
  dx = x_target - x_found
  dy = y_target - y_found
  dz = z_target - z_found
  final_distmin_squared = dx*dx + dy*dy + dz*dz

  !debug
  if (DEBUG) print *,'neighbors: best guess ',ispec_selected,xi_found,eta_found,gamma_found, &
                     'distance',sngl(sqrt(final_distmin_squared)*R_PLANET_KM)

  ! fill neighbors arrays
  !
  ! note: we add direct neighbors plus neighbors of neighbors.
  !       for very coarse meshes, the initial location guesses especially around doubling layers can be poor such that we need
  !       to enlarge the search of neighboring elements.
  index_neighbors(:) = 0
  num_neighbors = 0
  do ii = 1,xadj(ispec_selected+1)-xadj(ispec_selected)
    ! get neighbor
    ientry = xadj(ispec_selected) + ii
    ispec_ref = adjncy(ientry)

    ! checks
    if (ispec_ref < 1 .or. ispec_ref > nspec) stop 'Invalid ispec index in locate point search'

    ! checks if exists already in list
    do_neighbor = .true.
    do i_n = 1,num_neighbors
      if (index_neighbors(i_n) == ispec_ref) then
        do_neighbor = .false.
        exit
      endif
    enddo

    ! adds to search elements
    if (do_neighbor) then
      num_neighbors = num_neighbors + 1
      index_neighbors(num_neighbors) = ispec_ref
    endif

    ! adds neighbors of neighbor
    do jj = 1,xadj(ispec_ref+1)-xadj(ispec_ref)
      ! get neighbor
      ientry = xadj(ispec_ref) + jj
      ispec = adjncy(ientry)

      ! checks
      if (ispec < 1 .or. ispec > nspec) stop 'Invalid ispec index in locate point search'

      ! checks if exists already in list
      do_neighbor = .true.
      do i_n = 1,num_neighbors
        if (index_neighbors(i_n) == ispec) then
          do_neighbor = .false.
          exit
        endif
      enddo

      ! adds to search elements
      if (do_neighbor) then
        num_neighbors = num_neighbors + 1
        index_neighbors(num_neighbors) = ispec
      endif
    enddo
  enddo

  ! loops over neighboring elements
  do i_n = 1,num_neighbors
    ispec = index_neighbors(i_n)

    ! note: the final position location can be still off if we start too far away.
    !       here we guess the best "inner" GLL point inside this search neighbor element
    !       to be the starting initial guess for finding the local coordinates.

    ! gets first guess as starting point
    ix_initial_guess = MIDX
    iy_initial_guess = MIDY
    iz_initial_guess = MIDZ
    iglob = ibool(MIDX,MIDY,MIDZ,ispec)

    dx = x_target - dble(xstore(iglob))
    dy = y_target - dble(ystore(iglob))
    dz = z_target - dble(zstore(iglob))
    distmin_squared_guess = dx*dx + dy*dy + dz*dz

    ! loop only on points inside the element
    ! exclude edges to ensure this point is not shared with other elements
    do k = 2,NGLLZ-1
      do j = 2,NGLLY-1
        do i = 2,NGLLX-1
          iglob = ibool(i,j,k,ispec)
          dx = x_target - dble(xstore(iglob))
          dy = y_target - dble(ystore(iglob))
          dz = z_target - dble(zstore(iglob))
          dist_squared = dx*dx + dy*dy + dz*dz

          !  keep this point if it is closer to the receiver
          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          if (dist_squared < distmin_squared_guess) then
            distmin_squared_guess = dist_squared
            ix_initial_guess = i
            iy_initial_guess = j
            iz_initial_guess = k
          endif
        enddo
      enddo
    enddo

    ! gets xi/eta/gamma and corresponding x/y/z coordinates
    call find_local_coordinates(x_target,y_target,z_target,xi_n,eta_n,gamma_n,x_n,y_n,z_n, &
                                ispec,ix_initial_guess,iy_initial_guess,iz_initial_guess,POINT_CAN_BE_BURIED)

    ! final distance to target
    dx = x_target - x_n
    dy = y_target - y_n
    dz = z_target - z_n
    dist_squared = dx*dx + dy*dy + dz*dz

    ! debug
    if (DEBUG) print *,'  neighbor ',ispec,i_n,ientry,'ispec = ',ispec_selected,sngl(xi_n),sngl(eta_n),sngl(gamma_n), &
                       'distance',sngl(sqrt(dist_squared)*R_PLANET_KM),sngl(sqrt(final_distmin_squared)*R_PLANET_KM)

    ! flag to determine if taking new position
    is_better_location = .false.

    ! if we found a close point
    if (dist_squared < TOL_POINT_DISTANCE) then
      ! checks if position within element
      if (abs(xi_n) <= 1.d0 .and. abs(eta_n) <= 1.d0 .and. abs(gamma_n) <= 1.d0) then
        ! takes position if old position has xi/eta/gamma outside of element
        if (abs(xi_found) > 1.d0 .or. abs(eta_found) > 1.d0 .or. abs(gamma_found) > 1.d0) then
          ! old position is out of element, new inside an element; let's take the new one
          is_better_location = .true.
        endif
      endif
    endif

    ! if we have found an element that gives a shorter distance
    ! (we compare squared distances instead of distances themselves to significantly speed up calculations)
    if (dist_squared < final_distmin_squared) then
      ! check if location is better in case almost similar distances
      ! since we allow xi/eta/gamma to go a bit outside of the element range [-1.01,1.01],
      ! we keep the element where xi/eta/gamma is still within element
      !
      ! note: this can avoid problems when points on interfaces, say between acoustic/elastic, have been choosen.
      !       the point search takes whatever element comes last with a best position, even if the point
      !       has been located slightly away from the element. thus, we try to check if the previous location
      !       was accurate enough with coordinates still within an element.
      is_better_location = .true.
      if (abs(dist_squared - final_distmin_squared) < TOL_RELATIVE_POINT_DISTANCE) then
        if (abs(xi_n) > 1.d0 .or. abs(eta_n) > 1.d0 .or. abs(gamma_n) > 1.d0) then
          if (abs(xi_found) > 1.d0 .or. abs(eta_found) > 1.d0 .or. abs(gamma_found) > 1.d0) then
            ! old position is also out of element, let's take the new one
            is_better_location = .true.
          else
            ! old position is still within element and new distance has only little improved, so keep old one
            is_better_location = .false.
          endif
        endif
      endif
    endif

    ! takes this point if it is better
    if (is_better_location) then
      final_distmin_squared = dist_squared
      ! uses this as new location
      ispec_selected = ispec
      xi_found = xi_n
      eta_found = eta_n
      gamma_found = gamma_n
      x_found = x_n
      y_found = y_n
      z_found = z_n
    endif

    ! checks if position lies inside element (which usually means that located position is accurate)
    if (abs(xi_found) <= 1.d0 .and. abs(eta_found) <= 1.d0 .and. abs(gamma_found) <= 1.d0) exit

  enddo ! num_neighbors

  !debug
  if (DEBUG) print *,'neighbors: final ',ispec_selected,xi_found,eta_found,gamma_found, &
                     'distance',sqrt(final_distmin_squared)*R_PLANET_KM

  end subroutine find_best_neighbor


!
!-------------------------------------------------------------------------------------------------
!

  subroutine locate_MPI_slice(npoints_subset,ipoin_already_done, &
                              ispec_selected_subset, &
                              xyz_found_subset, &
                              xi_subset,eta_subset,gamma_subset, &
                              final_distance_subset, &
                              npoints_total, ispec_selected, islice_selected, &
                              xyz_found, &
                              xi_point, eta_point, gamma_point, &
                              final_distance)

! locate MPI slice which contains the point and bcast to all
! (locates subset of points in all slices)

  use constants, only: HUGEVAL,NDIM
  use specfem_par, only: NPROCTOT_VAL,myrank

  integer, intent(in) :: npoints_subset,ipoin_already_done
  integer, dimension(npoints_subset), intent(in)  :: ispec_selected_subset
  double precision, dimension(NDIM,npoints_subset), intent(in) :: xyz_found_subset
  double precision, dimension(npoints_subset), intent(in) :: xi_subset, eta_subset, gamma_subset
  double precision, dimension(npoints_subset), intent(in) :: final_distance_subset

  integer, intent(in) :: npoints_total
  integer, dimension(npoints_total), intent(inout)  :: ispec_selected, islice_selected
  double precision, dimension(NDIM,npoints_total), intent(inout)  :: xyz_found
  double precision, dimension(npoints_total), intent(inout)  :: xi_point, eta_point, gamma_point
  double precision, dimension(npoints_total), intent(inout)  :: final_distance

  ! local parameters
  integer :: ipoin,ipoin_in_this_subset,iproc
  double precision :: distmin

  ! gather arrays
  integer, dimension(npoints_subset,0:NPROCTOT_VAL-1) :: ispec_selected_all
  double precision, dimension(npoints_subset,0:NPROCTOT_VAL-1) :: xi_all,eta_all,gamma_all
  double precision, dimension(NDIM,npoints_subset,0:NPROCTOT_VAL-1) :: xyz_found_all
  double precision, dimension(npoints_subset,0:NPROCTOT_VAL-1) :: final_distance_all

  double precision :: xi,eta,gamma

  ! point inside element
  logical :: found_point_inside

  ! prefer slice with point location inside element rather than just selecting slice with closest location
  logical, parameter :: DO_PREFER_INSIDE_ELEMENT = .true.

  !debug
  !integer,dimension(1) :: iproc_min

  ! initializes with dummy values
  ispec_selected_all(:,:) = -1

  xi_all(:,:) = 0.d0
  eta_all(:,:) = 0.d0
  gamma_all(:,:) = 0.d0
  xyz_found_all(:,:,:) = 0.d0

  final_distance_all(:,:) = HUGEVAL

  ! gather all (on main process)
  call gather_all_i(ispec_selected_subset,npoints_subset,ispec_selected_all,npoints_subset,NPROCTOT_VAL)

  ! this is executed by main process only
  if (myrank == 0) then
    ! check that the gather operation went well
    if (any(ispec_selected_all(:,:) == -1)) then
      print *,'Error ispec all: procs = ',NPROCTOT_VAL,' subset size = ',npoints_subset
      print *,ispec_selected_all(:,:)
      call exit_MPI(myrank,'gather operation failed for point locations')
    endif
  endif

  call gather_all_dp(xyz_found_subset,NDIM*npoints_subset,xyz_found_all,NDIM*npoints_subset,NPROCTOT_VAL)

  call gather_all_dp(xi_subset,npoints_subset,xi_all,npoints_subset,NPROCTOT_VAL)
  call gather_all_dp(eta_subset,npoints_subset,eta_all,npoints_subset,NPROCTOT_VAL)
  call gather_all_dp(gamma_subset,npoints_subset,gamma_all,npoints_subset,NPROCTOT_VAL)

  call gather_all_dp(final_distance_subset,npoints_subset,final_distance_all,npoints_subset,NPROCTOT_VAL)

  ! find the slice and element to put the source
  if (myrank == 0) then

    ! we prefer slices with point locations inside an element over slices with a closer point distance
    ! but a point location outside the (best) element
    if (DO_PREFER_INSIDE_ELEMENT) then

      ! loops over subset
      do ipoin_in_this_subset = 1,npoints_subset
        ! mapping from station/source number in current subset to real station/source number in all the subsets
        ipoin = ipoin_in_this_subset + ipoin_already_done

        ! checks first if we have a close point inside an element
        found_point_inside = .false.
        distmin = HUGEVAL
        do iproc = 0,NPROCTOT_VAL-1
          ! point position
          xi = xi_all(ipoin_in_this_subset,iproc)
          eta = eta_all(ipoin_in_this_subset,iproc)
          gamma = gamma_all(ipoin_in_this_subset,iproc)

          ! points inside an element have xi/eta/gamma <= 1.0
          if (abs(xi) <= 1.d0 .and. abs(eta) <= 1.d0 .and. abs(gamma) <= 1.d0) then
            ! takes point if closer
            if (final_distance_all(ipoin_in_this_subset,iproc) < distmin) then
              found_point_inside = .true.

              xi_point(ipoin)    = xi
              eta_point(ipoin)   = eta
              gamma_point(ipoin) = gamma

              distmin =  final_distance_all(ipoin_in_this_subset,iproc)

              islice_selected(ipoin) = iproc
              ispec_selected(ipoin) = ispec_selected_all(ipoin_in_this_subset,iproc)

              xyz_found(:,ipoin) = xyz_found_all(:,ipoin_in_this_subset,iproc)
            endif
          endif
        enddo

        ! if we haven't found a close point inside an element, then look for just the closest possible
        if (.not. found_point_inside) then
          do iproc = 0,NPROCTOT_VAL-1
            if (final_distance_all(ipoin_in_this_subset,iproc) < distmin) then
              distmin =  final_distance_all(ipoin_in_this_subset,iproc)

              islice_selected(ipoin) = iproc
              ispec_selected(ipoin) = ispec_selected_all(ipoin_in_this_subset,iproc)

              xi_point(ipoin)    = xi_all(ipoin_in_this_subset,iproc)
              eta_point(ipoin)   = eta_all(ipoin_in_this_subset,iproc)
              gamma_point(ipoin) = gamma_all(ipoin_in_this_subset,iproc)

              xyz_found(:,ipoin) = xyz_found_all(:,ipoin_in_this_subset,iproc)
            endif
          enddo
        endif

        final_distance(ipoin) = distmin
        if (final_distance(ipoin) == HUGEVAL) call exit_MPI(myrank,'Error locating receiver')

        !debug
        !iproc_min = minloc(final_distance_all(ipoin_in_this_subset,:)) - 1 ! -1 to start procs at 0
        !print *,'debug: locate MPI point',ipoin,ipoin_in_this_subset,'dist',final_distance_all(ipoin_in_this_subset,:), &
        !  'iproc min',iproc_min(1), &
        !  'xi min',xi_all(ipoin_in_this_subset,iproc_min(1)), &
        !           eta_all(ipoin_in_this_subset,iproc_min(1)), &
        !           gamma_all(ipoin_in_this_subset,iproc_min(1)), &
        !  'iproc found',islice_selected(ipoin),'dist found',final_distance(ipoin), &
        !  'xi found',xi_point(ipoin),eta_point(ipoin),gamma_point(ipoin)
      enddo

    else
      ! old version takes closest point

      ! loops over subset
      do ipoin_in_this_subset = 1,npoints_subset

        ! mapping from station/source number in current subset to real station/source number in all the subsets
        ipoin = ipoin_in_this_subset + ipoin_already_done

        distmin = HUGEVAL
        do iproc = 0,NPROCTOT_VAL-1
          if (final_distance_all(ipoin_in_this_subset,iproc) < distmin) then
            distmin =  final_distance_all(ipoin_in_this_subset,iproc)

            islice_selected(ipoin) = iproc
            ispec_selected(ipoin) = ispec_selected_all(ipoin_in_this_subset,iproc)

            xi_point(ipoin)    = xi_all(ipoin_in_this_subset,iproc)
            eta_point(ipoin)   = eta_all(ipoin_in_this_subset,iproc)
            gamma_point(ipoin) = gamma_all(ipoin_in_this_subset,iproc)

            xyz_found(:,ipoin) = xyz_found_all(:,ipoin_in_this_subset,iproc)
          endif
        enddo
        final_distance(ipoin) = distmin
        if (final_distance(ipoin) == HUGEVAL) call exit_MPI(myrank,'Error locating receiver')
      enddo

    endif

  endif ! end of section executed by main process only

  end subroutine locate_MPI_slice
