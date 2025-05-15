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

  subroutine movie_volume_init_hdf5()

#ifdef USE_HDF5
    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    integer :: ier

    ! nglobs
    allocate(offset_poin_vol(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_poin_vol')
    allocate(offset_poin_vol_cm(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_poin_vol_cm')
    allocate(offset_poin_vol_oc(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_poin_vol_oc')
    allocate(offset_poin_vol_ic(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_poin_vol_ic')
    allocate(offset_nspec_vol(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_vol')
    allocate(offset_nspec_vol_cm(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_vol_cm')
    allocate(offset_nspec_vol_oc(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_vol_oc')
    allocate(offset_nspec_vol_ic(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_vol_ic')

    npoints_vol_mov_all_proc = 0
    npoints_vol_mov_all_proc_cm = 0
    npoints_vol_mov_all_proc_oc = 0
    npoints_vol_mov_all_proc_ic = 0

    nspec_vol_mov_all_proc = 0
    nspec_vol_mov_all_proc_cm = 0
    nspec_vol_mov_all_proc_oc = 0
    nspec_vol_mov_all_proc_ic = 0

    nspec_vol_mov_all_proc_cm_conn = 0
    nspec_vol_mov_all_proc_oc_conn = 0
    nspec_vol_mov_all_proc_ic_conn = 0

    ! if strain or vector output is requested
     if ((MOVIE_VOLUME_TYPE == 1 .or. MOVIE_VOLUME_TYPE == 2 .or. MOVIE_VOLUME_TYPE == 3 &
                                 .or. MOVIE_VOLUME_TYPE == 5 .or. MOVIE_VOLUME_TYPE == 6 )) then
      output_sv = .true.
    endif
    ! check if crust mantle region is used for movie
    if (MOVIE_VOLUME_TYPE == 4 .and. OUTPUT_CRUST_MANTLE) then
      output_cm = .true.
    endif
    ! check if outer core region is used for movie
    if ( (MOVIE_VOLUME_TYPE == 4 .or. MOVIE_VOLUME_TYPE == 7 .or. MOVIE_VOLUME_TYPE == 8 .or. MOVIE_VOLUME_TYPE == 9) &
      .and. OUTPUT_OUTER_CORE) then
      output_oc = .true.
    endif
    ! check if inner core region is used for movie
    if ( (MOVIE_VOLUME_TYPE == 4 .or. MOVIE_VOLUME_TYPE == 7 .or. MOVIE_VOLUME_TYPE == 8 .or. MOVIE_VOLUME_TYPE == 9) &
      .and. OUTPUT_INNER_CORE) then
      output_ic = .true.
    endif

    ! force false if number of elements is zero
    if (NSPEC_CRUST_MANTLE == 0) output_cm = .false.
    if (NSPEC_OUTER_CORE == 0) output_oc = .false.
    if (NSPEC_INNER_CORE == 0) output_ic = .false.

    ! print
    !print *, 'output_sv = ',output_sv
    !print *, 'output_cm = ',output_cm
    !print *, 'output_oc = ',output_oc
    !print *, 'output_ic = ',output_ic
    !print *, 'NSPEC_CRUST_MANTLE = ',NSPEC_CRUST_MANTLE
    !print *, 'NSPEC_CRUST_MANTLE_STRAIN_ONLY = ',NSPEC_CRUST_MANTLE_STRAIN_ONLY
    !print *, 'NSPEC_CRUST_MANTLE_STR_OR_ATT = ',NSPEC_CRUST_MANTLE_STR_OR_ATT
    !print *, 'NSPEC_OUTER_CORE = ',NSPEC_OUTER_CORE
    !print *, 'NSPEC_OUTER_CORE_3DMOVIE = ',NSPEC_OUTER_CORE_3DMOVIE
    !print *, 'NSPEC_INNER_CORE = ',NSPEC_INNER_CORE

#else
    write (*,*) 'Error: HDF5 is not enabled in this version of the code.'
    write (*,*) 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop
#endif

  end subroutine movie_volume_init_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine movie_volume_finalize_hdf5()

#ifdef USE_HDF5
    use specfem_par_movie_hdf5

    deallocate(offset_poin_vol)
    deallocate(offset_poin_vol_cm)
    deallocate(offset_poin_vol_oc)
    deallocate(offset_poin_vol_ic)
    deallocate(offset_nspec_vol)
    deallocate(offset_nspec_vol_cm)
    deallocate(offset_nspec_vol_oc)
    deallocate(offset_nspec_vol_ic)

#else
    write (*,*) 'Error: HDF5 is not enabled in this version of the code.'
    write (*,*) 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop
#endif

  end subroutine movie_volume_finalize_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_volume_mesh_hdf5(nu_3dmovie,num_ibool_3dmovie,mask_3dmovie,mask_ibool_3dmovie, &
                                          muvstore_crust_mantle_3dmovie,npoints_3dmovie,nelems_3dmovie_in)

  use specfem_par

#ifdef USE_HDF5
  use specfem_par_crustmantle, only: ibool_crust_mantle,rstore_crust_mantle
  use specfem_par_outercore, only: ibool_outer_core,rstore_outer_core
  use specfem_par_innercore, only: ibool_inner_core,rstore_inner_core
  use specfem_par_movie_hdf5
#endif

  implicit none

  integer,intent(in) :: npoints_3dmovie,nelems_3dmovie_in
  integer, dimension(NGLOB_CRUST_MANTLE_3DMOVIE),intent(in) :: num_ibool_3dmovie

  real(kind=CUSTOM_REAL), dimension(3,3,npoints_3dmovie),intent(inout) :: nu_3dmovie

  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE),intent(in) :: mask_3dmovie
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE),intent(in) :: muvstore_crust_mantle_3dmovie
  logical, dimension(NGLOB_CRUST_MANTLE_3DMOVIE),intent(in) :: mask_ibool_3dmovie

#ifdef USE_HDF5
  ! local parameters
  integer :: ipoints_3dmovie,ispec,i,j,k,iNIT,nelems_3dmovie
  integer :: iglob,iglob_center
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval,st,ct,sp,cp
  real(kind=CUSTOM_REAL), dimension(npoints_3dmovie)     :: store_val3D_x,    store_val3D_y,    store_val3D_z
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE)  :: store_val3D_x_cm, store_val3D_y_cm, store_val3D_z_cm
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE)    :: store_val3D_x_oc, store_val3D_y_oc, store_val3D_z_oc
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE)    :: store_val3D_x_ic, store_val3D_y_ic, store_val3D_z_ic
  real(kind=CUSTOM_REAL), dimension(npoints_3dmovie)     :: store_val3D_mu
  ! dummy num_ibool_3dmovie for cm oc ic
  integer, dimension(NGLOB_CRUST_MANTLE) :: num_ibool_3dmovie_cm
  integer, dimension(NGLOB_OUTER_CORE)   :: num_ibool_3dmovie_oc
  integer, dimension(NGLOB_INNER_CORE)   :: num_ibool_3dmovie_ic
  ! dummy mask_ibool_3dmovie for cm oc ic
  logical, dimension(NGLOB_CRUST_MANTLE) :: mask_ibool_3dmovie_cm
  logical, dimension(NGLOB_OUTER_CORE)   :: mask_ibool_3dmovie_oc
  logical, dimension(NGLOB_INNER_CORE)   :: mask_ibool_3dmovie_ic

  integer, dimension(:,:), allocatable :: elm_conn, elm_conn_cm, elm_conn_oc, elm_conn_ic

  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_vol_cm_conn, offset_nspec_vol_oc_conn, offset_nspec_vol_ic_conn

  integer :: nelems_3dmovie_cm, nelems_3dmovie_oc, nelems_3dmovie_ic

  ! initialize arrays for hdf5 volume movie output
  call movie_volume_init_hdf5()

  ! safety check
  if (NDIM /= 3) stop 'movie volume output requires NDIM = 3'

  ! output resolution
  if (MOVIE_COARSE) then
    iNIT = NGLLX-1
    nelems_3dmovie = nelems_3dmovie_in
  else
    iNIT = 1
    nelems_3dmovie = nelems_3dmovie_in * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
  endif

  ! outer core and inner core is always fine mesh
  nelems_3dmovie_cm = NSPEC_CRUST_MANTLE * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
  nelems_3dmovie_oc = NSPEC_OUTER_CORE   * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
  nelems_3dmovie_ic = NSPEC_INNER_CORE   * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)

  ! allocate elm_conn
  allocate(elm_conn(9,nelems_3dmovie))
  allocate(elm_conn_cm(9,nelems_3dmovie_cm))
  allocate(elm_conn_oc(9,nelems_3dmovie_oc))
  allocate(elm_conn_ic(9,nelems_3dmovie_ic))

  ! prepare offset array
  call gather_all_all_singlei(npoints_3dmovie,    offset_poin_vol,    NPROCTOT_VAL)
  call gather_all_all_singlei(NGLOB_CRUST_MANTLE, offset_poin_vol_cm, NPROCTOT_VAL)
  call gather_all_all_singlei(NGLOB_OUTER_CORE,   offset_poin_vol_oc, NPROCTOT_VAL)
  call gather_all_all_singlei(NGLOB_INNER_CORE,   offset_poin_vol_ic, NPROCTOT_VAL)

  call gather_all_all_singlei(nelems_3dmovie,    offset_nspec_vol,    NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC_CRUST_MANTLE, offset_nspec_vol_cm, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC_OUTER_CORE,   offset_nspec_vol_oc, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC_INNER_CORE,   offset_nspec_vol_ic, NPROCTOT_VAL)

  ! offset arrays for element connectivity
  call gather_all_all_singlei(nelems_3dmovie_cm, offset_nspec_vol_cm_conn, NPROCTOT_VAL)
  call gather_all_all_singlei(nelems_3dmovie_oc, offset_nspec_vol_oc_conn, NPROCTOT_VAL)
  call gather_all_all_singlei(nelems_3dmovie_ic, offset_nspec_vol_ic_conn, NPROCTOT_VAL)

  npoints_vol_mov_all_proc    = sum(offset_poin_vol)
  npoints_vol_mov_all_proc_cm = sum(offset_poin_vol_cm)
  npoints_vol_mov_all_proc_oc = sum(offset_poin_vol_oc)
  npoints_vol_mov_all_proc_ic = sum(offset_poin_vol_ic)

  nspec_vol_mov_all_proc    = sum(offset_nspec_vol)
  nspec_vol_mov_all_proc_cm = sum(offset_nspec_vol_cm)
  nspec_vol_mov_all_proc_oc = sum(offset_nspec_vol_oc)
  nspec_vol_mov_all_proc_ic = sum(offset_nspec_vol_ic)

  nspec_vol_mov_all_proc_cm_conn = sum(offset_nspec_vol_cm_conn)
  nspec_vol_mov_all_proc_oc_conn = sum(offset_nspec_vol_oc_conn)
  nspec_vol_mov_all_proc_ic_conn = sum(offset_nspec_vol_ic_conn)

  !if (myrank == 0) then
  !  print *, 'npoints_vol_mov_all_proc = ',    npoints_vol_mov_all_proc
  !  print *, 'npoints_vol_mov_all_proc_cm = ', npoints_vol_mov_all_proc_cm
  !  print *, 'npoints_vol_mov_all_proc_oc = ', npoints_vol_mov_all_proc_oc
  !  print *, 'npoints_vol_mov_all_proc_ic = ', npoints_vol_mov_all_proc_ic

  !  print *, 'nspec_vol_mov_all_proc = ',    nspec_vol_mov_all_proc
  !  print *, 'nspec_vol_mov_all_proc_cm = ', nspec_vol_mov_all_proc_cm
  !  print *, 'nspec_vol_mov_all_proc_oc = ', nspec_vol_mov_all_proc_oc
  !  print *, 'nspec_vol_mov_all_proc_ic = ', nspec_vol_mov_all_proc_ic

  !  print *, 'nspec_vol_mov_all_proc_cm_conn = ', nspec_vol_mov_all_proc_cm_conn
  !  print *, 'nspec_vol_mov_all_proc_oc_conn = ', nspec_vol_mov_all_proc_oc_conn
  !  print *, 'nspec_vol_mov_all_proc_ic_conn = ', nspec_vol_mov_all_proc_ic_conn
  !endif

  !
  ! create the xyz arrays for crust and mantle and strain and vector output
  !
  if (output_sv) then
    ! loops over all elements
    ipoints_3dmovie = 0
    do ispec = 1,NSPEC_CRUST_MANTLE

      ! checks center of element for movie flag
      iglob_center = ibool_crust_mantle((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)

      ! checks if movie element
      if (mask_ibool_3dmovie(iglob_center)) then

        ! stores element coordinates
        do k = 1,NGLLZ,iNIT
          do j = 1,NGLLY,iNIT
            do i = 1,NGLLX,iNIT
              ! only store points once
              if (mask_3dmovie(i,j,k,ispec)) then
                ! point increment
                ipoints_3dmovie = ipoints_3dmovie + 1

                ! gets point position
                iglob = ibool_crust_mantle(i,j,k,ispec)

                rval     = rstore_crust_mantle(1,iglob)
                thetaval = rstore_crust_mantle(2,iglob)
                phival   = rstore_crust_mantle(3,iglob)

                !x,y,z store have been converted to r theta phi already, need to revert back for xyz output
                call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

                store_val3D_x(ipoints_3dmovie)  = xval
                store_val3D_y(ipoints_3dmovie)  = yval
                store_val3D_z(ipoints_3dmovie)  = zval
                store_val3D_mu(ipoints_3dmovie) = muvstore_crust_mantle_3dmovie(i,j,k,ispec)

                st = sin(thetaval)
                ct = cos(thetaval)
                sp = sin(phival)
                cp = cos(phival)

                nu_3dmovie(1,1,ipoints_3dmovie) = -ct*cp
                nu_3dmovie(1,2,ipoints_3dmovie) = -ct*sp
                nu_3dmovie(1,3,ipoints_3dmovie) = st
                nu_3dmovie(2,1,ipoints_3dmovie) = -sp
                nu_3dmovie(2,2,ipoints_3dmovie) = cp
                nu_3dmovie(2,3,ipoints_3dmovie) = 0.d0
                nu_3dmovie(3,1,ipoints_3dmovie) = st*cp
                nu_3dmovie(3,2,ipoints_3dmovie) = st*sp
                nu_3dmovie(3,3,ipoints_3dmovie) = ct
              endif !mask_3dmovie
            enddo  !i
          enddo  !j
        enddo  !k
      endif

    enddo !ispec

    ! check if counters are correct
    if (ipoints_3dmovie /= npoints_3dmovie) then
      print *, 'Error: did not find the right number of points for 3D movie'
      print *, 'ipoints_3dmovie = ',ipoints_3dmovie,' npoints_3dmovie = ',npoints_3dmovie
      stop
    endif

  endif

  !
  ! create the xyz arrays for crust and mantle (not strain or vector output)
  !
  if (output_cm) then

    do ispec = 1,NSPEC_CRUST_MANTLE
      do k = 1,NGLLZ,1
        do j = 1,NGLLY,1
          do i = 1,NGLLX,1
            iglob           = ibool_crust_mantle(i,j,k,ispec)
            rval            = rstore_crust_mantle(1,iglob)
            thetaval        = rstore_crust_mantle(2,iglob)
            phival          = rstore_crust_mantle(3,iglob)

            call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

            store_val3D_x_cm(iglob) = xval
            store_val3D_y_cm(iglob) = yval
            store_val3D_z_cm(iglob) = zval

            ! dummy num_ibool_3dmovie for cm
            num_ibool_3dmovie_cm(iglob) = iglob
            ! all the mask_ibool_3dmovie_cm are true
            mask_ibool_3dmovie_cm(iglob) = .true.
          enddo
        enddo
      enddo
    enddo

  endif

  !
  ! create the xyz arrays for outer core
  !
  if (output_oc) then

    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1,NGLLZ,1
        do j = 1,NGLLY,1
          do i = 1,NGLLX,1
            iglob           = ibool_outer_core(i,j,k,ispec)
            rval            = rstore_outer_core(1,iglob)
            thetaval        = rstore_outer_core(2,iglob)
            phival          = rstore_outer_core(3,iglob)

            call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

            store_val3D_x_oc(iglob) = xval
            store_val3D_y_oc(iglob) = yval
            store_val3D_z_oc(iglob) = zval

            ! dummy num_ibool_3dmovie for oc
            num_ibool_3dmovie_oc(iglob) = iglob
            ! all the mask_ibool_3dmovie_oc are true
            mask_ibool_3dmovie_oc(iglob) = .true.
          enddo
        enddo
      enddo
    enddo

  endif

  !
  ! create the xyz arrays for inner core
  !
  if (output_ic) then

    do ispec = 1, NSPEC_INNER_CORE
      do k = 1,NGLLZ,1
        do j = 1,NGLLY,1
          do i = 1,NGLLX,1
            iglob           = ibool_inner_core(i,j,k,ispec)
            rval            = rstore_inner_core(1,iglob)
            thetaval        = rstore_inner_core(2,iglob)
            phival          = rstore_inner_core(3,iglob)

            call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

            store_val3D_x_ic(iglob) = xval
            store_val3D_y_ic(iglob) = yval
            store_val3D_z_ic(iglob) = zval

            ! dummy num_ibool_3dmovie for ic
            num_ibool_3dmovie_ic(iglob) = iglob
            ! all the mask_ibool_3dmovie_ic are true
            mask_ibool_3dmovie_ic(iglob) = .true.
          enddo
        enddo
      enddo
    enddo

  endif

  ! create elm_conn for movie
  ! for crust and mantle (strain and vector output)
  if (output_sv) call get_conn_for_movie(elm_conn,    sum(offset_poin_vol(0:myrank-1)), iNIT, nelems_3dmovie, &
                        npoints_3dmovie, NSPEC_CRUST_MANTLE, num_ibool_3dmovie, mask_ibool_3dmovie, ibool_crust_mantle)
  ! for crust and mantle (not strain or vector output)
  if (output_cm) call get_conn_for_movie(elm_conn_cm, sum(offset_poin_vol_cm(0:myrank-1)), 1, nelems_3dmovie_cm, &
                        NGLOB_CRUST_MANTLE, NSPEC_CRUST_MANTLE, num_ibool_3dmovie_cm, mask_ibool_3dmovie_cm, ibool_crust_mantle)
  ! for outer core
  if (output_oc) call get_conn_for_movie(elm_conn_oc, sum(offset_poin_vol_oc(0:myrank-1)), 1, nelems_3dmovie_oc, &
                        NGLOB_OUTER_CORE, NSPEC_OUTER_CORE, num_ibool_3dmovie_oc, mask_ibool_3dmovie_oc, ibool_outer_core)
  ! for inner core
  if (output_ic) call get_conn_for_movie(elm_conn_ic, sum(offset_poin_vol_ic(0:myrank-1)), 1, nelems_3dmovie_ic, &
                        NGLOB_INNER_CORE, NSPEC_INNER_CORE, num_ibool_3dmovie_ic, mask_ibool_3dmovie_ic, ibool_inner_core)

  ! TODO ADD IOSERVER

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! file name and group name
  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'mesh'

  ! create the file, group and dataset
  if (myrank == 0) then
    call h5_create_file(file_name)
    call h5_open_or_create_group(group_name)

    ! create the dataset
    ! for crust and mantle (strain and vector output)
    if (output_sv) then
      call h5_create_dataset_gen_in_group('elm_conn', (/9, nspec_vol_mov_all_proc/), 2, 1)
      call h5_create_dataset_gen_in_group('x', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('y', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('z', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    endif
    ! for crust and mantle (not strain or vector output)
    if (output_cm) then
      call h5_create_dataset_gen_in_group('elm_conn_cm', (/9, nspec_vol_mov_all_proc_cm_conn/), 2, 1)
      call h5_create_dataset_gen_in_group('x_cm', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('y_cm', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('z_cm', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
    endif
    ! for outer core
    if (output_oc) then
      call h5_create_dataset_gen_in_group('elm_conn_oc', (/9, nspec_vol_mov_all_proc_oc_conn/), 2, 1)
      call h5_create_dataset_gen_in_group('x_oc', (/npoints_vol_mov_all_proc_oc/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('y_oc', (/npoints_vol_mov_all_proc_oc/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('z_oc', (/npoints_vol_mov_all_proc_oc/), 1, CUSTOM_REAL)
    endif
    ! for inner core
    if (output_ic) then
      call h5_create_dataset_gen_in_group('elm_conn_ic', (/9, nspec_vol_mov_all_proc_ic_conn/), 2, 1)
      call h5_create_dataset_gen_in_group('x_ic', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('y_ic', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('z_ic', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
    endif

    ! close the group
    call h5_close_group()
    ! close the file
    call h5_close_file()

  endif

  ! synchronize
  call synchronize_all()

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  ! write the data
  ! for crust and mantle (strain and vector output)
  if (output_sv) then
    call h5_write_dataset_collect_hyperslab_in_group('elm_conn', elm_conn, (/0, sum(offset_nspec_vol(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('x', store_val3D_x, (/sum(offset_poin_vol(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('y', store_val3D_y, (/sum(offset_poin_vol(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('z', store_val3D_z, (/sum(offset_poin_vol(0:myrank-1))/), H5_COL)
  endif
  ! for crust and mantle (not strain or vector output)
  if (output_cm) then
    call h5_write_dataset_collect_hyperslab_in_group('elm_conn_cm', elm_conn_cm, &
                                                     (/0, sum(offset_nspec_vol_cm_conn(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('x_cm', store_val3D_x_cm, (/sum(offset_poin_vol_cm(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('y_cm', store_val3D_y_cm, (/sum(offset_poin_vol_cm(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('z_cm', store_val3D_z_cm, (/sum(offset_poin_vol_cm(0:myrank-1))/), H5_COL)
  endif
  ! for outer core
  if (output_oc) then
    call h5_write_dataset_collect_hyperslab_in_group('elm_conn_oc', elm_conn_oc, &
                                                     (/0, sum(offset_nspec_vol_oc_conn(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('x_oc', store_val3D_x_oc, (/sum(offset_poin_vol_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('y_oc', store_val3D_y_oc, (/sum(offset_poin_vol_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('z_oc', store_val3D_z_oc, (/sum(offset_poin_vol_oc(0:myrank-1))/), H5_COL)
  endif
  ! for inner core
  if (output_ic) then
    call h5_write_dataset_collect_hyperslab_in_group('elm_conn_ic', elm_conn_ic, &
                                                     (/0, sum(offset_nspec_vol_ic_conn(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('x_ic', store_val3D_x_ic, (/sum(offset_poin_vol_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('y_ic', store_val3D_y_ic, (/sum(offset_poin_vol_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group('z_ic', store_val3D_z_ic, (/sum(offset_poin_vol_ic(0:myrank-1))/), H5_COL)
  endif

  ! close the group
  call h5_close_group()
  ! close the file
  call h5_close_file_p()

  ! deallocate
  deallocate(elm_conn)
  deallocate(elm_conn_cm)
  deallocate(elm_conn_oc)
  deallocate(elm_conn_ic)

  ! write xdmf for all timesteps
  call write_xdmf_vol_hdf5(npoints_vol_mov_all_proc,    nspec_vol_mov_all_proc, &
                           npoints_vol_mov_all_proc_cm, nspec_vol_mov_all_proc_cm, &
                           npoints_vol_mov_all_proc_oc, nspec_vol_mov_all_proc_oc, &
                           npoints_vol_mov_all_proc_ic, nspec_vol_mov_all_proc_ic)

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = nelems_3dmovie_in

  idummy = size(num_ibool_3dmovie,kind=4)
  idummy = size(nu_3dmovie,kind=4)
  idummy = size(mask_3dmovie,kind=4)
  idummy = size(mask_ibool_3dmovie,kind=4)
  idummy = size(muvstore_crust_mantle_3dmovie,kind=4)

  print *, 'Error: HDF5 is not enabled in this version of the code.'
  print *, 'Please recompile with the HDF5 option enabled with --with-hdf5'
  stop

#endif

  end subroutine write_movie_volume_mesh_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_volume_strains_hdf5(vnspec_eps_cm, &
                                        eps_trace_over_3_crust_mantle, &
                                        vnspec_cm, &
                                        epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                        epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle)

  use constants_solver
#ifdef USE_HDF5
  use shared_parameters, only: OUTPUT_FILES,MOVIE_VOLUME_TYPE,MOVIE_COARSE,H5_COL
  use specfem_par, only: it
  use specfem_par_movie, only: npoints_3dmovie,muvstore_crust_mantle_3dmovie,mask_3dmovie,nu_3dmovie
  use specfem_par_movie_hdf5
#endif

  implicit none

  ! input
  integer,intent(in) :: vnspec_eps_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_eps_cm),intent(in) :: eps_trace_over_3_crust_mantle

  integer,intent(in) :: vnspec_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_cm),intent(in) :: &
    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle

#ifdef USE_HDF5
  ! variables
  real(kind=CUSTOM_REAL) :: muv_3dmovie
  real(kind=CUSTOM_REAL),dimension(3,3) :: eps_loc,eps_loc_new
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: store_val3d_NN,store_val3d_EE,store_val3d_ZZ, &
                                                     store_val3d_NE,store_val3d_NZ,store_val3d_EZ
  integer :: ipoints_3dmovie,i,j,k,ispec,iNIT,ier
  character(len=1) :: movie_prefix

  ! check
  if (NDIM /= 3) call exit_MPI(myrank, 'write_movie_volume_strains() requires NDIM = 3')
  if (vnspec_cm /= NSPEC_CRUST_MANTLE) call exit_MPI(myrank,'Invalid vnspec_cm value for write_movie_volume_strains() routine')

  ! allocates arrays
  allocate(store_val3d_NN(npoints_3dmovie), &
           store_val3d_EE(npoints_3dmovie), &
           store_val3d_ZZ(npoints_3dmovie), &
           store_val3d_NE(npoints_3dmovie), &
           store_val3d_NZ(npoints_3dmovie), &
           store_val3d_EZ(npoints_3dmovie), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating store_val3d_ .. arrays')

  if (MOVIE_VOLUME_TYPE == 1) then
    movie_prefix='E' ! strain
  else if (MOVIE_VOLUME_TYPE == 2) then
    movie_prefix='S' ! time integral of strain
  else if (MOVIE_VOLUME_TYPE == 3) then
    movie_prefix='P' ! potency, or integral of strain x \mu
  endif

  ! stepping
  if (MOVIE_COARSE) then
   iNIT = NGLLX-1
  else
   iNIT = 1
  endif

  ipoints_3dmovie = 0
  do ispec = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ,iNIT
      do j = 1,NGLLY,iNIT
        do i = 1,NGLLX,iNIT
          if (mask_3dmovie(i,j,k,ispec)) then
            ipoints_3dmovie = ipoints_3dmovie + 1
            muv_3dmovie = muvstore_crust_mantle_3dmovie(i,j,k,ispec)

            eps_loc(1,1) = eps_trace_over_3_crust_mantle(i,j,k,ispec) + epsilondev_xx_crust_mantle(i,j,k,ispec)
            eps_loc(2,2) = eps_trace_over_3_crust_mantle(i,j,k,ispec) + epsilondev_yy_crust_mantle(i,j,k,ispec)
            eps_loc(3,3) = eps_trace_over_3_crust_mantle(i,j,k,ispec) &
                           - epsilondev_xx_crust_mantle(i,j,k,ispec) &
                           - epsilondev_yy_crust_mantle(i,j,k,ispec)

            eps_loc(1,2) = epsilondev_xy_crust_mantle(i,j,k,ispec)
            eps_loc(1,3) = epsilondev_xz_crust_mantle(i,j,k,ispec)
            eps_loc(2,3) = epsilondev_yz_crust_mantle(i,j,k,ispec)

            eps_loc(2,1) = eps_loc(1,2)
            eps_loc(3,1) = eps_loc(1,3)
            eps_loc(3,2) = eps_loc(2,3)

            ! rotate eps_loc to spherical coordinates
            eps_loc_new(:,:) = matmul(matmul(nu_3dmovie(:,:,ipoints_3dmovie),eps_loc(:,:)), &
                                      transpose(nu_3dmovie(:,:,ipoints_3dmovie)))
            if (MOVIE_VOLUME_TYPE == 3) eps_loc_new(:,:) = eps_loc(:,:)*muv_3dmovie

            store_val3d_NN(ipoints_3dmovie) = eps_loc_new(1,1)
            store_val3d_EE(ipoints_3dmovie) = eps_loc_new(2,2)
            store_val3d_ZZ(ipoints_3dmovie) = eps_loc_new(3,3)
            store_val3d_NE(ipoints_3dmovie) = eps_loc_new(1,2)
            store_val3d_NZ(ipoints_3dmovie) = eps_loc_new(1,3)
            store_val3d_EZ(ipoints_3dmovie) = eps_loc_new(2,3)
          endif
        enddo
      enddo
    enddo
  enddo
  if (ipoints_3dmovie /= npoints_3dmovie) stop 'did not find the right number of points for 3D movie'

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! create group and datasets
  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  if (myrank == 0) then

    ! create the file, group and dataset
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    ! create the dataset
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'NN', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'EE', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'ZZ', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'NE', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'NZ', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'EZ', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)

    ! close the group
    call h5_close_group()
    ! close the file
    call h5_close_file()

  endif

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'NN', store_val3d_NN, &
                                                   (/sum(offset_poin_vol(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'EE', store_val3d_EE, &
                                                   (/sum(offset_poin_vol(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'ZZ', store_val3d_ZZ, &
                                                   (/sum(offset_poin_vol(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'NE', store_val3d_NE, &
                                                   (/sum(offset_poin_vol(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'NZ', store_val3d_NZ, &
                                                   (/sum(offset_poin_vol(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'EZ', store_val3d_EZ, &
                                                   (/sum(offset_poin_vol(0:myrank-1))/), H5_COL)

  call h5_close_group()
  call h5_close_file_p()

  deallocate(store_val3d_NN,store_val3d_EE,store_val3d_ZZ, &
             store_val3d_NE,store_val3d_NZ,store_val3d_EZ)


#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = size(eps_trace_over_3_crust_mantle,kind=4)
  idummy = size(epsilondev_xx_crust_mantle,kind=4)
  idummy = size(epsilondev_yy_crust_mantle,kind=4)
  idummy = size(epsilondev_xy_crust_mantle,kind=4)
  idummy = size(epsilondev_xz_crust_mantle,kind=4)
  idummy = size(epsilondev_yz_crust_mantle,kind=4)

  print *, 'Error: HDF5 is not enabled in this version of the code.'
  print *, 'Please recompile with the HDF5 option enabled with --with-hdf5'
  stop

#endif

  end subroutine write_movie_volume_strains_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_volume_divcurl_hdf5(vnspec_eps_cm,eps_trace_over_3_crust_mantle, &
                                       div_displ_outer_core, &
                                       accel_outer_core,kappavstore_outer_core,rhostore_outer_core,ibool_outer_core, &
                                       vnspec_eps_ic,eps_trace_over_3_inner_core, &
                                       vnspec_cm,epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                       epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                       vnspec_ic,epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
                                       epsilondev_xz_inner_core,epsilondev_yz_inner_core)

! outputs divergence and curl: MOVIE_VOLUME_TYPE == 4

  use constants_solver

#ifdef USE_HDF5
  use shared_parameters, only: OUTPUT_FILES
  use specfem_par, only: it
  use specfem_par_crustmantle, only: ibool_crust_mantle
  use specfem_par_innercore, only: ibool_inner_core
  use specfem_par_movie_hdf5
#endif

  implicit none

  integer,intent(in) :: vnspec_eps_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_eps_cm) :: eps_trace_over_3_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_3DMOVIE) :: div_displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: accel_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: rhostore_outer_core,kappavstore_outer_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  integer,intent(in) :: vnspec_eps_ic
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_eps_ic) :: eps_trace_over_3_inner_core

  integer,intent(in) :: vnspec_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_cm) :: &
    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle

  integer,intent(in) :: vnspec_ic
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_ic) :: &
    epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
    epsilondev_xz_inner_core,epsilondev_yz_inner_core

#ifdef USE_HDF5

  ! local parameters
  real(kind=CUSTOM_REAL) :: rhol,kappal
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: div_s_outer_core
  integer :: ispec,iglob,i,j,k,ier
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: tmp_data

  ! checks
  if (vnspec_cm /= NSPEC_CRUST_MANTLE) call exit_MPI(myrank,'Invalid vnspec_cm value for write_movie_volume_divcurl() routine')

  ! create group and datasets
  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  if (myrank == 0) then
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    if (MOVIE_OUTPUT_DIV) then
      call h5_create_dataset_gen_in_group('reg1_div_displ', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('reg2_div_displ', (/npoints_vol_mov_all_proc_oc/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('reg3_div_displ', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
    endif

    if (MOVIE_OUTPUT_CURL) then
      call h5_create_dataset_gen_in_group('crust_mantle_epsdev_displ_xx', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('crust_mantle_epsdev_displ_yy', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('crust_mantle_epsdev_displ_xy', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('crust_mantle_epsdev_displ_xz', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('crust_mantle_epsdev_displ_yz', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)

      call h5_create_dataset_gen_in_group('inner_core_epsdev_displ_xx', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('inner_core_epsdev_displ_yy', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('inner_core_epsdev_displ_xy', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('inner_core_epsdev_displ_xz', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('inner_core_epsdev_displ_yz', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
    endif

    if (MOVIE_OUTPUT_CURLNORM) then
      call h5_create_dataset_gen_in_group('reg1_epsdev_displ_norm', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('reg3_epsdev_displ_norm', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
    endif

    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  if (MOVIE_OUTPUT_DIV) then
    call write_array3dspec_as_1d_hdf5('reg1_div_displ', offset_nspec_vol_cm(myrank), offset_poin_vol_cm(myrank), &
                                      eps_trace_over_3_crust_mantle, sum(offset_poin_vol_cm(0:myrank-1)), ibool_crust_mantle)

    if (NSPEC_OUTER_CORE_3DMOVIE > 1) then
      call write_array3dspec_as_1d_hdf5('reg2_div_displ', offset_nspec_vol_oc(myrank), offset_poin_vol_oc(myrank), &
                                        div_displ_outer_core, sum(offset_poin_vol_oc(0:myrank-1)), ibool_outer_core)
    else
      allocate(div_s_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array div_s_outer_core')
      do ispec = 1, NSPEC_OUTER_CORE
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool_outer_core(i,j,k,ispec)
              rhol = rhostore_outer_core(i,j,k,ispec)
              kappal = kappavstore_outer_core(i,j,k,ispec)
              div_s_outer_core(i,j,k,ispec) = rhol * accel_outer_core(iglob) / kappal
            enddo
          enddo
        enddo
      enddo
      call write_array3dspec_as_1d_hdf5('reg2_div_displ', offset_nspec_vol_oc(myrank), offset_poin_vol_oc(myrank), &
                                        div_s_outer_core, sum(offset_poin_vol_oc(0:myrank-1)), ibool_outer_core)
      deallocate(div_s_outer_core)
    endif

    call write_array3dspec_as_1d_hdf5('reg3_div_displ', offset_nspec_vol_ic(myrank), offset_poin_vol_ic(myrank), &
                                      eps_trace_over_3_inner_core, sum(offset_poin_vol_ic(0:myrank-1)), ibool_inner_core)
  endif

  if (MOVIE_OUTPUT_CURL) then
    call write_array3dspec_as_1d_hdf5('crust_mantle_epsdev_displ_xx', offset_nspec_vol_cm(myrank), offset_poin_vol_cm(myrank), &
                                      epsilondev_xx_crust_mantle, sum(offset_poin_vol_cm(0:myrank-1)), ibool_crust_mantle)
    call write_array3dspec_as_1d_hdf5('crust_mantle_epsdev_displ_yy', offset_nspec_vol_cm(myrank), offset_poin_vol_cm(myrank), &
                                      epsilondev_yy_crust_mantle, sum(offset_poin_vol_cm(0:myrank-1)), ibool_crust_mantle)
    call write_array3dspec_as_1d_hdf5('crust_mantle_epsdev_displ_xy', offset_nspec_vol_cm(myrank), offset_poin_vol_cm(myrank), &
                                      epsilondev_xy_crust_mantle, sum(offset_poin_vol_cm(0:myrank-1)), ibool_crust_mantle)
    call write_array3dspec_as_1d_hdf5('crust_mantle_epsdev_displ_xz', offset_nspec_vol_cm(myrank), offset_poin_vol_cm(myrank), &
                                      epsilondev_xz_crust_mantle, sum(offset_poin_vol_cm(0:myrank-1)), ibool_crust_mantle)
    call write_array3dspec_as_1d_hdf5('crust_mantle_epsdev_displ_yz', offset_nspec_vol_cm(myrank), offset_poin_vol_cm(myrank), &
                                      epsilondev_yz_crust_mantle, sum(offset_poin_vol_cm(0:myrank-1)), ibool_crust_mantle)

    call write_array3dspec_as_1d_hdf5('inner_core_epsdev_displ_xx', offset_nspec_vol_ic(myrank), offset_poin_vol_ic(myrank), &
                                      epsilondev_xx_inner_core, sum(offset_poin_vol_ic(0:myrank-1)), ibool_inner_core)
    call write_array3dspec_as_1d_hdf5('inner_core_epsdev_displ_yy', offset_nspec_vol_ic(myrank), offset_poin_vol_ic(myrank), &
                                      epsilondev_yy_inner_core, sum(offset_poin_vol_ic(0:myrank-1)), ibool_inner_core)
    call write_array3dspec_as_1d_hdf5('inner_core_epsdev_displ_xy', offset_nspec_vol_ic(myrank), offset_poin_vol_ic(myrank), &
                                      epsilondev_xy_inner_core, sum(offset_poin_vol_ic(0:myrank-1)), ibool_inner_core)
    call write_array3dspec_as_1d_hdf5('inner_core_epsdev_displ_xz', offset_nspec_vol_ic(myrank), offset_poin_vol_ic(myrank), &
                                      epsilondev_xz_inner_core, sum(offset_poin_vol_ic(0:myrank-1)), ibool_inner_core)
    call write_array3dspec_as_1d_hdf5('inner_core_epsdev_displ_yz', offset_nspec_vol_ic(myrank), offset_poin_vol_ic(myrank), &
                                      epsilondev_yz_inner_core, sum(offset_poin_vol_ic(0:myrank-1)), ibool_inner_core)
  endif

  if (MOVIE_OUTPUT_CURLNORM) then
    allocate(tmp_data(NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
    ! Frobenius norm
    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            tmp_data(iglob) = sqrt( epsilondev_xx_crust_mantle(i,j,k,ispec)**2 &
                                  + epsilondev_yy_crust_mantle(i,j,k,ispec)**2 &
                                  + epsilondev_xy_crust_mantle(i,j,k,ispec)**2 &
                                  + epsilondev_xz_crust_mantle(i,j,k,ispec)**2 &
                                  + epsilondev_yz_crust_mantle(i,j,k,ispec)**2)
          enddo
        enddo
      enddo
    enddo
    call write_array3dspec_as_1d_hdf5('reg1_epsdev_displ_norm', offset_nspec_vol_cm(myrank), offset_poin_vol_cm(myrank), &
                                      tmp_data, sum(offset_poin_vol_cm(0:myrank-1)), ibool_crust_mantle)
    deallocate(tmp_data)

    ! Frobenius norm
    allocate(tmp_data(NGLOB_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            tmp_data(iglob) = sqrt( epsilondev_xx_inner_core(i,j,k,ispec)**2 &
                                  + epsilondev_yy_inner_core(i,j,k,ispec)**2 &
                                  + epsilondev_xy_inner_core(i,j,k,ispec)**2 &
                                  + epsilondev_xz_inner_core(i,j,k,ispec)**2 &
                                  + epsilondev_yz_inner_core(i,j,k,ispec)**2)
          enddo
        enddo
      enddo
    enddo
    call write_array3dspec_as_1d_hdf5('reg3_epsdev_displ_norm', offset_nspec_vol_ic(myrank), offset_poin_vol_ic(myrank), &
                                      tmp_data, sum(offset_poin_vol_ic(0:myrank-1)), ibool_inner_core)
    deallocate(tmp_data)
  endif

  call h5_close_group()
  call h5_close_file()

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = size(eps_trace_over_3_crust_mantle,kind=4)
  idummy = size(div_displ_outer_core,kind=4)

  idummy = size(accel_outer_core,kind=4)
  idummy = size(kappavstore_outer_core,kind=4)
  idummy = size(rhostore_outer_core,kind=4)
  idummy = size(ibool_outer_core,kind=4)

  idummy = size(eps_trace_over_3_inner_core,kind=4)
  idummy = size(epsilondev_xx_crust_mantle,kind=4)
  idummy = size(epsilondev_yy_crust_mantle,kind=4)
  idummy = size(epsilondev_xy_crust_mantle,kind=4)
  idummy = size(epsilondev_xz_crust_mantle,kind=4)
  idummy = size(epsilondev_yz_crust_mantle,kind=4)

  idummy = size(epsilondev_xx_inner_core,kind=4)
  idummy = size(epsilondev_yy_inner_core,kind=4)
  idummy = size(epsilondev_xy_inner_core,kind=4)
  idummy = size(epsilondev_xz_inner_core,kind=4)
  idummy = size(epsilondev_yz_inner_core,kind=4)

  print *, 'Error: HDF5 is not enabled in this version of the code.'
  print *, 'Please recompile with the HDF5 option enabled with --with-hdf5'
  stop

#endif

  end subroutine write_movie_volume_divcurl_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_volume_vector_hdf5(npoints_3dmovie, &
                                       ibool_crust_mantle,vector_crust_mantle, &
                                       scalingval,mask_3dmovie,nu_3dmovie)

! outputs displacement/velocity: MOVIE_VOLUME_TYPE == 5 / 6

  use constants_solver
#ifdef USE_HDF5
  use shared_parameters, only: OUTPUT_FILES,MOVIE_VOLUME_TYPE,MOVIE_COARSE,H5_COL
  use specfem_par, only: it
  use specfem_par_movie_hdf5
#endif

  implicit none

  ! input
  integer :: npoints_3dmovie
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  ! displacement or velocity array
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: vector_crust_mantle
  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE) :: mask_3dmovie

  double precision :: scalingval
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,npoints_3dmovie) :: nu_3dmovie

#ifdef USE_HDF5

  ! local variables
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val3d_N,store_val3d_E,store_val3d_Z
  real(kind=CUSTOM_REAL), dimension(NDIM) :: vector_local,vector_local_new

  integer :: ipoints_3dmovie,i,j,k,ispec,iNIT,iglob,ier
  character(len=2) :: movie_prefix

  ! check
  if (NDIM /= 3) call exit_MPI(myrank,'write_movie_volume requires NDIM = 3')

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  ! allocates arrays
  allocate(store_val3d_N(npoints_3dmovie), &
           store_val3d_E(npoints_3dmovie), &
           store_val3d_Z(npoints_3dmovie), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating store_val3d_N,.. movie arrays')

  if (MOVIE_VOLUME_TYPE == 5) then
    movie_prefix='DI' ! displacement
  else if (MOVIE_VOLUME_TYPE == 6) then
    movie_prefix='VE' ! velocity
  endif

  if (MOVIE_COARSE) then
   iNIT = NGLLX-1
  else
   iNIT = 1
  endif

  ipoints_3dmovie = 0

  ! stores field in crust/mantle region
  do ispec = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ,iNIT
      do j = 1,NGLLY,iNIT
        do i = 1,NGLLX,iNIT
          if (mask_3dmovie(i,j,k,ispec)) then
            ipoints_3dmovie = ipoints_3dmovie + 1
            iglob = ibool_crust_mantle(i,j,k,ispec)

            ! dimensionalizes field by scaling
            vector_local(:) = vector_crust_mantle(:,iglob)*real(scalingval,kind=CUSTOM_REAL)

            ! rotate eps_loc to spherical coordinates
            vector_local_new(:) = matmul(nu_3dmovie(:,:,ipoints_3dmovie), vector_local(:))

            ! stores field
            store_val3d_N(ipoints_3dmovie) = vector_local_new(1)
            store_val3d_E(ipoints_3dmovie) = vector_local_new(2)
            store_val3d_Z(ipoints_3dmovie) = vector_local_new(3)
          endif
        enddo
      enddo
   enddo
  enddo
  close(IOUT)

  ! checks number of processed points
  if (ipoints_3dmovie /= npoints_3dmovie) stop 'did not find the right number of points for 3D movie'

  ! create group and datasets
  if (myrank == 0) then
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'N', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'E', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'Z', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)

    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'N', store_val3d_N(1:npoints_3dmovie), &
                                                   (/offset_poin_vol(0:myrank-1)/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'E', store_val3d_E(1:npoints_3dmovie), &
                                                   (/offset_poin_vol(0:myrank-1)/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'Z', store_val3d_Z(1:npoints_3dmovie), &
                                                   (/offset_poin_vol(0:myrank-1)/), H5_COL)

  call h5_close_group()
  call h5_close_file()

  deallocate(store_val3d_N,store_val3d_E,store_val3d_Z)

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy
  double precision :: d_dummy

  idummy = size(ibool_crust_mantle,kind=4)
  idummy = size(mask_3dmovie,kind=4)
  idummy = size(nu_3dmovie,kind=4)
  idummy = size(vector_crust_mantle,kind=4)

  d_dummy = scalingval

  print *, 'Error: HDF5 is not enabled in this version of the code.'
  print *, 'Please recompile with the HDF5 option enabled with --with-hdf5'
  stop

#endif

  end subroutine write_movie_volume_vector_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_volume_displnorm_hdf5(displ_crust_mantle,displ_inner_core,displ_outer_core, &
                                         ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

! outputs norm of displacement: MOVIE_VOLUME_TYPE == 7

  use constants_solver

#ifdef USE_HDF5
  use shared_parameters, only: OUTPUT_FILES
  use specfem_par, only: it, scale_displ
  use specfem_par_movie_hdf5
#endif

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: displ_outer_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

#ifdef USE_HDF5

  ! local parameters
  integer :: ispec,iglob,i,j,k,ier
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: tmp_data

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  ! create group and datasets
  if (myrank == 0) then
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    if (OUTPUT_CRUST_MANTLE) then
      call h5_create_dataset_gen_in_group('reg1_displ', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
    endif
    if (OUTPUT_OUTER_CORE) then
      call h5_create_dataset_gen_in_group('reg2_displ', (/npoints_vol_mov_all_proc_oc/), 1, CUSTOM_REAL)
    endif
    if (OUTPUT_INNER_CORE) then
      call h5_create_dataset_gen_in_group('reg3_displ', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
    endif

    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  ! outputs norm of displacement
  if (OUTPUT_CRUST_MANTLE) then
    ! crust mantle
    ! these binary arrays can be converted into mesh format using the utility ./bin/xcombine_vol_data
    allocate(tmp_data(NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! norm
            tmp_data(iglob) = real(scale_displ,kind=CUSTOM_REAL) * sqrt( displ_crust_mantle(1,iglob)**2 &
                                          + displ_crust_mantle(2,iglob)**2 &
                                          + displ_crust_mantle(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call write_array3dspec_as_1d_hdf5('reg1_displ', offset_nspec_vol_cm(myrank-1), offset_poin_vol_cm(myrank-1), &
                                      tmp_data, sum(offset_poin_vol_cm(0:myrank-1)), ibool_crust_mantle)

    deallocate(tmp_data)
  endif

  if (OUTPUT_OUTER_CORE) then
    ! outer core
    allocate(tmp_data(NGLOB_OUTER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)
            ! norm
            ! note: disp_outer_core is potential, this just outputs the potential,
            !          not the actual displacement u = grad(rho * Chi) / rho
            tmp_data(iglob) = abs(displ_outer_core(iglob))
          enddo
        enddo
      enddo
    enddo

    call write_array3dspec_as_1d_hdf5('reg2_displ', offset_nspec_vol_oc(myrank-1), offset_poin_vol_oc(myrank-1), &
                                      tmp_data, sum(offset_poin_vol_oc(0:myrank-1)), ibool_outer_core)

    deallocate(tmp_data)
  endif

  if (OUTPUT_INNER_CORE) then
    ! inner core
    allocate(tmp_data(NGLOB_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm
            tmp_data(iglob) = real(scale_displ,kind=CUSTOM_REAL) * sqrt( displ_inner_core(1,iglob)**2 &
                                          + displ_inner_core(2,iglob)**2 &
                                          + displ_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call write_array3dspec_as_1d_hdf5('reg3_displ', offset_nspec_vol_ic(myrank-1), offset_poin_vol_ic(myrank-1), &
                                      tmp_data, sum(offset_poin_vol_ic(0:myrank-1)), ibool_inner_core)

    deallocate(tmp_data)
  endif

  call h5_close_group()
  call h5_close_file_p()


#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = size(displ_crust_mantle,kind=4)
  idummy = size(displ_inner_core,kind=4)
  idummy = size(displ_outer_core,kind=4)
  idummy = size(ibool_crust_mantle,kind=4)
  idummy = size(ibool_inner_core,kind=4)
  idummy = size(ibool_outer_core,kind=4)

  print *, 'Error: HDF5 is not enabled in this version of the code.'
  print *, 'Please recompile with the HDF5 option enabled with --with-hdf5'
  stop

#endif

  end subroutine write_movie_volume_displnorm_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_volume_velnorm_hdf5(veloc_crust_mantle,veloc_inner_core,veloc_outer_core, &
                                       ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

! outputs norm of velocity: MOVIE_VOLUME_TYPE == 8

  use constants_solver

#ifdef USE_HDF5
  use shared_parameters, only: OUTPUT_FILES
  use specfem_par, only: it, scale_veloc
  use specfem_par_movie_hdf5
#endif

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: veloc_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: veloc_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: veloc_inner_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

#ifdef USE_HDF5

  ! local parameters
  integer :: ispec,iglob,i,j,k,ier
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: tmp_data

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  ! create group and datasets
  if (myrank == 0) then
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    if (OUTPUT_CRUST_MANTLE) then
      call h5_create_dataset_gen_in_group('reg1_veloc', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
    endif
    if (OUTPUT_OUTER_CORE) then
      call h5_create_dataset_gen_in_group('reg2_veloc', (/npoints_vol_mov_all_proc_oc/), 1, CUSTOM_REAL)
    endif
    if (OUTPUT_INNER_CORE) then
      call h5_create_dataset_gen_in_group('reg3_veloc', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
    endif

    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  ! outputs norm of velocity
  if (OUTPUT_CRUST_MANTLE) then
    ! crust mantle
    ! these binary arrays can be converted into mesh format using the utility ./bin/xcombine_vol_data
    allocate(tmp_data(NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! norm of velocity
            tmp_data(iglob) = real(scale_veloc,kind=CUSTOM_REAL) * sqrt( veloc_crust_mantle(1,iglob)**2 &
                                          + veloc_crust_mantle(2,iglob)**2 &
                                          + veloc_crust_mantle(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call write_array3dspec_as_1d_hdf5('reg1_veloc', offset_nspec_vol_cm(myrank-1), offset_poin_vol_cm(myrank-1), &
                                      tmp_data, sum(offset_poin_vol_cm(0:myrank-1)), ibool_crust_mantle)

    deallocate(tmp_data)
  endif

  if (OUTPUT_OUTER_CORE) then
    ! outer core
    allocate(tmp_data(NGLOB_OUTER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)
            ! norm of velocity
            ! note: this outputs only the first time derivative of the potential,
            !          not the actual velocity v = grad(Chi_dot)
            tmp_data(iglob) = abs(veloc_outer_core(iglob))
          enddo
        enddo
      enddo
    enddo

    call write_array3dspec_as_1d_hdf5('reg2_veloc', offset_nspec_vol_oc(myrank-1), offset_poin_vol_oc(myrank-1), &
                                      tmp_data, sum(offset_poin_vol_oc(0:myrank-1)), ibool_outer_core)

    deallocate(tmp_data)
  endif

  if (OUTPUT_INNER_CORE) then
    ! inner core
    allocate(tmp_data(NGLOB_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm of velocity
            tmp_data(iglob) = real(scale_veloc,kind=CUSTOM_REAL) * sqrt( veloc_inner_core(1,iglob)**2 &
                                                                       + veloc_inner_core(2,iglob)**2 &
                                                                       + veloc_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call write_array3dspec_as_1d_hdf5('reg3_veloc', offset_nspec_vol_ic(myrank-1), offset_poin_vol_ic(myrank-1), &
                                      tmp_data, sum(offset_poin_vol_ic(0:myrank-1)), ibool_inner_core)

    deallocate(tmp_data)
  endif

  call h5_close_group()
  call h5_close_file_p()

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = size(veloc_crust_mantle,kind=4)
  idummy = size(veloc_inner_core,kind=4)
  idummy = size(veloc_outer_core,kind=4)
  idummy = size(ibool_crust_mantle,kind=4)
  idummy = size(ibool_inner_core,kind=4)
  idummy = size(ibool_outer_core,kind=4)

  print *, 'Error: HDF5 is not enabled in this version of the code.'
  print *, 'Please recompile with the HDF5 option enabled with --with-hdf5'
  stop

#endif

  end subroutine write_movie_volume_velnorm_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_volume_accelnorm_hdf5(accel_crust_mantle,accel_inner_core,accel_outer_core, &
                                         ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

! outputs norm of acceleration: MOVIE_VOLUME_TYPE == 1

  use constants_solver

#ifdef USE_HDF5
  use shared_parameters, only: OUTPUT_FILES
  use specfem_par, only: it, scale_t_inv,scale_veloc
  use specfem_par_movie_hdf5
#endif

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: accel_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: accel_outer_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

#ifdef USE_HDF5
  ! local parameters
  integer :: ispec,iglob,i,j,k,ier
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: tmp_data
  real(kind=CUSTOM_REAL) :: scale_accel

  ! dimensionalized scaling
  scale_accel = real(scale_veloc * scale_t_inv,kind=CUSTOM_REAL)

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  ! create group and datasets
  if (myrank == 0) then
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    if (OUTPUT_CRUST_MANTLE) then
      call h5_create_dataset_gen_in_group('reg1_accel', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
    endif
    if (OUTPUT_OUTER_CORE) then
      call h5_create_dataset_gen_in_group('reg2_accel', (/npoints_vol_mov_all_proc_oc/), 1, CUSTOM_REAL)
    endif
    if (OUTPUT_INNER_CORE) then
      call h5_create_dataset_gen_in_group('reg3_accel', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
    endif

    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  ! outputs norm of acceleration
  if (OUTPUT_CRUST_MANTLE) then
    ! acceleration
    ! these binary arrays can be converted into mesh format using the utility ./bin/xcombine_vol_data
    allocate(tmp_data(NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! norm
            tmp_data(iglob) = scale_accel * sqrt( accel_crust_mantle(1,iglob)**2 &
                                                + accel_crust_mantle(2,iglob)**2 &
                                                + accel_crust_mantle(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call write_array3dspec_as_1d_hdf5('reg1_accel', offset_nspec_vol_cm(myrank-1), offset_poin_vol_cm(myrank-1), &
                                      tmp_data, sum(offset_poin_vol_cm(0:myrank-1)), ibool_crust_mantle)

    deallocate(tmp_data)
  endif

  if (OUTPUT_OUTER_CORE) then
    ! outer core acceleration
    allocate(tmp_data(NGLOB_OUTER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)
            ! norm
            ! note: this outputs only the second time derivative of the potential,
            !          not the actual acceleration or pressure p = - rho * Chi_dot_dot
            tmp_data(iglob) = abs(accel_outer_core(iglob))
          enddo
        enddo
      enddo
    enddo

    call write_array3dspec_as_1d_hdf5('reg2_accel', offset_nspec_vol_oc(myrank-1), offset_poin_vol_oc(myrank-1), &
                                      tmp_data, sum(offset_poin_vol_oc(0:myrank-1)), ibool_outer_core)

    deallocate(tmp_data)
  endif

  if (OUTPUT_INNER_CORE) then
    ! inner core
    allocate(tmp_data(NGLOB_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm of acceleration
            tmp_data(iglob) = scale_accel * sqrt( accel_inner_core(1,iglob)**2 &
                                                + accel_inner_core(2,iglob)**2 &
                                                + accel_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call write_array3dspec_as_1d_hdf5('reg3_accel', offset_nspec_vol_ic(myrank-1), offset_poin_vol_ic(myrank-1), &
                                      tmp_data, sum(offset_poin_vol_ic(0:myrank-1)), ibool_inner_core)

    deallocate(tmp_data)
  endif

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = size(accel_crust_mantle,kind=4)
  idummy = size(accel_inner_core,kind=4)
  idummy = size(accel_outer_core,kind=4)
  idummy = size(ibool_crust_mantle,kind=4)
  idummy = size(ibool_inner_core,kind=4)
  idummy = size(ibool_outer_core,kind=4)

  print *, 'Error: HDF5 is not enabled in this version of the code.'
  print *, 'Please recompile with the HDF5 option enabled with --with-hdf5'
  stop

#endif

  end subroutine write_movie_volume_accelnorm_hdf5

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5

  subroutine get_conn_for_movie(elm_conn,offset,iNIT,nelems_3dmovie,npoints_3dmovie, &
                                nelems_in_this_region, &
                                num_ibool_3dmovie,mask_ibool_3dmovie,ibool_of_the_section)

  use specfem_par
  !use specfem_par_crustmantle, only: ibool_crust_mantle
  use constants_solver


  implicit none

  integer, intent(in) :: offset ! node id offset (starting global element id of each proc)
  integer, intent(in) :: iNIT
  integer, intent(in) :: nelems_3dmovie, npoints_3dmovie, nelems_in_this_region
  integer, dimension(npoints_3dmovie),intent(in) :: num_ibool_3dmovie
  logical, dimension(npoints_3dmovie),intent(in) :: mask_ibool_3dmovie
  integer, dimension(NGLLX,NGLLY,NGLLZ,nelems_in_this_region),intent(in) :: ibool_of_the_section
  integer, dimension(9,nelems_3dmovie), intent(out) :: elm_conn
  ! local parameters
  integer :: ispec
  integer,parameter :: cell_type = 9
  integer :: iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer :: n1,n2,n3,n4,n5,n6,n7,n8
  integer :: i,j,k
  integer :: ispecele, iglob_center

  ispecele = 0
  do ispec = 1, nelems_in_this_region

    ! checks center of element for movie flag
    iglob_center = ibool_of_the_section((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)

    ! checks if movie element
    if (mask_ibool_3dmovie(iglob_center)) then

      do k = 1,NGLLZ-1,iNIT
        do j = 1,NGLLY-1,iNIT
          do i = 1,NGLLX-1,iNIT
            ! this element is in the movie region
            ispecele = ispecele+1

            ! defines corners of a vtk element
            iglob1 = ibool_of_the_section(i,j,k,ispec)
            iglob2 = ibool_of_the_section(i+iNIT,j,k,ispec)
            iglob3 = ibool_of_the_section(i+iNIT,j+iNIT,k,ispec)
            iglob4 = ibool_of_the_section(i,j+iNIT,k,ispec)
            iglob5 = ibool_of_the_section(i,j,k+iNIT,ispec)
            iglob6 = ibool_of_the_section(i+iNIT,j,k+iNIT,ispec)
            iglob7 = ibool_of_the_section(i+iNIT,j+iNIT,k+iNIT,ispec)
            iglob8 = ibool_of_the_section(i,j+iNIT,k+iNIT,ispec)

            ! vtk indexing starts at 0 -> adds minus 1
            n1 = num_ibool_3dmovie(iglob1)-1
            n2 = num_ibool_3dmovie(iglob2)-1
            n3 = num_ibool_3dmovie(iglob3)-1
            n4 = num_ibool_3dmovie(iglob4)-1
            n5 = num_ibool_3dmovie(iglob5)-1
            n6 = num_ibool_3dmovie(iglob6)-1
            n7 = num_ibool_3dmovie(iglob7)-1
            n8 = num_ibool_3dmovie(iglob8)-1

            elm_conn(1, ispecele)  = cell_type
            elm_conn(2, ispecele)  = n1 + offset   ! node id starts 0 in xdmf rule
            elm_conn(3, ispecele)  = n2 + offset
            elm_conn(4, ispecele)  = n3 + offset
            elm_conn(5, ispecele)  = n4 + offset
            elm_conn(6, ispecele)  = n5 + offset
            elm_conn(7, ispecele)  = n6 + offset
            elm_conn(8, ispecele)  = n7 + offset
            elm_conn(9, ispecele)  = n8 + offset

            ! checks indices
            if (n1 < 0 .or. n2 < 0 .or. n3 < 0 .or. n4 < 0 .or. n5 < 0 .or. n6 < 0 .or. n7 < 0 .or. n8 < 0) then
              print *,'Error: movie element ',ispec,ispecele,'has invalid node index:',n1,n2,n3,n4,n5,n6,n7,n8
              call exit_mpi(myrank,'Error invalid movie element node index')
            endif

          enddo !i
        enddo !j
      enddo !k
    endif

  enddo !ispec

  ! check if ispecele is consistent with nelems_3dmovie
  if (ispecele /= nelems_3dmovie) then
    print *,'Error: number of movie elements is not consistent with nelems_3dmovie'
    print *,'ispecele = ',ispecele,'nelems_3dmovie = ',nelems_3dmovie
    call exit_mpi(myrank,'Error number of movie elements is not consistent with nelems_3dmovie')
  endif

  end subroutine get_conn_for_movie

#endif

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5

  subroutine elm2node_base(array_3dspec, array_1dmovie, nelms, npoints, ibool_of_the_section)

  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  integer, intent(in) :: nelms, npoints
  integer, dimension(NGLLX,NGLLY,NGLLZ,nelms), intent(in) :: ibool_of_the_section
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nelms), intent(in) :: array_3dspec
  real(kind=CUSTOM_REAL), dimension(npoints), intent(inout) :: array_1dmovie

  ! local parameters
  integer :: i,j,k,ispec,iglob

  ! convert 3d array to 1d array
  do ispec = 1, nelms

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_of_the_section(i,j,k,ispec)
          array_1dmovie(iglob) = array_3dspec(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  end subroutine elm2node_base

#endif

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5

  subroutine write_array3dspec_as_1d_hdf5(dset_name, nelms, npoints, array_3dspec, offset1d, ibool_of_the_section)

  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  character(len=*), intent(in) :: dset_name
  integer, intent(in) :: nelms, npoints
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nelms), intent(in) :: array_3dspec
  integer, intent(in) :: offset1d
  integer, dimension(NGLLX,NGLLY,NGLLZ,nelms), intent(in) :: ibool_of_the_section

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(npoints) :: array_1dmovie

  ! convert 3d array to 1d array
  call elm2node_base(array_3dspec, array_1dmovie, nelms, npoints, ibool_of_the_section)

  ! write 1d array to hdf5
  call h5_write_dataset_collect_hyperslab_in_group(dset_name, array_1dmovie, (/offset1d/), H5_COL)

  end subroutine write_array3dspec_as_1d_hdf5

#endif

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5

  subroutine write_xdmf_vol_hdf5_one_data(fname_h5, attr_name, dset_name, len_data, target_unit, it_str, value_on_node)

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    character(len=*), intent(in) :: fname_h5
    character(len=*), intent(in) :: attr_name
    character(len=*), intent(in) :: dset_name
    integer, intent(in) :: len_data
    integer, intent(in) :: target_unit
    character(len=*), intent(in) :: it_str
    logical, intent(in) :: value_on_node ! values are defined on nodes or elements

    character(len=20) :: center_str

    if (value_on_node) then
      center_str = 'Node' ! len data should be the number of nodes
    else
      center_str = 'Cell' ! len data should be the number of elements
    endif

    ! write a header for single data
    write(target_unit,*) '<Attribute Name="'//trim(attr_name)//'" AttributeType="Scalar" Center="'//trim(center_str)//'">'
    write(target_unit,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(len_data))//'">'
    write(target_unit,*) '      '//trim(fname_h5)//':/it_'//trim(it_str)//'/'//trim(dset_name)
    write(target_unit,*) '</DataItem>'
    write(target_unit,*) '</Attribute>'

  end subroutine write_xdmf_vol_hdf5_one_data

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_xdmf_vol_hdf5_header(nelems, nglobs, fname_h5_data_vol_xdmf, target_unit, region_flag)

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    integer, intent(in) :: nelems, nglobs
    character(len=*), intent(in) :: fname_h5_data_vol_xdmf
    integer, intent(in) :: target_unit
    integer, intent(in) :: region_flag ! 1: crust mantle, 2: outer core, 3: inner core

    character(len=64) :: nelm_str, nglo_str, elemconn_str, x_str, y_str, z_str

    if (region_flag == 1) then
      elemconn_str = 'elm_conn'
      x_str = 'x'
      y_str = 'y'
      z_str = 'z'
    else if (region_flag == 2) then
      elemconn_str = 'elm_conn_cm'
      x_str = 'x_cm'
      y_str = 'y_cm'
      z_str = 'z_cm'
    else if (region_flag == 3) then
      elemconn_str = 'elm_conn_oc'
      x_str = 'x_oc'
      y_str = 'y_oc'
      z_str = 'z_oc'
    else if (region_flag == 4) then
      elemconn_str = 'elm_conn_ic'
      x_str = 'x_ic'
      y_str = 'y_ic'
      z_str = 'z_ic'
    else
      print *,'Error: invalid region_flag in write_xdmf_vol_hdf5_header'
      call exit_mpi(myrank,'Error invalid region_flag in write_xdmf_vol_hdf5_header')
    endif

    ! convert integer to string
    nelm_str = i2c(nelems)
    nglo_str = i2c(nglobs)

    ! definition of topology and geometry
    ! refer only control nodes (8 or 27) as a coarse output
    ! data array need to be extracted from full data array on GLL points
    write(target_unit,'(a)') '<?xml version="1.0" ?>'
    write(target_unit,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(target_unit,*) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">'
    write(target_unit,*) '<Domain name="mesh">'

    ! loop for writing information of mesh partitions
    write(target_unit,*) '<Topology TopologyType="Mixed" NumberOfElements="'//trim(nelm_str)//'">'
    write(target_unit,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Int" Precision="4" Dimensions="'&
                         //trim(nelm_str)//' 9">'
    write(target_unit,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/'//trim(elemconn_str)
    write(target_unit,*) '</DataItem>'
    write(target_unit,*) '</Topology>'
    write(target_unit,*) '<Geometry GeometryType="X_Y_Z">'
    write(target_unit,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
    write(target_unit,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/'//trim(x_str)
    write(target_unit,*) '</DataItem>'
    write(target_unit,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
    write(target_unit,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/'//trim(y_str)
    write(target_unit,*) '</DataItem>'
    write(target_unit,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
    write(target_unit,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/'//trim(z_str)
    write(target_unit,*) '</DataItem>'
    write(target_unit,*) '</Geometry>'

    write(target_unit,*) '<Grid Name="time_col" GridType="Collection" CollectionType="Temporal">'

  end subroutine write_xdmf_vol_hdf5_header

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_xdmf_vol_hdf5_footer(target_unit)

    implicit none

    integer, intent(in) :: target_unit

    ! file finish
    write(target_unit,*) '</Grid>'
    write(target_unit,*) '</Domain>'
    write(target_unit,*) '</Xdmf>'

  end subroutine write_xdmf_vol_hdf5_footer

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_xdmf_vol_hdf5(npoints_3dmovie,    nelems_3dmovie, &
                                 npoints_3dmovie_cm, nelems_3dmovie_cm, &
                                 npoints_3dmovie_oc, nelems_3dmovie_oc, &
                                 npoints_3dmovie_ic, nelems_3dmovie_ic)

  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  integer, intent(in) :: npoints_3dmovie,    nelems_3dmovie
  integer, intent(in) :: npoints_3dmovie_cm, nelems_3dmovie_cm
  integer, intent(in) :: npoints_3dmovie_oc, nelems_3dmovie_oc
  integer, intent(in) :: npoints_3dmovie_ic, nelems_3dmovie_ic

  ! local parameters
  integer                       :: i, ii
  character(len=20)             :: it_str, movie_prefix
  character(len=MAX_STRING_LEN) :: fname_xdmf_vol, fname_xdmf_vol_oc, fname_xdmf_vol_ic
  character(len=MAX_STRING_LEN) :: fname_h5_data_vol_xdmf

  ! checks if anything do, only main process writes out xdmf file
  if (myrank /= 0) return

  !
  ! write out the crust mantle xdmf file (for strain and vector)
  !
  if (output_sv) then

    fname_xdmf_vol = trim(OUTPUT_FILES) // "/movie_volume.xmf"
    fname_h5_data_vol_xdmf = "./movie_volume.h5"  ! relative to movie_volume.xmf file

    ! open xdmf file
    open(unit=xdmf_vol, file=trim(fname_xdmf_vol), recl=256)

    call write_xdmf_vol_hdf5_header(nspec_vol_mov_all_proc, npoints_vol_mov_all_proc, fname_h5_data_vol_xdmf, xdmf_vol, 1)


    do i = 1, int(NSTEP/NTSTEP_BETWEEN_FRAMES)

      ii = i*NTSTEP_BETWEEN_FRAMES
      !write(it_str, "(i6.6)") ii
      it_str = i2c(ii)

      write(xdmf_vol,*) '<Grid Name="vol_mov" GridType="Uniform">'
      write(xdmf_vol,*) '<Time Value="'//trim(r2c(sngl((ii-1)*DT-t0)))//'" />'
      write(xdmf_vol,*) '<Topology Reference="/Xdmf/Domain/Topology" />'
      write(xdmf_vol,*) '<Geometry Reference="/Xdmf/Domain/Geometry" />'

      ! write headers for each dataset

      ! volume strain
      if (MOVIE_VOLUME_TYPE == 1 .or. MOVIE_VOLUME_TYPE == 2 .or. MOVIE_VOLUME_TYPE == 3) then
        if (MOVIE_VOLUME_TYPE == 1) then
          movie_prefix = 'E' ! strain
        else if (MOVIE_VOLUME_TYPE == 2) then
          movie_prefix = 'S' ! time integral of strain
        else if (MOVIE_VOLUME_TYPE == 3) then
          movie_prefix = 'P' ! potency, or itegral of strain x \mu
        endif

        ! movie_prefix/NN,EE,ZZ,NE,NZ,EZ
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, trim(movie_prefix)//'_NN', trim(movie_prefix)//'NN', &
                                          npoints_vol_mov_all_proc, xdmf_vol, it_str, .true.)
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, trim(movie_prefix)//'_EE', trim(movie_prefix)//'EE', &
                                          npoints_vol_mov_all_proc, xdmf_vol, it_str, .true.)
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, trim(movie_prefix)//'_ZZ', trim(movie_prefix)//'ZZ', &
                                          npoints_vol_mov_all_proc, xdmf_vol, it_str, .true.)
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, trim(movie_prefix)//'_NE', trim(movie_prefix)//'NE', &
                                          npoints_vol_mov_all_proc, xdmf_vol, it_str, .true.)
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, trim(movie_prefix)//'_NZ', trim(movie_prefix)//'NZ', &
                                          npoints_vol_mov_all_proc, xdmf_vol, it_str, .true.)
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, trim(movie_prefix)//'_EZ', trim(movie_prefix)//'EZ', &
                                          npoints_vol_mov_all_proc, xdmf_vol, it_str, .true.)

      ! volume vector
      else if (MOVIE_VOLUME_TYPE == 5 .or. MOVIE_VOLUME_TYPE == 6) then
        if (MOVIE_VOLUME_TYPE == 5) then
          movie_prefix = 'DI' ! displacement
        else if (MOVIE_VOLUME_TYPE == 6) then
          movie_prefix = 'VE' ! velocity
        endif

        ! movie_prefix/N,E,Z
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, trim(movie_prefix)//'_N', trim(movie_prefix)//'N', &
                                          npoints_vol_mov_all_proc, xdmf_vol, it_str, .true.)
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, trim(movie_prefix)//'_E', trim(movie_prefix)//'E', &
                                          npoints_vol_mov_all_proc, xdmf_vol, it_str, .true.)
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, trim(movie_prefix)//'_Z', trim(movie_prefix)//'Z', &
                                          npoints_vol_mov_all_proc, xdmf_vol, it_str, .true.)

      endif

      write(xdmf_vol,*) '</Grid>'

    enddo

    call write_xdmf_vol_hdf5_footer(xdmf_vol)

    ! close xdmf file
    close(xdmf_vol)

  endif ! output_sv

  !
  ! write out the crust and mantle xdmf file (for strain and vector)
  !
  if (output_cm) then

    fname_xdmf_vol = trim(OUTPUT_FILES) // "/movie_volume_cm.xmf"
    fname_h5_data_vol_xdmf = "./movie_volume.h5"  ! relative to movie_volume_cm.xmf file

    ! open xdmf file
    open(unit=xdmf_vol, file=trim(fname_xdmf_vol), recl=256)

    call write_xdmf_vol_hdf5_header(nspec_vol_mov_all_proc_cm_conn, npoints_vol_mov_all_proc_cm, &
                                    fname_h5_data_vol_xdmf, xdmf_vol, 2)

    do i = 1, int(NSTEP/NTSTEP_BETWEEN_FRAMES)

      ii = i*NTSTEP_BETWEEN_FRAMES
      !write(it_str, "(i6.6)") ii
      it_str = i2c(ii)

      write(xdmf_vol,*) '<Grid Name="vol_mov" GridType="Uniform">'
      write(xdmf_vol,*) '<Time Value="'//trim(r2c(sngl((ii-1)*DT-t0)))//'" />'
      write(xdmf_vol,*) '<Topology Reference="/Xdmf/Domain/Topology" />'
      write(xdmf_vol,*) '<Geometry Reference="/Xdmf/Domain/Geometry" />'

      ! write headers for each dataset
      ! volume divcurl (div)
      if (MOVIE_VOLUME_TYPE == 4 .and. MOVIE_OUTPUT_DIV) then
        ! reg1_div_displ
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'reg1_div_displ', 'reg1_div_displ', &
                                          npoints_vol_mov_all_proc_cm, xdmf_vol, it_str, .true.) ! value on element
      else if (MOVIE_VOLUME_TYPE == 4 .and. MOVIE_OUTPUT_CURL) then
        ! curst_mantle_epsdev_disple_xx,yy,xy,xz,yz
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'crust_mantle_epsdev_displ_xx', 'crust_mantle_epsdev_displ_xx', &
                                          npoints_vol_mov_all_proc_cm, xdmf_vol, it_str, .true.) ! value on element
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'crust_mantle_epsdev_displ_yy', 'crust_mantle_epsdev_displ_yy', &
                                          npoints_vol_mov_all_proc_cm, xdmf_vol, it_str, .true.) ! value on element
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'crust_mantle_epsdev_displ_xy', 'crust_mantle_epsdev_displ_xy', &
                                          npoints_vol_mov_all_proc_cm, xdmf_vol, it_str, .true.) ! value on element
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'crust_mantle_epsdev_displ_xz', 'crust_mantle_epsdev_displ_xz', &
                                          npoints_vol_mov_all_proc_cm, xdmf_vol, it_str, .true.) ! value on element
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'crust_mantle_epsdev_displ_yz', 'crust_mantle_epsdev_displ_yz', &
                                          npoints_vol_mov_all_proc_cm, xdmf_vol, it_str, .true.) ! value on element
      else if (MOVIE_VOLUME_TYPE == 4 .and. MOVIE_OUTPUT_CURLNORM) then
        ! reg1_epsdev_displ_norm
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'reg1_epsdev_displ_norm', 'reg1_epsdev_displ_norm', &
                                          npoints_vol_mov_all_proc_cm, xdmf_vol, it_str, .true.) ! value on element
      endif

      write(xdmf_vol,*) '</Grid>'

    enddo

    call write_xdmf_vol_hdf5_footer(xdmf_vol)

    ! close xdmf file
    close(xdmf_vol)

  endif ! output_cm

  !
  ! write out the outer core xdmf file
  !
  if (output_oc) then

    fname_xdmf_vol_oc = trim(OUTPUT_FILES) // "/movie_volume_oc.xmf"
    fname_h5_data_vol_xdmf = "./movie_volume.h5"  ! relative to movie_volume_oc.xmf file

    ! open xdmf file
    open(unit=xdmf_vol, file=trim(fname_xdmf_vol_oc), recl=256)

    call write_xdmf_vol_hdf5_header(nspec_vol_mov_all_proc_oc_conn, npoints_vol_mov_all_proc_oc, &
                                    fname_h5_data_vol_xdmf, xdmf_vol, 3)

    do i = 1, int(NSTEP/NTSTEP_BETWEEN_FRAMES)

      ii = i*NTSTEP_BETWEEN_FRAMES
      !write(it_str, "(i6.6)") ii
      it_str = i2c(ii)

      write(xdmf_vol,*) '<Grid Name="vol_mov" GridType="Uniform">'
      write(xdmf_vol,*) '<Time Value="'//trim(r2c(sngl((ii-1)*DT-t0)))//'" />'
      write(xdmf_vol,*) '<Topology Reference="/Xdmf/Domain/Topology" />'
      write(xdmf_vol,*) '<Geometry Reference="/Xdmf/Domain/Geometry" />'

      ! write headers for each dataset
      if (MOVIE_VOLUME_TYPE == 4 .and. MOVIE_OUTPUT_DIV) then
        ! reg2_div_displ
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'reg2_div_displ', 'reg2_div_displ', &
                                          npoints_vol_mov_all_proc_oc, xdmf_vol, it_str, .true.)
      else if (MOVIE_VOLUME_TYPE == 4 .and. MOVIE_OUTPUT_CURL) then
        ! no output
      else if (MOVIE_VOLUME_TYPE == 4 .and. MOVIE_OUTPUT_CURLNORM) then
        ! no output
      else if (MOVIE_VOLUME_TYPE == 7 .and. OUTPUT_OUTER_CORE) then
        ! reg2_displ
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'reg2_displ', 'reg2_displ', &
                                          npoints_vol_mov_all_proc_oc, xdmf_vol, it_str, .true.)
      else if (MOVIE_VOLUME_TYPE == 8 .and. OUTPUT_OUTER_CORE) then
        ! reg2_veloc
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'reg2_veloc', 'reg2_veloc', &
                                          npoints_vol_mov_all_proc_oc, xdmf_vol, it_str, .true.)
      else if (MOVIE_VOLUME_TYPE == 9 .and. OUTPUT_OUTER_CORE) then
        ! reg2_accel
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'reg2_accel', 'reg2_accel', &
                                          npoints_vol_mov_all_proc_oc, xdmf_vol, it_str, .true.)
      endif

      write(xdmf_vol,*) '</Grid>'

    enddo

    call write_xdmf_vol_hdf5_footer(xdmf_vol)

    ! close xdmf file
    close(xdmf_vol)

  endif ! output_oc

  !
  ! write out the inner core xdmf file
  !
  if (output_ic) then

    fname_xdmf_vol_ic = trim(OUTPUT_FILES) // "/movie_volume_ic.xmf"
    fname_h5_data_vol_xdmf = "./movie_volume.h5"  ! relative to movie_volume_ic.xmf file

    ! open xdmf file
    open(unit=xdmf_vol, file=trim(fname_xdmf_vol_ic), recl=256)

    call write_xdmf_vol_hdf5_header(nspec_vol_mov_all_proc_ic_conn, npoints_vol_mov_all_proc_ic, &
                                    fname_h5_data_vol_xdmf, xdmf_vol, 4)

    do i = 1, int(NSTEP/NTSTEP_BETWEEN_FRAMES)

      ii = i*NTSTEP_BETWEEN_FRAMES
      !write(it_str, "(i6.6)") ii
      it_str = i2c(ii)

      write(xdmf_vol,*) '<Grid Name="vol_mov" GridType="Uniform">'
      write(xdmf_vol,*) '<Time Value="'//trim(r2c(sngl((ii-1)*DT-t0)))//'" />'
      write(xdmf_vol,*) '<Topology Reference="/Xdmf/Domain/Topology" />'
      write(xdmf_vol,*) '<Geometry Reference="/Xdmf/Domain/Geometry" />'

      ! write headers for each dataset
      if (MOVIE_VOLUME_TYPE == 4 .and. MOVIE_OUTPUT_DIV) then
        ! reg3_div_displ
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'reg3_div_displ', 'reg3_div_displ', &
                                          npoints_vol_mov_all_proc_ic, xdmf_vol, it_str, .true.)
      else if (MOVIE_VOLUME_TYPE == 4 .and. MOVIE_OUTPUT_CURL) then
        ! inner_core_epsdev_disple_xx,yy,xy,xz,yz
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'inner_core_epsdev_displ_xx', 'inner_core_epsdev_displ_xx', &
                                          npoints_vol_mov_all_proc_ic, xdmf_vol, it_str, .true.)
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'inner_core_epsdev_displ_yy', 'inner_core_epsdev_displ_yy', &
                                          npoints_vol_mov_all_proc_ic, xdmf_vol, it_str, .true.)
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'inner_core_epsdev_displ_xy', 'inner_core_epsdev_displ_xy', &
                                          npoints_vol_mov_all_proc_ic, xdmf_vol, it_str, .true.)
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'inner_core_epsdev_displ_xz', 'inner_core_epsdev_displ_xz', &
                                          npoints_vol_mov_all_proc_ic, xdmf_vol, it_str, .true.)
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'inner_core_epsdev_displ_yz', 'inner_core_epsdev_displ_yz', &
                                          npoints_vol_mov_all_proc_ic, xdmf_vol, it_str, .true.)
      else if (MOVIE_VOLUME_TYPE == 4 .and. MOVIE_OUTPUT_CURLNORM) then
        ! reg3_epsdev_displ_norm
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'reg3_epsdev_displ_norm', 'reg3_epsdev_displ_norm', &
                                          npoints_vol_mov_all_proc_ic, xdmf_vol, it_str, .true.)
      else if (MOVIE_VOLUME_TYPE == 7 .and. OUTPUT_INNER_CORE) then
        ! reg3_displ
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'reg3_displ', 'reg3_displ', &
                                          npoints_vol_mov_all_proc_ic, xdmf_vol, it_str, .true.)
      else if (MOVIE_VOLUME_TYPE == 8 .and. OUTPUT_INNER_CORE) then
        ! reg3_veloc
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'reg3_veloc', 'reg3_veloc', &
                                          npoints_vol_mov_all_proc_ic, xdmf_vol, it_str, .true.)
      else if (MOVIE_VOLUME_TYPE == 9 .and. OUTPUT_INNER_CORE) then
        ! reg3_accel
        call write_xdmf_vol_hdf5_one_data(fname_h5_data_vol_xdmf, 'reg3_accel', 'reg3_accel', &
                                          npoints_vol_mov_all_proc_ic, xdmf_vol, it_str, .true.)
      endif

      write(xdmf_vol,*) '</Grid>'

    enddo

    call write_xdmf_vol_hdf5_footer(xdmf_vol)

    ! close xdmf file
    close(xdmf_vol)

  endif ! output_ic

  ! to avoid compiler warnings
  i = npoints_3dmovie
  i = npoints_3dmovie_cm
  i = npoints_3dmovie_ic
  i = npoints_3dmovie_oc

  i = nelems_3dmovie
  i = nelems_3dmovie_cm
  i = nelems_3dmovie_ic
  i = nelems_3dmovie_oc

  end subroutine write_xdmf_vol_hdf5

#endif
