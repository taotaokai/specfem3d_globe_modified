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

subroutine movie_surface_init_hdf5()
#ifdef USE_HDF5
  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  integer :: ier

  allocate(offset_poin(0:NPROCTOT_VAL-1),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_poin array')

  npoints_surf_mov_all_proc = 0

#else

      write(*,*) 'Error: HDF5 is not enabled in this version of Specfem3D_Globe.'
      write(*,*) 'Please recompile with the HDF5 option enabled with --with-hdf5'
      stop
#endif
end subroutine movie_surface_init_hdf5


subroutine movie_surface_finalize_hdf5()
#ifdef USE_HDF5
  use specfem_par_movie_hdf5

  implicit none

  deallocate(offset_poin)

#else

      write(*,*) 'Error: HDF5 is not enabled in this version of Specfem3D_Globe.'
      write(*,*) 'Please recompile with the HDF5 option enabled with --with-hdf5'
      stop
#endif
end subroutine movie_surface_finalize_hdf5

subroutine write_movie_surface_mesh_hdf5()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie

#ifdef USE_HDF5
  use specfem_par_movie_hdf5

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_x,store_val_y,store_val_z
  integer :: ipoin,ispec2D,ispec,i,j,k,ier,iglob1,iglob2,iglob3,iglob4,npoin
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval

  call movie_surface_init_hdf5()

  ! gather npoints on each process
  call gather_all_all_singlei(nmovie_points,offset_poin,NPROCTOT_VAL)
  ! total number of points on all processes
  npoints_surf_mov_all_proc = sum(offset_poin)

  ! allocates movie surface arrays
  allocate(store_val_x(nmovie_points), &
           store_val_y(nmovie_points), &
           store_val_z(nmovie_points),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating movie surface location arrays')

  ! initialize h5 file for surface movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! create file and dataset
  file_name = trim(OUTPUT_FILES)//"/movie_surface.h5"

  if (myrank == 0) then
    call h5_create_file(file_name)
    ! create group surf_coord
    call h5_create_group("surf_coord")
    call h5_open_group("surf_coord")
    ! create datasets x, y, z
    call h5_create_dataset_gen_in_group("x", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group("y", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group("z", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)

    ! close group surf_coord
    call h5_close_group()

    ! close file
    call h5_close_file()
  endif

  ! gets coordinates of surface mesh
  ipoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)
    ! in case of global, NCHUNKS_VAL == 6 simulations, be aware that for
    ! the cubed sphere, the mapping changes for different chunks,
    ! i.e. e.g. x(1,1) and x(5,5) flip left and right sides of the elements in geographical coordinates.
    ! for future consideration, like in create_movie_GMT_global.f90 ...
    k = NGLLZ
    ! loop on all the points inside the element
    if (.not. MOVIE_COARSE) then
      do j = 1, NGLLY-1, 1
        do i = 1, NGLLX-1, 1
          ! stores values
          iglob1 = ibool_crust_mantle(i,j,k,ispec)
          iglob2 = ibool_crust_mantle(i+1,j,k,ispec)
          iglob3 = ibool_crust_mantle(i+1,j+1,k,ispec)
          iglob4 = ibool_crust_mantle(i,j+1,k,ispec)
          ! iglob1
          ipoin    = ipoin + 1
          rval     = rstore_crust_mantle(1,iglob1) ! radius r (normalized)
          thetaval = rstore_crust_mantle(2,iglob1) ! colatitude theta (in radian)
          phival   = rstore_crust_mantle(3,iglob1) ! longitude phi (in radian)
          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
          store_val_x(ipoin) = xval
          store_val_y(ipoin) = yval
          store_val_z(ipoin) = zval
          ! iglob2
          ipoin = ipoin + 1
          rval     = rstore_crust_mantle(1,iglob2) ! radius r (normalized)
          thetaval = rstore_crust_mantle(2,iglob2) ! colatitude theta (in radian)
          phival   = rstore_crust_mantle(3,iglob2) ! longitude phi (in radian)
          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
          store_val_x(ipoin) = xval
          store_val_y(ipoin) = yval
          store_val_z(ipoin) = zval
          ! iglob3
          ipoin = ipoin + 1
          rval     = rstore_crust_mantle(1,iglob3) ! radius r (normalized)
          thetaval = rstore_crust_mantle(2,iglob3) ! colatitude theta (in radian)
          phival   = rstore_crust_mantle(3,iglob3) ! longitude phi (in radian)
          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
          store_val_x(ipoin) = xval
          store_val_y(ipoin) = yval
          store_val_z(ipoin) = zval
          ! iglob4
          ipoin = ipoin + 1
          rval     = rstore_crust_mantle(1,iglob4) ! radius r (normalized)
          thetaval = rstore_crust_mantle(2,iglob4) ! colatitude theta (in radian)
          phival   = rstore_crust_mantle(3,iglob4) ! longitude phi (in radian)
          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
          store_val_x(ipoin) = xval
          store_val_y(ipoin) = yval
          store_val_z(ipoin) = zval
        enddo
      enddo
    else ! MOVIE_COARSE
      iglob1 = ibool_crust_mantle(1,1,k,ispec)
      iglob2 = ibool_crust_mantle(NGLLX,1,k,ispec)
      iglob3 = ibool_crust_mantle(NGLLX,NGLLY,k,ispec)
      iglob4 = ibool_crust_mantle(1,NGLLY,k,ispec)

      ipoin = ipoin + 1
      rval  = rstore_crust_mantle(1,iglob1) ! radius r (normalized)
      thetaval = rstore_crust_mantle(2,iglob1) ! colatitude theta (in radian)
      phival   = rstore_crust_mantle(3,iglob1) ! longitude phi (in radian)
      call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
      store_val_x(ipoin) = xval
      store_val_y(ipoin) = yval
      store_val_z(ipoin) = zval

      ipoin = ipoin + 1
      rval  = rstore_crust_mantle(1,iglob2) ! radius r (normalized)
      thetaval = rstore_crust_mantle(2,iglob2) ! colatitude theta (in radian)
      phival   = rstore_crust_mantle(3,iglob2) ! longitude phi (in radian)
      call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
      store_val_x(ipoin) = xval
      store_val_y(ipoin) = yval
      store_val_z(ipoin) = zval

      ipoin = ipoin + 1
      rval  = rstore_crust_mantle(1,iglob3) ! radius r (normalized)
      thetaval = rstore_crust_mantle(2,iglob3) ! colatitude theta (in radian)
      phival   = rstore_crust_mantle(3,iglob3) ! longitude phi (in radian)
      call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
      store_val_x(ipoin) = xval
      store_val_y(ipoin) = yval
      store_val_z(ipoin) = zval

      ipoin = ipoin + 1
      rval  = rstore_crust_mantle(1,iglob4) ! radius r (normalized)
      thetaval = rstore_crust_mantle(2,iglob4) ! colatitude theta (in radian)
      phival   = rstore_crust_mantle(3,iglob4) ! longitude phi (in radian)
      call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
      store_val_x(ipoin) = xval
      store_val_y(ipoin) = yval
      store_val_z(ipoin) = zval

    endif
  enddo
  npoin = ipoin
  if (npoin /= nmovie_points ) call exit_mpi(myrank,'Error number of movie points not equal to nmovie_points')

  call synchronize_all()

  ! write data to h5 file
  call h5_open_file_p_collect(file_name)
  call h5_open_group("surf_coord")

  ! write x, y, z
  call h5_write_dataset_collect_hyperslab_in_group("x", store_val_x, (/sum(offset_poin(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group("y", store_val_y, (/sum(offset_poin(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group("z", store_val_z, (/sum(offset_poin(0:myrank-1))/), H5_COL)

  ! close group and file
  call h5_close_group()
  call h5_close_file_p()

  deallocate(store_val_x,store_val_y,store_val_z)

  if (myrank == 0) then
    ! write xdmf header
    call write_xdmf_surface_header(npoints_surf_mov_all_proc)
  endif

#else

    write(*,*) 'Error: HDF5 is not enabled in this version of Specfem3D_Globe.'
    write(*,*) 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop

#endif

end subroutine write_movie_surface_mesh_hdf5


subroutine write_movie_surface_hdf5()

#ifdef USE_HDF5
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie
  use specfem_par_movie_hdf5

  implicit none

  ! local parameters
  integer :: ipoin,ispec2D,ispec,i,j,k,iglob1,iglob2,iglob3,iglob4

  ! by default: save velocity here to avoid static offset on displacement for movies

  ! gets coordinates of surface mesh and surface displacement
  ipoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)

    ! in case of global, NCHUNKS_VAL == 6 simulations, be aware that for
    ! the cubed sphere, the mapping changes for different chunks,
    ! i.e. e.g. x(1,1) and x(5,5) flip left and right sides of the elements in geographical coordinates.
    ! for future consideration, like in create_movie_GMT_global.f90 ...
    k = NGLLZ

    ! loop on all the points inside the element
    if (.not. MOVIE_COARSE) then
      do j = 1, NGLLY-1, 1
        do i = 1, NGLLX-1, 1
          ! stores values
          iglob1 = ibool_crust_mantle(i,j,k,ispec)
          iglob2 = ibool_crust_mantle(i+1,j,k,ispec)
          iglob3 = ibool_crust_mantle(i+1,j+1,k,ispec)
          iglob4 = ibool_crust_mantle(i,j+1,k,ispec)

          if (MOVIE_VOLUME_TYPE == 5) then
            ! stores displacement
            ! iglob1
            ipoin = ipoin + 1
            store_val_ux(ipoin) = displ_crust_mantle(1,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = displ_crust_mantle(2,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = displ_crust_mantle(3,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob2
            ipoin = ipoin + 1
            store_val_ux(ipoin) = displ_crust_mantle(1,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = displ_crust_mantle(2,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = displ_crust_mantle(3,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob3
            ipoin = ipoin + 1
            store_val_ux(ipoin) = displ_crust_mantle(1,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = displ_crust_mantle(2,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = displ_crust_mantle(3,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob4
            ipoin = ipoin + 1
            store_val_ux(ipoin) = displ_crust_mantle(1,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = displ_crust_mantle(2,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = displ_crust_mantle(3,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
          else
            ! stores velocity
            ! iglob1
            ipoin = ipoin + 1
            store_val_ux(ipoin) = veloc_crust_mantle(1,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = veloc_crust_mantle(2,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = veloc_crust_mantle(3,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob2
            ipoin = ipoin + 1
            store_val_ux(ipoin) = veloc_crust_mantle(1,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = veloc_crust_mantle(2,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = veloc_crust_mantle(3,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob3
            ipoin = ipoin + 1
            store_val_ux(ipoin) = veloc_crust_mantle(1,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = veloc_crust_mantle(2,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = veloc_crust_mantle(3,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob4
            ipoin = ipoin + 1
            store_val_ux(ipoin) = veloc_crust_mantle(1,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = veloc_crust_mantle(2,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = veloc_crust_mantle(3,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
          endif
        enddo
      enddo
    else ! MOVIE_COARSE
      iglob1 = ibool_crust_mantle(1,1,k,ispec)
      iglob2 = ibool_crust_mantle(NGLLX-1,1,k,ispec)
      iglob3 = ibool_crust_mantle(NGLLX-1,NGLLY-1,k,ispec)
      iglob4 = ibool_crust_mantle(1,NGLLY-1,k,ispec)

      if (MOVIE_VOLUME_TYPE == 5) then
        ! stores displacement
        ! iglob1
        ipoin = ipoin + 1
        store_val_ux(ipoin) = displ_crust_mantle(1,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = displ_crust_mantle(2,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = displ_crust_mantle(3,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob2
        ipoin = ipoin + 1
        store_val_ux(ipoin) = displ_crust_mantle(1,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = displ_crust_mantle(2,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = displ_crust_mantle(3,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob3
        ipoin = ipoin + 1
        store_val_ux(ipoin) = displ_crust_mantle(1,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = displ_crust_mantle(2,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = displ_crust_mantle(3,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob4
        ipoin = ipoin + 1
        store_val_ux(ipoin) = displ_crust_mantle(1,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = displ_crust_mantle(2,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = displ_crust_mantle(3,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
      else
        ! stores velocity
        ! iglob1
        ipoin = ipoin + 1
        store_val_ux(ipoin) = veloc_crust_mantle(1,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = veloc_crust_mantle(2,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = veloc_crust_mantle(3,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob2
        ipoin = ipoin + 1
        store_val_ux(ipoin) = veloc_crust_mantle(1,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = veloc_crust_mantle(2,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = veloc_crust_mantle(3,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob3
        ipoin = ipoin + 1
        store_val_ux(ipoin) = veloc_crust_mantle(1,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = veloc_crust_mantle(2,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = veloc_crust_mantle(3,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob4
        ipoin = ipoin + 1
        store_val_ux(ipoin) = veloc_crust_mantle(1,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = veloc_crust_mantle(2,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = veloc_crust_mantle(3,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
      endif
    endif
  enddo
  ! TODO ADD IOSERVER

  ! initialize h5 file for surface movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! create file and dataset
  file_name = trim(OUTPUT_FILES)//"/movie_surface.h5"
  group_name = "it_"//trim(i2c(it))

  ! create dataset
  if (myrank == 0) then
    call h5_open_file(file_name)
    call h5_create_group(group_name)
    call h5_open_group(group_name)

    ! create datasets ux, uy, uz
    call h5_create_dataset_gen_in_group("ux", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group("uy", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group("uz", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)

    ! close group
    call h5_close_group()
    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! write data to h5 file
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  ! write ux, uy, uz
  call h5_write_dataset_collect_hyperslab_in_group("ux", store_val_ux, (/sum(offset_poin(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group("uy", store_val_uy, (/sum(offset_poin(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group("uz", store_val_uz, (/sum(offset_poin(0:myrank-1))/), H5_COL)

  ! close group and file
  call h5_close_group()
  call h5_close_file_p()

  ! write xdmf body
  call write_xdmf_surface_body(it, npoints_surf_mov_all_proc)

#else

    write(*,*) 'Error: HDF5 is not enabled in this version of Specfem3D_Globe.'
    write(*,*) 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop

#endif


end subroutine write_movie_surface_hdf5



!
! xdmf output routines
!
#ifdef USE_HDF5

  subroutine write_xdmf_surface_header(num_nodes)

  use specfem_par
  use specfem_par_movie_hdf5

  implicit none
  integer, intent(in) :: num_nodes
  ! local parameters
  integer :: num_elm
  character(len=MAX_STRING_LEN) :: fname_xdmf_surf
  character(len=MAX_STRING_LEN) :: fname_h5_data_surf_xdmf

  ! checks if anything do, only main process writes out xdmf file
  if (myrank /= 0) return

  ! writeout xdmf file for surface movie
  fname_xdmf_surf = trim(OUTPUT_FILES) // "/movie_surface.xmf"
  fname_h5_data_surf_xdmf = "./movie_surface.h5"   ! relative to movie_surface.xmf file
  ! note: this seems not to work and point to a wrong directory:
  !         fname_h5_data_surf_xdmf = trim(OUTPUT_FILES) // "/movie_surface.h5"

  num_elm = num_nodes / 4

  open(unit=xdmf_surf, file=trim(fname_xdmf_surf), recl=256)

  write(xdmf_surf,'(a)') '<?xml version="1.0" ?>'
  write(xdmf_surf,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(xdmf_surf,*) '<Xdmf Version="3.0">'
  write(xdmf_surf,*) '<Domain Name="mesh">'
  write(xdmf_surf,*) '<Topology Name="topo" TopologyType="Quadrilateral" NumberOfElements="'//trim(i2c(num_elm))//'"/>'
  write(xdmf_surf,*) '<Geometry GeometryType="X_Y_Z">'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                                                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '        '//trim(fname_h5_data_surf_xdmf)//':/surf_coord/x'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '        '//trim(fname_h5_data_surf_xdmf)//':/surf_coord/y'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                       //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '        '//trim(fname_h5_data_surf_xdmf)//':/surf_coord/z'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '</Geometry>'

  write(xdmf_surf,*) '<Grid Name="fensap" GridType="Collection" CollectionType="Temporal">'
  ! 17 lines

  ! file finish
  write(xdmf_surf,*) '</Grid>'
  write(xdmf_surf,*) '</Domain>'
  write(xdmf_surf,*) '</Xdmf>'
  ! 20 lines

  ! position where the additional data will be inserted
  surf_xdmf_pos = 17

  close(xdmf_surf)

  end subroutine write_xdmf_surface_header

#endif


#ifdef USE_HDF5

  subroutine write_xdmf_surface_body(it_io, num_nodes)

  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  integer, intent(in)    :: it_io
  integer, intent(in)    :: num_nodes
  ! local parameters
  integer :: i
  character(len=20) :: it_str
  character(len=MAX_STRING_LEN) :: fname_xdmf_surf
  character(len=MAX_STRING_LEN) :: fname_h5_data_surf_xdmf

  ! checks if anything do, only main process writes out xdmf file
  if (myrank /= 0) return

  ! append data section to xdmf file for surface movie
  fname_xdmf_surf = trim(OUTPUT_FILES)//"/movie_surface.xmf"
  fname_h5_data_surf_xdmf = "./movie_surface.h5"    ! relative to movie_surface.xmf file
  ! this seems to point to a wrong directory:
  !   fname_h5_data_surf_xdmf = trim(OUTPUT_FILES) // "/movie_surface.h5"

  ! open xdmf file
  open(unit=xdmf_surf, file=trim(fname_xdmf_surf), status='old', recl=256)

  ! skip lines till the position where we want to write new information
  do i = 1, surf_xdmf_pos
    read(xdmf_surf, *)
  enddo

  !write(it_str, "(i6.6)") it_io
  it_str = i2c(it_io)

  write(xdmf_surf,*) '<Grid Name="surf_mov" GridType="Uniform">'
  write(xdmf_surf,*) '<Time Value="'//trim(r2c(sngl((it_io-1)*DT-t0)))//'" />'
  write(xdmf_surf,*) '<Topology Reference="/Xdmf/Domain/Topology" />'
  write(xdmf_surf,*) '<Geometry Reference="/Xdmf/Domain/Geometry" />'
  write(xdmf_surf,*) '<Attribute Name="ux" AttributeType="Scalar" Center="Node">'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '      '//trim(fname_h5_data_surf_xdmf)//':/it_'//trim(it_str)//'/ux'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '</Attribute>'
  write(xdmf_surf,*) '<Attribute Name="uy" AttributeType="Scalar" Center="Node">'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '      '//trim(fname_h5_data_surf_xdmf)//':/it_'//trim(it_str)//'/uy'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '</Attribute>'
  write(xdmf_surf,*) '<Attribute Name="uz" AttributeType="Scalar" Center="Node">'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '      '//trim(fname_h5_data_surf_xdmf)//':/it_'//trim(it_str)//'/uz'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '</Attribute>'
  write(xdmf_surf,*) '</Grid>'
  ! 20 lines

  ! file finish
  write(xdmf_surf,*) '</Grid>'
  write(xdmf_surf,*) '</Domain>'
  write(xdmf_surf,*) '</Xdmf>'

  close(xdmf_surf)

  ! updates file record position
  surf_xdmf_pos = surf_xdmf_pos + 20

  end subroutine write_xdmf_surface_body

#endif

