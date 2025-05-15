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


  subroutine read_arrays_solver_hdf5(iregion_code, &
                                nspec,nglob,nglob_xy, &
                                nspec_iso,nspec_tiso,nspec_ani, &
                                rho_vp,rho_vs, &
                                xstore,ystore,zstore, &
                                xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                rhostore, kappavstore,muvstore,kappahstore,muhstore,eta_anisostore, &
                                c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                                c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                                c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                                mu0store, &
                                ibool,idoubling,ispec_is_tiso, &
                                rmassx,rmassy,rmassz, &
                                nglob_oceans,rmass_ocean_load, &
                                b_rmassx,b_rmassy)


  use constants_solver

#ifdef USE_HDF5
  use shared_parameters, only: H5_COL
  use specfem_par, only: &
    ABSORBING_CONDITIONS, &
    LOCAL_PATH,ABSORBING_CONDITIONS
  use manager_hdf5
#endif

  implicit none

  integer,intent(in) :: iregion_code
  integer,intent(in) :: nspec,nglob,nglob_xy
  integer,intent(in) :: nspec_iso,nspec_tiso,nspec_ani

  ! Stacey
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: rho_vp,rho_vs

  real(kind=CUSTOM_REAL), dimension(nglob),intent(inout) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  ! material properties
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_iso),intent(inout) :: &
    rhostore,kappavstore,muvstore

  ! additional arrays for anisotropy stored only where needed to save memory
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_tiso),intent(inout) :: &
    kappahstore,muhstore,eta_anisostore

  ! additional arrays for full anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani),intent(inout) :: &
    c11store,c12store,c13store,c14store,c15store,c16store, &
    c22store,c23store,c24store,c25store,c26store,c33store,c34store, &
    c35store,c36store,c44store,c45store,c46store,c55store,c56store,c66store

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: mu0store

  ! global addressing
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: ibool
  integer, dimension(nspec),intent(inout) :: idoubling
  logical, dimension(nspec),intent(inout) :: ispec_is_tiso

  ! mass matrices and additional ocean load mass matrix
  real(kind=CUSTOM_REAL), dimension(nglob_xy),intent(inout) :: rmassx,rmassy
  real(kind=CUSTOM_REAL), dimension(nglob_xy),intent(inout) :: b_rmassx,b_rmassy

  real(kind=CUSTOM_REAL), dimension(nglob),intent(inout)    :: rmassz

  integer,intent(in) :: nglob_oceans
  real(kind=CUSTOM_REAL), dimension(nglob_oceans),intent(inout) :: rmass_ocean_load

#ifdef USE_HDF5

  ! local parameters
  integer :: lnspec,lnglob
  ! group, dataset name
  character(len=64) :: gname_region

  ! MPI variables
  integer :: info, comm

  ! offset arrays
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nnodes
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nnodes_xy
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nnodes_oceans
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nelems

  name_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/solver_data.h5'
  ! group name
  write(gname_region, "('reg',i1)") iregion_code

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! open the hdf5 file
  call h5_open_file_p_collect(name_database_hdf5)
  ! open the group
  call h5_open_group(gname_region)

  call h5_read_dataset_collect_hyperslab_in_group("offset_nnodes", offset_nnodes, (/0/), H5_COL)
  call h5_read_dataset_collect_hyperslab_in_group("offset_nnodes_xy", offset_nnodes_xy, (/0/), H5_COL)
  call h5_read_dataset_collect_hyperslab_in_group("offset_nnodes_oceans", offset_nnodes_oceans, (/0/), H5_COL)
  call h5_read_dataset_collect_hyperslab_in_group("offset_nelems", offset_nelems, (/0/), H5_COL)

  lnspec = offset_nelems(myrank)
  lnglob = offset_nnodes(myrank)

  ! checks dimensions
  if (lnspec /= nspec) then
    call h5_close_group()
    call h5_close_file()
    print *, 'Error at rank ', myrank
    print *,'Error file dimension: nspec in file = ',lnspec,' but nspec desired:',nspec
    print *,'please check file ', name_database_hdf5
    call exit_mpi(myrank,'Error dimensions in solver_data.h5')
  endif
  if (lnglob /= nglob) then
    close(IIN)
    print *, 'Error at rank ', myrank
    print *,'Error file dimension: nglob in file = ',lnglob,' but nglob desired:',nglob
    print *,'please check file ', name_database_hdf5
    call exit_mpi(myrank,'Error dimensions in solver_data.h5')
  endif

  ! read xstore
  call h5_read_dataset_collect_hyperslab_in_group("xstore", xstore, (/sum(offset_nnodes(0:myrank-1))/), H5_COL)
  ! read ystore
  call h5_read_dataset_collect_hyperslab_in_group("ystore", ystore, (/sum(offset_nnodes(0:myrank-1))/), H5_COL)
  ! read zstore
  call h5_read_dataset_collect_hyperslab_in_group("zstore", zstore, (/sum(offset_nnodes(0:myrank-1))/), H5_COL)
  ! ibool
  call h5_read_dataset_collect_hyperslab_in_group("ibool", ibool, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! idoubling
  call h5_read_dataset_collect_hyperslab_in_group("idoubling", idoubling, (/sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! ispec_is_tiso
  call h5_read_dataset_collect_hyperslab_in_group("ispec_is_tiso", ispec_is_tiso, (/sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! xix
  call h5_read_dataset_collect_hyperslab_in_group("xixstore", xix, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! xiy
  call h5_read_dataset_collect_hyperslab_in_group("xiystore", xiy, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! xiz
  call h5_read_dataset_collect_hyperslab_in_group("xizstore", xiz, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! etax
  call h5_read_dataset_collect_hyperslab_in_group("etaxstore", etax, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! etay
  call h5_read_dataset_collect_hyperslab_in_group("etaystore", etay, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! etaz
  call h5_read_dataset_collect_hyperslab_in_group("etazstore", etaz, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! gammax
  call h5_read_dataset_collect_hyperslab_in_group("gammaxstore", gammax, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! gammay
  call h5_read_dataset_collect_hyperslab_in_group("gammaystore", gammay, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! gammaz
  call h5_read_dataset_collect_hyperslab_in_group("gammazstore", gammaz, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)

  if (iregion_code == IREGION_TRINFINITE .or. iregion_code == IREGION_INFINITE) then
    call h5_close_group()
    call h5_close_file()
    return
  endif

  ! rhostore
  call h5_read_dataset_collect_hyperslab_in_group("rhostore", rhostore, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! kappavstore
  call h5_read_dataset_collect_hyperslab_in_group("kappavstore", kappavstore, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)

  if (iregion_code /= IREGION_OUTER_CORE) then
    ! muvstore
    call h5_read_dataset_collect_hyperslab_in_group("muvstore", muvstore, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  endif

  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust/mantle
    if (ANISOTROPIC_3D_MANTLE_VAL) then
      ! c11store
      call h5_read_dataset_collect_hyperslab_in_group("c11store", c11store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c12store
      call h5_read_dataset_collect_hyperslab_in_group("c12store", c12store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c13store
      call h5_read_dataset_collect_hyperslab_in_group("c13store", c13store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c14store
      call h5_read_dataset_collect_hyperslab_in_group("c14store", c14store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c15store
      call h5_read_dataset_collect_hyperslab_in_group("c15store", c15store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c16store
      call h5_read_dataset_collect_hyperslab_in_group("c16store", c16store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c22store
      call h5_read_dataset_collect_hyperslab_in_group("c22store", c22store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c23store
      call h5_read_dataset_collect_hyperslab_in_group("c23store", c23store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c24store
      call h5_read_dataset_collect_hyperslab_in_group("c24store", c24store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c25store
      call h5_read_dataset_collect_hyperslab_in_group("c25store", c25store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c26store
      call h5_read_dataset_collect_hyperslab_in_group("c26store", c26store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c33store
      call h5_read_dataset_collect_hyperslab_in_group("c33store", c33store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c34store
      call h5_read_dataset_collect_hyperslab_in_group("c34store", c34store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c35store
      call h5_read_dataset_collect_hyperslab_in_group("c35store", c35store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c36store
      call h5_read_dataset_collect_hyperslab_in_group("c36store", c36store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c44store
      call h5_read_dataset_collect_hyperslab_in_group("c44store", c44store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c45store
      call h5_read_dataset_collect_hyperslab_in_group("c45store", c45store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c46store
      call h5_read_dataset_collect_hyperslab_in_group("c46store", c46store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c55store
      call h5_read_dataset_collect_hyperslab_in_group("c55store", c55store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c56store
      call h5_read_dataset_collect_hyperslab_in_group("c56store", c56store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c66store
      call h5_read_dataset_collect_hyperslab_in_group("c66store", c66store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    else
      if (TRANSVERSE_ISOTROPY_VAL) then
        ! kappahstore
        call h5_read_dataset_collect_hyperslab_in_group("kappahstore", kappahstore, &
                                                        (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
        ! muhstore
        call h5_read_dataset_collect_hyperslab_in_group("muhstore", muhstore, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
        ! eta_anisostore
        call h5_read_dataset_collect_hyperslab_in_group("eta_anisostore", eta_anisostore, &
                                                        (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      endif
    endif

    ! mu0store
    call h5_read_dataset_collect_hyperslab_in_group("mu0store", mu0store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)

  case (IREGION_INNER_CORE)
    ! inner core
    if (ANISOTROPIC_INNER_CORE_VAL) then
      ! c11store
      call h5_read_dataset_collect_hyperslab_in_group("c11store", c11store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c12store
      call h5_read_dataset_collect_hyperslab_in_group("c12store", c12store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c13store
      call h5_read_dataset_collect_hyperslab_in_group("c13store", c13store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c33store
      call h5_read_dataset_collect_hyperslab_in_group("c33store", c33store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! c44store
      call h5_read_dataset_collect_hyperslab_in_group("c44store", c44store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    endif
  end select

  ! Stacey
  if (ABSORBING_CONDITIONS) then
    if (iregion_code == IREGION_CRUST_MANTLE) then
      ! rho_vp
      call h5_read_dataset_collect_hyperslab_in_group("rho_vp", rho_vp, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
      ! rho_vs
      call h5_read_dataset_collect_hyperslab_in_group("rho_vs", rho_vs, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    else if (iregion_code == IREGION_OUTER_CORE) then
      ! rho_vp
      call h5_read_dataset_collect_hyperslab_in_group("rho_vp", rho_vp, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    endif
  endif

  ! mass matrices
  if (((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) .and. iregion_code == IREGION_INNER_CORE)) then
    ! rmassx
    call h5_read_dataset_collect_hyperslab_in_group("rmassx", rmassx, (/sum(offset_nnodes_xy(0:myrank-1))/), H5_COL)
    ! rmassy
    call h5_read_dataset_collect_hyperslab_in_group("rmassy", rmassy, (/sum(offset_nnodes_xy(0:myrank-1))/), H5_COL)
  endif

  ! rmassz
  call h5_read_dataset_collect_hyperslab_in_group("rmassz", rmassz, (/sum(offset_nnodes(0:myrank-1))/), H5_COL)

  if (((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) .and. iregion_code == IREGION_INNER_CORE)) then
    ! b_rmassx
    call h5_read_dataset_collect_hyperslab_in_group("b_rmassx", b_rmassx, (/sum(offset_nnodes_xy(0:myrank-1))/), H5_COL)
    ! b_rmassy
    call h5_read_dataset_collect_hyperslab_in_group("b_rmassy", b_rmassy, (/sum(offset_nnodes_xy(0:myrank-1))/), H5_COL)
  endif

  ! read additional ocean load mass matrix
  if (OCEANS_VAL .and. iregion_code == IREGION_CRUST_MANTLE) then
    ! rmass_ocean_load
    call h5_read_dataset_collect_hyperslab_in_group("rmass_ocean_load", rmass_ocean_load, &
                                                    (/sum(offset_nnodes_oceans(0:myrank-1))/), H5_COL)
  endif

  ! close group
  call h5_close_group()
  ! close file
  call h5_close_file_p()

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = iregion_code

  idummy = size(ibool,kind=4)
  idummy = size(idoubling,kind=4)
  idummy = size(ispec_is_tiso,kind=4)

  idummy = size(xstore,kind=4)
  idummy = size(ystore,kind=4)
  idummy = size(zstore,kind=4)

  idummy = size(xix,kind=4)
  idummy = size(xiy,kind=4)
  idummy = size(xiz,kind=4)
  idummy = size(etax,kind=4)
  idummy = size(etay,kind=4)
  idummy = size(etaz,kind=4)
  idummy = size(gammax,kind=4)
  idummy = size(gammay,kind=4)
  idummy = size(gammaz,kind=4)

  idummy = size(c11store,kind=4)
  idummy = size(c12store,kind=4)
  idummy = size(c13store,kind=4)
  idummy = size(c14store,kind=4)
  idummy = size(c15store,kind=4)
  idummy = size(c16store,kind=4)
  idummy = size(c22store,kind=4)
  idummy = size(c23store,kind=4)
  idummy = size(c24store,kind=4)
  idummy = size(c25store,kind=4)
  idummy = size(c26store,kind=4)
  idummy = size(c33store,kind=4)
  idummy = size(c34store,kind=4)
  idummy = size(c35store,kind=4)
  idummy = size(c36store,kind=4)
  idummy = size(c44store,kind=4)
  idummy = size(c45store,kind=4)
  idummy = size(c46store,kind=4)
  idummy = size(c55store,kind=4)
  idummy = size(c56store,kind=4)
  idummy = size(c66store,kind=4)

  idummy = size(kappavstore,kind=4)
  idummy = size(kappahstore,kind=4)
  idummy = size(muvstore,kind=4)
  idummy = size(muhstore,kind=4)
  idummy = size(mu0store,kind=4)
  idummy = size(rhostore,kind=4)
  idummy = size(eta_anisostore,kind=4)
  idummy = size(rho_vp,kind=4)
  idummy = size(rho_vs,kind=4)

  idummy = size(rmassx,kind=4)
  idummy = size(rmassy,kind=4)
  idummy = size(rmassz,kind=4)

  idummy = size(b_rmassx,kind=4)
  idummy = size(b_rmassy,kind=4)

  idummy = size(rmass_ocean_load,kind=4)

  print *
  print *,'ERROR: HDF5 support not enabled'
  print *, 'Please recompile with HDF5 support with the --with-hdf5 option'
  print *
  stop
#endif

  end subroutine read_arrays_solver_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_MPI_hdf5(iregion_code)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_outercore
  use specfem_par_innercore

  use specfem_par_trinfinite
  use specfem_par_infinite
  use specfem_par_full_gravity

#ifdef USE_HDF5
  use manager_hdf5
#endif

  implicit none

  integer, intent(in) :: iregion_code
#ifdef USE_HDF5

  ! local parameters
  integer :: ierr
  integer :: num_interfaces,max_nibool_interfaces, &
             num_phase_ispec,num_colors_outer,num_colors_inner
  integer :: nspec_inner,nspec_outer

  ! group, dataset name
  character(len=64) :: gname_region
  ! MPI variables
  integer :: info, comm

  ! offset arrays
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_num_interfaces
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_max_nibool_interfaces
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_num_phase_ispec
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_inner
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_outer
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_num_colors_outer
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_num_colors_inner


  call h5_initialize()

  name_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/solver_data_mpi.h5'
  write(gname_region, "('reg',i1)") iregion_code

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! open the hdf5 file
  call h5_open_file_p_collect(name_database_hdf5)
  ! open the group
  call h5_open_group(gname_region)

  call h5_read_dataset_collect_hyperslab_in_group('offset_num_interfaces', offset_num_interfaces, (/0/), H5_COL)
  call h5_read_dataset_collect_hyperslab_in_group('offset_max_nibool_interfaces', offset_max_nibool_interfaces, (/0/), H5_COL)
  call h5_read_dataset_collect_hyperslab_in_group('offset_num_phase_ispec', offset_num_phase_ispec, (/0/), H5_COL)
  call h5_read_dataset_collect_hyperslab_in_group('offset_nspec_inner', offset_nspec_inner, (/0/), H5_COL)
  call h5_read_dataset_collect_hyperslab_in_group('offset_nspec_outer', offset_nspec_outer, (/0/), H5_COL)
  call h5_read_dataset_collect_hyperslab_in_group('offset_num_colors_outer', offset_num_colors_outer, (/0/), H5_COL)
  call h5_read_dataset_collect_hyperslab_in_group('offset_num_colors_inner', offset_num_colors_inner, (/0/), H5_COL)

  num_interfaces = offset_num_interfaces(myrank)
  max_nibool_interfaces = offset_max_nibool_interfaces(myrank)
  num_phase_ispec = offset_num_phase_ispec(myrank)
  nspec_inner = offset_nspec_inner(myrank)
  nspec_outer = offset_nspec_outer(myrank)
  num_colors_outer = offset_num_colors_outer(myrank)
  num_colors_inner = offset_num_colors_inner(myrank)

  select case(iregion_code)
  case (IREGION_CRUST_MANTLE)

    num_interfaces_crust_mantle = num_interfaces
    max_nibool_interfaces_cm = max_nibool_interfaces
    num_phase_ispec_crust_mantle = num_phase_ispec
    nspec_inner_crust_mantle = nspec_inner
    nspec_outer_crust_mantle = nspec_outer
    num_colors_outer_crust_mantle = num_colors_outer
    num_colors_inner_crust_mantle = num_colors_inner

    allocate(my_neighbors_crust_mantle(num_interfaces), &
             nibool_interfaces_crust_mantle(num_interfaces), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating my_neighbors_crust_mantle etc.')
    if (num_interfaces > 0) then
      call h5_read_dataset_collect_hyperslab_in_group("my_neighbors", &
                                    my_neighbors_crust_mantle(1:num_interfaces), &
                                    (/sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab_in_group("nibool_interfaces", &
                                    nibool_interfaces_crust_mantle(1:num_interfaces), &
                                    (/sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
    endif

    allocate(ibool_interfaces_crust_mantle(max_nibool_interfaces,num_interfaces), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating ibool_interfaces_crust_mantle etc.')
    if (num_interfaces > 0) then
      ibool_interfaces_crust_mantle(:,:) = 0
      call h5_read_dataset_collect_hyperslab_in_group("ibool_interfaces", &
                                    ibool_interfaces_crust_mantle(1:max_nibool_interfaces,1:num_interfaces), &
                                    (/0,sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
    endif

    allocate(phase_ispec_inner_crust_mantle(num_phase_ispec,2), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating phase_ispec_inner_crust_mantle etc.')
    if (num_phase_ispec > 0) then
      phase_ispec_inner_crust_mantle(:,:) = 0
      call h5_read_dataset_collect_hyperslab_in_group("phase_ispec_inner", &
                                    phase_ispec_inner_crust_mantle(1:num_phase_ispec,1:2), &
                                    (/sum(offset_num_phase_ispec(0:myrank-1)),0/), H5_COL)
    endif

    allocate(num_elem_colors_crust_mantle(num_colors_outer+num_colors_inner), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle etc.')
    if (USE_MESH_COLORING_GPU) then
      call h5_read_dataset_collect_hyperslab_in_group("num_elem_colors", &
                                    num_elem_colors_crust_mantle(1:(num_colors_outer+num_colors_inner)), &
                                    (/sum(offset_num_colors_outer(0:myrank-1))/), H5_COL)
    endif

  case (IREGION_OUTER_CORE)

    num_interfaces_outer_core = num_interfaces
    max_nibool_interfaces_oc = max_nibool_interfaces
    num_phase_ispec_outer_core = num_phase_ispec
    nspec_inner_outer_core = nspec_inner
    nspec_outer_outer_core = nspec_outer
    num_colors_outer_outer_core = num_colors_outer
    num_colors_inner_outer_core = num_colors_inner

    allocate(my_neighbors_outer_core(num_interfaces), &
             nibool_interfaces_outer_core(num_interfaces), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating my_neighbors_outer_core etc.')
    if (num_interfaces > 0) then
      call h5_read_dataset_collect_hyperslab_in_group("my_neighbors", &
                                    my_neighbors_outer_core(1:num_interfaces), &
                                    (/sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab_in_group("nibool_interfaces", &
                                    nibool_interfaces_outer_core(1:num_interfaces), &
                                    (/sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
    endif

    allocate(ibool_interfaces_outer_core(max_nibool_interfaces,num_interfaces), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating ibool_interfaces_outer_core etc.')
    if (num_interfaces > 0) then
    ibool_interfaces_outer_core(:,:) = 0
      call h5_read_dataset_collect_hyperslab_in_group("ibool_interfaces", &
                                    ibool_interfaces_outer_core(1:max_nibool_interfaces,1:num_interfaces), &
                                    (/0,sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
    endif

    allocate(phase_ispec_inner_outer_core(num_phase_ispec,2), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating phase_ispec_inner_outer_core etc.')
    if (num_phase_ispec > 0) then
      phase_ispec_inner_outer_core(:,:) = 0
      call h5_read_dataset_collect_hyperslab_in_group("phase_ispec_inner", &
                                    phase_ispec_inner_outer_core(1:num_phase_ispec,1:2), &
                                    (/sum(offset_num_phase_ispec(0:myrank-1)),0/), H5_COL)
    endif

    allocate(num_elem_colors_outer_core(num_colors_outer+num_colors_inner), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core etc.')
    if (USE_MESH_COLORING_GPU) then
      call h5_read_dataset_collect_hyperslab_in_group("num_elem_colors", &
                                    num_elem_colors_outer_core(1:(num_colors_outer+num_colors_inner)), &
                                    (/sum(offset_num_colors_outer(0:myrank-1))/), H5_COL)
    endif

  case (IREGION_INNER_CORE)

      num_interfaces_inner_core = num_interfaces
      max_nibool_interfaces_ic = max_nibool_interfaces
      num_phase_ispec_inner_core = num_phase_ispec
      nspec_inner_inner_core = nspec_inner
      nspec_outer_inner_core = nspec_outer
      num_colors_outer_inner_core = num_colors_outer
      num_colors_inner_inner_core = num_colors_inner

      allocate(my_neighbors_inner_core(num_interfaces), &
              nibool_interfaces_inner_core(num_interfaces), stat=ierr)
      if (ierr /= 0) call exit_mpi(myrank,'Error allocating my_neighbors_inner_core etc.')
      if (num_interfaces > 0) then
        call h5_read_dataset_collect_hyperslab_in_group("my_neighbors", &
                                      my_neighbors_inner_core(1:num_interfaces), &
                                      (/sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
        call h5_read_dataset_collect_hyperslab_in_group("nibool_interfaces", &
                                      nibool_interfaces_inner_core(1:num_interfaces), &
                                      (/sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
      endif

      allocate(ibool_interfaces_inner_core(max_nibool_interfaces,num_interfaces), stat=ierr)
      if (ierr /= 0) call exit_mpi(myrank,'Error allocating ibool_interfaces_inner_core etc.')
      if (num_interfaces > 0) then
        ibool_interfaces_inner_core(:,:) = 0
        call h5_read_dataset_collect_hyperslab_in_group("ibool_interfaces", &
                                      ibool_interfaces_inner_core(1:max_nibool_interfaces,1:num_interfaces), &
                                      (/0,sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
      endif

      allocate(phase_ispec_inner_inner_core(num_phase_ispec,2), stat=ierr)
      if (ierr /= 0) call exit_mpi(myrank,'Error allocating phase_ispec_inner_inner_core etc.')
      if (num_phase_ispec > 0) then
        phase_ispec_inner_inner_core(:,:) = 0
        call h5_read_dataset_collect_hyperslab_in_group("phase_ispec_inner", &
                                      phase_ispec_inner_inner_core(1:num_phase_ispec,1:2), &
                                      (/sum(offset_num_phase_ispec(0:myrank-1)),0/), H5_COL)
      endif

      allocate(num_elem_colors_inner_core(num_colors_outer+num_colors_inner), stat=ierr)
      if (ierr /= 0) call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core etc.')
      if (USE_MESH_COLORING_GPU) then
        call h5_read_dataset_collect_hyperslab_in_group("num_elem_colors", &
                                      num_elem_colors_inner_core(1:(num_colors_outer+num_colors_inner)), &
                                      (/sum(offset_num_colors_outer(0:myrank-1))/), H5_COL)
      endif

  case (IREGION_TRINFINITE)

    num_interfaces_trinfinite = num_interfaces
    max_nibool_interfaces_trinfinite = max_nibool_interfaces
    num_phase_ispec_trinfinite = num_phase_ispec
    nspec_inner_trinfinite = nspec_inner
    nspec_outer_trinfinite = nspec_outer
    num_colors_outer_trinfinite = num_colors_outer
    num_colors_inner_trinfinite = num_colors_inner

    allocate(my_neighbors_trinfinite(num_interfaces), &
             nibool_interfaces_trinfinite(num_interfaces), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating my_neighbors_trinfinite etc.')
    if (num_interfaces > 0) then
      call h5_read_dataset_collect_hyperslab_in_group("my_neighbors", &
                                    my_neighbors_trinfinite(1:num_interfaces), &
                                    (/sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab_in_group("nibool_interfaces", &
                                    nibool_interfaces_trinfinite(1:num_interfaces), &
                                    (/sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
    endif

    allocate(ibool_interfaces_trinfinite(max_nibool_interfaces,num_interfaces), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating ibool_interfaces_trinfinite etc.')
    if (num_interfaces > 0) then
      ibool_interfaces_trinfinite(:,:) = 0
      call h5_read_dataset_collect_hyperslab_in_group("ibool_interfaces", &
                                    ibool_interfaces_trinfinite(1:max_nibool_interfaces,1:num_interfaces), &
                                    (/0,sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
    endif

    allocate(phase_ispec_inner_trinfinite(num_phase_ispec,2), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating phase_ispec_inner_trinfinite etc.')
    if (num_phase_ispec > 0) then
      phase_ispec_inner_trinfinite(:,:) = 0
      call h5_read_dataset_collect_hyperslab_in_group("phase_ispec_inner", &
                                    phase_ispec_inner_trinfinite(1:num_phase_ispec,1:2), &
                                    (/sum(offset_num_phase_ispec(0:myrank-1)),0/), H5_COL)
    endif

    allocate(num_elem_colors_trinfinite(num_colors_outer+num_colors_inner), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating num_elem_colors_trinfinite etc.')
    if (USE_MESH_COLORING_GPU) then
      call h5_read_dataset_collect_hyperslab_in_group("num_elem_colors", &
                                    num_elem_colors_trinfinite(1:(num_colors_outer+num_colors_inner)), &
                                    (/sum(offset_num_colors_outer(0:myrank-1))/), H5_COL)
    endif

  case (IREGION_INFINITE)

    num_interfaces_infinite = num_interfaces
    max_nibool_interfaces_infinite = max_nibool_interfaces
    num_phase_ispec_infinite = num_phase_ispec
    nspec_inner_infinite = nspec_inner
    nspec_outer_infinite = nspec_outer
    num_colors_outer_infinite = num_colors_outer
    num_colors_inner_infinite = num_colors_inner

    allocate(my_neighbors_infinite(num_interfaces), &
             nibool_interfaces_infinite(num_interfaces), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating my_neighbors_infinite etc.')
    if (num_interfaces > 0) then
      call h5_read_dataset_collect_hyperslab_in_group("my_neighbors", &
                                    my_neighbors_infinite(1:num_interfaces), &
                                    (/sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab_in_group("nibool_interfaces", &
                                    nibool_interfaces_infinite(1:num_interfaces), &
                                    (/sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
    endif

    allocate(ibool_interfaces_infinite(max_nibool_interfaces,num_interfaces), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating ibool_interfaces_infinite etc.')
    if (num_interfaces > 0) then
      ibool_interfaces_infinite(:,:) = 0
      call h5_read_dataset_collect_hyperslab_in_group("ibool_interfaces", &
                                    ibool_interfaces_infinite(1:max_nibool_interfaces,1:num_interfaces), &
                                    (/0,sum(offset_num_interfaces(0:myrank-1))/), H5_COL)
    endif

    allocate(phase_ispec_inner_infinite(num_phase_ispec,2), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating phase_ispec_inner_infinite etc.')
    if (num_phase_ispec > 0) then
      phase_ispec_inner_infinite(:,:) = 0
      call h5_read_dataset_collect_hyperslab_in_group("phase_ispec_inner", &
                                    phase_ispec_inner_infinite(1:num_phase_ispec,1:2), &
                                    (/sum(offset_num_phase_ispec(0:myrank-1)),0/), H5_COL)
    endif

    allocate(num_elem_colors_infinite(num_colors_outer+num_colors_inner), stat=ierr)
    if (ierr /= 0) call exit_mpi(myrank,'Error allocating num_elem_colors_infinite etc.')
    if (USE_MESH_COLORING_GPU) then
      call h5_read_dataset_collect_hyperslab_in_group("num_elem_colors", &
                                    num_elem_colors_infinite(1:(num_colors_outer+num_colors_inner)), &
                                    (/sum(offset_num_colors_outer(0:myrank-1))/), H5_COL)
    endif

  case default
    print *, 'ERROR: unknown region code'
    stop
  end select

  ! close group
  call h5_close_group()
  ! close file
  call h5_close_file_p()

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = iregion_code

  print *
  print *,'ERROR: HDF5 support not enabled'
  print *, 'Please recompile with HDF5 support with the --with-hdf5 option'
  print *
  stop
#endif

  end subroutine read_mesh_databases_MPI_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_coupling_hdf5()

#ifdef USE_HDF5
  use constants

  use meshfem_par, only: &
    myrank, LOCAL_PATH

  !use meshfem_models_par, only: &
  !  HONOR_1D_SPHERICAL_MOHO
    !SAVE_BOUNDARY_MESH,HONOR_1D_SPHERICAL_MOHO,SUPPRESS_CRUSTAL_MESH

  !use regions_mesh_par2, only: &
  !  NSPEC2D_MOHO, NSPEC2D_400, NSPEC2D_670, &
  !  ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot, &
  !  ibelm_670_top,ibelm_670_bot,normal_moho,normal_400,normal_670, &
  !  ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot, &
  !  ispec2D_670_top,ispec2D_670_bot ! prname

  use shared_parameters, only: FULL_GRAVITY

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use specfem_par_trinfinite
  use specfem_par_infinite

  use manager_hdf5
#endif

  implicit none

#ifdef USE_HDF5

  character(len=64) :: gname_region

  ! MPI variables
  integer :: info, comm

  ! offset arrays
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_xmin
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_xmax
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_ymin
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_ymax
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_top
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_bottom
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_moho_top
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_moho_bottom
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_400_top
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_400_bottom
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_670_top
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2D_670_bottom
  !integer, dimension(0:NPROCTOT_VAL-1) :: tmp_array


  ! dump integers
  integer :: tmp_nspec2d_moho, tmp_nspec2d_400, tmp_nspec2d_670

  call world_get_comm(comm)
  call world_get_info_null(info)

  ! initialize hdf5
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! open the hdf5 file
  name_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/boundary.h5'
  call h5_open_file_p_collect(name_database_hdf5)


  if (NSPEC_CRUST_MANTLE > 0) then
    ! open the group
    write(gname_region, "('reg',i1)") IREGION_CRUST_MANTLE
    call h5_open_group(gname_region)

    ! read actual number of elements
    ! nspec2D_xmin
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_xmin", nspec2D_xmin_crust_mantle, (/myrank/), H5_COL)
    ! nspec2D_xmax
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_xmax", nspec2D_xmax_crust_mantle, (/myrank/), H5_COL)
    ! nspec2D_ymin
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_ymin", nspec2D_ymin_crust_mantle, (/myrank/), H5_COL)
    ! nspec2D_ymax
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_ymax", nspec2D_ymax_crust_mantle, (/myrank/), H5_COL)

    ! read offset arrays (stored length)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_xmin", offset_nspec2D_xmin, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_xmax", offset_nspec2D_xmax, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_ymin", offset_nspec2D_ymin, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_ymax", offset_nspec2D_ymax, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_top", offset_nspec2D_top, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_bottom", offset_nspec2D_bottom, (/0/), H5_COL)

    ! ibelm_xmin
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_xmin", ibelm_xmin_crust_mantle, &
                                                    (/sum(offset_nspec2D_xmin(0:myrank-1))/), H5_COL)
    ! ibelm_xmax
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_xmax", ibelm_xmax_crust_mantle, &
                                                    (/sum(offset_nspec2D_xmax(0:myrank-1))/), H5_COL)
    ! ibelm_ymin
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_ymin", ibelm_ymin_crust_mantle, &
                                                    (/sum(offset_nspec2D_ymin(0:myrank-1))/), H5_COL)
    ! ibelm_ymax
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_ymax", ibelm_ymax_crust_mantle, &
                                                    (/sum(offset_nspec2D_ymax(0:myrank-1))/), H5_COL)
    ! ibelm_top
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_top", ibelm_top_crust_mantle, &
                                                    (/sum(offset_nspec2D_top(0:myrank-1))/), H5_COL)
    ! ibelm_bottom
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_bottom", ibelm_bottom_crust_mantle, &
                                                    (/sum(offset_nspec2D_bottom(0:myrank-1))/), H5_COL)

    ! normal_xmin
    call h5_read_dataset_collect_hyperslab_in_group("normal_xmin", normal_xmin_crust_mantle, &
                                                    (/0,0,0,sum(offset_nspec2D_xmin(0:myrank-1))/), H5_COL)
    ! normal_xmax
    call h5_read_dataset_collect_hyperslab_in_group("normal_xmax", normal_xmax_crust_mantle, &
                                                    (/0,0,0,sum(offset_nspec2D_xmax(0:myrank-1))/), H5_COL)
    ! normal_ymin
    call h5_read_dataset_collect_hyperslab_in_group("normal_ymin", normal_ymin_crust_mantle, &
                                                    (/0,0,0,sum(offset_nspec2D_ymin(0:myrank-1))/), H5_COL)
    ! normal_ymax
    call h5_read_dataset_collect_hyperslab_in_group("normal_ymax", normal_ymax_crust_mantle, &
                                                    (/0,0,0,sum(offset_nspec2D_ymax(0:myrank-1))/), H5_COL)
    ! normal_top
    call h5_read_dataset_collect_hyperslab_in_group("normal_top", normal_top_crust_mantle, &
                                                    (/0,0,0,sum(offset_nspec2D_top(0:myrank-1))/), H5_COL)
    ! normal_bottom
    call h5_read_dataset_collect_hyperslab_in_group("normal_bottom", normal_bottom_crust_mantle, &
                                                    (/0,0,0,sum(offset_nspec2D_bottom(0:myrank-1))/), H5_COL)

    ! jacobian_xmin
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_xmin", jacobian2D_xmin_crust_mantle, &
                                                    (/0,0,sum(offset_nspec2D_xmin(0:myrank-1))/), H5_COL)
    ! jacobian_xmax
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_xmax", jacobian2D_xmax_crust_mantle, &
                                                    (/0,0,sum(offset_nspec2D_xmax(0:myrank-1))/), H5_COL)
    ! jacobian_ymin
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_ymin", jacobian2D_ymin_crust_mantle, &
                                                    (/0,0,sum(offset_nspec2D_ymin(0:myrank-1))/), H5_COL)
    ! jacobian_ymax
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_ymax", jacobian2D_ymax_crust_mantle, &
                                                    (/0,0,sum(offset_nspec2D_ymax(0:myrank-1))/), H5_COL)
    ! jacobian_top
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_top", jacobian2D_top_crust_mantle, &
                                                    (/0,0,sum(offset_nspec2D_top(0:myrank-1))/), H5_COL)
    ! jacobian_bottom
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_bottom", jacobian2D_bottom_crust_mantle, &
                                                    (/0,0,sum(offset_nspec2D_bottom(0:myrank-1))/), H5_COL)

    ! close group
    call h5_close_group()
  endif ! NSPEC_CRUST_MANTLE > 0

  if (NSPEC_OUTER_CORE > 0) then
    ! change group name
    write(gname_region, "('reg',i1)") IREGION_OUTER_CORE
    call h5_open_group(gname_region)

    ! read actual number of elements
    ! nspec2D_xmin
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_xmin", nspec2D_xmin_outer_core, (/myrank/), H5_COL)
    ! nspec2D_xmax
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_xmax", nspec2D_xmax_outer_core, (/myrank/), H5_COL)
    ! nspec2D_ymin
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_ymin", nspec2D_ymin_outer_core, (/myrank/), H5_COL)
    ! nspec2D_ymax
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_ymax", nspec2D_ymax_outer_core, (/myrank/), H5_COL)

    ! read offset arrays (stored length)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_xmin", offset_nspec2D_xmin, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_xmax", offset_nspec2D_xmax, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_ymin", offset_nspec2D_ymin, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_ymax", offset_nspec2D_ymax, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_top", offset_nspec2D_top, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_bottom", offset_nspec2D_bottom, (/0/), H5_COL)

    ! ibelm_xmin
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_xmin", ibelm_xmin_outer_core, &
                                                    (/sum(offset_nspec2D_xmin(0:myrank-1))/), H5_COL)
    ! ibelm_xmax
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_xmax", ibelm_xmax_outer_core, &
                                                    (/sum(offset_nspec2D_xmax(0:myrank-1))/), H5_COL)
    ! ibelm_ymin
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_ymin", ibelm_ymin_outer_core, &
                                                    (/sum(offset_nspec2D_ymin(0:myrank-1))/), H5_COL)
    ! ibelm_ymax
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_ymax", ibelm_ymax_outer_core, &
                                                    (/sum(offset_nspec2D_ymax(0:myrank-1))/), H5_COL)
    ! ibelm_top
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_top", ibelm_top_outer_core, &
                                                    (/sum(offset_nspec2D_top(0:myrank-1))/), H5_COL)
    ! ibelm_bottom
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_bottom", ibelm_bottom_outer_core, &
                                                    (/sum(offset_nspec2D_bottom(0:myrank-1))/), H5_COL)

    ! normal_xmin
    call h5_read_dataset_collect_hyperslab_in_group("normal_xmin", normal_xmin_outer_core, &
                                                    (/0,0,0,sum(offset_nspec2D_xmin(0:myrank-1))/), H5_COL)
    ! normal_xmax
    call h5_read_dataset_collect_hyperslab_in_group("normal_xmax", normal_xmax_outer_core, &
                                                    (/0,0,0,sum(offset_nspec2D_xmax(0:myrank-1))/), H5_COL)
    ! normal_ymin
    call h5_read_dataset_collect_hyperslab_in_group("normal_ymin", normal_ymin_outer_core, &
                                                    (/0,0,0,sum(offset_nspec2D_ymin(0:myrank-1))/), H5_COL)
    ! normal_ymax
    call h5_read_dataset_collect_hyperslab_in_group("normal_ymax", normal_ymax_outer_core, &
                                                    (/0,0,0,sum(offset_nspec2D_ymax(0:myrank-1))/), H5_COL)
    ! normal_top
    call h5_read_dataset_collect_hyperslab_in_group("normal_top", normal_top_outer_core, &
                                                    (/0,0,0,sum(offset_nspec2D_top(0:myrank-1))/), H5_COL)
    ! normal_bot
    call h5_read_dataset_collect_hyperslab_in_group("normal_bottom", normal_bottom_outer_core, &
                                                    (/0,0,0,sum(offset_nspec2D_bottom(0:myrank-1))/), H5_COL)

    ! jacobian_xmin
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_xmin", jacobian2D_xmin_outer_core, &
                                                    (/0,0,sum(offset_nspec2D_xmin(0:myrank-1))/), H5_COL)
    ! jacobian_xmax
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_xmax", jacobian2D_xmax_outer_core, &
                                                    (/0,0,sum(offset_nspec2D_xmax(0:myrank-1))/), H5_COL)
    ! jacobian_ymin
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_ymin", jacobian2D_ymin_outer_core, &
                                                    (/0,0,sum(offset_nspec2D_ymin(0:myrank-1))/), H5_COL)
    ! jacobian_ymax
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_ymax", jacobian2D_ymax_outer_core, &
                                                    (/0,0,sum(offset_nspec2D_ymax(0:myrank-1))/), H5_COL)
    ! jacobian_top
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_top", jacobian2D_top_outer_core, &
                                                    (/0,0,sum(offset_nspec2D_top(0:myrank-1))/), H5_COL)
    ! jacobian_bottom
    call h5_read_dataset_collect_hyperslab_in_group("jacobian2D_bottom", jacobian2D_bottom_outer_core, &
                                                    (/0,0,sum(offset_nspec2D_bottom(0:myrank-1))/), H5_COL)

    ! close group
    call h5_close_group()
  endif ! NSPEC_OUTER_CORE > 0

  if (NSPEC_INNER_CORE > 0) then

    ! change group name
    write(gname_region, "('reg',i1)") IREGION_INNER_CORE
    call h5_open_group(gname_region)

    ! read actual number of elements
    ! nspec2D_xmin
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_xmin", nspec2D_xmin_inner_core, (/myrank/), H5_COL)
    ! nspec2D_xmax
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_xmax", nspec2D_xmax_inner_core, (/myrank/), H5_COL)
    ! nspec2D_ymin
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_ymin", nspec2D_ymin_inner_core, (/myrank/), H5_COL)
    ! nspec2D_ymax
    call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_ymax", nspec2D_ymax_inner_core, (/myrank/), H5_COL)

    ! read offset arrays (stored length)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_xmin", offset_nspec2D_xmin, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_xmax", offset_nspec2D_xmax, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_ymin", offset_nspec2D_ymin, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_ymax", offset_nspec2D_ymax, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_top", offset_nspec2D_top, (/0/), H5_COL)
    call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_bottom", offset_nspec2D_bottom, (/0/), H5_COL)

    ! ibelm_xmin
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_xmin", ibelm_xmin_inner_core, &
                                                    (/sum(offset_nspec2D_xmin(0:myrank-1))/), H5_COL)
    ! ibelm_xmax
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_xmax", ibelm_xmax_inner_core, &
                                                    (/sum(offset_nspec2D_xmax(0:myrank-1))/), H5_COL)
    ! ibelm_ymin
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_ymin", ibelm_ymin_inner_core, &
                                                    (/sum(offset_nspec2D_ymin(0:myrank-1))/), H5_COL)
    ! ibelm_ymax
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_ymax", ibelm_ymax_inner_core, &
                                                    (/sum(offset_nspec2D_ymax(0:myrank-1))/), H5_COL)
    ! ibelm_top
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_top", ibelm_top_inner_core, &
                                                    (/sum(offset_nspec2D_top(0:myrank-1))/), H5_COL)
    ! ibelm_bottom
    call h5_read_dataset_collect_hyperslab_in_group("ibelm_bottom", ibelm_bottom_inner_core, &
                                                    (/sum(offset_nspec2D_bottom(0:myrank-1))/), H5_COL)

    ! close group
    call h5_close_group()
  endif ! NSPEC_INNER_CORE > 0

  if (FULL_GRAVITY) then
    if (ADD_TRINF) then
      if (NSPEC_TRINFINITE > 0) then
        ! change group name
        write(gname_region, "('reg',i1)") IREGION_TRINFINITE
        call h5_open_group(gname_region)

        ! read actual number of elements
        ! nspec2D_xmin
        call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_xmin", nspec2D_xmin_trinfinite, (/myrank/), H5_COL)
        ! nspec2D_xmax
        call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_xmax", nspec2D_xmax_trinfinite, (/myrank/), H5_COL)
        ! nspec2D_ymin
        call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_ymin", nspec2D_ymin_trinfinite, (/myrank/), H5_COL)
        ! nspec2D_ymax
        call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_ymax", nspec2D_ymax_trinfinite, (/myrank/), H5_COL)

        ! read offset arrays (stored length)
        call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_xmin", offset_nspec2D_xmin, (/0/), H5_COL)
        call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_xmax", offset_nspec2D_xmax, (/0/), H5_COL)
        call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_ymin", offset_nspec2D_ymin, (/0/), H5_COL)
        call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_ymax", offset_nspec2D_ymax, (/0/), H5_COL)
        call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_top", offset_nspec2D_top, (/0/), H5_COL)
        call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_bottom", offset_nspec2D_bottom, (/0/), H5_COL)

        ! ibelm_xmin
        call h5_read_dataset_collect_hyperslab_in_group("ibelm_xmin", ibelm_xmin_trinfinite, &
                                                        (/sum(offset_nspec2D_xmin(0:myrank-1))/), H5_COL)
        ! ibelm_xmax
        call h5_read_dataset_collect_hyperslab_in_group("ibelm_xmax", ibelm_xmax_trinfinite, &
                                                        (/sum(offset_nspec2D_xmax(0:myrank-1))/), H5_COL)
        ! ibelm_ymin
        call h5_read_dataset_collect_hyperslab_in_group("ibelm_ymin", ibelm_ymin_trinfinite, &
                                                        (/sum(offset_nspec2D_ymin(0:myrank-1))/), H5_COL)
        ! ibelm_ymax
        call h5_read_dataset_collect_hyperslab_in_group("ibelm_ymax", ibelm_ymax_trinfinite, &
                                                        (/sum(offset_nspec2D_ymax(0:myrank-1))/), H5_COL)
        ! ibelm_top
        call h5_read_dataset_collect_hyperslab_in_group("ibelm_top", ibelm_top_trinfinite, &
                                                        (/sum(offset_nspec2D_top(0:myrank-1))/), H5_COL)
        ! ibelm_bottom
        call h5_read_dataset_collect_hyperslab_in_group("ibelm_bottom", ibelm_bottom_trinfinite, &
                                                        (/sum(offset_nspec2D_bottom(0:myrank-1))/), H5_COL)

        ! close group
        call h5_close_group()
      endif ! NSPEC_TRINFINITE > 0
    endif ! ADD_TRINF

    if (NSPEC_INFINITE > 0) then
      ! change group name
      write(gname_region, "('reg',i1)") IREGION_INFINITE
      call h5_open_group(gname_region)

      ! read actual number of elements
      ! nspec2D_xmin
      call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_xmin", nspec2D_xmin_infinite, (/myrank/), H5_COL)
      ! nspec2D_xmax
      call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_xmax", nspec2D_xmax_infinite, (/myrank/), H5_COL)
      ! nspec2D_ymin
      call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_ymin", nspec2D_ymin_infinite, (/myrank/), H5_COL)
      ! nspec2D_ymax
      call h5_read_dataset_scalar_collect_hyperslab_in_group("nspec2D_ymax", nspec2D_ymax_infinite, (/myrank/), H5_COL)

      ! read offset arrays (stored length)
      call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_xmin", offset_nspec2D_xmin, (/0/), H5_COL)
      call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_xmax", offset_nspec2D_xmax, (/0/), H5_COL)
      call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_ymin", offset_nspec2D_ymin, (/0/), H5_COL)
      call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_ymax", offset_nspec2D_ymax, (/0/), H5_COL)
      call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_top", offset_nspec2D_top, (/0/), H5_COL)
      call h5_read_dataset_collect_hyperslab_in_group("sub_nspec2D_bottom", offset_nspec2D_bottom, (/0/), H5_COL)

      ! ibelm_xmin
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_xmin", ibelm_xmin_infinite, &
                                                      (/sum(offset_nspec2D_xmin(0:myrank-1))/), H5_COL)
      ! ibelm_xmax
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_xmax", ibelm_xmax_infinite, &
                                                      (/sum(offset_nspec2D_xmax(0:myrank-1))/), H5_COL)
      ! ibelm_ymin
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_ymin", ibelm_ymin_infinite, &
                                                      (/sum(offset_nspec2D_ymin(0:myrank-1))/), H5_COL)
      ! ibelm_ymax
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_ymax", ibelm_ymax_infinite, &
                                                      (/sum(offset_nspec2D_ymax(0:myrank-1))/), H5_COL)
      ! ibelm_top
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_top", ibelm_top_infinite, &
                                                      (/sum(offset_nspec2D_top(0:myrank-1))/), H5_COL)
      ! ibelm_bottom
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_bottom", ibelm_bottom_infinite, &
                                                      (/sum(offset_nspec2D_bottom(0:myrank-1))/), H5_COL)

      ! close group
      call h5_close_group()
    endif ! NSPEC_INFINITE > 0
  endif ! FULL_GRAVITY

  ! close file
  call h5_close_file_p()

    ! boundary mesh for crust and mantle
  if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then
    ! open boundary_disc.h5 file
    name_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/boundary_disc.h5'
    call h5_open_file_p_collect(name_database_hdf5)

    if (NSPEC_CRUST_MANTLE > 0) then
      write(gname_region, "('reg',i1)") IREGION_CRUST_MANTLE
      call h5_open_group(gname_region)

      ! read NSPEC2D_MOHO
      call h5_read_dataset_scalar_collect_hyperslab_in_group("NSPEC2D_MOHO", tmp_nspec2d_moho, (/myrank/), H5_COL)
      ! read NSPEC2D_400
      call h5_read_dataset_scalar_collect_hyperslab_in_group("NSPEC2D_400", tmp_nspec2d_400, (/myrank/), H5_COL)
      ! read NSPEC2D_670
      call h5_read_dataset_scalar_collect_hyperslab_in_group("NSPEC2D_670", tmp_nspec2d_670, (/myrank/), H5_COL)

      ! checks setup
      if (tmp_nspec2d_moho /= NSPEC2D_MOHO .or. tmp_nspec2d_400 /= NSPEC2D_400 .or. tmp_nspec2d_670 /= NSPEC2D_670) then
        print *,'Error: invalid NSPEC2D values read in for solver: ', &
                 tmp_nspec2d_moho,tmp_nspec2d_400,tmp_nspec2d_670,'(boundary_disc.h5)'
        print *,'       should be MOHO/400/670 : ',NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,'(mesh_parameters.h5)'
        call exit_mpi(myrank, 'Error reading boundary_disc.h5 file')
      endif

      ! read the actual number of elements
      ! sub_NSPEC2D_MOHO_top
      call h5_read_dataset_collect_hyperslab_in_group("sub_NSPEC2D_MOHO_top", offset_nspec2D_moho_top, (/0/), H5_COL)
      ! sub_NSPEC2D_MOHO_bottom
      call h5_read_dataset_collect_hyperslab_in_group("sub_NSPEC2D_MOHO_bottom", offset_nspec2D_moho_bottom, (/0/), H5_COL)
      ! sub_NSPEC2D_400_top
      call h5_read_dataset_collect_hyperslab_in_group("sub_NSPEC2D_400_top", offset_nspec2D_400_top, (/0/), H5_COL)
      ! sub_NSPEC2D_400_bottom
      call h5_read_dataset_collect_hyperslab_in_group("sub_NSPEC2D_400_bottom", offset_nspec2D_400_bottom, (/0/), H5_COL)
      ! sub_NSPEC2D_670_top
      call h5_read_dataset_collect_hyperslab_in_group("sub_NSPEC2D_670_top", offset_nspec2D_670_top, (/0/), H5_COL)
      ! sub_NSPEC2D_670_bottom
      call h5_read_dataset_collect_hyperslab_in_group("sub_NSPEC2D_670_bot", offset_nspec2D_670_bottom, (/0/), H5_COL)

      ! ibelm_moho_top
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_moho_top", ibelm_moho_top, &
                                                      (/sum(offset_nspec2D_moho_top(0:myrank-1))/), H5_COL)
      ! ibelm_moho_bot
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_moho_bot", ibelm_moho_bot, &
                                                      (/sum(offset_nspec2D_moho_bottom(0:myrank-1))/), H5_COL)
      ! ibelm_400_top
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_400_top", ibelm_400_top, &
                                                      (/sum(offset_nspec2D_400_top(0:myrank-1))/), H5_COL)
      ! ibelm_400_bot
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_400_bot", ibelm_400_bot, &
                                                      (/sum(offset_nspec2D_400_bottom(0:myrank-1))/), H5_COL)
      ! ibelm_670_top
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_670_top", ibelm_670_top, &
                                                      (/sum(offset_nspec2D_670_top(0:myrank-1))/), H5_COL)
      ! ibelm_670_bot
      call h5_read_dataset_collect_hyperslab_in_group("ibelm_670_bot", ibelm_670_bot, &
                                                      (/sum(offset_nspec2D_670_bottom(0:myrank-1))/), H5_COL)
      ! normal_moho
      call h5_read_dataset_collect_hyperslab_in_group("normal_moho", normal_moho, &
                                                      (/0,0,0,sum(offset_nspec2D_moho_top(0:myrank-1))/), H5_COL)
      ! normal_400
      call h5_read_dataset_collect_hyperslab_in_group("normal_400", normal_400, &
                                                      (/0,0,0,sum(offset_nspec2D_400_top(0:myrank-1))/), H5_COL)
      ! normal_670
      call h5_read_dataset_collect_hyperslab_in_group("normal_670", normal_670, &
                                                      (/0,0,0,sum(offset_nspec2D_670_top(0:myrank-1))/), H5_COL)

      ! close group
      call h5_close_group()
    endif ! NSPEC_CRUST_MANTLE > 0
  endif ! SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3

#else
  print *
  print *,'ERROR: HDF5 support not enabled'
  print *, 'Please recompile with HDF5 support with the --with-hdf5 option'
  print *
  stop
#endif

  end subroutine read_mesh_databases_coupling_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_stacey_hdf5

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

#ifdef USE_HDF5
  use manager_hdf5
#endif

  implicit none

#ifdef USE_HDF5

  ! offset
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_num_abs_boundary_faces

  ! hdf5 variables
  integer :: comm, info, ier
  character(len=64) :: gname_region


  call world_get_comm(comm)
  call world_get_info_null(info)

  ! initialize hdf5
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! name of the file and group
  name_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/stacey.h5'

  ! open the hdf5 file
  call h5_open_file_p_collect(name_database_hdf5)

  if (NSPEC_CRUST_MANTLE > 0) then

    ! open the group
    write(gname_region, "('reg',i1)") IREGION_CRUST_MANTLE
    call h5_open_group(gname_region)

    ! read the offset_num_abs_boundary_faces
    call h5_read_dataset_collect_hyperslab_in_group("num_abs_boundary_faces", offset_num_abs_boundary_faces, (/0/), H5_COL)
    ! for this process
    num_abs_boundary_faces_crust_mantle = offset_num_abs_boundary_faces(myrank)

    if (num_abs_boundary_faces_crust_mantle > 0) then
      ! allocates absorbing boundary arrays
      allocate(abs_boundary_ispec_crust_mantle(num_abs_boundary_faces_crust_mantle),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
      allocate(abs_boundary_ijk_crust_mantle(3,NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
      allocate(abs_boundary_jacobian2Dw_crust_mantle(NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
      allocate(abs_boundary_normal_crust_mantle(NDIM,NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_normal')
      allocate(abs_boundary_npoin_crust_mantle(num_abs_boundary_faces_crust_mantle),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
      if (ier /= 0) stop 'Error allocating array abs_boundary_ispec etc.'

      ! abs_boundary_ispec
      call h5_read_dataset_collect_hyperslab_in_group("abs_boundary_ispec", abs_boundary_ispec_crust_mantle, &
                                                      (/sum(offset_num_abs_boundary_faces(0:myrank-1))/), H5_COL)
      ! abs_boundary_npoin
      call h5_read_dataset_collect_hyperslab_in_group("abs_boundary_npoin", abs_boundary_npoin_crust_mantle, &
                                                      (/sum(offset_num_abs_boundary_faces(0:myrank-1))/), H5_COL)
      ! abs_boundary_ijk
      call h5_read_dataset_collect_hyperslab_in_group("abs_boundary_ijk", abs_boundary_ijk_crust_mantle, &
                                                      (/0,0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), H5_COL)
      ! abs_boundary_jacobian2Dw
      call h5_read_dataset_collect_hyperslab_in_group("abs_boundary_jacobian2Dw", abs_boundary_jacobian2Dw_crust_mantle, &
                                                      (/0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), H5_COL)
      ! abs_boundary_normal
      call h5_read_dataset_collect_hyperslab_in_group("abs_boundary_normal", abs_boundary_normal_crust_mantle, &
                                                      (/0,0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), H5_COL)

    else ! dummy
      ! dummy arrays
      allocate(abs_boundary_ispec_crust_mantle(1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
      allocate(abs_boundary_ijk_crust_mantle(1,1,1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
      allocate(abs_boundary_jacobian2Dw_crust_mantle(1,1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
      allocate(abs_boundary_normal_crust_mantle(1,1,1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_normal')
      allocate(abs_boundary_npoin_crust_mantle(1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
      abs_boundary_ispec_crust_mantle(:) = 0; abs_boundary_npoin_crust_mantle(:) = 0
      abs_boundary_ijk_crust_mantle(:,:,:) = 0
      abs_boundary_jacobian2Dw_crust_mantle(:,:) = 0.0; abs_boundary_normal_crust_mantle(:,:,:) = 0.0

    endif  ! num_abs_boundary_faces_crust_mantle > 0

    ! close group
    call h5_close_group()
  endif ! NSPEC_CRUST_MANTLE > 0

  if (NSPEC_OUTER_CORE > 0) then
    ! open the group
    write(gname_region, "('reg',i1)") IREGION_OUTER_CORE
    call h5_open_group(gname_region)

    ! read the offset_num_abs_boundary_faces
    call h5_read_dataset_collect_hyperslab_in_group("num_abs_boundary_faces", offset_num_abs_boundary_faces, (/0/), H5_COL)
    ! for this process
    num_abs_boundary_faces_outer_core = offset_num_abs_boundary_faces(myrank)

    if (num_abs_boundary_faces_outer_core > 0) then
      ! allocates absorbing boundary arrays
      allocate(abs_boundary_ispec_outer_core(num_abs_boundary_faces_outer_core),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
      allocate(abs_boundary_ijk_outer_core(3,NGLLSQUARE,num_abs_boundary_faces_outer_core),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
      allocate(abs_boundary_jacobian2Dw_outer_core(NGLLSQUARE,num_abs_boundary_faces_outer_core),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
      !allocate(abs_boundary_normal_outer_core(NDIM,NGLLSQUARE,num_abs_boundary_faces_outer_core),stat=ier)
      !if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_normal')
      allocate(abs_boundary_npoin_outer_core(num_abs_boundary_faces_outer_core),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
      if (ier /= 0) stop 'Error allocating array abs_boundary_ispec etc.'

      ! abs_boundary_ispec
      call h5_read_dataset_collect_hyperslab_in_group("abs_boundary_ispec", abs_boundary_ispec_outer_core, &
                                                      (/sum(offset_num_abs_boundary_faces(0:myrank-1))/), H5_COL)
      ! abs_boundary_npoin
      call h5_read_dataset_collect_hyperslab_in_group("abs_boundary_npoin", abs_boundary_npoin_outer_core, &
                                                      (/sum(offset_num_abs_boundary_faces(0:myrank-1))/), H5_COL)
      ! abs_boundary_ijk
      call h5_read_dataset_collect_hyperslab_in_group("abs_boundary_ijk", abs_boundary_ijk_outer_core, &
                                                      (/0,0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), H5_COL)
      ! abs_boundary_jacobian2Dw
      call h5_read_dataset_collect_hyperslab_in_group("abs_boundary_jacobian2Dw", abs_boundary_jacobian2Dw_outer_core, &
                                                      (/0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), H5_COL)
      ! abs_boundary_normal
      !call h5_read_dataset_collect_hyperslab_in_group("abs_boundary_normal", abs_boundary_normal_outer_core, &
      !                                                (/0,0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), H5_COL)

    else ! dummy
      ! dummy arrays
      allocate(abs_boundary_ispec_outer_core(1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
      allocate(abs_boundary_ijk_outer_core(1,1,1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
      allocate(abs_boundary_jacobian2Dw_outer_core(1,1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
      !allocate(abs_boundary_normal_outer_core(1,1,1),stat=ier)
      !if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_normal')
      allocate(abs_boundary_npoin_outer_core(1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
      abs_boundary_ispec_outer_core(:) = 0
      abs_boundary_npoin_outer_core(:) = 0
      abs_boundary_ijk_outer_core(:,:,:) = 0
      abs_boundary_jacobian2Dw_outer_core(:,:) = 0.0
      !abs_boundary_normal_outer_core(:,:,:) = 0.0

    endif  ! num_abs_boundary_faces_outer_core > 0

    ! close group
    call h5_close_group()
  endif ! NSPEC_OUTER_CORE > 0

  ! close file
  call h5_close_file_p()


#else
  print *
  print *,'ERROR: HDF5 support not enabled'
  print *, 'Please recompile with HDF5 support with the --with-hdf5 option'
  print *
  stop
#endif

  end subroutine read_mesh_databases_stacey_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_attenuation_hdf5(iregion_code, factor_common, scale_factor, tau_s, vnspec, f_c_source)

  use constants_solver

#ifdef USE_HDF5
  use shared_parameters, only: H5_COL
  use specfem_par, only: ATTENUATION_VAL,LOCAL_PATH
  use manager_hdf5
#endif

  implicit none

  integer,intent(in) :: iregion_code

  integer,intent(in) :: vnspec
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec),intent(inout) :: scale_factor
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec),intent(inout) :: factor_common
  double precision, dimension(N_SLS),intent(inout) :: tau_s
  double precision,intent(inout) :: f_c_source

#ifdef USE_HDF5
  ! offset
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nelems
  ! hdf5 variables
  character(len=64) :: gname_region
  integer :: comm, info
  double precision, dimension(0:NPROCTOT_VAL-1) :: tmp_dp_arr

  if (.not. ATTENUATION_VAL) return

  call world_get_comm(comm)
  call world_get_info_null(info)

  ! initialize hdf5
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! get offset_nelems from solver_data.h5
  name_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/solver_data.h5'
  write(gname_region, "('reg',i1)") iregion_code
  ! open the hdf5 file
  call h5_open_file_p_collect(name_database_hdf5)
  ! open the group
  call h5_open_group(gname_region)
  ! read the offset_nelems
  call h5_read_dataset_collect_hyperslab_in_group("offset_nelems", offset_nelems, (/0/), H5_COL)
  ! close group and file
  call h5_close_group()
  call h5_close_file_p()

  ! open attenuation.h5 file
  name_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/attenuation.h5'
  ! open the hdf5 file
  call h5_open_file_p_collect(name_database_hdf5)
  ! open the group
  call h5_open_group(gname_region)

  ! tau_s
  call h5_read_dataset_collect_hyperslab_in_group("tau_s_store", tau_s, (/myrank*N_SLS/), H5_COL)
  ! factor_common (tau_e_store)
  call h5_read_dataset_collect_hyperslab_in_group("tau_e_store", factor_common, (/0,0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! scale_factor (Qmu_store)
  call h5_read_dataset_collect_hyperslab_in_group("Qmu_store", scale_factor, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  ! f_c_source
  call h5_read_dataset_collect_hyperslab_in_group("att_f_c_source", tmp_dp_arr, (/0/), H5_COL)
  f_c_source = tmp_dp_arr(myrank)

  ! close group
  call h5_close_group()
  ! close file
  call h5_close_file_p()

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy
  double precision :: d_dummy

  idummy = iregion_code
  d_dummy = f_c_source

  idummy = size(factor_common,kind=4)
  idummy = size(scale_factor,kind=4)
  idummy = size(tau_s,kind=4)

  print *
  print *,'ERROR: HDF5 support not enabled'
  print *, 'Please recompile with HDF5 support with the --with-hdf5 option'
  print *
  stop
#endif

  end subroutine read_attenuation_hdf5
