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


  subroutine write_kernels_strength_noise_hdf5()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_noise

#ifdef USE_HDF5
  use manager_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! offset array
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm_adj

  ! local parameters
  character(len=MAX_STRING_LEN) :: file_name
  integer :: info, comm

  ! gather the number of elements in each region
  call gather_all_all_singlei(size(sigma_kl_crust_mantle,4), offset_nspec_cm_adj, NPROCTOT_VAL)

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)
    ! create dataset
    call h5_create_dataset_gen('sigma_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  call h5_open_file_p_collect(file_name)

  ! write data
  call h5_write_dataset_collect_hyperslab('sigma_kernel', sigma_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

  ! close hdf5
  call h5_close_file_p()

#else

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_strength_noise_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_cm_ani_hdf5(alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
                                   betav_kl_crust_mantle,betah_kl_crust_mantle, &
                                   eta_kl_crust_mantle, &
                                   bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
                                   bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle, &
                                   Gc_prime_kl_crust_mantle, Gs_prime_kl_crust_mantle)


  use specfem_par
  use specfem_par_crustmantle

#ifdef USE_HDF5
  use manager_hdf5
#endif

  implicit none

  ! input Parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
      betav_kl_crust_mantle,betah_kl_crust_mantle, &
      eta_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
      bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      Gc_prime_kl_crust_mantle, Gs_prime_kl_crust_mantle

#ifdef USE_HDF5
  ! local parameters
  character(len=MAX_STRING_LEN) :: file_name
  integer :: info, comm

  ! offset array
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm_adj

  ! check if anything to do
  if (.not. ANISOTROPIC_KL) return

  ! gather the number of elements in each region
  !call gather_all_all_singlei(size(alphav_kl_crust_mantle,4), offset_nspec_cm_adj, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC_CRUST_MANTLE_ADJOINT, offset_nspec_cm_adj, NPROCTOT_VAL)

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! create or open file
    call h5_create_or_open_file(file_name)
    ! create dataset
    if (SAVE_TRANSVERSE_KL_ONLY) then
      call h5_create_dataset_gen('alphav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('alphah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('betav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('betah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('eta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)

      call h5_create_dataset_gen('bulk_c_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('bulk_betav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('bulk_betah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)

      call h5_create_dataset_gen('alpha_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('bulk_beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)

    else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
      call h5_create_dataset_gen('alphav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('alphah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('betav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('betah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)

      call h5_create_dataset_gen('bulk_c_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('bulk_betav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('bulk_betah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)

      call h5_create_dataset_gen('eta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)

      call h5_create_dataset_gen('Gc_prime_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('Gs_prime_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)

      ! check isotropic kernel
      if (.false.) then
        call h5_create_dataset_gen('alpha_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('bulk_beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      endif
      ! check anisotropic kernels
      if (.false.) then
        call h5_create_dataset_gen('A_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('C_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('L_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('N_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('F_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('Gc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('Gs_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('Jc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('Kc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('Mc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('Bc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('Hc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('Ec_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('Dc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      endif

    else
      ! fully anisotropic kernels
      call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('cijkl_kernel', (/21,NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 5, CUSTOM_REAL)

    endif

    ! close file
    call h5_close_file()

  endif ! myrank == 0

  ! synchronize all
  call synchronize_all()

  ! write data from all ranks
  call h5_open_file_p_collect(file_name)

  ! write data
  if (SAVE_TRANSVERSE_KL_ONLY) then
    call h5_write_dataset_collect_hyperslab('alphav_kernel', alphav_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('alphah_kernel', alphah_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('betav_kernel', betav_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('betah_kernel', betah_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('eta_kernel', eta_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

    call h5_write_dataset_collect_hyperslab('bulk_c_kernel', bulk_c_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('bulk_betav_kernel', bulk_betav_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('bulk_betah_kernel', bulk_betah_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

    call h5_write_dataset_collect_hyperslab('alpha_kernel', alpha_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('beta_kernel', beta_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('bulk_beta_kernel', bulk_beta_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

  else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
    call h5_write_dataset_collect_hyperslab('alphav_kernel', alphav_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('alphah_kernel', alphah_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('betav_kernel', betav_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('betah_kernel', betah_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

    call h5_write_dataset_collect_hyperslab('bulk_c_kernel', bulk_c_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('bulk_betav_kernel', bulk_betav_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('bulk_betah_kernel', bulk_betah_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

    call h5_write_dataset_collect_hyperslab('eta_kernel', eta_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

    call h5_write_dataset_collect_hyperslab('Gc_prime_kernel', Gc_prime_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('Gs_prime_kernel', Gs_prime_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

    ! check isotropic kernel
    if (.false.) then
      call h5_write_dataset_collect_hyperslab('alpha_kernel', alpha_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('beta_kernel', beta_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('bulk_beta_kernel', bulk_beta_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    endif

    ! check anisotropic kernels
    !if (.false.) then
    !  call h5_write_dataset_collect_hyperslab('A_kernel', A_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('C_kernel', C_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('L_kernel', L_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('N_kernel', N_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('F_kernel', F_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('Gc_kernel', Gc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('Gs_kernel', Gs_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('Jc_kernel', Jc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('Kc_kernel', Kc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('Mc_kernel', Mc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('Bc_kernel', Bc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('Hc_kernel', Hc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('Ec_kernel', Ec_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !  call h5_write_dataset_collect_hyperslab('Dc_kernel', Dc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    !endif

  else

    ! fully anisotropic kernels
    call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('cijkl_kernel', cijkl_kl_crust_mantle, &
                                            (/0,0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

  endif

  ! close hdf5
  call h5_close_file_p()

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = size(alphav_kl_crust_mantle,kind=4)
  idummy = size(alphah_kl_crust_mantle,kind=4)
  idummy = size(betav_kl_crust_mantle,kind=4)
  idummy = size(betah_kl_crust_mantle,kind=4)
  idummy = size(eta_kl_crust_mantle,kind=4)
  idummy = size(bulk_c_kl_crust_mantle,kind=4)
  idummy = size(bulk_beta_kl_crust_mantle,kind=4)
  idummy = size(bulk_betav_kl_crust_mantle,kind=4)
  idummy = size(bulk_betah_kl_crust_mantle,kind=4)
  idummy = size(Gc_prime_kl_crust_mantle,kind=4)
  idummy = size(Gs_prime_kl_crust_mantle,kind=4)

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_cm_ani_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_cm_iso_hdf5(mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle, &
                                    bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle)

  use specfem_par
  use specfem_par_crustmantle

#ifdef USE_HDF5
  use manager_hdf5
#endif

  implicit none

  ! Parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle, &
      bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle

#ifdef USE_HDF5
  ! offset array
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm_adj

  ! local parameters
  character(len=MAX_STRING_LEN) :: file_name
  integer :: info, comm

  ! checks if anything to do
  if (ANISOTROPIC_KL) return

  ! gather the number of elements in each region
  call gather_all_all_singlei(NSPEC_CRUST_MANTLE_ADJOINT, offset_nspec_cm_adj, NPROCTOT_VAL)

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)
    ! create dataset
    call h5_create_dataset_gen('rhonotprime_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
    call h5_create_dataset_gen('kappa_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
    call h5_create_dataset_gen('mu_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)

    call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
    call h5_create_dataset_gen('alpha_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
    call h5_create_dataset_gen('beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)

    call h5_create_dataset_gen('bulk_c_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
    call h5_create_dataset_gen('bulk_beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)

    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  call h5_open_file_p_collect(file_name)

  ! write data
  call h5_write_dataset_collect_hyperslab('rhonotprime_kernel', rhonotprime_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('kappa_kernel', kappa_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('mu_kernel', mu_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

  call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('alpha_kernel', alpha_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('beta_kernel', beta_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

  call h5_write_dataset_collect_hyperslab('bulk_c_kernel', bulk_c_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('bulk_beta_kernel', bulk_beta_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

  ! close hdf5
  call h5_close_file_p()

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = size(mu_kl_crust_mantle,kind=4)
  idummy = size(kappa_kl_crust_mantle,kind=4)
  idummy = size(rhonotprime_kl_crust_mantle,kind=4)
  idummy = size(bulk_c_kl_crust_mantle,kind=4)
  idummy = size(bulk_beta_kl_crust_mantle,kind=4)

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_cm_iso_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_oc_hdf5()

  use specfem_par
  use specfem_par_outercore

#ifdef USE_HDF5
  use manager_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! offset array
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_oc_adj

  ! local parameters
  character(len=MAX_STRING_LEN) :: file_name
  integer :: info, comm

  ! gather the number of elements in each region
  call gather_all_all_singlei(NSPEC_OUTER_CORE_ADJOINT, offset_nspec_oc_adj, NPROCTOT_VAL)

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)
    ! create dataset
    call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_oc_adj)/), 4, CUSTOM_REAL)
    call h5_create_dataset_gen('alpha_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_oc_adj)/), 4, CUSTOM_REAL)
    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  call h5_open_file_p_collect(file_name)

  ! write data
  call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_outer_core, &
                                          (/0,0,0,sum(offset_nspec_oc_adj(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('alpha_kernel', alpha_kl_outer_core, &
                                          (/0,0,0,sum(offset_nspec_oc_adj(0:myrank-1))/), H5_COL)

  ! close hdf5
  call h5_close_file_p()

#else

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_oc_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_ic_hdf5()

  use specfem_par
  use specfem_par_innercore

#ifdef USE_HDF5
  use manager_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! offset array
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_ic_adj

  ! local parameters
  character(len=MAX_STRING_LEN) :: file_name
  integer :: info, comm

  ! gather the number of elements in each region
  call gather_all_all_singlei(NSPEC_INNER_CORE_ADJOINT, offset_nspec_ic_adj, NPROCTOT_VAL)

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)
    ! create dataset
    call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_ic_adj)/), 4, CUSTOM_REAL)
    call h5_create_dataset_gen('alpha_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_ic_adj)/), 4, CUSTOM_REAL)
    call h5_create_dataset_gen('beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_ic_adj)/), 4, CUSTOM_REAL)
    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  call h5_open_file_p_collect(file_name)

  ! write data
  call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_inner_core, &
                                          (/0,0,0,sum(offset_nspec_ic_adj(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('alpha_kernel', alpha_kl_inner_core, &
                                          (/0,0,0,sum(offset_nspec_ic_adj(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('beta_kernel', beta_kl_inner_core, &
                                          (/0,0,0,sum(offset_nspec_ic_adj(0:myrank-1))/), H5_COL)

  ! close hdf5
  call h5_close_file_p()

#else

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_ic_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_boundary_kl_hdf5()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

#ifdef USE_HDF5
  use manager_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! offset array
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2d_moho
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2d_400
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2d_670
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2d_cmb
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2d_icb


  ! local parameters
  character(len=MAX_STRING_LEN) :: file_name
  integer :: info, comm

  if (.not. SAVE_KERNELS_BOUNDARY) return

  ! gather the number of elements in each region
  call gather_all_all_singlei(NSPEC2D_MOHO, offset_nspec2d_moho, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC2D_400, offset_nspec2d_400, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC2D_670, offset_nspec2d_670, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC2D_CMB, offset_nspec2d_cmb, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC2D_ICB, offset_nspec2d_icb, NPROCTOT_VAL)

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)
    ! create dataset
    if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
      call h5_create_dataset_gen('moho_kernel', (/0,0,sum(offset_nspec2d_moho(0:myrank-1))/), 3, CUSTOM_REAL)
    endif
    call h5_create_dataset_gen('d400_kernel', (/0,0,sum(offset_nspec2d_400(0:myrank-1))/), 3, CUSTOM_REAL)
    call h5_create_dataset_gen('d670_kernel', (/0,0,sum(offset_nspec2d_670(0:myrank-1))/), 3, CUSTOM_REAL)
    call h5_create_dataset_gen('CMB_kernel',  (/0,0,sum(offset_nspec2d_cmb(0:myrank-1))/), 3, CUSTOM_REAL)
    call h5_create_dataset_gen('ICB_kernel',  (/0,0,sum(offset_nspec2d_icb(0:myrank-1))/), 3, CUSTOM_REAL)

    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  call h5_open_file_p_collect(file_name)

  ! write data
  if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
    call h5_write_dataset_collect_hyperslab('moho_kernel', moho_kl, (/0,0,sum(offset_nspec2d_moho(0:myrank-1))/), H5_COL)
  endif
  call h5_write_dataset_collect_hyperslab('d400_kernel', d400_kl, (/0,0,sum(offset_nspec2d_400(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('d670_kernel', d670_kl, (/0,0,sum(offset_nspec2d_670(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('CMB_kernel', CMB_kl, (/0,0,sum(offset_nspec2d_CMB(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('ICB_kernel', ICB_kl, (/0,0,sum(offset_nspec2d_ICB(0:myrank-1))/), H5_COL)

  ! close hdf5
  call h5_close_file_p()

#else

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_boundary_kl_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_Hessian_hdf5()

  use specfem_par
  use specfem_par_crustmantle

#ifdef USE_HDF5
  use manager_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! offset array
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm_adj

  ! local parameters
  character(len=MAX_STRING_LEN) :: file_name
  integer :: info, comm

  ! gather the number of elements in each region
  call gather_all_all_singlei(NSPEC_CRUST_MANTLE_ADJOINT, offset_nspec_cm_adj, NPROCTOT_VAL)

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)
    ! create dataset
    call h5_create_dataset_gen('hess_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
    call h5_create_dataset_gen('hess_rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
    call h5_create_dataset_gen('hess_kappa_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
    call h5_create_dataset_gen('hess_mu_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm_adj)/), 4, CUSTOM_REAL)
    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  call h5_open_file_p_collect(file_name)

  ! write data
  call h5_write_dataset_collect_hyperslab('hess_kernel', hess_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('hess_rho_kernel', hess_rho_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('hess_kappa_kernel', hess_kappa_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('hess_mu_kernel', hess_mu_kl_crust_mantle, &
                                          (/0,0,0,sum(offset_nspec_cm_adj(0:myrank-1))/), H5_COL)

  ! close hdf5
  call h5_close_file_p()

#else

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_Hessian_hdf5

