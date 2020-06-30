!>\file gsd_chem_plume_wrapper.F90
!! This file is GSD Chemistry plume wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 05/2020

 module gsd_chem_plume_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use gsd_chem_config
!   use dep_dry_mod
!   use dep_wet_ls_mod
!   use gocart_settling_mod
!   use gocart_aerosols_mod
!   use gocart_dmsemis_mod
!   use gocart_chem_mod
!   use opt_mod
!   use seas_mod,        only : gocart_seasalt_driver
   use plume_data_mod
!   use dust_gocart_mod, only : gocart_dust_driver
!   use dust_afwa_mod,   only : gocart_dust_afwa_driver
!   use dust_fengsha_mod,only : gocart_dust_fengsha_driver
!   use dust_data_mod
   use plume_rise_mod
!   use vash_settling_mod

   implicit none

   private

   public :: gsd_chem_plume_wrapper_init, gsd_chem_plume_wrapper_run, gsd_chem_plume_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine gsd_chem_plume_wrapper_init()
      end subroutine gsd_chem_plume_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_gsd_chem_plume_wrapper_finalize Argument Table
!!
      subroutine gsd_chem_plume_wrapper_finalize()
      end subroutine gsd_chem_plume_wrapper_finalize

!> \defgroup gsd_chem_plume_group GSD Chem seas wrapper Module
!! This is the gsd chemistry
!>\defgroup gsd_chem_plume_wrapper GSD Chem seas wrapper Module  
!> \ingroup gsd_chem_plume_group
!! This is the GSD Chem seas wrapper Module
!! \section arg_table_gsd_chem_plume_wrapper_run Argument Table
!! \htmlinclude gsd_chem_plume_wrapper_run.html
!!
!>\section gsd_chem_plume_wrapper GSD Chemistry Scheme General Algorithm
!> @{
    subroutine gsd_chem_plume_wrapper_run(im, kte, kme, ktau, dt, garea, land, &
                   u10m, v10m, ustar, rlat, rlon, tskin,julian,xcosz,   &
                   rain_cpl, rainc_cpl, hf2d, pb2d,                     &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum,  &
                   w, exch, dqdt,                                       &
                   nsoil, smc, vegtype, soiltyp, sigmaf,jdate,idat,     & 
                   dswsfc, zorl,snow_cpl,                               &
                   dust_in,emi_in,emi2_in,fire_GBBEPx,fire_MODIS,       &
                   nseasalt,ntrac,                                      &
                   ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                     &
                   ntbc1,ntbc2,ntoc1,ntoc2,                             &
                   ntss1,ntss2,ntss3,ntss4,ntss5,                       &
                   ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,      &
                   gq0,ebu,tile_num,                                    &
                   errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau,nsoil,jdate(8),idat(8),tile_num
    integer,        intent(in) :: nseasalt,ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer,        intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    real(kind_phys),intent(in) :: dt,julian

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    integer, dimension(im), intent(in) :: land, vegtype, soiltyp        
    real(kind_phys), dimension(im,nsoil), intent(in) :: smc
    real(kind_phys), dimension(im,    5), intent(in) :: dust_in
    real(kind_phys), dimension(im,   10), intent(in) :: emi_in
    real(kind_phys), dimension(im,64, 3), intent(in) :: emi2_in
    real(kind_phys), dimension(im,    5), intent(in) :: fire_GBBEPx
    real(kind_phys), dimension(im,   13), intent(in) :: fire_MODIS
    real(kind_phys), dimension(im), intent(in) :: u10m, v10m, ustar,              &
                garea, rlat,rlon, tskin, rain_cpl, rainc_cpl,                     &
                hf2d, pb2d, sigmaf, dswsfc, zorl, snow_cpl, xcosz
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d,        &
                us3d, vs3d, spechum, w, exch, dqdt
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_ebu), intent(inout) :: ebu
!    real(kind_phys), dimension(im,nseasalt), intent(inout) :: ssem
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

!   integer, parameter :: num_moist=3, num_chem=20, num_emis_seas=5, num_emis_dust=5
    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy, u_phy, v_phy,       &
                     p_phy, z_at_w, dz8w, p8w, t8w, rho_phy, vvel, zmid,        &
                     exch_h, dqdti

    real(kind_phys), dimension(ims:im, jms:jme) :: u10, v10, ust, tsk,            &
                     xland, xlat, xlong, dxy, rcav, rnav, hfx, pbl

!>- sea salt & chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem )  ::                  &
                     var_rmv, dry_fall, tr_fall, sedim
    real(kind_phys), dimension(ims:im, 1, jms:jme, 1:num_emis_seas  ) :: emis_seas
    real(kind_phys), dimension(ims:im, jms:jme) :: seashelp

    integer :: ide, ime, ite, kde, julday

!   integer, parameter :: SEAS_OPT_DEFAULT = 1
    integer, parameter :: chem_in_opt = 0  ! 0 for coldstart, 1 for restart
    logical, parameter :: readrestart = .false.
    integer, parameter :: nvl_gocart  = 64  ! number of input levels from gocart file
   
!>- dust & chemistry variables
    real(kind_phys), dimension(ims:im, jms:jme, 3) ::    erod ! read from input?
    real(kind_phys), dimension(ims:im, jms:jme) :: ssm, rdrag, uthr, snowh  ! fengsha dust
    real(kind_phys), dimension(ims:im, jms:jme) :: vegfrac, rmol, gsw, znt, clayf, sandf, dms_0
    real(kind_phys), dimension(ims:im, nsoil, jms:jme) :: smois
    real(kind_phys), dimension(ims:im, 1:1, jms:jme, 1:num_emis_dust) :: emis_dust
    real(kind_phys), dimension(ims:im, 1:1, jms:jme, 1:5)             :: srce_dust
    real(kind_phys), dimension(ims:im, jms:jme) :: dusthelp
    integer,         dimension(ims:im, jms:jme) :: isltyp, ivgtyp

    integer :: current_month

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: pm10, pm2_5_dry, pm2_5_dry_ec
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: ac3, ahno3, anh3, asulf, cor3, h2oai, h2oaj, nu3
    real(kind_phys), dimension(ims:im, jms:jme) :: dep_vel_o3, e_co, ash_fall

!>- chemical background variables

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: oh_t, h2o2_t, no3_t
    real(kind_phys), dimension(ims:im, jms:jme) :: ttday, tcosz

!>- optical variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: relhum
    real(kind_phys), dimension(ims:im,         jms:jme) :: aod
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:nbands) :: extt
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:nbands) :: ssca
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:nbands) :: asympar
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:4) ::                  &
        tauaersw, gaersw, waersw, bscoefsw,                                        &
        l2aer,  l3aer, l4aer, l5aer, l6aer, l7aer           
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:16) :: tauaerlw
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_ext_coef) :: ext_coeff
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_bscat_coef) :: bscat_coeff
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_asym_par)   :: asym_par
!>- optical variables
    real(kind_phys), dimension(im) :: aod2d
    real(kind_phys), dimension(im, kte, 1:nbands) :: ext_cof, sscal, asymp

!>- plume variables
    ! -- buffers
   !real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_ebu) :: ebu
    real(kind_phys), dimension(ims:im, jms:jme, num_ebu_in) :: ebu_in
    real(kind_phys), dimension(ims:im, jms:jme) ::                             &
         mean_fct_agef, mean_fct_aggr, mean_fct_agsv, mean_fct_agtf,            &
         firesize_agef, firesize_aggr, firesize_agsv, firesize_agtf
    real(kind_phys), dimension(ims:im, jms:jme, num_frp_plume ) :: plume_frp
    real(kind_phys), dimension(ims:im, kms:kemit, jms:jme, 1:num_emis_ant) :: emis_ant
    real(kind_phys) :: dtstep
    integer,parameter :: plumerise_flag = 2  ! 1=MODIS, 2=GBBEPx
    logical :: call_plume, scale_fire_emiss
    logical, save :: firstfire = .true.
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

    ! -- output tracers
    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem) :: wet_dep
!*  real(kind_phys), dimension(ims:im, jms:jme, 1:nvl, 1:nbands) :: ext_cof, sscal, asymp
    real(kind_phys), dimension(ims:im, jms:jme, 1:kme) :: p10, pm25, ebu_oc 
    real(kind_phys), dimension(ims:im, jms:jme, 1:kme) :: oh_bg, h2o2_bg, no3_bg


!>-- local variables
    real(kind_phys) :: curr_secs
    real(kind_phys) :: factor, factor2, factor3
    logical :: call_gocart, call_radiation
    logical :: store_arrays
    integer :: nbegin, nv, nvv
    integer :: i, j, jp, k, kp, n
  

    errmsg = ''
    errflg = 0

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    ! -- set control flags
!    call_plume = (biomass_burn_opt == BURN_OPT_ENABLE) .and. (plumerisefire_frq
!    > 0)
!    if (call_plume) &
!       call_plume    = (mod(ktau, max(1, int(60*plumerisefire_frq/dt))) == 0) &
!                        .or. (ktau == 1) .or. firstfire
!    call_gocart      = (mod(ktau, call_chemistry) == 0) .or. (ktau == 1)
!    call_radiation   = (mod(int(curr_secs), max(1, 60*aer_ra_frq)) == 0) .or.
!    (ktau == 1)
    call_plume       = (biomass_burn_opt == BURN_OPT_ENABLE) .and. (plumerisefire_frq > 0)
    if (call_plume) &
       call_plume    = (mod(int(curr_secs), max(1, 60*plumerisefire_frq)) == 0)      &
                        .or. (ktau == 1)
    call_gocart      = (mod(ktau, call_chemistry) == 0) .or. (ktau == 1)
    call_radiation   = (mod(int(curr_secs), max(1, 60*aer_ra_frq)) == 0) .or. (ktau == 1)
    scale_fire_emiss = .false.

    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
     ! -- initialize buffers
!      ebu   = 0._kind_phys
    end if

!>- get ready for chemistry run
!   print*,'GSD vegtype,soiltyp,sigmaf',vegtype,soiltyp,sigmaf
    call gsd_chem_prep_plume(                                             &
        readrestart,chem_in_opt,ktau,dtstep,xcosz,                      &
        u10m,v10m,ustar,land,garea,rlat,rlon,tskin,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                 &
        exch,dqdt,                                                      &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                   &
        snow_cpl,dust_in,emi_in,emi2_in,                                &
        fire_GBBEPx,fire_MODIS,hf2d,pb2d,                               &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        t8w,exch_h,dqdti,                                               &
        z_at_w,vvel,zmid,                                               &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                                &
        ntbc1,ntbc2,ntoc1,ntoc2,                                        &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                  &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                 &
        ntrac,gq0,                                                      &
        num_chem, num_moist,num_ebu_in,                                 &
        call_gocart,nvl_gocart,                                         &
        ttday,tcosz,julday,                                         &
        plumerise_flag,num_plume_data,num_emis_ant,                     &
        emis_ant,ppm2ugkg,                                              &
        mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,        &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,        &
        moist,chem,plume_frp,ebu_in,                                    &
        smois,ivgtyp,isltyp,vegfrac,rmol,gsw,znt,hfx,pbl,               &
        relhum,snowh,clayf,rdrag,sandf,ssm,uthr,dms_0,erod,             &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)


    ! compute wild-fire plumes
    if (call_plume) then
      call plumerise_driver (ktau,dtstep,num_chem,num_ebu,num_ebu_in,   &
        ebu,ebu_in,                                                     &
        mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,        &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,        &
        'GOCART','BIOMASSB', t_phy,moist(:,:,:,p_qv),                   &
        rho_phy,vvel,u_phy,v_phy,p_phy,                                 &
        z_at_w,scale_fire_emiss,plume_frp,plumerise_flag,               &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte                                )
    end if

    ! -- add biomass burning emissions at every timestep
    if (biomass_burn_opt == BURN_OPT_ENABLE) then
      jp = jte
      factor3 = 0._kind_phys
      select case (plumerise_flag)
        case (FIRE_OPT_MODIS)
          factor3 = 4.828e-04_kind_phys/60.
          kp = kte    ! full column
        case (FIRE_OPT_GBBEPx)
          factor3 = 1.e-03_kind_phys * mwdry / mw_so2_aer
          if (plumerisefire_frq > 0) then
            kp = kte  ! full column
          else
            kp = kts  ! surface only
          end if
        case default
          ! -- no further options available, skip this step
          jp = jts - 1
      end select

      if (kp == kts) then
        ! -- only include surface emissions
        k = kts
        do j = jts, jp
          do i = its, ite
            ! -- factor for pm emissions, factor2 for burn emissions
            factor  = dt*rri(i,k,j)/dz8w(i,k,j)
            factor2 = factor * factor3
            chem(i,k,j,p_oc1) = chem(i,k,j,p_oc1) + factor  * ebu_in(i,j,p_ebu_in_oc  )
            chem(i,k,j,p_bc1) = chem(i,k,j,p_bc1) + factor  * ebu_in(i,j,p_ebu_in_bc  )
            chem(i,k,j,p_p25) = chem(i,k,j,p_p25) + factor  * ebu_in(i,j,p_ebu_in_pm25)
            chem(i,k,j,p_p10) = chem(i,k,j,p_p10) + factor  * ebu_in(i,j,p_ebu_in_pm10)
            chem(i,k,j,p_so2) = chem(i,k,j,p_so2) + factor2 * ebu_in(i,j,p_ebu_in_so2 )
          end do
        end do

      else
        ! -- use full-column emissions
        do j = jts, jp
          do k = kts, kp
            do i = its, ite
              ! -- factor for pm emissions, factor2 for burn emissions
              factor  = dt*rri(i,k,j)/dz8w(i,k,j)
              factor2 = factor * factor3
              chem(i,k,j,p_oc1) = chem(i,k,j,p_oc1) + factor  * ebu(i,k,j,p_ebu_oc  )
              chem(i,k,j,p_bc1) = chem(i,k,j,p_bc1) + factor  * ebu(i,k,j,p_ebu_bc  )
              chem(i,k,j,p_p25) = chem(i,k,j,p_p25) + factor  * ebu(i,k,j,p_ebu_pm25)
              chem(i,k,j,p_p10) = chem(i,k,j,p_p10) + factor  * ebu(i,k,j,p_ebu_pm10)
              chem(i,k,j,p_so2) = chem(i,k,j,p_so2) + factor2 * ebu(i,k,j,p_ebu_so2 )
            end do
          end do
        end do
      end if

    end if


    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntso2  )=ppm2ugkg(p_so2   ) * max(epsilc,chem(i,k,1,p_so2))
       gq0(i,k,ntsulf )=ppm2ugkg(p_sulf  ) * max(epsilc,chem(i,k,1,p_sulf))
       gq0(i,k,ntdms  )=ppm2ugkg(p_dms   ) * max(epsilc,chem(i,k,1,p_dms))
       gq0(i,k,ntmsa  )=ppm2ugkg(p_msa   ) * max(epsilc,chem(i,k,1,p_msa))
       gq0(i,k,ntpp25 )=ppm2ugkg(p_p25   ) * max(epsilc,chem(i,k,1,p_p25))
       gq0(i,k,ntbc1  )=ppm2ugkg(p_bc1   ) * max(epsilc,chem(i,k,1,p_bc1))
       gq0(i,k,ntbc2  )=ppm2ugkg(p_bc2   ) * max(epsilc,chem(i,k,1,p_bc2))
       gq0(i,k,ntoc1  )=ppm2ugkg(p_oc1   ) * max(epsilc,chem(i,k,1,p_oc1))
       gq0(i,k,ntoc2  )=ppm2ugkg(p_oc2   ) * max(epsilc,chem(i,k,1,p_oc2))
       gq0(i,k,ntdust1)=ppm2ugkg(p_dust_1) * max(epsilc,chem(i,k,1,p_dust_1))
       gq0(i,k,ntdust2)=ppm2ugkg(p_dust_2) * max(epsilc,chem(i,k,1,p_dust_2))
       gq0(i,k,ntdust3)=ppm2ugkg(p_dust_3) * max(epsilc,chem(i,k,1,p_dust_3))
       gq0(i,k,ntdust4)=ppm2ugkg(p_dust_4) * max(epsilc,chem(i,k,1,p_dust_4))
       gq0(i,k,ntdust5)=ppm2ugkg(p_dust_5) * max(epsilc,chem(i,k,1,p_dust_5))
       gq0(i,k,ntss1  )=ppm2ugkg(p_seas_1) * max(epsilc,chem(i,k,1,p_seas_1))
       gq0(i,k,ntss2  )=ppm2ugkg(p_seas_2) * max(epsilc,chem(i,k,1,p_seas_2))
       gq0(i,k,ntss3  )=ppm2ugkg(p_seas_3) * max(epsilc,chem(i,k,1,p_seas_3))
       gq0(i,k,ntss4  )=ppm2ugkg(p_seas_4) * max(epsilc,chem(i,k,1,p_seas_4))
       gq0(i,k,ntss5  )=ppm2ugkg(p_seas_5) * max(epsilc,chem(i,k,1,p_seas_5))
       gq0(i,k,ntpp10 )=ppm2ugkg(p_p10   ) * max(epsilc,chem(i,k,1,p_p10))
     enddo
    enddo


!   call gsd_chem_post() ! postprocessing for diagnostics

!
   end subroutine gsd_chem_plume_wrapper_run
!> @}

  subroutine gsd_chem_prep_plume(                                        &
        readrestart,chem_in_opt,ktau,dtstep,xcosz,                     &
        u10m,v10m,ustar,land,garea,rlat,rlon,ts2d,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                &
        exch,dqdt,                                                     &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                  &
        snow_cpl,dust_in,emi_in,emi2_in,                               &
        fire_GBBEPx,fire_MODIS,hf2d,pb2d,                              &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                          &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                  &
        t8w,exch_h,dqdti,                                              &
        z_at_w,vvel,zmid,                                              &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                               &
        ntbc1,ntbc2,ntoc1,ntoc2,                                       &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                 &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                &
        ntrac,gq0,                                                     &
        num_chem, num_moist,num_ebu_in,                                &
        call_gocart,nvl_gocart,                                        &
        ttday,tcosz,julday,                                        &
        plumerise_flag,num_plume_data,num_emis_ant,                    &
        emis_ant,ppm2ugkg,                                             &
        mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,       &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,       &
        moist,chem,plumedist,ebu_in,                                   &
        smois,ivgtyp,isltyp,vegfrac,rmol,gsw,znt,hfx,pbl,              &
        relhum,snowh,clayf,rdrag,sandf,ssm,uthr,dms_0,erod,            &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    logical, intent(in) :: readrestart
    integer, intent(in) :: chem_in_opt, ktau, julday
    real(kind=kind_phys), intent(in) :: dtstep

    !FV3 input variables
    integer, intent(in) :: nsoil
    integer, dimension(ims:ime), intent(in) :: land, vegtype, soiltyp
    integer, intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer, intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         u10m, v10m, ustar, garea, rlat, rlon, ts2d, sigmaf, dswsfc,       &
         zorl, snow_cpl, hf2d, pb2d,xcosz
    real(kind=kind_phys), dimension(ims:ime, nsoil),   intent(in) :: smc 
    real(kind=kind_phys), dimension(ims:ime,     5),   intent(in) :: dust_in
    real(kind=kind_phys), dimension(ims:ime,    10),   intent(in) :: emi_in
    real(kind=kind_phys), dimension(ims:ime, 64, 3),   intent(in) :: emi2_in
    real(kind=kind_phys), dimension(ims:ime,     5),   intent(in) :: fire_GBBEPx
    real(kind=kind_phys), dimension(ims:ime,    13),   intent(in) :: fire_MODIS
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) ::     &
         pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,exch,dqdt
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_moist, num_ebu_in,                &
                           plumerise_flag, num_plume_data, num_emis_ant,   &
                           nvl_gocart
    logical,intent(in) ::  call_gocart
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg

    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in),intent(out) :: ebu_in
    real(kind_phys), dimension(ims:ime, kms:kemit, jms:jme, num_emis_ant), intent(inout) :: emis_ant
    
    integer,dimension(ims:ime, jms:jme), intent(out) :: isltyp, ivgtyp
    real(kind_phys), dimension(ims:ime, jms:jme, 3), intent(inout) :: erod
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, u_phy, v_phy, p_phy, rho_phy, dz8w, p8w, t8w, vvel, zmid,         &
         exch_h,dqdti
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         u10, v10, ust, tsk, xland, xlat, xlong, dxy, vegfrac, rmol, gsw, znt, hfx,    &
         pbl, snowh, clayf, rdrag, sandf, ssm, uthr, dms_0, ttday, tcosz
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w, relhum
    real(kind_phys), dimension(ims:ime, nsoil, jms:jme), intent(out) :: smois
    real(kind_phys), dimension(ims:ime, jms:jme, num_frp_plume), intent(out) :: plumedist
    real(kind_phys), dimension(ims:ime, jms:jme   ), intent(out) ::                    &
                   mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,            &
                   firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr       
    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in) :: emiss_ab
    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in) :: emiss_abu
    real(kind_phys), dimension(ims:ime, jms:jme, num_plume_data) :: plume
    real(kind_phys), parameter :: frac_so2_ant = 0.5_kind_phys  ! antropogenic so2 fraction
    real(kind_phys), parameter :: frp2plume = 1.e+06_kind_phys  ! FRP-to-plume conversion factor
    real(kind_phys), parameter :: frpc  = 1.e+09_kind_phys      ! FRP conversion factor
    real(kind_phys), dimension(nvl_gocart) :: p_gocart

    ! -- local variables
!   real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(kind_phys) ::  factor,factor2,pu,pl,aln,pwant
    real(kind_phys) ::  xhour,xmin,gmtp,xlonn,xtime,real_time
    real(kind_phys), DIMENSION (1,1) :: sza,cosszax
    integer i,ip,j,jp,k,kp,kk,kkp,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour

    p_gocart = (/ 1000., 992.5, 985., 977.5, 970., 955., 940., 925., 910.,               &
          895., 880., 865., 850., 825., 800., 775., 750., 712.5,  675., 637.5, 600.,     &
          562.5, 525., 487.5, 450., 412.5, 375., 337.5, 288.08, 244.88, 208.15, 176.93,  &
          150.39, 127.84, 108.66, 92.37, 78.51, 66.6, 56.39, 47.64, 40.18, 33.81, 28.37, &
          23.73, 19.79,  16.46, 13.64, 11.28, 9.29, 7.62, 6.22, 5.05, 4.08, 3.28, 2.62,  &
          2.08, 1.65, 1.3, 1.02, 0.8, 0.62, 0.48, 0.37, 0.28 /)

    ! -- initialize output arrays
    ebu_in         = 0._kind_phys
    isltyp         = 0._kind_phys
    ivgtyp         = 0._kind_phys
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    u_phy          = 0._kind_phys
    v_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    t8w            = 0._kind_phys
    vvel           = 0._kind_phys
    zmid           = 0._kind_phys
    exch_h         = 0._kind_phys
    dqdti          = 0._kind_phys
    u10            = 0._kind_phys
    v10            = 0._kind_phys
    ust            = 0._kind_phys
    tsk            = 0._kind_phys
    xland          = 0._kind_phys
    xlat           = 0._kind_phys
    xlong          = 0._kind_phys
    dxy            = 0._kind_phys
    vegfrac        = 0._kind_phys
    rmol           = 0._kind_phys
    gsw            = 0._kind_phys
    znt            = 0._kind_phys
    hfx            = 0._kind_phys
    pbl            = 0._kind_phys
    snowh          = 0._kind_phys
    clayf          = 0._kind_phys
    rdrag          = 0._kind_phys
    sandf          = 0._kind_phys 
    ssm            = 0._kind_phys
    uthr           = 0._kind_phys
    dms_0          = 0._kind_phys
    ttday          = 0._kind_phys
    tcosz          = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys
    relhum         = 0._kind_phys

    ! -- initialize fire emissions
    plume          = 0._kind_phys
    plumedist      = 0._kind_phys
    mean_fct_agtf  = 0._kind_phys
    mean_fct_agef  = 0._kind_phys
    mean_fct_agsv  = 0._kind_phys
    mean_fct_aggr  = 0._kind_phys
    firesize_agtf  = 0._kind_phys
    firesize_agef  = 0._kind_phys
    firesize_agsv  = 0._kind_phys
    firesize_aggr  = 0._kind_phys


    do i=its,ite
     u10  (i,1)=u10m (i)
     v10  (i,1)=v10m (i)
     tsk  (i,1)=ts2d (i)
     ust  (i,1)=ustar(i)
     dxy  (i,1)=garea(i)
     xland(i,1)=real(land(i))
     xlat (i,1)=rlat(i)*180./pi
     xlong(i,1)=rlon(i)*180./pi
     gsw  (i,1)=dswsfc(i)
     znt  (i,1)=zorl(i)*0.01
     hfx  (i,1)=hf2d(i)
     pbl  (i,1)=pb2d(i)
     snowh(i,1)=snow_cpl(i)*0.001
     clayf(i,1)=dust_in(i,1)
     rdrag(i,1)=dust_in(i,2)
     sandf(i,1)=dust_in(i,3)
     ssm  (i,1)=dust_in(i,4)
     uthr (i,1)=dust_in(i,5)
     ivgtyp (i,1)=vegtype(i)
     isltyp (i,1)=soiltyp(i)
     vegfrac(i,1)=sigmaf (i)
     dms_0(i,1  )=emi_in(i,7) ! --dm0
     erod (i,1,1)=emi_in(i,8) ! --ero1
     erod (i,1,2)=emi_in(i,9) ! --ero2
     erod (i,1,3)=emi_in(i,10)! --ero3
    enddo
   
    rmol=0.

    do k=1,nsoil
     do j=jts,jte
      do i=its,ite
       smois(i,k,j)=smc(i,k)
      enddo
     enddo
    enddo

    if (ktau <= 1) then
      emis_ant = 0.
     !emis_vol = 0.
    end if

    do j=jts,jte
      jp = j - jts + 1
      do i=its,ite
         ip = i - its + 1
         z_at_w(i,kts,j)=max(0.,ph3d(ip,1)/g)
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=abs(ph3d(ip,kp+1)-ph3d(ip,kp))/g
          z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          p8w(i,k,j)=pr3d(ip,kp)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kk=min(k,kte)
        kkp = kk - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
          t_phy(i,k,j)=tk3d(ip,kkp)
          p_phy(i,k,j)=prl3d(ip,kkp)
          u_phy(i,k,j)=us3d(ip,kkp)
          dqdti(i,k,j)=dqdt(ip,kkp)
          v_phy(i,k,j)=vs3d(ip,kkp)
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          vvel(i,k,j)=-w(ip,kkp)*rri(i,k,j)/g 
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=gq0(ip,kkp,p_atm_shum)
          if (t_phy(i,k,j) > 265.) then
            moist(i,k,j,2)=gq0(ip,kkp,p_atm_cldq)
            moist(i,k,j,3)=0.
            if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          else
            moist(i,k,j,2)=0.
            moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldq)
            if(moist(i,k,j,3) < 1.e-8)moist(i,k,j,3)=0.
          endif
          relhum(i,k,j) = .95
          relhum(i,k,j) = MIN( .95, moist(i,k,j,1) / &
            (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
            (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
          relhum(i,k,j)=max(0.1,relhum(i,k,j))
          !--
          zmid(i,k,j)=phl3d(ip,kkp)/g
        enddo
      enddo
    enddo

    ! -- the imported atmospheric heat diffusivity is only available up to kte-1
    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte-1
        kkp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          exch_h(i,k,j)=exch(ip,kkp)
        enddo
      enddo
    enddo

    do j=jts,jte
      do k=2,kte
        do i=its,ite
          t8w(i,k,j)=.5*(t_phy(i,k,j)+t_phy(i,k-1,j))
        enddo
      enddo
    enddo

    ! -- only used in phtolysis....
    do j=jts,jte
      do i=its,ite
        t8w(i,1,j)=t_phy(i,1,j)
        t8w(i,kte+1,j)=t_phy(i,kte,j)
      enddo
    enddo

    ! -- fire
    emiss_ab  = 0.   ! background
    emiss_abu = 0.   ! fire emission
    do j=jts,jte
     do i=its,ite
      emiss_ab(i,j,p_e_bc)   =emi_in(i,1)
      emiss_ab(i,j,p_e_oc)   =emi_in(i,2)
      emiss_ab(i,j,p_e_sulf) =emi_in(i,3)
      emiss_ab(i,j,p_e_pm_25)=emi_in(i,4)
      emiss_ab(i,j,p_e_so2)  =emi_in(i,5)
      emiss_ab(i,j,p_e_pm_10)=emi_in(i,6)
     enddo
    enddo

    !print*,'hli ',plumerise_flag,FIRE_OPT_MODIS,FIRE_OPT_GBBEPx
    select case (plumerise_flag)
      case (FIRE_OPT_MODIS)
        do j=jts,jte
         do i=its,ite
          emiss_abu(i,j,p_e_bc)   =fire_MODIS(i,1)
          emiss_abu(i,j,p_e_oc)   =fire_MODIS(i,2)
          emiss_abu(i,j,p_e_pm_25)=fire_MODIS(i,3)
          emiss_abu(i,j,p_e_so2)  =fire_MODIS(i,4)
          emiss_abu(i,j,p_e_pm_10)=fire_MODIS(i,5)
          plume(i,j,1)            =fire_MODIS(i,6)
          plume(i,j,2)            =fire_MODIS(i,7)
          plume(i,j,3)            =fire_MODIS(i,8)
          plume(i,j,4)            =fire_MODIS(i,9)
          plume(i,j,5)            =fire_MODIS(i,10)
          plume(i,j,6)            =fire_MODIS(i,11)
          plume(i,j,7)            =fire_MODIS(i,12)
          plume(i,j,8)            =fire_MODIS(i,13)
         enddo
        enddo
      case (FIRE_OPT_GBBEPx)
        do j=jts,jte
         do i=its,ite
          emiss_abu(i,j,p_e_bc)   =fire_GBBEPx(i,1)
          emiss_abu(i,j,p_e_oc)   =fire_GBBEPx(i,2)
          emiss_abu(i,j,p_e_pm_25)=fire_GBBEPx(i,3)
          emiss_abu(i,j,p_e_so2)  =fire_GBBEPx(i,4)
          plume(i,j,1)            =fire_GBBEPx(i,5)
         enddo
        enddo
!        print*,'hli GBBEPx plume',maxval(plume(:,:,1))
      case default
          ! -- no further option available
    end select


    factor=0.
    jmax=0
    jmaxi=0
    k=kts
    if (p_bc2 > 1) then
      do j=jts,jte
        do i=its,ite
          emis_ant(i,k,j,p_e_bc)=emiss_ab(i,j,p_e_bc)
          emis_ant(i,k,j,p_e_oc)=emiss_ab(i,j,p_e_oc) + emiss_ab(i,j,p_e_pm_25)
          emis_ant(i,k,j,p_e_sulf)=emiss_ab(i,j,p_e_sulf)
          emis_ant(i,k,j,p_e_so2)=frac_so2_ant * emiss_ab(i,j,p_e_so2)
          emis_ant(i,k,j,p_e_dms)= 0. !emiss_ab(j,p_e_dms)
          emis_ant(i,k,j,p_e_pm_25)=emiss_ab(i,j,p_e_pm_25)
          emis_ant(i,k,j,p_e_pm_10)=emiss_ab(i,j,p_e_pm_10)


          ebu_in(i,j,p_ebu_in_pm10)=emiss_abu(i,j,p_e_pm_10)
          ebu_in(i,j,p_ebu_in_dms)= 0._kind_phys

          select case (plumerise_flag)
            case (FIRE_OPT_MODIS)
              ebu_in(i,j,p_ebu_in_oc)   = emiss_abu(i,j,p_e_oc)
              ebu_in(i,j,p_ebu_in_bc)   = emiss_abu(i,j,p_e_bc)
              ebu_in(i,j,p_ebu_in_pm25) = emiss_abu(i,j,p_e_pm_25)
              ebu_in(i,j,p_ebu_in_so2)  = emiss_abu(i,j,p_e_so2)
              mean_fct_agtf(i,j)=plume(i,j,1)
              mean_fct_agef(i,j)=plume(i,j,2)
              mean_fct_agsv(i,j)=plume(i,j,3)
              mean_fct_aggr(i,j)=plume(i,j,4)
              firesize_agtf(i,j)=plume(i,j,5)
              firesize_agef(i,j)=plume(i,j,6)
              firesize_agsv(i,j)=plume(i,j,7)
              firesize_aggr(i,j)=plume(i,j,8)
            case (FIRE_OPT_GBBEPx)
              ebu_in(i,j,p_ebu_in_oc)   = frpc * (emiss_abu(i,j,p_e_pm_25) - emiss_abu(i,j,p_e_bc))
              ebu_in(i,j,p_ebu_in_bc)   = frpc * emiss_abu(i,j,p_e_bc)
              ebu_in(i,j,p_ebu_in_pm25) = frpc * (emiss_abu(i,j,p_e_pm_25) - emiss_abu(i,j,p_e_bc) - emiss_abu(i,j,p_e_oc))
              ebu_in(i,j,p_ebu_in_so2)  = frpc * emiss_abu(i,j,p_e_so2)
              plumedist(i,j,p_frp_flam_frac) = flaming(catb(ivgtyp(i,j)))
              plumedist(i,j,p_frp_mean     ) = frp2plume * plume(i,j,1)
              plumedist(i,j,p_frp_std      ) = 0.3_kind_phys   * frp2plume * plume(i,j,1)
              plumedist(i,j,p_frp_mean_size) = msize(ivgtyp(i,j)) * frp2plume * plume(i,j,1)
              plumedist(i,j,p_frp_std_size ) = 0.5_kind_phys * plumedist(i,j,p_frp_mean_size)
            case default
              ! -- no further option available
          end select
        enddo
      enddo
    endif
!    print*,'hli plumedist(:,:,p_frp_mean)',maxval(plumedist(:,:,p_frp_mean))

 
    do k=kms,kte
     do i=ims,ime
       chem(i,k,jts,p_so2   )=max(epsilc,gq0(i,k,ntso2  )/ppm2ugkg(p_so2))
       chem(i,k,jts,p_sulf  )=max(epsilc,gq0(i,k,ntsulf )/ppm2ugkg(p_sulf))
       chem(i,k,jts,p_dms   )=max(epsilc,gq0(i,k,ntdms  )/ppm2ugkg(p_dms))
       chem(i,k,jts,p_msa   )=max(epsilc,gq0(i,k,ntmsa  )/ppm2ugkg(p_msa))
       chem(i,k,jts,p_p25   )=max(epsilc,gq0(i,k,ntpp25 )/ppm2ugkg(p_p25))
       chem(i,k,jts,p_bc1   )=max(epsilc,gq0(i,k,ntbc1  )/ppm2ugkg(p_bc1))
       chem(i,k,jts,p_bc2   )=max(epsilc,gq0(i,k,ntbc2  )/ppm2ugkg(p_bc2))
       chem(i,k,jts,p_oc1   )=max(epsilc,gq0(i,k,ntoc1  )/ppm2ugkg(p_oc1))
       chem(i,k,jts,p_oc2   )=max(epsilc,gq0(i,k,ntoc2  )/ppm2ugkg(p_oc2))
       chem(i,k,jts,p_dust_1)=max(epsilc,gq0(i,k,ntdust1)/ppm2ugkg(p_dust_1))
       chem(i,k,jts,p_dust_2)=max(epsilc,gq0(i,k,ntdust2)/ppm2ugkg(p_dust_2))
       chem(i,k,jts,p_dust_3)=max(epsilc,gq0(i,k,ntdust3)/ppm2ugkg(p_dust_3))
       chem(i,k,jts,p_dust_4)=max(epsilc,gq0(i,k,ntdust4)/ppm2ugkg(p_dust_4))
       chem(i,k,jts,p_dust_5)=max(epsilc,gq0(i,k,ntdust5)/ppm2ugkg(p_dust_5))
       chem(i,k,jts,p_seas_1)=max(epsilc,gq0(i,k,ntss1  )/ppm2ugkg(p_seas_1))
       chem(i,k,jts,p_seas_2)=max(epsilc,gq0(i,k,ntss2  )/ppm2ugkg(p_seas_2))
       chem(i,k,jts,p_seas_3)=max(epsilc,gq0(i,k,ntss3  )/ppm2ugkg(p_seas_3))
       chem(i,k,jts,p_seas_4)=max(epsilc,gq0(i,k,ntss4  )/ppm2ugkg(p_seas_4))
       chem(i,k,jts,p_seas_5)=max(epsilc,gq0(i,k,ntss5  )/ppm2ugkg(p_seas_5))
       chem(i,k,jts,p_p10   )=max(epsilc,gq0(i,k,ntpp10 )/ppm2ugkg(p_p10))
     enddo
    enddo

    !
    ! -- gocart background fields only if gocart is called
    !

!   emis_ant=0.
    nv=1
    k=kts
    factor2=0.
    factor=0.
    if (p_bc2 > 1)then
      if (chem_opt == CHEM_OPT_GOCART) then
        do j=jts,jte
          do i=its,ite
            factor=dtstep*rri(i,k,j)/dz8w(i,k,j)
            factor2=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
            chem(i,k,j,p_bc1)=chem(i,k,j,p_bc1)+emis_ant(i,k,j,p_e_bc)*factor
            chem(i,k,j,p_oc1)=chem(i,k,j,p_oc1)+emis_ant(i,k,j,p_e_oc)*factor
            chem(i,k,j,p_p25)=chem(i,k,j,p_p25)+emis_ant(i,k,j,p_e_pm_25)*factor
            chem(i,k,j,p_p10)=chem(i,k,j,p_p10)+emis_ant(i,k,j,p_e_pm_10)*factor
            chem(i,k,j,p_sulf)=chem(i,k,j,p_sulf)+emis_ant(i,k,j,p_e_sulf)*factor
            chem(i,k,j,p_so2)=chem(i,k,j,p_so2)+emis_ant(i,k,j,p_e_so2)*factor2
          enddo
        enddo
      endif
    else if (p_tr2 > 1)then    !co2 here
      do j=jts,jte
        do i=its,ite
!           factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
          !factor2=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
          !chem(i,k,j,p_tr1)=chem(i,k,j,p_tr1)+emis_ant(i,k,j,p_e_tr1)*factor2
          !chem(i,k,j,p_tr2)=chem(i,k,j,p_tr2)+emis_ant(i,k,j,p_e_tr2)*factor2
        enddo
      enddo
    else if ((p_tr2 > 1) .and. (p_bc2 > 1))then
      !call chem_rc_set(CHEM_RC_FAILURE, msg="Inconsistent options detected.", &
      !  file=__FILE__, line=__LINE__, rc=rc)
      return
    endif

    ! -- real-time application, keeping eruption constant
!
!    if (ktau <= 2) then
!      ! -- volcanic emissions
!      if (num_emis_vol > 0) then
!         !------------------
!      end if
!     end if

  end subroutine gsd_chem_prep_plume

!> @}
  end module gsd_chem_plume_wrapper
