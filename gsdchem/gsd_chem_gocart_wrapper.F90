!>\file gsd_chem_gocart_wrapper.F90
!! This file is GSD Chemistry gocart wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 06/2020

 module gsd_chem_gocart_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use gsd_chem_config
!   use dep_dry_mod
!   use dep_wet_ls_mod
!   use gocart_settling_mod
   use gocart_aerosols_mod
!   use gocart_dmsemis_mod
   use gocart_chem_mod
!   use opt_mod
!   use seas_mod,        only : gocart_seasalt_driver
!   use dust_gocart_mod, only : gocart_dust_driver
!   use dust_afwa_mod,   only : gocart_dust_afwa_driver
!   use dust_fengsha_mod,only : gocart_dust_fengsha_driver
   use dust_data_mod
!   use plume_rise_mod
!   use vash_settling_mod

   implicit none

   private

   public :: gsd_chem_gocart_wrapper_init, gsd_chem_gocart_wrapper_run, gsd_chem_gocart_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine gsd_chem_gocart_wrapper_init()
      end subroutine gsd_chem_gocart_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_gsd_chem_gocart_wrapper_finalize Argument Table
!!
      subroutine gsd_chem_gocart_wrapper_finalize()
      end subroutine gsd_chem_gocart_wrapper_finalize

!> \defgroup gsd_chem_group GSD Chem driver Module
!! This is the gsd chemistry
!>\defgroup gsd_chem_gocart_wrapper GSD Chem driver Module  
!> \ingroup gsd_chem_group
!! This is the GSD Chem driver Module
!! \section arg_table_gsd_chem_gocart_wrapper_run Argument Table
!! \htmlinclude gsd_chem_gocart_wrapper_run.html
!!
!>\section gsd_chem_gocart_wrapper GSD Chemistry Scheme General Algorithm
!> @{
    subroutine gsd_chem_gocart_wrapper_run(im, kte, kme, ktau, dt, garea, land, &
                   u10m, v10m, ustar, rlat, rlon, tskin,julian,xcosz,   &
                   rain_cpl, rainc_cpl, hf2d, pb2d,                     &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum,  &
                   w, exch, dqdt,                                       &
                   nsoil, smc, vegtype, soiltyp, sigmaf,jdate,idat,     & 
                   dswsfc, zorl,snow_cpl,                               &
                   emi2_in,       &
                   nseasalt,ntrac,                                      &
                   ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                     &
                   ntbc1,ntbc2,ntoc1,ntoc2,                             &
                   ntss1,ntss2,ntss3,ntss4,ntss5,                       &
                   ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,      &
                   gq0,tile_num,                                    &
                   cplchm_rad_opt,lmk,faersw_cpl,                       &
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
    real(kind_phys), dimension(im,64, 3), intent(in) :: emi2_in
    real(kind_phys), dimension(im), intent(in) :: u10m, v10m, ustar,              &
                garea, rlat,rlon, tskin, rain_cpl, rainc_cpl,                     &
                hf2d, pb2d, sigmaf, dswsfc, zorl, snow_cpl, xcosz
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d,        &
                us3d, vs3d, spechum, w, exch, dqdt
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0
    integer,        intent(in) :: lmk
    real(kind_phys), dimension(im, lmk, 14, 3),intent(inout) :: faersw_cpl
    logical, intent(in) :: cplchm_rad_opt
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
    real(kind_phys), dimension(ims:im, jms:jme) :: vegfrac, rmol, gsw, znt
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
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: backg_oh,backg_h2o2,backg_no3

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
    real(kind_phys) :: dtstep, gmt
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

    ! -- output tracers
    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem) :: wet_dep
!*  real(kind_phys), dimension(ims:im, jms:jme, 1:nvl, 1:nbands) :: ext_cof, sscal, asymp
    real(kind_phys), dimension(ims:im, jms:jme, 1:kme) :: p10, pm25!, ebu_oc 
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

!    print*,'hli test1 ktau',ktau
    h2oai = 0.
    h2oaj = 0.
    nu3   = 0.
    ac3   = 0.
    cor3  = 0.
    asulf = 0.
    ahno3 = 0.
    anh3  = 0.
    e_co  = 0.
    dep_vel_o3 = 0.
    ash_fall   = 0.
    extt =0.
    ssca   =0.
    asympar=0.


    gmt = real(idat(5))
    julday = real(julian)                                       

    current_month=jdate(2)
    curr_secs = ktau * dt

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
    call_gocart      = (mod(ktau, call_chemistry) == 0) .or. (ktau == 1)
    call_radiation   = (mod(int(curr_secs), max(1, 60*aer_ra_frq)) == 0) .or. (ktau == 1)

    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
     ! -- initialize buffers
!      ebu   = 0._kind_phys
    end if

    ! -- compute incremental convective and large-scale rainfall
    do i=its,ite
     rcav(i,1)=max(rainc_cpl(i)*1000.              , 0.) ! meter to mm
     rnav(i,1)=max((rain_cpl(i)-rainc_cpl(i))*1000., 0.) ! meter to mm
    enddo

!!!

!>- get ready for chemistry run
!   print*,'GSD vegtype,soiltyp,sigmaf',vegtype,soiltyp,sigmaf
    call gsd_chem_gocart_prep(                                             &
        readrestart,chem_in_opt,ktau,dtstep,xcosz,                      &
        u10m,v10m,ustar,land,garea,rlat,rlon,tskin,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                 &
        exch,dqdt,                                                      &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                   &
        snow_cpl,emi2_in,                                &
        hf2d,pb2d,                               &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        t8w,exch_h,dqdti,                                               &
        z_at_w,vvel,zmid,                                               &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                                &
        ntbc1,ntbc2,ntoc1,ntoc2,                                        &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                  &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                 &
        ntrac,gq0,                                                      &
        num_chem, num_moist,                                 &
        call_gocart,nvl_gocart,                                         &
        ttday,tcosz,gmt,julday,                                         &
        backg_oh,backg_h2o2,backg_no3,                                  &
        ppm2ugkg,                                              &
        moist,chem,                                   &
        smois,ivgtyp,isltyp,vegfrac,rmol,gsw,znt,hfx,pbl,               &
        relhum,             &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)
!    print*,'hli test2 ktau',ktau

    if (call_gocart) then
      call gocart_chem_driver(ktau,dt,dtstep,gmt,julday,xcosz,          &
           t_phy,moist,chem,rho_phy,dz8w,p8w,backg_oh,oh_t,             &
           backg_h2o2,h2o2_t,backg_no3,no3_t,                           &
           dxy,g,xlat,xlong,ttday,tcosz,                             &
           chem_opt,num_chem,num_moist,                                 &
           ids,ide, jds,jde, kds,kde,                                   &
           ims,ime, jms,jme, kms,kme,                                   &
           its,ite, jts,jte, kts,kte                        )
      call gocart_aerosols_driver(ktau,dtstep,t_phy,moist,              &
           chem,rho_phy,dz8w,p8w,dxy,g,                              &
           chem_opt,num_chem,num_moist,                                 &
           ids,ide, jds,jde, kds,kde,                                   &
           ims,ime, jms,jme, kms,kme,                                   &
           its,ite, jts,jte, kts,kte                        )
    endif

    call sum_pm_gocart (                                                &
         rri, chem,pm2_5_dry, pm2_5_dry_ec, pm10,                       &
         num_chem,chem_opt,                                             &
         ids,ide, jds,jde, kds,kde,                                     &
         ims,ime, jms,jme, kms,kme,                                     &
         its,ite, jts,jte, kts,kte)

    ! -- pm25 and pm10 for output , not for tracer options
    do j = jts, jte
      do k = kts, kte
        do i = its, ite
          pm25  (i,j,k) = pm2_5_dry(i,k,j)
          p10   (i,j,k) = pm10     (i,k,j)
          !ebu_oc(i,j,k) = ebu      (i,k,j,p_ebu_oc)
        end do
      end do
    end do

    if (call_gocart) then
      do j = jts, jte
        do k = kts, kte
          do i = its, ite
            oh_bg  (i,j,k) = max(0., oh_t  (i,k,j))
            h2o2_bg(i,j,k) = max(0., h2o2_t(i,k,j))
            no3_bg (i,j,k) = max(0., no3_t (i,k,j))
          end do
        end do
      end do
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

!    print*,'hli test3 ktau',ktau
!   call gsd_chem_post() ! postprocessing for diagnostics

!
   end subroutine gsd_chem_gocart_wrapper_run
!> @}

  subroutine gsd_chem_gocart_prep(                                        &
        readrestart,chem_in_opt,ktau,dtstep,xcosz,                     &
        u10m,v10m,ustar,land,garea,rlat,rlon,ts2d,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                &
        exch,dqdt,                                                     &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                  &
        snow_cpl,emi2_in,                               &
        hf2d,pb2d,                              &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                          &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                  &
        t8w,exch_h,dqdti,                                              &
        z_at_w,vvel,zmid,                                              &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                               &
        ntbc1,ntbc2,ntoc1,ntoc2,                                       &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                 &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                &
        ntrac,gq0,                                                     &
        num_chem, num_moist,                                &
        call_gocart,nvl_gocart,                                        &
        ttday,tcosz,gmt,julday,                                        &
        backg_oh,backg_h2o2,backg_no3,                                 &
        ppm2ugkg,                                             &
        moist,chem,                                   &
        smois,ivgtyp,isltyp,vegfrac,rmol,gsw,znt,hfx,pbl,              &
        relhum,            &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    logical, intent(in) :: readrestart
    integer, intent(in) :: chem_in_opt, ktau, julday
    real(kind=kind_phys), intent(in) :: dtstep, gmt

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
    real(kind=kind_phys), dimension(ims:ime, 64, 3),   intent(in) :: emi2_in
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) ::     &
         pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,exch,dqdt
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_moist,                 &
                           nvl_gocart
    logical,intent(in) ::  call_gocart
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme),    intent(out) ::          &
                           backg_oh,backg_h2o2,backg_no3

    
    integer,dimension(ims:ime, jms:jme), intent(out) :: isltyp, ivgtyp
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, u_phy, v_phy, p_phy, rho_phy, dz8w, p8w, t8w, vvel, zmid,         &
         exch_h,dqdti
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         u10, v10, ust, tsk, xland, xlat, xlong, dxy, vegfrac, rmol, gsw, znt, hfx,    &
         pbl, ttday, tcosz
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w, relhum
    real(kind_phys), dimension(ims:ime, nsoil, jms:jme), intent(out) :: smois
    real(kind_phys), parameter :: frac_so2_ant = 0.5_kind_phys  ! antropogenic so2 fraction
    real(kind_phys), parameter :: frp2plume = 1.e+06_kind_phys  ! FRP-to-plume conversion factor
    real(kind_phys), parameter :: frpc  = 1.e+09_kind_phys      ! FRP conversion factor
    real(kind_phys), dimension(nvl_gocart) :: p_gocart

    ! -- local variables
!   real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(kind_phys), dimension(ims:ime, jms:jme, nvl_gocart) :: oh_backgd,h2o2_backgd,no3_backgd
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
    backg_oh       = 0._kind_phys
    backg_h2o2     = 0._kind_phys 
    backg_no3      = 0._kind_phys
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
    ttday          = 0._kind_phys
    tcosz          = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys
    relhum         = 0._kind_phys


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
     ivgtyp (i,1)=vegtype(i)
     isltyp (i,1)=soiltyp(i)
     vegfrac(i,1)=sigmaf (i)
    enddo
   
    rmol=0.

    do k=1,nsoil
     do j=jts,jte
      do i=its,ite
       smois(i,k,j)=smc(i,k)
      enddo
     enddo
    enddo

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

    if (.NOT. readrestart) then
      if (chem_in_opt == 0 ) then
        if(ktau.le.1)then
!           if(chem_opt > 0 ) then
          do j=jts,jte
            jp = j - jts + 1
            do k=kts,kte
              do i=its,ite
                ip = i - its + 1
                if (chem_opt == CHEM_OPT_GOCART) then
                  do n=1,num_chem
                    chem(i,k,j,n)=1.e-12
                  enddo
                endif  ! chem_opt==300
                chem(i,k,j,p_so2)=5.e-6
                chem(i,k,j,p_sulf)=3.e-6
                if ((chem_opt >= CHEM_OPT_GOCART) .and. (chem_opt < CHEM_OPT_MAX)) then
                  chem(i,k,j,p_msa)=0.1e-6
                  chem(i,k,j,p_dms)=0.1e-6
                  chem(i,k,j,p_bc1)=0.1e-3
                  chem(i,k,j,p_bc2)=0.1e-3
                  chem(i,k,j,p_oc1)=0.1e-3
                  chem(i,k,j,p_oc2)=0.1e-3
                  chem(i,k,j,p_p25)=0.1e-3 !lzhang
                  chem(i,k,j,p_p10)=0.1e-3 !lzhang
                endif !chem_opt >= 300 .and. chem_opt <  500

!                if ((chem_opt == CHEM_OPT_GOCART_RACM) .or. (chem_opt == CHEM_OPT_RACM_SOA_VBS)) then  !added o3 background !lzhang
!                  kk=min(k,kte)
!                  kkp = kk - kts + 1
!                  ! -- add initial constant into O3,CH4 and CO ect.
!                  chem(i,k,j,p_o3)=epsilc
!                  ! -- this section needs to be revisited before enabling the
!                  ! corresponding chem_opt options
!                  ! maxth=min(400.,th_pvsrf(i,j))
!                  ! if (tr3d(ip,jp,kkp,p_atm_ptem) > maxth) then
!                  !   chem(i,k,j,p_o3)=(airmw/48.)*tr3d(ip,jp,kkp,p_atm_o3mr)*1e6
!                  !   !convert kg/kg to ppm
!                  ! else
!                  !   chem(i,k,j,p_o3)=0.03 !ppm
!                  ! endif
!                  chem(i,k,j,p_ch4)=1.85 !ppm
!                  chem(i,k,j,p_co)=0.06 !ppm
!                  chem(i,k,j,p_co2)=380.
!                  chem(i,k,j,p_ete)=epsilc
!                  chem(i,k,j,p_udd)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_hket)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_api)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_lim)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_dien)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_macr)=chem(i,k,j,p_ete)
!                endif !( (chem_opt == 301.or.chem_opt==108))
              enddo
            enddo
          enddo
        endif !(ktau<=1)

      else !(chem_in_opt == 0 )

        if ((ktau<=1).and.((chem_opt == CHEM_OPT_GOCART_RACM).or.(chem_opt == CHEM_OPT_RACM_SOA_VBS))) then  !added GFS o3 background above 380K!lzhang
          do j=jts,jte
            jp = j - jts + 1
            do k=kts,kte+1
              kk=min(k,kte)
              kkp = kk - kts + 1
              do i=its,ite
                ip = i - its + 1
                ! -- this section needs to be revisited before enabling the
                ! corresponding chem_opt options
                ! maxth=min(400.,th_pvsrf(i,j))
                ! if (tr3d(ip,jp,kkp,p_atm_ptem) >= maxth) then
                !   chem(i,k,j,p_o3)=(airmw/48.)*tr3d(ip,jp,kkp,p_atm_o3mr)*1e6 !convert kg/kg to ppm
                ! endif !380K
              enddo
            enddo
          enddo
        endif ! chem_opt == 301.or.chem_opt==108

      endif !(chem_in_opt == 1 )
     endif ! readrestart

     !-- assgin read in 3D background chemical species
     do i=its,ite
       do k=1,nvl_gocart
          h2o2_backgd(i,1,k)=emi2_in(i,k,1)
          no3_backgd (i,1,k)=emi2_in(i,k,2)
          oh_backgd  (i,1,k)=emi2_in(i,k,3)
       enddo
     enddo

    !
    ! -- gocart background fields only if gocart is called
    !
    !if (.NOT. readrestart) then
    if (call_gocart .and. (chem_opt == CHEM_OPT_GOCART))then
      do j=jts,jte
        do i=its,ite
          do k=kts,kte
            do ll=2,nvl_gocart
              l=ll
              if (p_gocart(l) < .01*p_phy(i,k,j)) exit
            enddo
            pu=alog(p_gocart(l))
            pl=alog(p_gocart(l-1))
            pwant=alog(.01*p_phy(i,k,j))
            if (pwant > pl)then
              backg_oh(i,k,j)=oh_backgd(i,j,l)
              backg_h2o2(i,k,j)=h2o2_backgd(i,j,l)
              backg_no3(i,k,j)=no3_backgd(i,j,l)
            else
              aln=(oh_backgd(i,j,l)*(pwant-pl)+            &
                oh_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
              backg_oh(i,k,j)=aln
              aln=(h2o2_backgd(i,j,l)*(pwant-pl)+            &
                h2o2_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
              backg_h2o2(i,k,j)=aln
              aln=(no3_backgd(i,j,l)*(pwant-pl)+            &
                no3_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
              backg_no3(i,k,j)=aln
            endif
          enddo
        enddo
      enddo
    endif   ! end gocart stuff
    !endif !restart


    if ((chem_opt == CHEM_OPT_RACM_SOA_VBS) .or. (chem_opt >= CHEM_OPT_GOCART .and. chem_opt < CHEM_OPT_MAX)) then
      !ndystep=86400/ifix(dtstepc)
      ndystep=86400/ifix(dtstep)
      do j=jts,jte
        do i=its,ite
          tcosz(i,j)=0.
          ttday(i,j)=0.
!         rlat=xlat(i,j)*3.1415926535590/180.
          xlonn=xlong(i,j)
          do n=1,ndystep
            xtime=n*dtstep/60.
            ixhour=ifix(gmt+.01)+ifix(xtime/60.)
            xhour=float(ixhour)
            xmin=60.*gmt+(xtime-xhour*60.)
            gmtp=mod(xhour,24.)
            gmtp=gmtp+xmin/60.
            CALL szangle(1, 1, julday, gmtp, sza, cosszax,xlonn,rlat(i))
            TCOSZ(i,j)=TCOSZ(I,J)+cosszax(1,1)
            if (cosszax(1,1) > 0.) ttday(i,j)=ttday(i,j)+dtstep
            !--use physics inst cosine zenith -- hli 03/06/2020
!            TCOSZ(i,j)=TCOSZ(I,J)+xcosz(i)
!            if (xcosz(i) > 0.) ttday(i,j)=ttday(i,j)+dtstep
          enddo
        enddo
      enddo
    endif !chem_opt >= 300 .and. chem_opt <  500


  end subroutine gsd_chem_gocart_prep

!> @}
  end module gsd_chem_gocart_wrapper
