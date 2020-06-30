!>\file gsd_chem_seas_wrapper.F90
!! This file is GSD Chemistry sea salt wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 05/2020

 module gsd_chem_seas_wrapper

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
   use seas_mod,        only : gocart_seasalt_driver
!   use plume_data_mod
!   use dust_gocart_mod, only : gocart_dust_driver
!   use dust_afwa_mod,   only : gocart_dust_afwa_driver
!   use dust_fengsha_mod,only : gocart_dust_fengsha_driver
!   use dust_data_mod
!   use plume_rise_mod
!   use vash_settling_mod

   implicit none

   private

   public :: gsd_chem_seas_wrapper_init, gsd_chem_seas_wrapper_run, gsd_chem_seas_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine gsd_chem_seas_wrapper_init()
      end subroutine gsd_chem_seas_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_gsd_chem_seas_wrapper_finalize Argument Table
!!
      subroutine gsd_chem_seas_wrapper_finalize()
      end subroutine gsd_chem_seas_wrapper_finalize

!> \defgroup gsd_chem_seas_group GSD Chem seas wrapper Module
!! This is the gsd chemistry
!>\defgroup gsd_chem_seas_wrapper GSD Chem seas wrapper Module  
!> \ingroup gsd_chem_seas_group
!! This is the GSD Chem seas wrapper Module
!! \section arg_table_gsd_chem_seas_wrapper_run Argument Table
!! \htmlinclude gsd_chem_seas_wrapper_run.html
!!
!>\section gsd_chem_seas_wrapper GSD Chemistry Scheme General Algorithm
!> @{
    subroutine gsd_chem_seas_wrapper_run(im, kte, kme, ktau, dt, garea, land, &
                   u10m, v10m, ustar, rlat, rlon, tskin,julian,xcosz,   &
                   rain_cpl, rainc_cpl, hf2d, pb2d,                     &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum,  &
                   w, exch, dqdt,                                       &
                   nsoil, smc, vegtype, soiltyp, sigmaf,jdate,idat,     & 
                   dswsfc, zorl,snow_cpl,                               &
                   nseasalt,ntrac,                                      &
                   ntss1,ntss2,ntss3,ntss4,ntss5,                       &
                   gq0,tile_num,                                        &
                   errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau,nsoil,jdate(8),idat(8),tile_num
    integer,        intent(in) :: nseasalt,ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    real(kind_phys),intent(in) :: dt,julian

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    integer, dimension(im), intent(in) :: land, vegtype, soiltyp        
    real(kind_phys), dimension(im,nsoil), intent(in) :: smc
    real(kind_phys), dimension(im), intent(in) :: u10m, v10m, ustar,              &
                garea, rlat,rlon, tskin, rain_cpl, rainc_cpl,                     &
                hf2d, pb2d, sigmaf, dswsfc, zorl, snow_cpl, xcosz
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d,        &
                us3d, vs3d, spechum, w, exch, dqdt
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0
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


!>- plume variables
    ! -- buffers
    real(kind_phys) :: dtstep, gmt
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg


!>-- local variables
    real(kind_phys) :: curr_secs
    logical :: call_radiation
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
    call gsd_chem_prep_seas(                                            &
        readrestart,chem_in_opt,ktau,dtstep,xcosz,                      &
        u10m,v10m,ustar,land,garea,rlat,rlon,tskin,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                 &
        exch,dqdt,                                                      &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                   &
        snow_cpl,                                                       &
        hf2d,pb2d,                               &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        t8w,exch_h,dqdti,                                               &
        z_at_w,vvel,zmid,                                               &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                  &
        ntrac,gq0,                                                      &
        num_chem, num_moist,                                 &
        ppm2ugkg,                                              &
        moist,chem,                                                     &
        smois,ivgtyp,isltyp,vegfrac,rmol,gsw,znt,hfx,pbl,               &
        snowh,clayf,rdrag,sandf,ssm,uthr,dms_0,erod,             &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)


    ! -- compute sea salt
    if (seas_opt >= SEAS_OPT_DEFAULT) then
    call gocart_seasalt_driver(ktau,dt,rri,t_phy,moist,                 &
        u_phy,v_phy,chem,rho_phy,dz8w,u10,v10,ust,p8w,tsk,              &
        xland,xlat,xlong,dxy,g,emis_seas,                               &
        seashelp,num_emis_seas,num_moist,num_chem,seas_opt,             &
        ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde,                                      &
        ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme,                                      &
        its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte)
    endif



    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntss1  )=ppm2ugkg(p_seas_1) * max(epsilc,chem(i,k,1,p_seas_1))
       gq0(i,k,ntss2  )=ppm2ugkg(p_seas_2) * max(epsilc,chem(i,k,1,p_seas_2))
       gq0(i,k,ntss3  )=ppm2ugkg(p_seas_3) * max(epsilc,chem(i,k,1,p_seas_3))
       gq0(i,k,ntss4  )=ppm2ugkg(p_seas_4) * max(epsilc,chem(i,k,1,p_seas_4))
       gq0(i,k,ntss5  )=ppm2ugkg(p_seas_5) * max(epsilc,chem(i,k,1,p_seas_5))
     enddo
    enddo

!   call gsd_chem_post() ! postprocessing for diagnostics

!
   end subroutine gsd_chem_seas_wrapper_run
!> @}

  subroutine gsd_chem_prep_seas(                                       &
        readrestart,chem_in_opt,ktau,dtstep,xcosz,                     &
        u10m,v10m,ustar,land,garea,rlat,rlon,ts2d,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                &
        exch,dqdt,                                                     &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                  &
        snow_cpl,                                                      &
        hf2d,pb2d,                              &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                          &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                  &
        t8w,exch_h,dqdti,                                              &
        z_at_w,vvel,zmid,                                              &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                 &
        ntrac,gq0,                                                     &
        num_chem, num_moist,                                &
        ppm2ugkg,                                             &
        moist,chem,                                                     &
        smois,ivgtyp,isltyp,vegfrac,rmol,gsw,znt,hfx,pbl,              &
        snowh,clayf,rdrag,sandf,ssm,uthr,dms_0,erod,            &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    logical, intent(in) :: readrestart
    integer, intent(in) :: chem_in_opt, ktau!, julday
    real(kind=kind_phys), intent(in) :: dtstep

    !FV3 input variables
    integer, intent(in) :: nsoil
    integer, dimension(ims:ime), intent(in) :: land, vegtype, soiltyp
    integer, intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         u10m, v10m, ustar, garea, rlat, rlon, ts2d, sigmaf, dswsfc,       &
         zorl, snow_cpl, hf2d, pb2d,xcosz
    real(kind=kind_phys), dimension(ims:ime, nsoil),   intent(in) :: smc 
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) ::     &
         pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,exch,dqdt
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_moist
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                       &
                           ims,ime, jms,jme, kms,kme,                       &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg

    integer,dimension(ims:ime, jms:jme), intent(out) :: isltyp, ivgtyp
    real(kind_phys), dimension(ims:ime, jms:jme, 3), intent(inout) :: erod
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, u_phy, v_phy, p_phy, rho_phy, dz8w, p8w, t8w, vvel, zmid,         &
         exch_h,dqdti
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         u10, v10, ust, tsk, xland, xlat, xlong, dxy, vegfrac, rmol, gsw, znt, hfx,    &
         pbl, snowh, clayf, rdrag, sandf, ssm, uthr, dms_0
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w
    real(kind_phys), dimension(ims:ime, nsoil, jms:jme), intent(out) :: smois

    ! -- local variables
!   real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(kind_phys) ::  factor,factor2,pu,pl,aln,pwant
    real(kind_phys) ::  xhour,xmin,gmtp,xlonn,xtime,real_time
    integer i,ip,j,jp,k,kp,kk,kkp,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour

    ! -- initialize output arrays
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
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys


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
       chem(i,k,jts,p_seas_1)=max(epsilc,gq0(i,k,ntss1  )/ppm2ugkg(p_seas_1))
       chem(i,k,jts,p_seas_2)=max(epsilc,gq0(i,k,ntss2  )/ppm2ugkg(p_seas_2))
       chem(i,k,jts,p_seas_3)=max(epsilc,gq0(i,k,ntss3  )/ppm2ugkg(p_seas_3))
       chem(i,k,jts,p_seas_4)=max(epsilc,gq0(i,k,ntss4  )/ppm2ugkg(p_seas_4))
       chem(i,k,jts,p_seas_5)=max(epsilc,gq0(i,k,ntss5  )/ppm2ugkg(p_seas_5))
     enddo
    enddo


  end subroutine gsd_chem_prep_seas
!> @}
  end module gsd_chem_seas_wrapper
