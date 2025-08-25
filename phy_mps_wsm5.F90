!
#undef VD_RAIN
#undef ICE_PARAM
#undef EFFIC
!-------------------------------------------------------------------------------
   module module_mp_wsm5
!-------------------------------------------------------------------------------!
!
!  constant   :
!
!    dtcldcr       - maximum time step for minor loops
!    n0r           - intercept parameter rain
!    n0s           - intercept parameter snow
!    n0smax        - maximum n0s (t=-90C unlimited)
!    alpha         - 0.122 exponen factor for n0s
!    avtr, bvtr    - a constant for terminal velocity of rain
!    avts, bvts    - a constant for terminal velocity of snow
!    lamdarmax     - limited maximum value for slope parameter of rain
!    lamdasmax     - limited maximum value for slope parameter of snow
!    r0            - 8 microm  in contrast to 10 micro m
!    peaut         - collection efficiency
!    xncr          - maritime cloud in contrast to 3.e8 in tc80
!    xmyu          - the dynamic viscosity kg m-1 s-1 
!    dicon         - constant for the cloud-ice diamter
!    dimax         - limited maximum value for the cloud-ice diamter
!    pfrz1, pfrz2  - constant in Biggs freezing
!    qcrmin        - minimum values for qr and qs
!    eacrc         - now/cloud-water collection efficiency
!
!-------------------------------------------------------------------------------
   use module_mp_radar
   use control_phy, only : l_refl
!
   real, parameter, private :: dtcldcr = 120.
   real, parameter, private :: n0r     = 8.e6
   real, parameter, private :: n0s     = 2.e6
   real, parameter, private :: n0g     = 4.e6
   real, parameter, private :: n0smax  =  1.e11
   real, parameter, private :: alpha   = .12
#ifdef VD_RAIN
   real, parameter, private :: avtr    = 5881
   real, parameter, private :: bvtr    = 1.03
   real, parameter, private :: fr      = 202.4
#else
   real, parameter, private :: avtr    = 841.9
   real, parameter, private :: bvtr    = 0.8
#endif   
   real, parameter, private :: avts    = 11.72
   real, parameter, private :: bvts    = .41
   real, parameter, private :: avti    = 14900.
   real, parameter, private :: bvti    = 1.31
#ifdef ICE_PARAM
   real, parameter, private :: dmi     = 2.09
   real,            private :: rdmi    = 1./dmi
   real, parameter, private :: evq     = 3.29
   real, parameter, private :: fvq     = 0.16
   real, parameter, private :: alpha1  = 1.0
#endif
   real, parameter, private :: lamdarmax = 8.e4
   real, parameter, private :: lamdasmax = 1.e5
   real, parameter, private :: r0      = .8e-5
   real, parameter, private :: peaut   = .55
   real, parameter, private :: xncr    = 3.e8
   real, parameter, private :: xmyu    = 1.718e-5
   real, parameter, private :: dicon   = 11.9
   real, parameter, private :: dimax   = 500.e-6
   real, parameter, private :: pfrz1   = 100.
   real, parameter, private :: pfrz2   = 0.66
   real, parameter, private :: qcrmin  = 1.e-9
   real, parameter, private :: eacrc   = 1.0
   real, parameter, private :: decfl   = 8.0
   real, save ::                                                               &
                     qc0, qck1,pidnc,bvtr1,bvtr2,bvtr3,bvtr4,g1pbr,            &
                     g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr,                     &
                     precr1,precr2,xmmax,roqimax,bvts1,                        &
                     bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs,                      &
                     g5pbso2,pvts,pacrs,precs1,precs2,pidn0r,                  &
                     pidn0s,xlv1,pacrc,                                        &
#ifdef ICE_PARAM
                     ani,bni,cmi,                                              &
#endif
                     rslopermax,rslopesmax,rslopegmax,                         &
                     rsloperbmax,rslopesbmax,rslopegbmax,                      &
                     rsloper2max,rslopes2max,rslopeg2max,                      &
                     rsloper3max,rslopes3max,rslopeg3max
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine phy_mps_wsm5(t1, q1,                                             &
                            q2, ncloud, p1, delz1,                             &
                            deltim, g, cpd, cpv, rd, rv, t0c, ttp, pi,         &
                            tv,                                                &
                            ep1, ep2, qmin,                                    &
                            XLS, XLV0, XLF0, den0, denr,                       &
                            cliq,cice,psat,                                    &
                            dxmeter1,                                          &
                            rain,rainncv,                                      &
                            sr,                                                &
                            qrscps,qrsscv,taucps,tauscv,                       &
                            refl_10cm,                                         &
                            ids, ide, jds, jde, kds, kde,                      &
                            ims, ime, jms, jme, kms, kme,                      &
                            its, ite, jts, jte, kts, kte,                      &
                            snow, snowncv                                      &
                            )
!-------------------------------------------------------------------------------
!
!  abstract : 
!    - this code is a 5-class mixed ice microphyiscs scheme (wsm5) of the
!    wrf single-moment microphyiscs (wsmmp). the wsmmp assumes that ice nuclei
!    number concentration is a function of temperature, and seperate assumption
!    is developed, in which ice crystal number concentration is a function
!    of ice amount. a theoretical background of the ice-microphysics and related
!    processes in the wsmmps are described in hong et al. (2004).
!    production terms in the wsm6 scheme are described in hong and lim (2006).
!    all units are in m.k.s. and source/sink terms in kgkg-1s-1.
!    - wsm5 cloud scheme
!    - coded by song-you hong (yonsei univ.)
!               jimy dudhia (ncar) and shu-hua chen (uc davis)
!               summer 2002
!    - implemented by song-you hong (yonsei univ.) and jimy dudhia (ncar)
!                     summer 2003
!    - further modifications :
!      semi-lagrangian sedimentation (jh,2010),hong, aug 2009
!      ==> higher accuracy and efficient at lower resolutions
!      reflectivity computation from greg thompson, lim, jun 2011
!      ==> only dignostic, but with removal of too large drops
!      effetive radius of hydrometeors, bae from kiaps, jan 2015
!      ==> consistency in solar insolation of rrtmg radiation
!
!  history log :
!    2015-07-01    soo ya bae     define top layer for computing processes
!                                 add subroutine effective radius
!    2015-10-01    soo ya bae     bug fix in melting term with a density factor 
!                                 move subroutine effecive radius to radiation
!    2016-10-01    soo ya bae     fvps function for specific humidity
!    2022-07-01    soo ya bae     1-d
!
!  references :
!    hong, dudhia and chen (hdc, 2004, mwr)
!    rutledge and hobbs (rh83, 1983, jas)
!    hong and lim (hl, 2006, j. korean meteor. soc.)
!
!  structure :
!    phy_main_solver --- phy_mps_wsm5 --- slope_wsm5
!                                      |
!                                      |- nislfv_precp
!                                      |
!                                      |- refl10cm_wsm5
! 
!  input :
!    deltim                       - timestep
!    g, cpd, cpv, t0c, den0,      - constant 
!    rd, rv, ep1, ep2, qmin, 
!    xls, xlv0, xlf0, denr,
!    cliq, cice, psat
!    ids, ide, jds, jde, kds, kde - dimension
!    ims, ime, jms, jme, kms, kme 
!    its, ite, jts, jte, kts, kte 
!    ncloud                       - number of hydrometeor
!    p                            - pressure, Pa
!    delz                         - depth of model layer, m
!    
!  inout : 
!    t1                           - temperautre
!    q1                           - specific humidity
!    q2                           - mixing ratio of cloud, rain, ice, and snow
!                                   qc, qr, qi, qc
!    rain, rainncv                - precipitation
!    sr                           - ratio of snow to rain
!
!  all units are in m.k.s. and source/sink terms in kg kg-1 s-1.
!-------------------------------------------------------------------------------
   use comfcst,           only: kmaxin
   use phy_funct,         only: mpp1
!
   implicit none
!
   integer                                , intent(in   ) :: ncloud
   integer                                , intent(in   ) ::                  &
                                                   ids,ide, jds,jde, kds,kde, &
                                                   ims,ime, jms,jme, kms,kme, &
                                                   its,ite, jts,jte, kts,kte
   real                                   , intent(in   ) :: deltim
   real                                   , intent(in   ) :: g, cpd, cpv, t0c, &
                                                             den0, rd, rv,     &
                                                             ep1, ep2, qmin,   &
                                                             xls, xlv0, xlf0,  &
                                                             cliq, cice, psat, &
                                                             denr,ttp,pi
   real, dimension(ims:ime)               , intent(in   ) :: dxmeter1
   real, dimension(ims:ime)               , intent(in   ) :: taucps
   real, dimension(ims:ime)               , intent(in   ) :: tauscv
   real, dimension(ims:ime,kms:kme)       , intent(in   ) :: qrscps
   real, dimension(ims:ime,kms:kme)       , intent(in   ) :: qrsscv
   real, dimension(ims:ime,kms:kme)       , intent(in   ) :: p1
   real, dimension(ims:ime,kms:kme)       , intent(in   ) :: delz1
   real, dimension(ims:ime,kms:kme)       , intent(in   ) :: tv
   real, dimension(ims:ime,kms:kme)       , intent(inout) :: t1
   real, dimension(ims:ime,kms:kme)       , intent(inout) :: q1
   real, dimension(ims:ime,kms:kme,ncloud), intent(inout) :: q2
   real, dimension(ims:ime,kms:kme)       , intent(inout) :: refl_10cm
   real, dimension(ims:ime)               , intent(inout) :: rain
   real, dimension(ims:ime)               , intent(inout) :: rainncv
   real, dimension(ims:ime)               , intent(inout) :: sr
   real, dimension(ims:ime), optional     , intent(inout) :: snow
   real, dimension(ims:ime), optional     , intent(inout) :: snowncv
!
! local variables
!
   real, dimension(kts:kte,its:ite)   :: den
   real, dimension(kts:kte,its:ite)   :: dend
   real, dimension(kts:kte,its:ite)   :: delz
   real, dimension(kts:kte,its:ite)   :: q
   real, dimension(kts:kte,its:ite)   :: t
   real, dimension(kts:kte,its:ite)   :: p
   real, dimension(kts:kte,its:ite,2) :: qrs, qci
   real, dimension(kts:kte,its:ite)   :: denfac
   real, dimension(kts:kte,its:ite)   :: xl
   real, dimension(kts:kte,its:ite)   :: cpm
   real, dimension(kts:kte,2)         :: rslope, rslope2, rslope3, rslopeb
   real, dimension(kts:kte,2)         :: vt
#ifdef VD_RAIN
   real, dimension(kts:kte)           :: rslopef, rslopef2, rslopef3, rslopefb
   real, dimension(kts:kte)           :: rslopehf, rslopehf2, rslopehf3
   real, dimension(kts:kte)           :: rslopehfb
#endif
   real, dimension(kts:kte)           :: rh_wat
   real, dimension(kts:kte)           :: rh_ice
   real, dimension(kts:kte)           :: qsat_wat
   real, dimension(kts:kte)           :: qsat_ice
   real, dimension(kts:kte)           :: ab_wat
   real, dimension(kts:kte)           :: ab_ice
   real, dimension(kts:kte)           :: venfac
   real, dimension(kts:kte)           :: xlf
   real, dimension(kts:kte)           :: denqr
   real, dimension(kts:kte)           :: denqs
   real, dimension(kts:kte)           :: denqi
   real, dimension(kts:kte)           :: vtr
   real, dimension(kts:kte)           :: vts
   real, dimension(kts:kte)           :: vti
   real, dimension(kts:kte)           :: vt2i
   real, dimension(kts:kte)           :: vt2s
   real, dimension(kts:kte)           :: factor1
   real, dimension(kts:kte)           :: cldf
   real, dimension(kts:kte)           :: xal, xbl
   real, dimension(kts:kte)           :: tr, logtr
   real, dimension(kts:kte)           :: ni, mi, diameter
   real, dimension(kts:kte)           :: supcol, supcolt, supsat
   real, dimension(kts:kte)           :: satdt, supice
   real, dimension(kts:kte)           :: total, value1, xlvalue
   real, dimension(kts:kte)           :: n0sfac
   real, dimension(kts:kte)           :: ni0, roqi0, qimax
   real, dimension(kts:kte)           :: coeres, eacrs, acrfac
   real, dimension(kts:kte)           :: pigen, pidep, psdep, praut, psaut,    &
                                         prevp, psevp, pracw, psacw, psaci,    &
                                         pcond, psmlt, pihtf, psfrz
   integer, dimension(kts:kte)        :: ifsat
   real, dimension(its:ite,kts:kte)   :: refl
   real, dimension(its:ite,kts:kte)   :: qrs1
   real, dimension(its:ite,kts:kte)   :: qrs2
   real, dimension(its:ite,kts:kte)   :: qrs3
   real, dimension(its:ite,kts:kte)   :: qrsub
   real, dimension(its:ite,kts:kte)   :: qssub
   real, dimension(its:ite)           :: tstepsnow
   real, dimension(its:ite)           :: dxmeter
   real                               :: rdtcld
   real                               :: delt
   real                               :: zfac
!
! precipitation
!
   real                               :: precp_sum, precp_ice
   real                               :: precp_r, precp_s, precp_i
   real                               :: qrpath, qspath, qipath
   real                               :: dtcfl
   integer                            :: niter
   integer                            :: nstep
!
! find cloud top
!
   integer                            :: ktopqc, ktopqr, ktopqi, ktopqs
   integer                            :: ktoprh
   logical                            :: lqc, lqr, lqi, lqs
   logical                            :: flgcld
!
! top level for computing of cloud microphysical process
!
   integer                            :: ktop
   real                               :: cpmcal, xlcal, diffus,                &
                                         viscos, xka, conden, diffac,          &
                                         x, y, z, a, b, c, d, e,               &
                                         qdt, holdrr, holdrs, pvt, dtcld,      &
                                         factor, source, value, holdc, holdci
!
! variables for optimization
!
   real, dimension(its:ite)           :: tvec1
   real                               :: temp
   integer                            :: i, j, k, mstepmax, iprt, latd, lond,  &
                                         loop, loops, n, idim, kdim
!
! temporaries used for inlining fsvp function
!
   real                               :: dldti, xb, xai, xbi, xa,              &
                                         hvap, cvap, hsub, dldt
!------------------------------------------------------------------------------
!
! compute internal functions
!
   cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
   xlcal(x)  = xlv0-xlv1*(x-t0c)
!
   idim = ite-its+1
   kdim = kte-kts+1
!
! define top layer for computing of cloud microphysical process
! sigma level>0.05 -> kmaxin = k+1 calculated 
!
   ktop = kmaxin
!
! padding 0 for negative values generated by dynamics
!
   do k = kts,ktop
     do i = its,ite
       qci(k,i,1) = max(q2(i,k,1),0.0)
       qrs(k,i,1) = max(q2(i,k,2),0.0)
       qci(k,i,2) = max(q2(i,k,3),0.0)
       qrs(k,i,2) = max(q2(i,k,4),0.0)
     enddo
   enddo
!
   do i = its,ite
     do k = ktop+1,kte
       qci(k,i,1) = 0.0
       qrs(k,i,1) = 0.0
       qci(k,i,2) = 0.0
       qrs(k,i,2) = 0.0
     enddo
   enddo
!
! initialize the surface rain and snow
!
   do i = its,ite
     rain(i) = 0.
     rainncv(i) = 0.
     snow(i) = 0.
     if (present(snowncv).and.present(snow)) snowncv(i) = 0.
     sr(i) = 0.
     tstepsnow(i) = 0.
     dxmeter(i) = dxmeter1(i)
   enddo
!
   delt = deltim
!
   do k =  kts,kte
     do i = its,ite
       t(k,i) = t1(i,k)
       q(k,i) = q1(i,k)
       p(k,i) = p1(i,k)
       delz(k,i) = delz1(i,k)
       den(k,i)  = p(k,i)/(rd*tv(i,k))
       dend(k,i) = den(k,i)*(1.+q(k,i)*rv/rd)/(1.+q(k,i))
       denfac(k,i) = sqrt(den0/den(k,i))
       zfac = mpp1(t(k,i))
       qrsub(i,k) = (1.-zfac)*(qrscps(i,k)*taucps(i)+qrsscv(i,k)*tauscv(i))    &
                   /delt
       qssub(i,k) = zfac*(qrscps(i,k)*taucps(i)+qrsscv(i,k)*tauscv(i))/delt
     enddo
   enddo
!
! latent heat for phase changes and heat capacity. neglect the
! changes during microphysical process calculation
! emanuel(1994)
!
   do i = its,ite
     do k = kts,kte
       cpm(k,i) = cpmcal(q(k,i))
       xl(k,i) = xlcal(t(k,i))
     enddo
   enddo
!
   refl(:,:) = -35.
   qrs1(:,:) = 0.    ;  qrs2(:,:) = 0.    ;  qrs3(:,:) = 0.
!
   hsub = xls
   hvap = xlv0
   cvap = cpv
   dldt = cvap-cliq
   xa   = -dldt/rv
   xb   = xa+hvap/(rv*ttp)
   dldti= cvap-cice
   xai  = -dldti/rv
   xbi  = xai+hsub/(rv*ttp)
!
! compute the minor time steps.
!
   loops = max(nint(delt/dtcldcr),1)
   dtcld = delt/loops
   if (delt.le.dtcldcr) dtcld = delt       
!===============================================================================
!
! inner loop with the time step of dtcldcr (default = 120 sec)
!
!-------------------------------------------------------------------------------
   inner_loop : do loop = 1,loops
!
! i-loop for one-dimension code
!
     i_loop : do i = its,ite
!-------------------------------------------------------------------------------
!
! initialization for fsvp
!
       xal(:) = 0.       ; xbl(:) = 0.      ; tr(:) = 0.  ; logtr(:) = 0.
       qsat_wat(:) = 0.  ; qsat_ice(:) = 0.
       rh_wat(:) = 0.    ; rh_ice(:) = 0.
!
! Inline expansion for fsvp
!    qsat_wat(k) = fsvp(t(k,i),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!    qsat_ice(k) = fsvp(t(k,i),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!
       ktop = kmaxin
!
       do k = kts,ktop
         if (t(k,i).lt.ttp) then
           xal(k) = xai
           xbl(k) = xbi
         else
           xal(k) = xa
           xbl(k) = xb
         endif
       enddo
!
       do k = kts,ktop
         tr(k) = ttp/t(k,i)
         logtr(k) = log(tr(k))
         qsat_wat(k) = psat*exp(logtr(k)*(xa)+xb*(1.-tr(k)))
         qsat_wat(k) = min(qsat_wat(k),0.99*p(k,i))
         qsat_wat(k) = ep2*qsat_wat(k)/(p(k,i)-qsat_wat(k))
         qsat_wat(k) = max(qsat_wat(k),qmin)
         rh_wat(k)   = max(q(k,i)/qsat_wat(k),qmin)
         qsat_ice(k) = psat*exp(logtr(k)*(xal(k))+xbl(k)*(1.-tr(k)))
         qsat_ice(k) = min(qsat_ice(k),0.99*p(k,i))
         qsat_ice(k) = ep2*qsat_ice(k)/(p(k,i)-qsat_ice(k))
         qsat_ice(k) = max(qsat_ice(k),qmin)
         rh_ice(k)   = max(q(k,i)/qsat_ice(k),qmin)
       enddo
!
! initialize the variables for microphysical physics
!
       prevp(:) = 0.     ; psdep(:) = 0.     ; praut(:) = 0.     
       psaut(:) = 0.     ; pracw(:) = 0.     ; psaci(:) = 0.     
       psacw(:) = 0.     ; pigen(:) = 0.     ; pidep(:) = 0.     
       pcond(:) = 0.     ; psmlt(:) = 0.     ; psevp(:) = 0.
       pihtf(:) = 0.     ; psfrz(:) = 0.
       factor1(:) = 1.   ; cldf(:) = 0.
       ni(:)   = 1.e3    ; mi(:) = 0.        ; diameter(:) = 0.
       ni0(:)  = 1.e3    ; roqi0(:) = 0.     ; qimax(:) = 0.
       denqr(:) = 0.     ; denqs(:) = 0.     ; denqi(:) = 0.
       qrpath   = 0.     ; qspath   = 0.     ; qipath   = 0.
       precp_r  = 0.     ; precp_s  = 0.     ; precp_i  = 0.
       precp_sum = 0.    ; precp_ice = 0.
       vtr(:)   = 0.     ; vts(:) = 0.       ; vti(:) = 0.
       vt2i(:)  = 0.     ; vt2s(:) = 0.     
       rslope(:,:) = 0.  ; rslopeb(:,:) = 0. ; rslope2(:,:) = 0. 
       rslope3(:,:) = 0.
#ifdef VD_RAIN
       rslopef(:)  = 0.  ; rslopefb(:) = 0.  ; rslopef2(:) = 0.
       rslopef3(:) = 0.  ; rslopehf(:) = 0.  ; rslopehfb(:) = 0.
       rslopehf2(:) = 0. ; rslopehf3(:) = 0.
#endif
       ab_wat(:) = 0.    ; ab_ice(:) = 0.    ; venfac(:) = 0.
       supcol(:) = 0.    ; supcolt(:)= 0.    ; supsat(:) = 0.
       satdt(:)  = 0.    ; supice(:) = 0.    ; ifsat(:) = 0
       xlf(:)    = 0.    ; coeres(:) = 0.    ; n0sfac(:) = 0.
       eacrs(:)  = 0.    ; acrfac(:) = 0.    ; total(:) = 0.
!
! redefine cloud top for numerical efficiencies
!
       lqc    = .false.
       lqi    = .false.
       lqr    = .false.
       lqs    = .false.
       flgcld = .false.
!
       call find_cloud_top(ktop,qci(kts:kte,i,1),0.0,kts,kte,ktopqc)
       call find_cloud_top(ktop,qci(kts:kte,i,2),0.0,kts,kte,ktopqi)
       call find_cloud_top(ktop,qrs(kts:kte,i,1),0.0,kts,kte,ktopqr)
       call find_cloud_top(ktop,qrs(kts:kte,i,2),0.0,kts,kte,ktopqs)
       call find_cloud_top(ktop,rh_ice(kts:kte), 1.0,kts,kte,ktoprh)
!
       if (ktopqc.gt.0) lqc = .true.
       if (ktopqi.gt.0) lqi = .true.
       if (ktopqr.gt.0) lqr = .true.
       if (ktopqs.gt.0) lqs = .true.
!
! early checkout
!
       if (lqc .or. lqi .or. lqr .or. lqs) flgcld = .true.
       if ((.not.flgcld) .and. ktoprh.gt.0) flgcld = .true.
       if (.not.flgcld) cycle i_loop
!
!===============================================================================
!
! compute the fallout
!
!===============================================================================
!
! for rain
!
       ktop = max(ktopqr,ktopqs,ktopqi)
       if (lqr) then
!
! vertical velocity 
!
         call slope_rain(qrs(kts:kte,i,1),den(kts:kte,i),denfac(kts:kte,i),    &
                         rslope(kts:kte,1),rslopeb(kts:kte,1),                 &
                         rslope2(kts:kte,1),rslope3(kts:kte,1),                &
#ifdef VD_RAIN
                         rslopef(kts:kte),rslopefb(kts:kte),                   &
                         rslopef2(kts:kte),rslopef3(kts:kte),                  &
                         rslopehf(kts:kte),rslopehfb(kts:kte),                 &
                         rslopehf2(kts:kte),rslopehf3(kts:kte),                &
#endif
                         vtr(kts:kte),kts,kte,ktop)
!
         niter = 1 
         dtcfl = dtcld/niter
!
! subloop
!
         do n = 1,niter
           if (n.ge.2) then
             call slope_rain(qrs(kts:kte,i,1),den(kts:kte,i),denfac(kts:kte,i),&
                             rslope(kts:kte,1),rslopeb(kts:kte,1),             &
                             rslope2(kts:kte,1),rslope3(kts:kte,1),            &
#ifdef VD_RAIN
                             rslopef(kts:kte),rslopefb(kts:kte),               &
                             rslopef2(kts:kte),rslopef3(kts:kte),              &
                             rslopehf(kts:kte),rslopehfb(kts:kte),             &
                             rslopehf2(kts:kte),rslopehf3(kts:kte),            &
#endif
                             vtr(kts:kte),kts,kte,ktop)
           endif
!
           do k = ktop,kts,-1
             denqr(k) = dend(k,i)*qrs(k,i,1)
           enddo
!
           call nislfv_precp(1,kdim,dend(kts:kte,i),denfac(kts:kte,i),         &
                             delz(kts:kte,i),vtr(kts:kte),denqr(kts:kte),      &
                             qrpath,dtcfl,ktop)
!
           do k = kts,ktop
             qrs(k,i,1) = max(denqr(k)/dend(k,i),0.)
           enddo
!
           precp_r = precp_r+qrpath/dtcfl     ! kg m-2 s-1
         enddo     ! n loop
       endif       ! lqr
!-------------------------------------------------------------------------------
!
! for snow
!
       if (lqs) then
!
! vertical velocity
!
         call slope_snow(qrs(kts:kte,i,2),den(kts:kte,i),denfac(kts:kte,i),    &
                         t(kts:kte,i),rslope(kts:kte,2),rslopeb(kts:kte,2),    &
                         rslope2(kts:kte,2),rslope3(kts:kte,2),vts(kts:kte),   &
                         t0c,kts,kte,ktop)
!
         niter = 1
         dtcfl = dtcld/niter
!
! subloop
!
         do n = 1,niter
           if (n.ge.2) then
             call slope_snow(qrs(kts:kte,i,2),den(kts:kte,i),denfac(kts:kte,i),&
                            t(kts:kte,i),rslope(kts:kte,2),rslopeb(kts:kte,2), &
                            rslope2(kts:kte,2),rslope3(kts:kte,2),vts(kts:kte),&
                            t0c,kts,kte,ktop)
           endif
!
           do k = ktop,kts,-1
             denqs(k) = dend(k,i)*qrs(k,i,2)
           enddo
!
           call nislfv_precp(1,kdim,dend(kts:kte,i),denfac(kts:kte,i),         &
                             delz(kts:kte,i),vts(kts:kte),denqs(kts:kte),      &
                             qspath,dtcfl,ktop)
!
           do k = kts,ktop
             qrs(k,i,2) = max(denqs(k)/dend(k,i),0.)
           enddo
!
           precp_s = precp_s+qspath/dtcfl      ! kg m-2 s-1
         enddo       ! n loop
       endif         ! lqs
!-------------------------------------------------------------------------------
!
! for cloud ice
!
       if (lqi) then
         temp = 0.
!
! vertical velocity
!
         do k = ktop,kts,-1
           if (qci(k,i,2).le.0.) then
             vti(k) = 0.
           else
             temp = (den(k,i)*max(qci(k,i,2),qmin))
#ifdef ICE_PARAM
             ni(k) = min(max(ani*(temp**bni),1.e3),1.e6)
#else
             temp  = sqrt(sqrt(temp*temp*temp))
             ni(k) = min(max(5.38e7*temp,1.e3),1.e6)
#endif             
             mi(k) = den(k,i)*qci(k,i,2)/ni(k)
#ifdef ICE_PARAM
             diameter(k) = max(min(((1./cmi)**rdmi)*(mi(k)**rdmi),dimax)   &
                               ,1.e-25)
#else
             diameter(k) = max(min(dicon*sqrt(mi(k)),dimax),1.e-25)
#endif             
             vti(k) = avti*exp(log(diameter(k))*(bvti))
           endif
         enddo
!
         niter = 1
         dtcfl = dtcld/niter
!
! subloop
!
         do n = 1,niter
           if (n.ge.2) then
             do k = ktop,kts,-1
               if (qci(k,i,2).le.0.) then
                 vti(k) = 0.
               else
#ifdef ICE_PARAM
                 ni(k) = min(max(ani*(temp**bni),1.e3),1.e6)
#else

                 temp   = (den(k,i)*max(qci(k,i,2),qmin))
                 temp  = sqrt(sqrt(temp*temp*temp))
                 ni(k) = min(max(5.38e7*temp,1.e3),1.e6)
#endif                 
                 mi(k) = den(k,i)*qci(k,i,2)/ni(k)
#ifdef ICE_PARAM
                 diameter(k) = max(min(((1./cmi)**rdmi)*(mi(k)**rdmi),dimax)   &
                                   ,1.e-25)
#else
                 diameter(k) = max(min(dicon*sqrt(mi(k)),dimax),1.e-25)
#endif                 
                 vti(k) = avti*exp(log(diameter(k))*(bvti))
               endif ! qi<=0.
             enddo   ! k loop
           endif     ! n>=2
!
           do k = ktop,kts,-1
             denqi(k) = dend(k,i)*qci(k,i,2)
           enddo
!
           call nislfv_precp(1,kdim,dend(kts:kte,i),denfac(kts:kte,i),         &
                             delz(kts:kte,i),vti(kts:kte),denqi(kts:kte),      &
                             qipath,dtcfl,ktop)
!
           do k = kts,ktop
             qci(k,i,2) = max(denqi(k)/dend(k,i),0.)
           enddo
!
           precp_i = precp_i+qipath/dtcfl
         enddo     ! n loop
       endif       ! lqi
!-------------------------------------------------------------------------------
!
! total precp (den*qrsi*dz/dt; kg m-2 s-1) ===> rain (precp/denr*dt; m)
!
       precp_sum = precp_r+precp_s+precp_i
       precp_ice = precp_s+precp_i
       if (precp_sum.gt.0.) then
         rainncv(i) = precp_sum/denr*dtcld+rainncv(i)
         rain(i)    = precp_sum/denr*dtcld+rain(i)
       endif
       if (precp_ice.gt.0.) then
         tstepsnow(i) = precp_ice/denr*dtcld+tstepsnow(i)
         if (present(snowncv) .and. present(snow)) then
           snowncv(i) = precp_ice/denr*dtcld+snowncv(i)
           snow(i)    = precp_ice/denr*dtcld+snow(i)
         endif
       endif
!
       if (precp_sum.gt.0.) sr(i) = tstepsnow(i)/(rainncv(i)+1.e-12)
!
!-------------------------------------------------------------------------------
!
! change condensate variables to in-cloud variables
!
       factor1 = 1.
       ktop = max(ktopqc,ktopqi)
       call cldf_mps_diag(t(kts:kte,i),p(kts:kte,i),q(kts:kte,i),              &
                          qci(kts:kte,i,:),dxmeter(i),cldf(kts:kte),           &
                          kts,kte,ktop)
!
       do k = kts,ktop
         if(cldf(k).gt.0.) then
           qci(k,i,1) = qci(k,i,1)/cldf(k)
           qci(k,i,2) = qci(k,i,2)/cldf(k)
           qrs(k,i,1) = qrs(k,i,1)/cldf(k)
           qrs(k,i,2) = qrs(k,i,2)/cldf(k)
           factor1(k) = cldf(k)
         endif
       enddo
!
! update the slope parameters for microphysics computation
!
       ktop = ktopqr
       call slope_rain(qrs(kts:kte,i,1),den(kts:kte,i),denfac(kts:kte,i),      &
                       rslope(kts:kte,1),rslopeb(kts:kte,1),                   &
                       rslope2(kts:kte,1),rslope3(kts:kte,1),                  &
#ifdef VD_RAIN
                       rslopef(kts:kte),rslopefb(kts:kte),                     &
                       rslopef2(kts:kte),rslopef3(kts:kte),                    &
                       rslopehf(kts:kte),rslopehfb(kts:kte),                   &
                       rslopehf2(kts:kte),rslopehf3(kts:kte),                  &
#endif

                       vtr(kts:kte),kts,kte,ktop)
!   
       ktop = ktopqs
       call slope_snow(qrs(kts:kte,i,2),den(kts:kte,i),denfac(kts:kte,i),      &
                       t(kts:kte,i),rslope(kts:kte,2),rslopeb(kts:kte,2),      &
                       rslope2(kts:kte,2),rslope3(kts:kte,2),vts(kts:kte),     &
                       t0c,kts,kte,ktop)
!
!===============================================================================
!
! melting and freezing
!
!===============================================================================
!
       if (lqs) then
         ktop = ktopqs
         do k = ktop,kts,-1
           supcol(k) = t0c-t(k,i)
           n0sfac(k) = max(min(exp(alpha*supcol(k)),n0smax/n0s),1.)
           if (t(k,i).gt.t0c .and. qrs(k,i,2).gt.0.) then
!
! psmlt: melting of snow [HL A33] [RH83 A25]
!       (T>T0: S->R)
!
             xlf(k) = xlf0
             venfac(k)= (exp(log(((1.496e-6*((t(k,i))*sqrt(t(k,i)))            &
                        /((t(k,i))+120.)/(den(k,i)))/(8.794e-5                 &
                        *exp(log(t(k,i))*(1.81))/p(k,i))))                     &
                        *((.3333333)))/sqrt((1.496e-6*((t(k,i))                &
                        *sqrt(t(k,i)))/((t(k,i))+120.)/(den(k,i))))            &
                        *sqrt(sqrt(den0/(den(k,i)))))
             coeres(k) = rslope2(k,2)*sqrt(rslope(k,2)*rslopeb(k,2))
             psmlt(k) = (1.414e3*(1.496e-6*((t(k,i))*sqrt(t(k,i)))             &
                        /((t(k,i))+120.)/(den(k,i)) )*(den(k,i)))              &
                        /xlf(k)*(t0c-t(k,i))*pi/2.                             &
                        *n0sfac(k)*(precs1*rslope2(k,2)+precs2                 &
                        *venfac(k)*coeres(k))/den(k,i)
             psmlt(k) = min(max(psmlt(k)*dtcld,-qrs(k,i,2)),0.)
             qrs(k,i,2) = qrs(k,i,2)+psmlt(k)
             qrs(k,i,1) = qrs(k,i,1)-psmlt(k)
             t(k,i) = t(k,i)+xlf(k)/cpm(k,i)*psmlt(k)*factor1(k)
           endif   ! t>t0c & qs>0
         enddo     ! k loop 
       endif       ! lqs
!-------------------------------------------------------------------------------
!
       ktop = max(ktopqc,ktopqi,ktopqs)
       do k = kts,ktop
         supcol(k) = t0c-t(k,i)
         xlf(k) = xls-xl(k,i)
         if (supcol(k).lt.0.) xlf(k) = xlf0
!
! pimlt: instantaneous melting of cloud ice [HL A47] [RH83 A28]
!       (T>T0: I->C)
!
         if (lqi) then
           if (supcol(k).lt.0 .and. qci(k,i,2).gt.0.) then
             qci(k,i,1) = qci(k,i,1)+qci(k,i,2)
             t(k,i)     = t(k,i)-xlf(k)/cpm(k,i)*qci(k,i,2)*factor1(k)
             qci(k,i,2) = 0.
           endif   ! t>0 & qi>0
         endif     ! lqi
!
! pihmf: homogeneous freezing of cloud water below -40c [HL A45]
!        (T<-40C: C->I)
!
         if (lqc) then
           if (supcol(k).gt.40. .and. qci(k,i,1).gt.0.) then
             qci(k,i,2) = qci(k,i,2)+qci(k,i,1)
             t(k,i)     = t(k,i)+xlf(k)/cpm(k,i)*qci(k,i,1)*factor1(k)
             qci(k,i,1) = 0.
           endif   ! t<-40 & qc>0
!
! pihtf: heterogeneous freezing of cloud water [HL A44]
!        (T0>T>-40C: C->I)
!
           if (supcol(k).gt.0. .and. qci(k,i,1).gt.0.) then
             supcolt(k) = min(supcol(k),50.)
             pihtf(k)   = min(pfrz1*(exp(pfrz2*supcolt(k))-1.)                 &
                              *den(k,i)/denr/xncr*qci(k,i,1)*qci(k,i,1)*dtcld, &  
                              qci(k,i,1))
             qci(k,i,2) = qci(k,i,2)+pihtf(k)
             qci(k,i,1) = qci(k,i,1)-pihtf(k)
             t(k,i) = t(k,i)+xlf(k)/cpm(k,i)*pihtf(k)*factor1(k)
           endif   ! t<0 & qc>0
         endif     ! lqc
!
! psfrz: freezing of rain water [HL A20] [LFO 45]
!        (T<T0, R->S)
!
         temp = 0.
         if (lqs) then
           if (supcol(k).gt.0. .and. qrs(k,i,1).gt.0.) then
             supcolt(k) = min(supcol(k),50.)
             temp = rslope(k,1)
             temp = temp*temp*temp*temp*temp*temp*temp
             psfrz(k) = min(20.*(pi*pi)*pfrz1*n0r*denr/den(k,i)                &
                           *(exp(pfrz2*supcolt(k))-1.)*temp*dtcld,             &
                             qrs(k,i,1))
             qrs(k,i,2) = qrs(k,i,2)+psfrz(k)
             qrs(k,i,1) = qrs(k,i,1)-psfrz(k)
             t(k,i) = t(k,i)+xlf(k)/cpm(k,i)*psfrz(k)*factor1(k)
           endif   ! t<0 & qs>0
         endif     ! lqs
       enddo       ! k loop
!
!----------------------------------------------------------------------------
!
! check hydrometeors
!
       ktop = kmaxin
       lqc = .false.
       lqi = .false.
       lqr = .false.
       lqs = .false.
!
       call find_cloud_top(ktop,qrs(kts:kte,i,1),0.0,kts,kte,ktopqr)
       call find_cloud_top(ktop,qrs(kts:kte,i,2),0.0,kts,kte,ktopqs)
       call find_cloud_top(ktop,qci(kts:kte,i,1),0.0,kts,kte,ktopqc)
       call find_cloud_top(ktop,qci(kts:kte,i,2),0.0,kts,kte,ktopqi)
       call find_cloud_top(ktop,rh_ice(kts:kte), 1.0,kts,kte,ktoprh)
!
       if (ktopqc.gt.0) lqc = .true.
       if (ktopqi.gt.0) lqi = .true.
       if (ktopqr.gt.0) lqr = .true.
       if (ktopqs.gt.0) lqs = .true.
!
!-------------------------------------------------------------------------------
!
! change in-cloud condensate variables to grid-mean variables
!
        ktop = max(ktopqc,ktopqi)
        do k = kts,ktop
          if(cldf(k).gt.0.) then
            qci(k,i,1) = qci(k,i,1)*cldf(k)
            qci(k,i,2) = qci(k,i,2)*cldf(k)
            qrs(k,i,1) = qrs(k,i,1)*cldf(k)
            qrs(k,i,2) = qrs(k,i,2)*cldf(k)
          endif
        enddo
!
! change condensate variables to in-cloud variables
!
       call cldf_mps_diag(t(kts:kte,i),p(kts:kte,i),q(kts:kte,i),              &
                          qci(kts:kte,i,:),dxmeter(i),cldf(kts:kte),           &
                          kts,kte,ktop)
!
       do k = kts,ktop
          if(cldf(k).gt.0.) then
            qci(k,i,1) = qci(k,i,1)/cldf(k)
            qci(k,i,2) = qci(k,i,2)/cldf(k)
            qrs(k,i,1) = qrs(k,i,1)/cldf(k)
            qrs(k,i,2) = qrs(k,i,2)/cldf(k)
          endif
       enddo
!
! update the slope parameters for microphysics computation
!
       ktop = ktopqr
       call slope_rain(qrs(kts:kte,i,1),den(kts:kte,i),denfac(kts:kte,i),      &
                       rslope(kts:kte,1),rslopeb(kts:kte,1),                   &
                       rslope2(kts:kte,1),rslope3(kts:kte,1),                  &
#ifdef VD_RAIN
                       rslopef(kts:kte),rslopefb(kts:kte),                     &
                       rslopef2(kts:kte),rslopef3(kts:kte),                    &
                       rslopehf(kts:kte),rslopehfb(kts:kte),                   &
                       rslopehf2(kts:kte),rslopehf3(kts:kte),                  &
#endif

                       vtr(kts:kte),kts,kte,ktop)
!
       ktop = ktopqs
       call slope_snow(qrs(kts:kte,i,2),den(kts:kte,i),denfac(kts:kte,i),      &
                       t(kts:kte,i),rslope(kts:kte,2),rslopeb(kts:kte,2),      &
                       rslope2(kts:kte,2),rslope3(kts:kte,2),vts(kts:kte),     &
                       t0c,kts,kte,ktop)
!
! ab  :  the thermodynamic term in the denominator associated with
!         heat conduction and vapor diffusion
!         (ry88, y93, h85)
!         ab_wat(k) = diffac(xl(k,i),p(k,i),t(k,i),den(k,i),qsat_wat(k))
!         ab_ice(k) = diffac(xls(k,i),p(k,i),t(k,i),den(k,i),qsat_ice(k))
! venfac: parameter associated with the ventilation effects(y93)
!        venfac(k) = venfac(p(k,i),t(k,i),den(k,i))
!
       ktop = kmaxin
       do k = kts,ktop
         ab_wat(k) = ((((den(k,i))*(xl(k,i))*(xl(k,i)))*((t(k,i))+120.)     &
                       *(den(k,i)))/(1.414e3*(1.496e-6*((t(k,i))*sqrt(t(k,i))))&
                       *(den(k,i))*(rv*(t(k,i))*(t(k,i)))))                    &
                      +  p(k,i)/((qsat_wat(k))*(8.794e-5*exp(log(t(k,i))*(1.81))))
         ab_ice(k) = ((((den(k,i))*(xls)*(xls))*((t(k,i))+120.)*(den(k,i))) &
                       /(1.414e3*(1.496e-6*((t(k,i))*sqrt(t(k,i))))*(den(k,i)) &
                       *(rv*(t(k,i))*(t(k,i))))                                &
                      + p(k,i)/(qsat_ice(k)*(8.794e-5*exp(log(t(k,i))*(1.81)))))
         venfac(k) = (exp(.3333333*log(((1.496e-6 * ((t(k,i))*sqrt(t(k,i))))  &
                     *p(k,i))/(((t(k,i))+120.)*den(k,i)*(8.794e-5              &
                     *exp(log(t(k,i))*(1.81))))))*sqrt(sqrt(den0/(den(k,i))))) &
                      /sqrt((1.496e-6*((t(k,i))*sqrt(t(k,i))))                 &
                      /(((t(k,i))+120.)*den(k,i)))
       enddo 
!
!===============================================================================
!
! warm rain processes
!
! - follows the processes in RH83 and LFO except for autoconcersion
!
!===============================================================================
!
       ktop = max(ktopqc,ktopqr)
       do k = kts,ktop
         supsat(k) = max(q(k,i),qmin)-qsat_wat(k)
         satdt(k) = supsat(k)/dtcld
!
! praut: auto conversion rate from cloud to rain [HDC 16]
!        (C->R)
!
         if (lqc) then
           if (qci(k,i,1).gt.qc0) then
             praut(k) = qck1*exp(log(qci(k,i,1))*((7./3.)))
             praut(k) = min(praut(k),qci(k,i,1)/dtcld)
           endif   ! qc>qc0
         endif     ! lqc
!
! pracw: accretion of cloud water by rain [HL A40] [LFO 51]
!        (C->R)
!
         if (lqr) then
           if (qrs(k,i,1).gt.qcrmin.and.qci(k,i,1).gt.qmin) then
#ifdef VD_RAIN
             pracw(k) = min(pacrr*rslopef3(k)*rslopefb(k)                      &
#else
             pracw(k) = min(pacrr*rslope3(k,1)*rslopeb(k,1)                    &
#endif             
                           *qci(k,i,1)*denfac(k,i),qci(k,i,1)/dtcld)
           endif   ! qr>qcrmin & qc>qmin
!
! prevp: evaporation/condensation rate of rain [HDC 14]
!        (V->R or R->V)
!
           if (qrs(k,i,1).gt.0.) then
#ifdef VD_RAIN
             coeres(k) = rslopehf2(k)*sqrt(rslopehf(k)*rslopehfb(k))
#else
             coeres(k) = rslope2(k,1)*sqrt(rslope(k,1)*rslopeb(k,1))
#endif             
             prevp(k) = (rh_wat(k)-1.)*(precr1*rslope2(k,1)                      &
                        +precr2*venfac(k)*coeres(k))/ab_wat(k)
             prevp(k) = min(prevp(k),0.0)
             if (prevp(k).lt.0.) then
               prevp(k) = max(prevp(k),-qrs(k,i,1)/dtcld)
               prevp(k) = max(prevp(k),satdt(k)/2)
             else
               prevp(k) = min(prevp(k),satdt(k)/2)
             endif ! prevp<0
           endif   ! qr>0
         endif     ! lqr
       enddo       ! k loop
!
!===============================================================================
!
! cold rain processes
!
! - follows the revised ice microphysics processes in HDC
! - the processes same as in RH83 and RH84  and LFO behave
!   following ice crystal hapits defined in HDC, inclduing
!   intercept parameter for snow (n0s), ice crystal number
!   concentration (ni), ice nuclei number concentration
!   (n0i), ice diameter (d)
!
!===============================================================================
!
       rdtcld = 1./dtcld
       ktop = max(ktopqs,ktopqi,ktopqc,ktoprh)
!
       do k = kts,ktop
         supcol(k) = t0c-t(k,i)
         n0sfac(k) = max(min(exp(alpha*supcol(k)),n0smax/n0s),1.)
         supsat(k) = max(q(k,i),qmin)-qsat_ice(k)
         satdt(k)  = supsat(k)/dtcld
         ifsat(k)  = 0
!
! ni: ice crystal number concentraiton   [HDC 5c]
!
         temp = 0.
         temp = (den(k,i)*max(qci(k,i,2),qmin))
#ifdef ICE_PARAM
         ni(k) = min(max(ani*(temp**bni),1.e3),1.e6)
#else
         temp = sqrt(sqrt(temp*temp*temp))
         ni(k) = min(max(5.38e7*temp,1.e3),1.e6)
#endif         
         eacrs(k) = exp(0.07*(-supcol(k)))
!
         if (supcol(k).gt.0) then
           if (qrs(k,i,2).gt.qcrmin.and.qci(k,i,2).gt.qmin) then
             mi(k) = den(k,i)*qci(k,i,2)/ni(k)
#ifdef ICE_PARAM
             diameter(k) = min((1./cmi)**rdmi*(mi(k)**rdmi),dimax)
#else             
             diameter(k)  = min(dicon*sqrt(mi(k)),dimax)
#endif             
             vt2i(k) = avti*diameter(k)**bvti
             vt2s(k) = pvts*rslopeb(k,2)*denfac(k,i)
!
! psaci: Accretion of cloud ice by rain [HDC 10]
!        (T<T0: I->S)
!
             acrfac(k) = 2.*rslope3(k,2)+2.*diameter(k)*rslope2(k,2)           &
                       +diameter(k)**2*rslope(k,2)
             psaci(k) = pi*qci(k,i,2)*eacrs(k)*n0s*n0sfac(k)                   &
                       *abs(vt2s(k)-vt2i(k))*acrfac(k)/4.
           endif   ! qs>qcrmin, qi>qmin
         endif     ! t<t0c
!
! psacw: Accretion of cloud water by snow  [HL A7] [LFO 24]
!        (T<T0: C->S, and T>=T0: C->R)
!
         if (qrs(k,i,2).gt.qcrmin .and. qci(k,i,1).gt.qmin) then
           psacw(k) = min(pacrc*n0sfac(k)*rslope3(k,2)*rslopeb(k,2)            &
#ifdef EFFIC
                          *min(max(0.0,qrs(k,i,2)/qci(k,i,1)),1.)**2.          &
#endif           
                          *qci(k,i,1)*denfac(k,i),qci(k,i,1)*rdtcld)
         endif
!
         if (supcol(k).gt.0) then
!
! pidep: Deposition/Sublimation rate of ice [HDC 9]
!       (T<T0: V->I or I->V)
!
           if (qci(k,i,2).gt.0 .and. ifsat(k).ne.1) then
             mi(k)       = den(k,i)*qci(k,i,2)/ni(k)
#ifdef ICE_PARAM
             diameter(k) = min((1./cmi)**rdmi*(mi(k)**rdmi),dimax)
#else
             diameter(k) = dicon*sqrt(mi(k))
#endif             
             pidep(k)    = 4.*diameter(k)*ni(k)*(rh_ice(k)-1.)/ab_ice(k)
             supice(k)   = satdt(k)-prevp(k)
             if (pidep(k).lt.0.) then
               pidep(k) = max(max(pidep(k),satdt(k)*.5),supice(k))
               pidep(k) = max(pidep(k),-qci(k,i,2)*rdtcld)
             else
               pidep(k) = min(min(pidep(k),satdt(k)*.5),supice(k))
             endif
             if (abs(prevp(k)+pidep(k)).ge.abs(satdt(k))) ifsat(k) = 1
           endif   ! qi>0 & ifsat=0
!
! psdep: deposition/sublimation rate of snow [HDC 14]
!        (V->S or S->V)
!
           if (qrs(k,i,2).gt.0..and.ifsat(k).ne.1) then
             coeres(k) = rslope2(k,2)*sqrt(rslope(k,2)*rslopeb(k,2))
             psdep(k)  = (rh_ice(k)-1.)*n0sfac(k)*(precs1*rslope2(k,2)         &
                          +precs2*venfac(k)*coeres(k))/ab_ice(k)
             supice(k) = satdt(k)-prevp(k)-pidep(k)
             if (psdep(k).lt.0.) then
               psdep(k) = max(psdep(k),-qrs(k,i,2)*rdtcld)
               psdep(k) = max(max(psdep(k),satdt(k)*.5),supice(k))
             else
               psdep(k) = min(min(psdep(k),satdt(k)*.5),supice(k))
             endif
             if (abs(prevp(k)+pidep(k)+psdep(k)).ge.abs(satdt(k))) ifsat(k) = 1
           endif   ! qs>0 & ifsat=0
!
! pigen: generation(nucleation) of ice from vapor [HL A50] [HDC 7-8]
!       (T<T0: V->I)
!
           if (supsat(k).gt.0.and.ifsat(k).ne.1) then
             supice(k) = satdt(k)-prevp(k)-pidep(k)-psdep(k)
             ni0(k)    = 1.e3*exp(0.1*supcol(k))
#ifdef ICE_PARAM
             temp      = bvti/(bvti-dmi*fvq)
             roqi0(k)  = ((cmi*(evq/avti)**(dmi/bvti))**temp)*(ni0(k)**temp)
#else
             roqi0(k)  = 4.92e-11*(ni0(k)**1.33)
#endif            
            ! roqi0(k,i)  = 4.92e-11*exp(log(xni0(k,i))*(1.33))
             pigen(k)  = max(0.,(roqi0(k)/den(k,i)-max(qci(k,i,2),0.))     &
                           *rdtcld)
             pigen(k) = min(min(pigen(k),satdt(k)),supice(k))
           endif   ! subsat>0 & ifsat=0
!
! psaut: conversion(aggregation) of ice to snow [HDC 12]
!       (T<T0: I->S)
!
           if (qci(k,i,2).gt.0.) then
#ifdef ICE_PARAM
             qimax(k) = (avti/evq*dimax**bvti)**(1./fvq)/den(k,i)
             psaut(k) = max(0.,alpha1*(qci(i,k,2)-qimax(k))*rdtcld)
#else
             qimax(k) = roqimax/den(k,i)
             psaut(k) = max(0.,(qci(k,i,2)-qimax(k))*rdtcld)
#endif             
           endif   ! qi>0
         endif     ! t<t0c (pidep, psdep, pigen, psaut)
!
! psevp: Evaporation of melting snow [HL A35] [RH83 A27]
!       (T>T0: S->V)
!
         if (supcol(k).lt.0.) then
           if (qrs(k,i,2).gt.0..and.rh_wat(k).lt.1.) then                          
             psevp(k) = psdep(k)*ab_ice(k)/ab_wat(k)
           endif
           psevp(k) = min(max(psevp(k),-qrs(k,i,2)*rdtcld),0.)
         endif     ! t>t0c
       enddo       ! k loop
!-------------------------------------------------------------------------------
!
! change in-cloud condensate variables to grid-mean variables
! get grid-mean rates
!
       ktop = max(ktopqc,ktopqi)
       do k = kts,ktop
         if(cldf(k).gt.0.) then
           qci(k,i,1) = qci(k,i,1)*cldf(k)
           qci(k,i,2) = qci(k,i,2)*cldf(k)
           qrs(k,i,1) = qrs(k,i,1)*cldf(k)
           qrs(k,i,2) = qrs(k,i,2)*cldf(k)
!
           praut(k) = praut(k)*cldf(k)
           pracw(k) = pracw(k)*cldf(k)
           prevp(k) = prevp(k)*cldf(k)
           psaci(k) = psaci(k)*cldf(k)
           psacw(k) = psacw(k)*cldf(k)
           pidep(k) = pidep(k)*cldf(k)
           psdep(k) = psdep(k)*cldf(k)
           pigen(k) = pigen(k)*cldf(k)
           psaut(k) = psaut(k)*cldf(k)
           psevp(k) = psevp(k)*cldf(k)
         endif
       enddo
!
!===============================================================================
!
! update
!   - check mass conservation of generation terms and feedback to the
!     large scale
!
!===============================================================================
!
       ktop = kmaxin
       do k = kts,ktop
         if (t(k,i).le.t0c) then
!
! 1) t<=t0c
!
! cloud water
!
           value = max(qmin,qci(k,i,1))
           source = (praut(k)+pracw(k)+psacw(k))*dtcld
           if (source.gt.value) then
             factor = value/source
             praut(k) = praut(k)*factor
             pracw(k) = pracw(k)*factor
             psacw(k) = psacw(k)*factor
           endif
!
! cloud ice
!
           value = max(qmin,qci(k,i,2))
           source = (psaut(k)+psaci(k)-pigen(k)-pidep(k))*dtcld
           if (source.gt.value) then
             factor = value/source
             psaut(k) = psaut(k)*factor
             psaci(k) = psaci(k)*factor
             pigen(k) = pigen(k)*factor
             pidep(k) = pidep(k)*factor
           endif
!
! rain
!
           value = max(qmin,qrs(k,i,1))
           source = (-praut(k)-pracw(k)-prevp(k))*dtcld
           if (source.gt.value) then
             factor = value/source
             praut(k) = praut(k)*factor
             pracw(k) = pracw(k)*factor
             prevp(k) = prevp(k)*factor
           endif
!
! snow
!
           value = max(qmin,qrs(k,i,2))
           source = (-psdep(k)-psaut(k)-psaci(k)-psacw(k))*dtcld  
           if (source.gt.value) then
             factor = value/source
             psdep(k) = psdep(k)*factor
             psaut(k) = psaut(k)*factor
             psaci(k) = psaci(k)*factor
             psacw(k) = psacw(k)*factor
           endif
!
! update
!
           q(k,i)     = q(k,i)-(prevp(k)+psdep(k)+pigen(k)+pidep(k))*dtcld
           qci(k,i,1) = max(qci(k,i,1)-(praut(k)+pracw(k)+psacw(k))*dtcld,0.)
           qrs(k,i,1) = max(qrs(k,i,1)+(praut(k)+pracw(k)+prevp(k))*dtcld,0.)
           qci(k,i,2) = max(qci(k,i,2)-(psaut(k)+psaci(k)-pigen(k)-pidep(k))   &
                            *dtcld,0.)
           qrs(k,i,2) = max(qrs(k,i,2)+(psdep(k)+psaut(k)+psaci(k)+psacw(k))   &
                            *dtcld,0.)
           xlf(k)     = xls-xl(k,i)
           xlvalue(k) = -xls*(psdep(k)+pidep(k)+pigen(k))-xl(k,i)*prevp(k)     &
                         -xlf(k)*psacw(k)
           t(k,i)     = t(k,i)-xlvalue(k)/cpm(k,i)*dtcld
         else
!-------------------------------------------------------------------------------
!
! 2) t>t0c
!
! cloud water
!
           value = max(qmin,qci(k,i,1))
           source=(praut(k)+pracw(k)+psacw(k))*dtcld
           if (source.gt.value) then
             factor = value/source
             praut(k) = praut(k)*factor
             pracw(k) = pracw(k)*factor
             psacw(k) = psacw(k)*factor
           endif
!
! rain
!
           value = max(qmin,qrs(k,i,1))
           source = (-praut(k)-pracw(k)-prevp(k)-psacw(k))*dtcld
           if (source.gt.value) then
             factor = value/source
             praut(k) = praut(k)*factor
             pracw(k) = pracw(k)*factor
             prevp(k) = prevp(k)*factor
             psacw(k) = psacw(k)*factor
           endif  
!
! snow
!
           value = max(qcrmin,qrs(k,i,2))
           source=(-psevp(k))*dtcld
           if (source.gt.value) then
             factor = value/source
             psevp(k) = psevp(k)*factor
           endif
!
! update
!
           q(k,i)     = q(k,i)-(prevp(k)+psevp(k))*dtcld
           qci(k,i,1) = max(qci(k,i,1)-(praut(k)+pracw(k)+psacw(k))*dtcld,0.)
           qrs(k,i,1) = max(qrs(k,i,1)+(praut(k)+pracw(k)+prevp(k)+psacw(k))   &
                            *dtcld,0.)
           qrs(k,i,2) = max(qrs(k,i,2)+psevp(k)*dtcld,0.)
           xlvalue(k) = -xl(k,i)*(prevp(k)+psevp(k))
           t(k,i)     = t(k,i)-xlvalue(k)/cpm(k,i)*dtcld
         endif
       enddo
!
!-------------------------------------------------------------------------------
!
! Inline expansion for fsvp
!         qsat_wat(k) = fsvp(t(k,i),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!
! constant for fsvp
!
       hsub = xls
       hvap = xlv0
       cvap = cpv
       dldt = cvap-cliq
       xa = -dldt/rv
       xb = xa+hvap/(rv*ttp)
!
       ktop = kmaxin
       do k = kts,ktop
         tr(k) = ttp/t(k,i)
         logtr(k) = log(tr(k))
         qsat_wat(k) = psat*exp(logtr(k)*(xa)+xb*(1.-tr(k)))
         qsat_wat(k) = min(qsat_wat(k),0.99*p(k,i))
         qsat_wat(k) = ep2*qsat_wat(k)/(p(k,i)-qsat_wat(k))
         qsat_wat(k) = max(qsat_wat(k),qmin)
       enddo
!
!  pcond: condensational/evaporational rate of cloud water [HL A46] [RH83 A6]
!     if there exists additional water vapor condensated/if
!     evaporation of cloud water is not enough to remove subsaturation
!
!     value1(k) = conden(t(k,i),q(k,i),qsat_wat(k),xl(k,i),cpm(k,i))
!
       do k = kts,ktop
         value1(k) = ((max(q(k,i),qmin)-(qsat_wat(k)))/(1.+(xl(k,i))           &
                        *(xl(k,i))/(rv*(cpm(k,i)))*(qsat_wat(k))               &
                        /((t(k,i))*(t(k,i)))))
         pcond(k) = min(max(value1(k)/dtcld,0.),max(q(k,i),0.)/dtcld)
!
         if (qci(k,i,1).gt.0. .and. value1(k).lt.0.) then
           pcond(k) = max(value1(k),-qci(k,i,1))/dtcld
         endif
         q(k,i) = q(k,i)-pcond(k)*dtcld
         qci(k,i,1) = max(qci(k,i,1)+pcond(k)*dtcld,0.)
         t(k,i) = t(k,i)+pcond(k)*xl(k,i)/cpm(k,i)*dtcld
       enddo
!-------------------------------------------------------------------------------
!
! padding for small values
!
       do k = kts,ktop
         t1(i,k) = t(k,i)
         q1(i,k) = q(k,i)
         q2(i,k,1) = max(qci(k,i,1),0.)
         q2(i,k,2) = max(qrs(k,i,1),0.)
         q2(i,k,3) = max(qci(k,i,2),0.)
         q2(i,k,4) = max(qrs(k,i,2),0.)
       enddo
     enddo i_loop   ! i_loop
   enddo inner_loop ! inner_loop
!
! compute reflectivity
!
     if(l_refl) then
!
       do k = kts,ktop
         do i = its,ite
           qrs1(i,k) = (qrs(i,k,1)+qrsub(i,k))*dend(i,k)
           qrs2(i,k) = (qrs(i,k,2)+qssub(i,k))*dend(i,k)
         enddo
       enddo
!
       call refl10cm(idim, kdim, qrs1, qrs2, qrs3, t, n0r, n0s, n0smax, n0g,   &
                     t0c, alpha, refl, ktop)
!
       do k = kts,ktop
         do i = its,ite
           refl_10cm(i,k) = refl(i,k)
         enddo
       enddo
     endif
!
!   
   end subroutine phy_mps_wsm5
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine wsm5init(den0,denr,dens,pi,cl,cpv,allowed_to_read)
!-------------------------------------------------------------------------------
!
!  abstract :
!
!  input :
!    den0             - density
!    denr             - density of water (rain)
!    dens             - density of snow
!    cl               -
!    cpv              - 
!    pi               - 
!    allowed_to_read
!
!-------------------------------------------------------------------------------
   use phy_funct,  only : rgmma
!   
   implicit none
!
! constants which may not be tunable
!
   real   , intent(in   ) :: den0
   real   , intent(in   ) :: denr
   real   , intent(in   ) :: dens
   real   , intent(in   ) :: pi
   real   , intent(in   ) :: cl
   real   , intent(in   ) :: cpv
   logical, intent(in   ) :: allowed_to_read
!-------------------------------------------------------------------------------
   xlv1 = cl-cpv
!
   qc0  = 4./3.*pi*denr*r0**3*xncr/den0  ! 0.419e-3 -- .61e-3
   qck1 = .104*9.8*peaut/(xncr*denr)**(1./3.)/xmyu*den0**(4./3.) ! 7.03
   pidnc = pi*denr/6.
!
   bvtr1 = 1.+bvtr
   bvtr2 = 2.5+.5*bvtr
   bvtr3 = 3.+bvtr
   bvtr4 = 4.+bvtr
   g1pbr = rgmma(bvtr1)
   g3pbr = rgmma(bvtr3)
   g4pbr = rgmma(bvtr4)            ! 17.837825
   g5pbro2 = rgmma(bvtr2)          ! 1.8273
   pvtr = avtr*g4pbr/6.
   eacrr = 1.0
   pacrr = pi*n0r*avtr*g3pbr*.25*eacrr
   precr1 = 2.*pi*n0r*.78
   precr2 = 2.*pi*n0r*.31*avtr**.5*g5pbro2
   xmmax = (dimax/dicon)**2
   roqimax = 2.08e22*dimax**8
!
   bvts1 = 1.+bvts
   bvts2 = 2.5+.5*bvts
   bvts3 = 3.+bvts
   bvts4 = 4.+bvts
   g1pbs = rgmma(bvts1)    !.8875
   g3pbs = rgmma(bvts3)
   g4pbs = rgmma(bvts4)    ! 12.0786
   g5pbso2 = rgmma(bvts2)
   pvts  = avts*g4pbs/6.
   pacrs = pi*n0s*avts*g3pbs*.25
   precs1 = 4.*n0s*.65
   precs2 = 4.*n0s*.44*avts**.5*g5pbso2
   pidn0r = pi*denr*n0r
   pidn0s = pi*dens*n0s
   pacrc  = pi*n0s*avts*g3pbs*.25*eacrc
!
#ifdef ICE_PARAM
   cmi = (1./9.05138)**2.08
   ani = 1./cmi*(avti/evq)**(dmi/bvti)
   bni = 1.-dmi*fvq/bvti
!
#endif
   rslopermax = 1./lamdarmax
   rslopesmax = 1./lamdasmax
   rsloperbmax = rslopermax ** bvtr
   rslopesbmax = rslopesmax ** bvts
   rsloper2max = rslopermax * rslopermax
   rslopes2max = rslopesmax * rslopesmax
   rsloper3max = rsloper2max * rslopermax
   rslopes3max = rslopes2max * rslopesmax
!
   end subroutine wsm5init
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine slope_wsm5(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,      &
                         vt,kts,kte,ktop)
!-------------------------------------------------------------------------------
!
!  abstract : calculate intercept parameters and settling velocity
!
!  input :
!    qrs                      - mixing ratio of rain and snow
!    den                      - air density
!    denfac                   - sqrt(den0/den)
!    t                        - temperature
!    t0c                      - ice/water mix temperature
!    pi                       - 
!    its, ite, kts, kte, ktop - dimension
!
!  output :
!    rslope, rslopeb, rslope2, rslope3  - related intercept parameters
!    vt                                 - settling velocity
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer                   , intent(in   ) :: kts,kte
   integer                   , intent(in   ) :: ktop
   real, dimension(kts:kte)  , intent(in   ) :: den
   real, dimension(kts:kte)  , intent(in   ) :: denfac
   real, dimension(kts:kte)  , intent(in   ) :: t
   real, dimension(kts:kte,2), intent(in   ) :: qrs
   real, dimension(kts:kte,2), intent(  out) :: rslope
   real, dimension(kts:kte,2), intent(  out) :: rslopeb
   real, dimension(kts:kte,2), intent(  out) :: rslope2
   real, dimension(kts:kte,2), intent(  out) :: rslope3
   real, dimension(kts:kte,2), intent(  out) :: vt
   real, parameter                           :: t0c = 273.15  ! input
! 
! local variables
!
   real, dimension(kts:kte) :: n0sfac
   real, dimension(kts:kte) :: supcol
   real                     :: lamdar, lamdas
   real                     :: x, y, z
   integer                  :: i, j, k
!-------------------------------------------------------------------------------
!
! size distributions: (x=mixing ratio, y=air density):
! valid for mixing ratio > 1.e-9 kg/kg.
!
   lamdar(x,y)  = sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
   lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
!
   do k = kts,ktop
     supcol(k) = t0c-t(k)
!
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!
     n0sfac(k) = max(min(exp(alpha*supcol(k)),n0smax/n0s),1.)
!
     if (qrs(k,1).le.qcrmin)then
       rslope(k,1)  = rslopermax
       rslopeb(k,1) = rsloperbmax
       rslope2(k,1) = rsloper2max
       rslope3(k,1) = rsloper3max
     else
       rslope(k,1)  = 1./lamdar(qrs(k,1),den(k))
       rslopeb(k,1) = exp(log(rslope(k,1))*(bvtr))
       rslope2(k,1) = rslope(k,1)*rslope(k,1)
       rslope3(k,1) = rslope2(k,1)*rslope(k,1)
     endif
!
     if (qrs(k,2).le.qcrmin)then
       rslope(k,2)  = rslopesmax
       rslopeb(k,2) = rslopesbmax
       rslope2(k,2) = rslopes2max
       rslope3(k,2) = rslopes3max
     else
       rslope(k,2)  = 1./lamdas(qrs(k,2),den(k),n0sfac(k))
       rslopeb(k,2) = exp(log(rslope(k,2))*(bvts))
       rslope2(k,2) = rslope(k,2)*rslope(k,2)
       rslope3(k,2) = rslope2(k,2)*rslope(k,2)
     endif
!
     vt(k,1) = pvtr*rslopeb(k,1)*denfac(k)
     vt(k,2) = pvts*rslopeb(k,2)*denfac(k)
!
     if (qrs(k,1).le.0.0) vt(k,1) = 0.0
     if (qrs(k,2).le.0.0) vt(k,2) = 0.0
   enddo
!
   end subroutine slope_wsm5
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine slope_rain(qr,den,denfac,rslope,rslopeb,rslope2,rslope3,         &
#ifdef VD_RAIN
                         rslopef,rslopefb,rslopef2,rslopef3,                   &
                         rslopehf,rslopehfb,rslopehf2,rslopehf3,               &
#endif
                         vt,kts,kte,ktop)
!-------------------------------------------------------------------------------
!
!  abstract : calculate intercept parameters and settling velocity of rain
!
!  input :
!    qr                       - mixing ratio of rain 
!    den                      - air density
!    denfac                   -
!    t                        - temperature
!    t0c                      - 273.15 K
!    pi                       - 
!    its, ite, kts, kte, ktop - dimension
!
!  output :
!    rslope, rslopeb, rslope2, rslope3  - related intercept parameters
!    vt                                 - settling velocity
!
!-------------------------------------------------------------------------------
   implicit none
!
! parameters
!
   real, parameter                         :: t0c = 273.15  ! input
!
! passing variables
!
   integer                 , intent(in   ) :: kts,kte
   integer                 , intent(in   ) :: ktop
   real, dimension(kts:kte), intent(in   ) :: den
   real, dimension(kts:kte), intent(in   ) :: denfac
   real, dimension(kts:kte), intent(in   ) :: qr
   real, dimension(kts:kte), intent(  out) :: rslope
   real, dimension(kts:kte), intent(  out) :: rslopeb
   real, dimension(kts:kte), intent(  out) :: rslope2
   real, dimension(kts:kte), intent(  out) :: rslope3
#ifdef VD_RAIN
   real, dimension(kts:kte), intent(  out) :: rslopef
   real, dimension(kts:kte), intent(  out) :: rslopefb
   real, dimension(kts:kte), intent(  out) :: rslopef2
   real, dimension(kts:kte), intent(  out) :: rslopef3
   real, dimension(kts:kte), intent(  out) :: rslopehf
   real, dimension(kts:kte), intent(  out) :: rslopehfb
   real, dimension(kts:kte), intent(  out) :: rslopehf2
   real, dimension(kts:kte), intent(  out) :: rslopehf3
#endif
   real, dimension(kts:kte), intent(  out) :: vt
! 
! local variables
!
   real       :: lamdar
   real       :: x, y, z
   integer    :: i, j, k
!-------------------------------------------------------------------------------
!
! size distributions: (x=mixing ratio, y=air density):
! valid for mixing ratio > 1.e-9 kg/kg.
!
   lamdar(x,y) = sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
!
   do k = kts,ktop
     if (qr(k).le.qcrmin) then
       rslope(k)  = rslopermax
       rslopeb(k) = rsloperbmax
       rslope2(k) = rsloper2max
       rslope3(k) = rsloper3max
#ifdef VD_RAIN
       rslopef(k)   = rslopermax
       rslopefb(k)  = rsloperbmax
       rslopef2(k)  = rsloper2max
       rslopef3(k)  = rsloper3max
       rslopehf(k)  = rslopermax
       rslopehfb(k) = rsloperbmax
       rslopehf2(k) = rsloper2max
       rslopehf3(k) = rsloper3max
#endif

     else
       rslope(k)  = 1./lamdar(qr(k),den(k))
       rslopeb(k) = exp(log(rslope(k))*(bvtr))
       rslope2(k) = rslope(k)*rslope(k)
       rslope3(k) = rslope2(k)*rslope(k)
#ifdef VD_RAIN
       rslopef(k)   = 1./(lamdar(qr(k),den(k))+fr)
       rslopefb(k)  = exp(log(rslopef(k))*(bvtr))
       rslopef2(k)  = rslopef(k)*rslopef(k)
       rslopef3(k)  = rslopef2(k)*rslopef(k)
       rslopehf(k)  = 1./(lamdar(qr(k),den(k))+fr*0.5)
       rslopehfb(k) = exp(log(rslopehf(k))*(bvtr))
       rslopehf2(k) = rslopehf(k)*rslopehf(k)
       rslopehf3(k) = rslopehf2(k)*rslopehf(k)
#endif
     endif
!
#ifdef VD_RAIN
     vt(k) = pvtr*rslopefb(k)*denfac(k)*(rslopef2(k)/rslope2(k))**2.
#else
     vt(k) = pvtr*rslopeb(k)*denfac(k)
#endif     
     if (qr(k).le.0.0) vt(k) = 0.0
   enddo
!
   end subroutine slope_rain
!------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine slope_snow(qs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,      &
                         vt,t0c,kts,kte,ktop)
!-------------------------------------------------------------------------------
!
!  abstract : calculate intercept parameters and settling velocity of snow
!
!  input :
!    qs                       - mixing ratio of snow
!    den                      - air density
!    denfac                   - sqrt(den0/den)
!    t                        - temperature
!    t0c                      - 273.15 K
!    pi                       - 
!    its, ite, kts, kte, ktop - dimension
!
!  output :
!    rslope, rslopeb, rslope2, rslope3  - related intercept parameters
!    vt                                 - settling velocity
!
!-------------------------------------------------------------------------------
   implicit none
!
! passing variables
!
   integer                 , intent(in   ) :: kts,kte
   integer                 , intent(in   ) :: ktop
   real                    , intent(in   ) :: t0c
   real, dimension(kts:kte), intent(in   ) :: den
   real, dimension(kts:kte), intent(in   ) :: denfac
   real, dimension(kts:kte), intent(in   ) :: t
   real, dimension(kts:kte), intent(in   ) :: qs
   real, dimension(kts:kte), intent(  out) :: rslope
   real, dimension(kts:kte), intent(  out) :: rslopeb
   real, dimension(kts:kte), intent(  out) :: rslope2
   real, dimension(kts:kte), intent(  out) :: rslope3
   real, dimension(kts:kte), intent(  out) :: vt
! 
! local variables
!
   real, dimension(kts:kte) :: n0sfac
   real, dimension(kts:kte) :: supcol
   real                     :: lamdas
   real                     :: x, y, z
   integer                  :: i, j, k
!-------------------------------------------------------------------------------
!
! size distributions: (x=mixing ratio, y=air density):
! valid for mixing ratio > 1.e-9 kg/kg.
!
   lamdas(x,y,z) = sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
!
   do k = kts,ktop
     supcol(k) = t0c-t(k)
!
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!
     n0sfac(k) = max(min(exp(alpha*supcol(k)),n0smax/n0s),1.)
     if (qs(k).le.qcrmin) then
       rslope(k)  = rslopesmax
       rslopeb(k) = rslopesbmax
       rslope2(k) = rslopes2max
       rslope3(k) = rslopes3max
     else
       rslope(k)  = 1./lamdas(qs(k),den(k),n0sfac(k))
       rslopeb(k) = exp(log(rslope(k))*(bvts))
       rslope2(k) = rslope(k)*rslope(k)
       rslope3(k) = rslope2(k)*rslope(k)
     endif
!
     vt(k) = pvts*rslopeb(k)*denfac(k)
     if (qs(k).le.0.0) vt(k) = 0.0
   enddo
!
   end subroutine slope_snow
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nislfv_precp(im,km,den,denfac,dz,vt,rq,precip,dt,ktop)
!-------------------------------------------------------------------------------
!
!  abstract : non-iteration semi-Lagrangain forward advection for cloud
!             with mass conservation and positive definite advection
!             2nd order interpolation with monotonic piecewise linear method
!             this routine is under assumption of decfl < 1 for semi_Lagrangian
!
!  input :
!    im, km  - dimension
!    den     - dry air density
!    denfac  - sqrt(den0/den)
!    t       - temperature in K
!    dz      - depth of model layer in meter
!    vt      - terminal velocity at model layer in m/s
!    dt      - time step
!    ktop    - top layer for computing
!
!  inout :
!    precip   - precipitation
!    rql      - air density*mixing ratio
!     
!-------------------------------------------------------------------------------
   implicit none
!
! parameters
!
   real, parameter  :: fa1 = 9./16.
   real, parameter  :: fa2 = 1./16.
!
! passing variables
!
   integer            , intent(in   ) :: im, km
   integer            , intent(in   ) :: ktop
   real               , intent(in   ) :: dt
!   real, dimension(km,im), intent(in   ) :: dzl
!   real, dimension(km,im), intent(in   ) :: denl
!   real, dimension(km,im), intent(in   ) :: denfacl
!   real, dimension(km,im), intent(in   ) :: tkl
!   real, dimension(km),    intent(inout) :: vtl
!   real, dimension(km),    intent(inout) :: rq
!
   real, dimension(km), intent(in   ) :: dz
   real, dimension(km), intent(in   ) :: den
   real, dimension(km), intent(in   ) :: denfac
!   real, dimension(km), intent(in   ) :: t
   real, dimension(km), intent(in   ) :: vt
   real, dimension(km), intent(inout) :: rq
   real,                intent(inout) :: precip
!
!  local variables
!
 !  real, dimension(km)   :: dz
 !  real, dimension(km)   :: vt
 !  real, dimension(km)   :: wa
 !  real, dimension(km)   :: was
 !  real, dimension(km)   :: den
 !  real, dimension(km)   :: denfac
 !  real, dimension(km)   :: tk
!
   real, dimension(km)   :: qq
   real, dimension(km)   :: qn
   real, dimension(km)   :: qr
   real, dimension(km)   :: vtd
   real, dimension(km+1) :: vti
   real, dimension(km+1) :: zi
   real, dimension(km+2) :: za
   real, dimension(km+1) :: dza
   real, dimension(km+1) :: qa
   real, dimension(km+1) :: qmi
   real, dimension(km+1) :: qpi
   integer               ::  i, k, n, m, kk, kb, kt, iter
   real                  :: tl, tl2, qql, dql, qqd
   real                  :: th, th2, qqh, dqh
   real                  :: zsum, qsum, dim, dip, c1, con1
   real                  :: zsumt, qsumt, zsumb, qsumb
   real                  :: allold, allnew, zz, dzamin, cflmax, decfl
!-------------------------------------------------------------------------------
   precip = 0.0
   zi  = 0.   ; za   = 0.  ; dza = 0. 
   qq  = 0.   ; qn   = 0.  ; qr  =  0.   ; qa  = 0.   
   vti = 0.   ; vtd = 0.   
   qmi = 0.   ; qpi = 0.
!
   i_loop : do i = 1,im
     qq(:) = rq(:)
!     dz(:) = dzl(:,i)
!     vt(:) = vtl(:)
!     den(:) = denl(:,i)
!     denfac(:) = denfacl(:,i)
!     tk(:) = tkl(:,i)
!
! skip for no precipitation for all layers
!
     allold = 0.0
     do k = 1,ktop
       allold = allold + qq(k)
     enddo
     if (allold.le.0.0) then
       cycle i_loop
     endif
!
! compute interface values
!
     zi(1) = 0.0
     do k = 1,km
       zi(k+1) = zi(k)+dz(k)
     enddo
!
! 3rd order interpolation to get vti
!
     vti(1) = vt(1)
     vti(2) = 0.5*(vt(2)+vt(1))
     do k = 3,km-1
       vti(k) = fa1*(vt(k)+vt(k-1))-fa2*(vt(k+1)+vt(k-2))
     enddo
     vti(km) = 0.5*(vt(km)+vt(km-1))
     vti(km+1) = vt(km)
!
! terminate of top of raingroup
!
     do k = 2,ktop
       if (vt(k).eq.0.0) vti(k) = vt(k-1)
     enddo
!
! diffusivity of vti
!
     con1 = 0.05
     do k = km,1,-1
       decfl = (vti(k+1)-vti(k))*dt/dz(k)
       if (decfl.gt.con1 ) then
         vti(k) = vti(k+1)-con1*dz(k)/dt
       endif
     enddo
!
! compute arrival point (za)
!
     do k = 1,km+1
       za(k) = zi(k)-vti(k)*dt
     enddo
     za(km+2) = zi(km+1)
!
     do k = 1,km+1
       dza(k) = za(k+1)-za(k)
     enddo
!
! computer deformation at arrival point
!
     do k = 1,ktop
       qa(k) = qq(k)*dz(k)/dza(k)
     enddo
     qa(ktop+1) = 0.0
!
! estimate values at arrival cell interface with monotone
!
     do k = 2,ktop
       dip = (qa(k+1)-qa(k))/(dza(k+1)+dza(k))
       dim = (qa(k)-qa(k-1))/(dza(k-1)+dza(k))
       if (dip*dim.le.0.0) then
         qmi(k) = qa(k)
         qpi(k) = qa(k)
       else
         qpi(k) = qa(k)+0.5*(dip+dim)*dza(k)
         qmi(k) = 2.0*qa(k)-qpi(k)
         if (qpi(k).lt.0.0 .or. qmi(k).lt.0.0) then
           qpi(k) = qa(k)
           qmi(k) = qa(k)
         endif
       endif
     enddo
!
     qpi(1) = qa(1)
     qmi(1) = qa(1)
     qmi(ktop+1) = qa(ktop+1)
     qpi(ktop+1) = qa(ktop+1)
!
! interpolation to regular point
!
     qn = 0.0
     kb = 1
     kt = 1
     intp : do k = 1,ktop
       kb = max(kb-1,1)
       kt = max(kt-1,1)
!
! find kb and kt
!
       if (zi(k).ge.za(km+1)) then
         exit intp
       else
         find_kb : do kk = kb,ktop
           if (zi(k).le.za(kk+1)) then
             kb = kk
             exit find_kb
           else
             cycle find_kb
           endif
         enddo find_kb
         find_kt : do kk = kt,ktop+2
           if (zi(k+1).le.za(kk)) then
             kt = kk
             exit find_kt
           else
             cycle find_kt
           endif
         enddo find_kt
         kt = kt-1
!
! compute q with piecewise constant method
!
         if (kt.eq.kb) then
           tl = (zi(k)-za(kb))/dza(kb)
           th = (zi(k+1)-za(kb))/dza(kb)
           tl2 = tl*tl
           th2 = th*th
           qqd = 0.5*(qpi(kb)-qmi(kb))
           qqh = qqd*th2+qmi(kb)*th
           qql = qqd*tl2+qmi(kb)*tl
           qn(k) = (qqh-qql)/(th-tl)
         else if (kt.gt.kb) then
           tl = (zi(k)-za(kb))/dza(kb)
           tl2 = tl*tl
           qqd = 0.5*(qpi(kb)-qmi(kb))
           qql = qqd*tl2+qmi(kb)*tl
           dql = qa(kb)-qql
           zsum  = (1.-tl)*dza(kb)
           qsum  = dql*dza(kb)
           if (kt-kb.gt.1) then
             do m = kb+1,kt-1
               zsum = zsum+dza(m)
               qsum = qsum+qa(m)*dza(m)
             enddo
           endif
           th = (zi(k+1)-za(kt))/dza(kt)
           th2 = th*th
           qqd = 0.5*(qpi(kt)-qmi(kt))
           dqh = qqd*th2+qmi(kt)*th
           zsum  = zsum+th*dza(kt)
           qsum  = qsum+dqh*dza(kt)
           qn(k) = qsum/zsum
         endif
         cycle intp
       endif
!
     enddo intp
!
! rain out
!
     sum_precip: do k = 1,ktop
       if (za(k).lt.0.0 .and. za(k+1).le.0.0) then
         precip = precip+qa(k)*dza(k)
         cycle sum_precip
       else if ( za(k).lt.0.0 .and. za(k+1).gt.0.0 ) then    
         th = (0.0-za(k))/dza(k)               
         th2 = th*th                           
         qqd = 0.5*(qpi(k)-qmi(k))             
         qqh = qqd*th2+qmi(k)*th               
         precip = precip + qqh*dza(k)    
         exit sum_precip
       endif
       exit sum_precip
!
     enddo sum_precip
!
! replace the new values
!
     do k = 1,ktop
       rq(k) = qn(k)
     enddo
!
!-------------------------------------------------------------------------------
   enddo i_loop
!
   end subroutine nislfv_precp
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cldf_mps_diag(t, p, q, qci, dx, cldf, kts, kte, ktop)
!-------------------------------------------------------------------------------
   implicit none
!
! parameters
!
   real, parameter :: cldmin = 10.
   real, parameter :: cldmax = 100.
   real, parameter :: clddiff = cldmax - cldmin
   real, parameter :: cldf_min = 0.1
!
! passing variables
!
   integer,                    intent(in   ) :: kts, kte, ktop
   real,                       intent(in   ) :: dx
   real, dimension(kts:kte),   intent(in   ) :: t, p, q
   real, dimension(kts:kte,2), intent(in   ) :: qci
   real, dimension(kts:kte),   intent(  out) :: cldf
!
! local variables
!
   integer :: i, k, kk
   real    :: mps10, mps100
   real    :: dxkm
!----------------------------------------------------------------------------------
   do k = kts,kte
     cldf(k) = 0.
   enddo
!
! diagnostic method (kiaps)
!
   do k = kts,ktop
     mps10  = 5.57*(max(0.,qci(k,1)+qci(k,2))*1000.)**0.78
     mps100 = 4.37*(max(0.,qci(k,1)+qci(k,2))*1000.)**0.77
     dxkm  = dx/1000.
     cldf(k) = (mps10*(cldmax-dxkm)+mps100*(dxkm-cldmin))/clddiff
   enddo
!
   do k = kts,ktop
     cldf(k) = min(1.,max(0.,cldf(k)))
     if (qci(k,1)+qci(k,2).lt.1.e-6) cldf(k) = 0.
     if (cldf(k).lt.0.01) cldf(k) = 0.
     if (cldf(k).gt.0.99) cldf(k) = 1.
     if (cldf(k).ge.0.01 .and. cldf(k).lt.cldf_min) cldf(k) = cldf_min ! min
   enddo
!
   end subroutine cldf_mps_diag
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine find_cloud_top(ktop, qq, value, kts, kte, kout)
!-------------------------------------------------------------------------------
!
!  abstract :
!
!  input :
!    kts,kte  - dimension
!    q        - mixing ratio of hydrometeor, kg/kg
!    ktop     - top layer for computing
!
!  out :
!    kout     -
!
!-------------------------------------------------------------------------------
   implicit none
!
! passing variables
!
   integer,                     intent(in   ) :: kts, kte
   integer,                     intent(in   ) :: ktop
   real,                        intent(in   ) :: value
   real, dimension(kts:kte),    intent(in   ) :: qq
   integer,                     intent(  out) :: kout
!
! local variables
!
   integer  :: k
!-------------------------------------------------------------------------------
   kout = 0
!
   find_ktop : do k = ktop,kts,-1
     if (qq(k).gt.value) then
       kout = k
       exit find_ktop
     else
       cycle find_ktop
     endif
   enddo find_ktop
!
   return
!
   end subroutine find_cloud_top
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module module_mp_wsm5
!-------------------------------------------------------------------------------
