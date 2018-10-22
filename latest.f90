 MODULE WANG    
   implicit none
   integer, parameter:: in=10,im=2,igen=5           !! 'in' is the number of species in each patch, 'im' is the number of patches in the metacommunity
   integer, parameter:: imn=im*in, npar=imn*(imn+3)+im**2,ifact=1  !! npar is the number of parameters, ifact gives number of community pairs
   INTEGER, PARAMETER :: mlyap=0,NEQ=(mlyap+1)*imn,niter=50000,ntrans=10000,ninic=2000,net_Tscl=niter+10  !! ntrans is transient time, niter is iteration beyond transients
   real(kind=4), parameter:: d_max=14000.d0, a_max=1000.d0, xL=1.d0
   real(kind=4):: par(npar)

    CONTAINS

      SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER, intent(in):: NEQ
        real(kind=4), intent(in):: T, Y(NEQ) 
        real(kind=4), intent(out):: YDOT(NEQ)
        integer(kind=4)::i,j,iti,itj,kc,icptch,iflag
        real(kind=4):: cptsum,cdsum
        real(kind=4):: ar(imn),akc(imn),acpt(imn,imn),ad(imn),xadj(im,im) !! parameters

        !*************************
        kc=0; do i=1,imn                                       !!       
                kc=kc+1; ar(i)=par(kc)                         !!
              enddo                                            !! maximal growth rates
!             write(*,24),ar; pause 
              do i=1,imn                                       !!
                kc=kc+1; akc(i)=par(kc)                        !!
              enddo                                            !! carrying capacities
!            write(*,24) akc; pause  
              acpt=0.0                                         !! competition initialization
              do i=1,imn; do j=1,imn                           !! 
                kc=kc+1; acpt(i,j)=par(kc)                     !! competition terms
              enddo; enddo                                     !!
!              write(*,24) acpt; pause
              do i=1,imn                                       !!
                kc=kc+1; ad(i)=par(kc)                         !!
              enddo                                            !! dispersal terms
!           write(*,24),ad; pause
              do i=1,im
               do j=1,im
                kc=kc+1; xadj(i,j)=par(kc)
               enddo
              enddo 
        !*************************
        do iti=1,imn
        
           cptsum=0.0                                                   !!
           do itj=1,imn                                                 !! competition summation
            cptsum=cptsum+acpt(iti,itj)*Y(itj)                          !!
           enddo                                                        !!   

           cdsum=0.0                                                                            !!
           icptch=int(float(iti-1)/float(in))+1                                                 !!
           iflag=mod(iti,in); if(iflag==0) iflag=in                                             !!
           kc=0
           do itj=1,im                                                                          !! dispersal
            if(itj/=icptch) then
              kc=kc+1; cdsum=cdsum+ad(iti)*xadj(icptch,itj)*Y(iflag+(itj-1)*in)                 !!
            endif 
           !if(itj/=icptch) write(*,24) dfloat(iti),dfloat(icptch),dfloat(iflag+(itj-1)*in),ad(iti); pause
           enddo
            YDOT(iti)=ar(iti)*Y(iti)*(1.0-(Y(iti)+cptsum)/akc(iti))+(-ad(iti)*Y(iti)+1.d0*cdsum/float(im-1))
            
       enddo
       24  Format(30F15.7)
      RETURN
      END SUBROUTINE DERIVS
    END MODULE WANG

!******************************************************************

    PROGRAM DEMOWANG

      USE WANG
      USE random
!     Type declarations:
      IMPLICIT NONE
      
!!! REGRESSION ROUTINE
      EXTERNAL func
        
        INTERFACE
            REAL FUNCTION pxpand (order, xv, coeff)
                IMPLICIT NONE
                INTEGER,INTENT(IN)           :: order
                REAL,DIMENSION(1),INTENT(IN) :: coeff   ! was dim(order)
                REAL,INTENT(IN)              :: xv
            END FUNCTION pxpand
        END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=4):: h,T,TOUT, xrho1, xrho2, xvar1, xvar2, kppa
      real(kind=4):: ald,aldi,gmd,gmdi,btd, reg_avg
      real(kind=4):: alcv,gmcv,btcv
      real(kind=4):: sum_patch_l,svsum,spsum,tot_mean
      real(kind=4):: p_waxmn, p_mc
      integer(kind=4):: inic,IOUT,i,ii,j,jj,l,kcout,kij
      real(kind=4),allocatable:: Y(:),YDER(:),ctmp(:,:), noi_spec(:), noi_ptch(:), ypopl(:,:)
      real(kind=4),allocatable:: yavg(:,:),yprel(:,:),ysave(:,:,:),lo_simp(:), omg(:), sysame(:)
      real(kind=4),allocatable:: ymean(:), ysd(:), ycov(:,:)
      real(kind=4),allocatable:: adjm(:,:), beta(:), dist_mat(:,:)
      real(kind=4),allocatable:: simp_conc(:),xp_i(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! variables for the random number generator
    INTEGER, ALLOCATABLE  :: ix(:), seed(:) 

    REAL, ALLOCATABLE     :: xep(:), xptch(:), xnew(:), xmeanep(:), xmeanptch(:), xmeannew(:), &
                             covep(:), covptch(:), covnew(:), chol_fep(:), chol_fptch(:), chol_fnew(:), &
                             xseed(:)

    INTEGER               :: nep, nptch, nnew, ndf, n_binom, k, kep, kptch, knew, ier, pos1, pos2, which, iindx  

    REAL                  :: average, stdvn, shape, a, b, mu, p, pop_mean, pop_sd,   &
                            one = 1.0, zero = 0.0, sk, two = 2.0, pcntile(9), &
                            middle, xmax, xmin, xtseed
    LOGICAL               :: first

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR REGRESSION
         REAL, DIMENSION(20)     :: aR, wR
         REAL, DIMENSION(200,20) :: uR
         REAL, DIMENSION(20,20)  :: vR
         INTEGER                 :: jR
         REAL                    :: xxR, yyR, slopeaa,slopeab,slopebb,slopegg 
         INTEGER                 :: n1R, n2R, norderR
!         
!         !   raw data to be fit
         INTEGER, PARAMETER      :: nvalsR  = ninic
         REAL, DIMENSION(nvalsR)  :: xR !     = (/ 1.0, 2.0, 3.0, 4.0, 5.0 /)
         REAL, DIMENSION(nvalsR)  :: yR !     = (/ 1.0, 16.0, 80.9, 256.5, 625.0 /)
         real, dimension(nvalsR) :: aldR, btdR, gmdR, alcvR, btcvR, gmcvR

         norderR = 2; n1R=1; n2R=nvalsR     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   allocate(Y(NEQ),YDER(NEQ),ctmp(imn,imn),noi_spec(imn),noi_ptch(imn)); par=0.0
   allocate(yavg(im,in), yprel(im,in),lo_simp(im),sysame(in),ysave(im,in,niter))
   allocate(ymean(im),ysd(im),ycov(im,im))
   allocate(simp_conc(im),omg(im),xp_i(in))
   allocate(adjm(im,im),beta(im),dist_mat(im,im))    !!!! NETWORK PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RANDOM NUMBER GENERATOR PREPARATION
!     Set the random number seed.
 CALL RANDOM_SEED(size=k)
 ALLOCATE (seed(k),xseed(k))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FIXING THE SEED
!do i=1,k; call cpu_time(xtseed); xseed(i)=xtseed; enddo
 call random_number(xseed)
 seed=int(50000*xseed)
 CALL RANDOM_SEED(put=seed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!! STOCHASTIC PARAMETERS FROM PAPER
xrho1=0.0  !!print*, "Between species correlation";
!xrho2=0.0
read*,xrho2; print*, xrho2, "started" !! between patch correlation
xvar1=0.04  !! species specific variance
xvar2=0.04  !! patch specific variance
kppa=3.0    !! kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! allocation for the random number for species specific environmental response
    kep=in                    !! species/patch 
    nep = kep*(kep+1)/2
    ALLOCATE( xep(kep), xmeanep(kep), covep(nep), chol_fep(nep))

!!!!!! MEAN VECTOR
    xmeanep=0.0
! ! 
!!!!!! COVARIANCE MATRIX FOR SPECIES SPECIFIC ENVIRONMENTAL RESPONSE
        pos2 = 0
        DO i = 1, kep
         pos1 = pos2 + 1
         pos2 = pos2 + i
         covep(pos1:pos2)=xrho1*xvar1; covep(pos2)=xvar1
        END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! allocation for the random number for patch specific environmental response
    kptch=im
    nptch = kptch*(kptch+1)/2
    ALLOCATE ( xptch(kptch), xmeanptch(kptch), covptch(nptch), chol_fptch(nptch) )
!!!!!! MEAN VECTOR FOR PATCHES
    xmeanptch=0.0
! ! 
!!!!!! COVARIANCE MATRIX FOR BETWEEN PATCH ENVIRONMENTAL RESPONSE
    pos2 = 0
    DO i = 1, kptch
        pos1 = pos2 + 1
        pos2 = pos2 + i
        covptch(pos1:pos2)=xrho2*xvar2; covptch(pos2)=xvar2
    END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! allocation for the random number for combined noise response
!!! knew,nnew,allocate:: xnew,xmeannew,covnew,chol_fnew
    knew=im*in
    nnew = knew*(knew+1)/2
    ALLOCATE ( xnew(knew), xmeannew(knew), covnew(nnew), chol_fnew(nnew) )
!!!!!! MEAN VECTOR FOR PATCHES
    xmeannew=0.0
! ! 
!!!!!! COVARIANCE MATRIX FOR COMBINED RESPONSE
    pos2 = 0
    DO i = 1, knew
        pos1 = pos2 + 1
        pos2 = pos2 + i
        covnew(pos1:pos2)=xrho1*xrho2*xvar2; covnew(pos2)=xvar2
    END DO
    
!!!!!!!!!!!!!!!!! NETWORK STUFF

       do i=1,im
         beta(i)=0.5d0!0.3d0+0.4d0*dfloat(i)/dfloat(NEQ)  ! AREA PROXIES for Waxman
       enddo
       
     dist_mat=0.d0
      do i=1,im
       do j=i+1,im
        call random_number(dist_mat(j,i)); dist_mat(i,j)=dist_mat(j,i)  !! dist_mat is the distance matrix for the KHs
        enddo
      enddo   
      adjm=0.d0; kij=0; write(1,*)'digraph{'
        do i=1,im                                      ! 
          do j=1,im
            kij=kij+1                                          !
            !p_waxmn=(beta(i)*beta(j))*exp(-dist_mat(i,j)/xL)!(beta(i)*beta(j))*dexp(-dist_mat(i,j)/xL)  !
            !call random_number(p_mc)                         ! Adjacency matrix
            !if(p_mc<p_waxmn) adjm(i,j)=1.d0
            adjm(i,j)=1.d0
            if(i==j) adjm(i,j)=0.d0                           !
            if(adjm(i,j)==1) write(1,*) i,'->',j!
          enddo                                              ! 
         enddo
         write(1,*)'}'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RANDOM NUMBER GEN. PREPARATION DONE
      do inic=1,ninic!; print*,inic
      par=0.0
!     Set  the problem parameters:
      call random_number(par(1:imn)); 
      par(1:imn)=0.5+par(1:imn)*1.0              !! M*N MAXIMAL GROWTH RATES
      
      call random_number(par(imn+1:2*imn))
      par(imn+1:2*imn)=0.5+par(imn+1:2*imn)*1.0  !! M*N CARRYING CAPACITIES
      ctmp=0.0
          do i=1,im
           call random_number(ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in))
           ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in)=0.25+ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in)*0.5
          if (i>1) ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in)=ctmp(1:in,1:in) !! weird addition following the Shaopeng supplied code - ask?
          enddo
          kcout=0
          do i=1,imn
           do j=1,imn
            kcout=kcout+1
            if(i==j) ctmp(i,j)=0.0
            par(2*imn+kcout)=ctmp(i,j)            !! (M*N)^2 COMPETITION COEFFICIENTS FOR M PATCHES
           enddo
          enddo 
          !write(*,24)ctmp; pause
      do i=1,im
        call random_number(par(2*imn+imn**2+1+i*in-igen:2*imn+imn**2+i*in))                                         !! DISPERSAL RATES
       ! print*,par(2*imn+imn**2+1+i*in-igen:2*imn+imn**2+i*in)
        par(2*imn+imn**2+1+i*in-igen:2*imn+imn**2+i*in)=0.1+par(2*imn+imn**2+1+i*in-igen:2*imn+imn**2+i*in)*0.05         !!
        !print*,par(2*imn+imn**2+1+i*in-igen:2*imn+imn**2+i*in)
      enddo
      kcout=2*imn+imn**2+imn
      do i=1,im
        do j=1,im
         kcout=kcout+1; par(kcout)=adjm(i,j)
        enddo
      enddo
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !     Set the initial conditions:
       call random_number(Y); YDER=0.0
! !     Output points:
       T = 0.0
       h=0.01
! !     Perform the integration:
       DO IOUT = 1, ntrans
         call DERIVS(NEQ,T,Y,YDER)
         Y=Y+h*YDER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
        if(mod(IOUT,net_Tscl)==0) then  !! NETWORK RESET
         adjm=0.d0; kij=0 !! INITIALIZE ADJACENCY MATRIX
         do i=1,im                                                             ! 
          do j=1,im
            kij=kij+1                                                             !
            p_waxmn=(beta(i)*beta(j))*exp(-dist_mat(i,j)/xL)                   !
            call random_number(p_mc)                         ! Adjacency matrix
            if(p_mc<p_waxmn) adjm(i,j)=1.d0                   !
            !adjm(i,j)=1.d0
            if(i==j) adjm(i,j)=0.d0                           !
          enddo                                              ! 
         enddo
             kcout=2*imn+imn**2+imn
              do i=1,im
               do j=1,im
                kcout=kcout+1; par(kcout)=adjm(i,j)
               enddo
              enddo
       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       T=T+h
       END DO ! IOUT
       Y=Y+0.1
  !!! TO TAKE CARE OF ZERO POPULATIONS
       first = .true.
       
       DO IOUT = 1, niter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EVOLUTION BEGINS            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! between species random number
        CALL random_mvnorm(kep, xmeanep, covep, chol_fep, first, xep, ier)
        IF (ier .NE. 0) THEN
            WRITE(*, *) '** Covariance matrix1 is not +ve definite **'
            PAUSE; !EXIT replaced since the do loop was removed
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! between patch random number
        CALL random_mvnorm(kptch, xmeanptch, covptch, chol_fptch, first, xptch, ier)
        IF (ier .NE. 0) THEN
        WRITE(*, *) '** Covariance matrix2 is not +ve definite **'
        PAUSE; !EXIT replaced since the do loop was removed
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! between patch random number
!         CALL random_mvnorm(knew, xmeannew, covnew, chol_fnew, first, xnew, ier)
!         IF (ier .NE. 0) THEN
!         WRITE(*, *) '** Covariance matrix2 is not +ve definite **'
!         PAUSE; !EXIT replaced since the do loop was removed
!         END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,im
           noi_spec((i-1)*in+1:i*in)=xep
        enddo         !!!!!!!!!!!!!noise arrangement for integration
             
        kcout=0
        do i=1,im
          do j=1,in
           kcout=kcout+1
           noi_ptch(kcout)=xptch(i)
          enddo
        enddo
         call DERIVS(NEQ,T,Y,YDER)
         do i=1,imn
           Y(i)=Y(i)+h*YDER(i)+sqrt(h)*Y(i)*(noi_spec(i)+noi_ptch(i))! &
          ! & +kppa*xnew(i)) !!  noise added; T=T+h
         enddo
        kcout=0
        do l=1,im
         do i=1,in
           kcout=kcout+1
           ysave(l,i,IOUT)=Y(kcout)
          enddo
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
        if(mod(IOUT,net_Tscl)==0) then  !! NETWORK RESET
         adjm=0.d0; kij=0 !! INITIALIZE ADJACENCY MATRIX
         do i=1,im                                                             ! 
          do j=1,im
            kij=kij+1                                                             !
            p_waxmn=(beta(i)*beta(j))*exp(-dist_mat(i,j)/xL)                   !
            call random_number(p_mc)                         ! Adjacency matrix
            if(p_mc<p_waxmn) adjm(i,j)=1.d0                   !
            !adjm(i,j)=1.d0
            if(i==j) adjm(i,j)=0.d0                           !
          enddo                                              ! 
         enddo
             kcout=2*imn+imn**2+imn
              do i=1,im
               do j=1,im
                kcout=kcout+1; par(kcout)=adjm(i,j)
               enddo
              enddo
       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        T=T+h
       END DO ! IOUT
!!!!!!!!!!!!!!!!!!!!!!!! ENS based
   do l=1,im
     sum_patch_l=sum(ysave(l,:,niter))
    do i=1,in
     yprel(l,i)=ysave(l,i,niter)/sum_patch_l
    enddo
   enddo
 
  aldi=0.0
    do l=1,im
    simp_conc(l)=sum(yprel(l,:)**2)
    omg(l)=sum(ysave(l,:,niter))/sum(ysave(:,:,niter))
      aldi=aldi+omg(l)*simp_conc(l)
    enddo

    ald=1.0/aldi                              !! alpha_d
    !ald, btd, gmd, alcv, btcv, gmcv
    do i=1,in
     xp_i(i)=sum(ysave(:,i,niter))/sum(ysave(:,:,niter))
    enddo
    
    gmdi=sum(xp_i**2) !! gamma_d inverse
    gmd=1.0/gmdi                                !! gamma_d
    btd=gmd/ald                                  !! beta_d
 
! !!!!! VARIABILITY BASED MEASURES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      svsum=0.0
        do i=1,im        !!!! temporal mean for each patch
        spsum=0.0
         do IOUT=1,niter
          spsum=spsum+sum(ysave(i,:,IOUT))
         enddo
         ymean(i)=spsum/float(niter) !! temporal mean of community biomass in patch i
         svsum=svsum+ymean(i)
        enddo          !                                       QC done
        tot_mean=svsum!/float(niter)   ! temporal mean of metacommunity biomass
        reg_avg=tot_mean/float(im)    !regional average of temporal means of local community biomass
        
         ycov=0.d0
            do IOUT=1,niter
             do i=1,im
              do j=1,im
               ycov(i,j)=ycov(i,j)+(sum(ysave(i,:,IOUT))-ymean(i))*(sum(ysave(j,:,IOUT))-ymean(j))
              enddo
             enddo
            enddo
              ycov=ycov/float(niter)  !! covariance matrix between community biomasses
            do i=1,im
             ysd(i)=sqrt(ycov(i,i))
            enddo 
      !      print*,ycov; pause
        alcv=(sum(ysd)/tot_mean)**2                  !! alpha_cv  
        gmcv=sum(ycov)/(tot_mean)**2                 !! gamma_cv
        btcv=alcv/gmcv  !! is it correct?            !! beta_cv
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(11,24) ald, btd, gmd, alcv, btcv, gmcv 
        write(*,24) ald, btd, gmd, alcv, btcv, gmcv
        aldR(inic)=ald; btdR(inic)=btd; gmdR(inic)=gmd; alcvR(inic)=alcv; btcvR(inic)=btcv; gmcvR(inic)=gmcv
   enddo !!!! INIC
   
   !!!!! Regression start
   aR=0.d0; wR=0.d0; uR=0.d0; vR=0.d0;xr=0.d0; yR=0.d0

  xR=aldR; yR=alcvR
  CALL svdfit (xR(n1R), yR(n1R), n2R-n1R+1, aR, norderR, uR, vR, wR, nvalsR+5, 20, func) !! fit alcv and ald
  slopeaa=aR(2)
  
  yR=btcvR
  CALL svdfit (xR(n1R), yR(n1R), n2R-n1R+1, aR, norderR, uR, vR, wR, nvalsR+5, 20, func) !! fit btdcv and ald
  slopeab=aR(2)

  xR=btdR; yR=btcvR
  CALL svdfit (xR(n1R), yR(n1R), n2R-n1R+1, aR, norderR, uR, vR, wR, nvalsR+5, 20, func) !! fit btdcv and btd
  slopebb=aR(2)
  
  xR=gmdR; yR=gmcvR
  CALL svdfit (xR(n1R), yR(n1R), n2R-n1R+1, aR, norderR, uR, vR, wR, nvalsR+5, 20, func) !! fit gmcvcv and gmdd
  slopegg=aR(2)
  
  !!!!! Regression enddo
  write(21,24) xrho2,slopeaa,slopeab,slopebb,slopegg
!! OUTPUT FILE FORMAT
    
      24  Format(20F20.7)

END PROGRAM DEMOWANG


!!! POLYFIT

!     PROGRAM polyfit
!         IMPLICIT NONE
!         EXTERNAL func
!         
!         INTERFACE
!             REAL FUNCTION pxpand (order, xv, coeff)
!                 IMPLICIT NONE
!                 INTEGER,INTENT(IN)           :: order
!                 REAL,DIMENSION(1),INTENT(IN) :: coeff   ! was dim(order)
!                 REAL,INTENT(IN)              :: xv
!             END FUNCTION pxpand
!         END INTERFACE
!         
!         REAL, DIMENSION(20)     :: a, w
!         REAL, DIMENSION(200,20) :: u
!         REAL, DIMENSION(20,20)  :: v
!         INTEGER                 :: j
!         REAL                    :: xx, yy 
!         INTEGER                 :: n1, n2, norder
!         
!         !   raw data to be fit
!         INTEGER, PARAMETER      :: nvals  = 5
!         REAL, DIMENSION(nvals)  :: x      = (/ 1.0, 2.0, 3.0, 4.0, 5.0 /)1
!         REAL, DIMENSION(nvals)  :: y      = (/ 1.0, 16.0, 80.9, 256.5, 625.0 /) 
! 
!         !   polynomial order for fitting
!         norder = 4
!  
!         !   fit a polynomial of order norder to the 
!         !   [n1,n2] subset of arrays x() and y();
!         !   the array a() is returned with the best-fit coefficients:
!         !   y(x) = a(1) + a(2)*x + a(3)*x**2 +...+ a(n)*x**(n-1)
!         n1 = 1
!         n2 = nvals
!         CALL svdfit (x(n1), y(n1), n2-n1+1, a, norder,  &
!                      u, v, w, nvals+5, 20, func)
!         
!         OPEN (3, FILE='c:\tcfit\ircon\polyfit\POLYFIT.DAT', STATUS = 'APPEND')
!         WRITE (3,'(/5X,"Fit order is ",I2//)') norder
!         
!         DO j = 1, norder
!             WRITE (3,*) j-1, a(j)
!         END DO
! 
!         DO j = n1, n2
!             xx = x(j)
!             yy = pxpand (norder, xx, a)
!             WRITE (3,*) xx, y(j), yy, yy-y(j)
!         END DO
!         
!         CLOSE (3)
! 
!     END PROGRAM polyfit
! 

    SUBROUTINE svdfit (x, y, ndata, a, ma, u, v, w, mp, np, funcs)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                  :: ndata, ma, mp, np
       
        REAL, DIMENSION(ndata),  INTENT(IN)  :: x, y
        REAL, DIMENSION(ma),     INTENT(OUT) :: a
        REAL, DIMENSION(mp, np), INTENT(OUT) :: u
        REAL, DIMENSION(np, np), INTENT(OUT) :: v
        REAL, DIMENSION(np),     INTENT(OUT) :: w
        EXTERNAL funcs
        
        INTEGER, PARAMETER                   :: nmax = 1000, mmax = 50
        REAL, PARAMETER                      :: tol = 1.E-5
        
        REAL, DIMENSION(nmax)                :: b
        REAL, DIMENSION(mmax)                :: afunc
        INTEGER                              :: i, j
        REAL                                 :: wmax, thresh
        
        DO i = 1, ndata
            CALL funcs (x(i), afunc, ma)
            DO j = 1, ma
                u(i,j) = afunc(j)
            END DO                                  ! 11
            b(i) = y(i)
        END DO                                      ! 12
        
        CALL svdcmp (u, ndata, ma, mp, np, w, v)
        
        wmax = 0.
        DO j = 1, ma
            IF (w(j) > wmax) wmax = w(j)
        END DO
        
        thresh = tol*wmax
        
        DO j = 1, ma
            IF (w(j) < thresh) w(j) = 0.
        END DO
        
        CALL svbksb (u, w, v, ndata, ma, mp, np, b, a)
        
    END SUBROUTINE svdfit
        

    SUBROUTINE svdcmp (a, m, n, mp, np, w, v)
        IMPLICIT NONE
        INTEGER, PARAMETER                      :: nmax = 100
        INTEGER, INTENT(IN)                     :: m, n, mp, np
        REAL, DIMENSION(mp,np), INTENT(INOUT)   :: a
        REAL, DIMENSION(np),    INTENT(INOUT)   :: w
        REAL, DIMENSION(np,np), INTENT(INOUT)   :: v
        REAL, DIMENSION(nmax)                   :: rv1
        REAL                                    :: c, f, g, h, s, scale, anorm, x, y, z
        INTEGER                                 :: i, j, k, l, nm, its
        LOGICAL                                 :: block1
        
        g     = 0.
        scale = 0.
        anorm = 0.
        
        DO i = 1, n                 ! 25
            l = i+1
            rv1(i) = scale*g
            g      = 0.
            s      = 0.
            scale  = 0.
            IF (i <= m) THEN
                DO k = i, m
                    scale = scale + ABS(a(k,i))
                END DO                                      ! 11
                IF (scale /= 0.) THEN
                    
                    DO k = i, m
                        a(k,i) = a(k,i)/scale
                        s = s + a(k,i)*a(k,i)
                    END DO                                  ! 12
                    
                    f      = a(i,i)
                    g      = -SIGN(SQRT(s), f)
                    h      = f*g - s
                    a(i,i) = f - g
                    
                    IF (i /= n) THEN
                        DO j = l, n         ! 15
                            s = 0.
                            DO k = i, m
                                s = s + a(k,i)*a(k,j)
                            END DO                          ! 13
                            f = s/h
                            DO k = i, m
                                a(k,j) = a(k,j) + f*a(k,i)
                            END DO                          ! 14
                        END DO                              ! 15
                    END IF
                    
                    DO k = i, m
                        a(k,i) = scale*a(k,i)
                    END DO                                  ! 16
                END IF
            END IF
            
            w(i)  = scale*g
            g     = 0.
            s     = 0.
            scale = 0.
            
            IF (i <= m .AND. i /= n) THEN
                DO k = l, n
                    scale = scale + ABS(a(i,k))
                END DO                                      ! 17
                
                IF (scale /= 0.) THEN
                    DO k = l, n
                        a(i,k) = a(i,k)/scale
                        s = s + a(i,k)*a(i,k)
                    END DO                                  ! 18
                    
                    f      = a(i,l)
                    g      = -SIGN(SQRT(s), f)
                    h      = f*g - s
                    a(i,l) = f - g
                    
                    DO k = l, n
                        rv1(k) = a(i,k)/h
                    END DO                                  ! 19
                    
                    IF (i /= m) THEN
                        DO j = l, m
                            s = 0.
                            DO k = l, n
                                s = s + a(j,k)*a(i,k)
                            END DO                          ! 21
                            DO k = l, n
                                a(j,k) = a(j,k) + s*rv1(k)
                            END DO                          ! 22
                        END DO                              ! 23
                    END IF
                    
                    DO k = l, n
                        a(i,k) = scale*a(i,k)
                    END DO                                  ! 24
                END IF
            END IF
            
            anorm = AMAX1(anorm, (ABS(w(i)) + ABS(rv1(i)))) 
        END DO                                              ! 25
        
        DO i = n, 1, -1             ! 32
            IF (i < n) THEN
                IF (g /= 0.) THEN
                    DO j = l, n
                        v(j,i) = (a(i,j)/a(i,l))/g
                    END DO                                  ! 26
                    DO j = l, n
                        s = 0.
                        DO k = l, n
                            s = s + a(i,k)*v(k,j)
                        END DO                              ! 27
                        DO k = l, n
                            v(k,j) = v(k,j) + s*v(k,i)
                        END DO                              ! 28
                    END DO                                  ! 29
                END IF
                
                DO j = l, n
                    v(i,j) = 0.
                    v(j,i) = 0.
                END DO                                      ! 31
            END IF
            
            v(i,i) = 1.
            g      = rv1(i)
            l      = i
        END DO                                              ! 32
        
        DO i = n, 1, -1
            l = i + 1
            g = w(i)
            
            IF (i < n) THEN
                DO j = l, n
                    a(i,j) = 0.
                END DO                                      ! 33
            END IF
            
            IF (g /= 0.) THEN
                g = 1./g
                IF (i /= n) THEN
                    DO j = l, n
                        s = 0.
                        DO k = l, m
                            s = s + a(k,i)*a(k,j)
                        END DO                              ! 34
                        f = (s/a(i,i))*g
                        DO k = i, m
                            a(k,j) = a(k,j) + f*a(k,i)
                        END DO                              ! 35
                    END DO                                  ! 36
                END IF
                
                DO j = i, m
                    a(j,i) = a(j,i)*g
                END DO                                      ! 37
            
            ELSE
                DO j = i, m
                    a(j,i) = 0.
                END DO                                      ! 38
            END IF
            
            a(i,i) = a(i,i) + 1.
        END DO                                              ! 39
        
        loop49: DO k = n, 1, -1     ! 49
            DO its = 1, 30          ! 48
                
                block1 = .TRUE.
                DO l = k, 1, -1
                    nm = l - 1
                    IF      ((ABS(rv1(l)) + anorm) == anorm) THEN
                        block1 = .FALSE.
                        EXIT
                    ELSE IF ((ABS(w(nm)) + anorm)  == anorm) THEN
                        EXIT
                    END IF
                END DO                                      ! 41
                
                IF (block1) THEN
                    c = 0.                                  ! 1
                    s = 1.
                    DO i = l, k         ! 43
                        f = s*rv1(i)
                        IF ((ABS(f) + anorm) /= anorm) THEN
                            g    = w(i)
                            h    = SQRT(f*f + g*g)
                            w(i) = h
                            h    = 1./h
                            c    = g*h
                            s    = -(f*h)
                            DO j = 1, m
                                y       = a(j,nm)
                                z       = a(j,i)
                                a(j,nm) =   y*c  + z*s
                                a(j,i)  = -(y*s) + z*c
                            END DO                          ! 42
                        END IF
                    END DO                                  ! 43
                END IF
                
                !   block 2
                z = w(k)
                IF (l == k) THEN
                    IF (z < 0.) THEN
                        w(k) = -z
                        DO j = 1, n
                            v(j,k) = - v(j,k)
                        END DO                              ! 44
                    END IF
                    CYCLE loop49
                END IF
                
                IF (its == 30) PAUSE 'No convergence in 30 iterations'
                
                x  = w(l)
                nm = k - 1
                y  = w(nm)
                g  = rv1(nm)
                h  = rv1(k)
                f  = ((y-z)*(y+z) + (g-h)*(g+h)) / (2.*h*y)
                g  = SQRT(f*f + 1.)
                f  = ((x-z)*(x+z) + h*((y/(f + SIGN(g,f))) - h))/x
                c  = 1.
                s  = 1.

                DO j = l, nm    ! 47
                    i = j + 1
                    g = rv1(i)
                    y = w(i)
                    h = s*g
                    g = c*g
                    z = SQRT(f*f + h*h)
                    rv1(j) = z
                    c = f/z
                    s = h/z
                    f =   x*c  + g*s
                    g = -(x*s) + g*c
                    h = y*s
                    y = y*c
                    DO nm = 1, n        ! 45
                        x = v(nm,j)
                        z = v(nm,i)
                        v(nm,j) =   x*c  + z*s
                        v(nm,i) = -(x*s) + z*c
                    END DO                          ! 45
                    z = SQRT(f*f + h*h)
                    w(j) = z
                    IF (z /= 0.) THEN
                        z = 1./z
                        c = f*z
                        s = h*z
                    END IF
                    f =   c*g  + s*y
                    x = -(s*g) + c*y
                    DO nm = 1, m        ! 46
                        y = a(nm,j)
                        z = a(nm,i)
                        a(nm,j) =   y*c  + z*s
                        a(nm,i) = -(y*s) + z*c
                    END DO                              ! 46
                END DO                                  ! 47
                
                rv1(l) = 0.
                rv1(k) = f
                w(k)   = x
            END DO                                      ! 48
           
        END DO loop49                                   ! 49
    
    END SUBROUTINE svdcmp


    SUBROUTINE svbksb (u, w, v, m, n, mp, np, b, x)
        IMPLICIT NONE
        INTEGER,                INTENT(IN)   :: m, n, mp, np
        REAL, DIMENSION(mp,np), INTENT(IN)   :: u
        REAL, DIMENSION(np,np), INTENT(IN)   :: v
        REAL, DIMENSION(np)   , INTENT(IN)   :: w
        REAL, DIMENSION(np)   , INTENT(OUT)  :: x
        REAL, DIMENSION(mp)   , INTENT(IN)   :: b
        
        INTEGER, PARAMETER                   :: nmax = 100
        REAL, DIMENSION(nmax)                :: tmp
        INTEGER                              :: i, j, jj
        REAL                                 :: s
        
        DO j = 1, n
            s = 0.
            IF (w(j) /= 0.) THEN
                DO i = 1, m
                    s = s + u(i,j)*b(i)
                END DO                          ! 11
                s = s/w(j)
            END IF
            tmp(j) = s
        END DO                                  ! 12
            
        DO j = 1, n
            s = 0.
            DO jj = 1, n
                s = s + v(j,jj)*tmp(jj)
            END DO                              ! 13
            x(j) = s
        END DO                                  ! 14
 
    END SUBROUTINE svbksb
     
     
    SUBROUTINE func (x, afunc, ma)
        IMPLICIT NONE
        REAL, INTENT(IN)                    :: x
        INTEGER, INTENT(IN)                 :: ma
        REAL, DIMENSION(ma), INTENT(OUT)    :: afunc
        INTEGER                             :: j
        afunc(1) = 1.
        DO j = 2, ma
            afunc(j) = x*afunc(j-1)
        END DO
    END SUBROUTINE func
                  
        
    REAL FUNCTION pxpand (order, xv, coeff) RESULT (polyval)
        IMPLICIT NONE
        INTEGER,INTENT(IN)           :: order
        REAL,DIMENSION(1),INTENT(IN) :: coeff   ! was dim(order)
        REAL,INTENT(IN)              :: xv
        INTEGER                      :: j
        polyval = 0.
        DO j = order, 2, -1
            polyval = xv*(polyval + coeff(j))
        END DO
        polyval = polyval + coeff(1)
    END FUNCTION pxpand

