 MODULE WANG    
   implicit none
   integer, parameter:: in=10,im=2,igen=5           !! 'in' is the number of species in each patch, 'im' is the number of patches in the metacommunity
   integer, parameter:: imn=im*in, npar=imn*(imn+3)+im**2,ifact=1  !! npar is the number of parameters, ifact gives number of community pairs
   INTEGER, PARAMETER :: mlyap=0,NEQ=(mlyap+1)*imn,niter=30000,ntrans=20000,ninic=2000,net_Tscl=1000000  !! ntrans is transient time, niter is iteration beyond transients
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
           !write(*,24),ad; pause
              do i=1,im
               do j=1,im
               kc=kc+1; xadj(i,j)=par(kc)
               enddo
              enddo 
            !  print*,xadj
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
           do itj=1,im                               !! dispersal
            if(itj/=icptch) then
              cdsum=cdsum+xadj(itj,icptch)*Y(iflag+(itj-1)*in)                 !!
            endif
           !if(itj/=icptch) write(*,24) dfloat(iti),dfloat(icptch),dfloat(iflag+(itj-1)*in),ad(iti); pause
           enddo
           
           if (sum(xadj(:,icptch))>0.00001d0) then 
             cdsum=cdsum/sum(xadj(:,icptch))
            else
             cdsum=0.d0
            endif
            
            YDOT(iti)=ar(iti)*Y(iti)*(1.0-(Y(iti)+cptsum)/akc(iti))+ad(iti)*(cdsum-Y(iti))
            
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

    REAL, ALLOCATABLE     :: xep(:), xptch(:), xmeanep(:), xmeanptch(:), covep(:), covptch(:), chol_fep(:), &
                            chol_fptch(:), xseed(:)

    INTEGER               :: nep, nptch, ndf, n_binom, k, kep, kptch, ier, pos1, pos2, which, iindx

    REAL                  :: average, stdvn, shape, a, b, mu, p, pop_mean, pop_sd,   &
                            one = 1.0, zero = 0.0, sk, two = 2.0, pcntile(9), &
                            middle, xmax, xmin, xtseed
    LOGICAL               :: first

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
kppa=8.0    !! kappa

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
      adjm=0.d0; kij=0; write(1,*) 'digraph {'
        do i=1,im                                      ! 
          do j=1,im
            kij=kij+1                                          !
            p_waxmn=(beta(i)*beta(j))*exp(-dist_mat(i,j)/xL)!(beta(i)*beta(j))*dexp(-dist_mat(i,j)/xL)  !
            call random_number(p_mc)                         ! Adjacency matrix
            !if(p_mc<p_waxmn) adjm(i,j)=1.d0
            adjm(i,j)=1.d0
            if(i==j) adjm(i,j)=0.d0                           !
            if(adjm(i,j)==1) write(1,*) i,'->',j!
          enddo                                              ! 
         enddo
         write(1,*) '}'
         
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
            !if(p_mc<p_waxmn) adjm(i,j)=1.d0                   !
            adjm(i,j)=1.d0
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
           Y(i)=Y(i)+h*YDER(i)+sqrt(h)*Y(i)*(noi_spec(i)+noi_ptch(i) &
           & +kppa*noi_ptch(i)*noi_spec(i)) !!  noise added; T=T+h
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
            !if(p_mc<p_waxmn) adjm(i,j)=1.d0                   !
            adjm(i,j)=1.d0
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
   enddo !!!! INIC     
!! OUTPUT FILE FORMAT
      24  Format(20F20.7)

END PROGRAM DEMOWANG
