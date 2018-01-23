MODULE WANG    
   implicit none
   integer, parameter:: in=10,im=2,igen=5            !! 'in' is the number of species in each patch, 'im' is the number of patches in the metacommunity
   integer, parameter:: imn=im*in, npar=imn*(imn+3)  !! npar is the number of parameters
   INTEGER, PARAMETER :: mlyap=0,NEQ=(mlyap+1)*imn,niter=400000,ntrans=400000  !! ntrans is transient time, niter is iteration beyond transients
   real(kind=8):: par(npar)

    CONTAINS

      SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER, intent(in):: NEQ
        real(kind=8), intent(in):: T, Y(NEQ) 
        real(kind=8), intent(out):: YDOT(NEQ)
        integer(kind=8)::i,j,iti,itj,kc,icptch,iflag
        real(kind=8):: cptsum,cdsum
        real(kind=8):: ar(imn),akc(imn),acpt(imn,imn),ad(imn) !! parameters

        !*************************
        kc=0; do i=1,imn                                       !!       
                kc=kc+1; ar(i)=par(kc)                         !!
              enddo                                            !! maximal growth rates
             !write(*,24),ar; pause 
              do i=1,imn                                       !!
                kc=kc+1; akc(i)=par(kc)                        !!
              enddo                                            !! carrying capacities

              acpt=0.0                                         !! competition initialization
              do i=1,imn; do j=1,imn                           !! 
                kc=kc+1; acpt(i,j)=par(kc)                     !! competition terms
              enddo; enddo                                     !!

              do i=1,imn                                       !!
                kc=kc+1; ad(i)=par(kc)                         !!
              enddo                                            !! dispersal terms
           ! write(*,24),ad; pause

        !*************************
        do iti=1,imn

           cptsum=0.0                                                   !!
           do itj=1,imn                                                 !! competition summation
            cptsum=cptsum+acpt(iti,itj)*Y(itj)                          !!
           enddo                                                        !!   

           cdsum=0.0                                                           !!
           icptch=int(float(iti-1)/float(in))+1                                !!
           iflag=mod(iti,in); if(iflag==0) iflag=in                            !!  
           do itj=1,im                                                         !! dispersal
            if(itj/=icptch) cdsum=cdsum+ad(iti)*Y(iflag+(itj-1)*in)            !!
          ! if(itj/=icptch) print*,iti,icptch,iflag+(itj-1)*in; pause
           enddo
            YDOT(iti)=ar(iti)*Y(iti)*(1.0-(Y(iti)+cptsum)/akc(iti))+(-ad(iti)*Y(iti)+cdsum/float(im-1))
            
       enddo
      RETURN
      END SUBROUTINE DERIVS

    END MODULE WANG

!******************************************************************

    PROGRAM DEMOWANG

      USE WANG
      USE random
!     Type declarations:
      IMPLICIT NONE
      real(kind=8):: h,T,TOUT, xrho1, xrho2, xvar1, xvar2, kppa
      integer(kind=8):: IOUT,i,j,l,kcout
      real(kind=8),allocatable:: Y(:),YDER(:),ctmp(:,:), noi_spec(:), noi_ptch(:)
      
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

   allocate(Y(NEQ),YDER(NEQ),ctmp(imn,imn),noi_spec(imn),noi_ptch(imn)); par=0.d0
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
!xrho2=0.0  !!print*, "Between patch correlation"; 
read*,xrho2; print*, xrho2, "started"
xrho1=0.0   !!print*, "Between species correlation";
xvar1=0.04  !! species specific variance
xvar2=0.04  !! patch specific variance
kppa=8.0    !! kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! allocation for the random number for species specific environmental response
    kep=imn                    !! total species
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RANDOM NUMBER GEN. PREPARATION DONE
!     Set  the problem parameters:
      call random_number(par(1:imn)); 
      par(1:imn)=0.5+par(1:imn)*1.0              !! M*N MAXIMAL GROWTH RATES
      
      call random_number(par(imn+1:2*imn))
      par(imn+1:2*imn)=0.5+par(imn+1:2*imn)*1.0  !! M*N CARRYING CAPACITIES

      ctmp=0.d0
          do i=1,im
           call random_number(ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in))
           ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in)=0.25d0+ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in)*0.5d0
!            do j=(i-1)*in+1,i*in-1; do k=j+1,i*in  !!!!!!!
!             ctmp(k,j)=ctmp(j,k)                   !!!!! REMOVE IF SYMMETRY IN COMPETITION PARAMETERS IS NOT REQUIRED
!            enddo; enddo                           !!!!!!!
          enddo
          kcout=0
          do i=1,imn
           do j=1,imn
            kcout=kcout+1
            if(i==j) ctmp(i,j)=0.d0
            par(2*imn+kcout)=ctmp(i,j)            !! (M*N)^2 COMPETITION COEFFICIENTS FOR M PATCHES
           enddo
          enddo 
          
      do i=1,im
        call random_number(par(2*imn+imn**2+1+i*in-igen:2*imn+imn**2+i*in))                                         !! DISPERSAL RATES
        par(2*imn+imn**2+1+i*in-igen:2*imn+imn**2+i*in)=0.05+par(2*imn+imn**2+1+i*in:2*imn+imn**2+i*in)*0.1         !!
      enddo
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !     Set the initial conditions:
       call random_number(Y); YDER=0.0
! !     Output points:
       T = 0.0d0
       h=0.01d0
! !     Perform the integration:
       DO IOUT = 1, ntrans
         call DERIVS(NEQ,T,Y,YDER)
         Y=Y+h*YDER; T=T+h
       END DO ! IOUT
       
       write(*,24) Y(1:10) !! check output after transients
       
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
        do i=1,imn
           noi_spec(i)=xep(i)
        enddo         !!!!!!!!!!!!!noise arrangement for integration
        kcout=0
        do i=1,im
          do j=1,in
           kcout=kcout+1
           noi_ptch(kcout)=xptch(i)
          enddo
        enddo

         call DERIVS(NEQ,T,Y,YDER)
         Y=Y+h*YDER+sqrt(h)*Y*(noi_spec+noi_ptch+kppa*noi_ptch*noi_spec) !!  noise added; T=T+h
         write(11,24) T,Y(1:10)
         write(12,24) T,Y(11:20)
         T=T+h
         !write(21,24)noi_ptch
         !write(22,24)noi_spec
       END DO ! IOUT
       write(*,24) Y(1:10) !! check output
!! OUTPUT FILE FORMAT
      24  Format(11F15.7)

END PROGRAM DEMOWANG
