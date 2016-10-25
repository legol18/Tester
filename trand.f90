   MODULE TESTER 
   implicit none
   integer, parameter:: in=5,im=2 !! 'in' is the number of species in each patch, 'im' is the number of patches in the metacommunity
   integer, parameter:: imn=in*(im+1), npar=imn*(in+3)
   real(kind=4):: par(npar)

    CONTAINS

      SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        REAL T, Y, YDOT
        DIMENSION Y(NEQ), YDOT(NEQ)
        integer(kind=4)::i,j,k,iti,iti1,itj,itj1,kc,itc,iflag,idix
        real(kind=4):: cptsum,cdsum
        real(kind=4):: ar(im,in),akc(im,in),acpt(im,in,in),ad(im,in) !! parameters

        !*************************
        kc=0; do i=1,im; do j=1,in         !!       
                kc=kc+1; ar(i,j)=par(kc) !!
              enddo; enddo               !! maximal growth rates
              
              do i=1,im; do j=1,in         !!
                kc=kc+1; akc(i,j)=par(kc)!!
              enddo; enddo               !! carrying capacities

              acpt=0.0                      !! competition initialization
              do i=1,im; do j=1,in; do k=1,in   !!
                kc=kc+1; acpt(i,j,k)=par(kc) !! 
                if(j==k) acpt(i,j,k)=0.0    !! competition terms
              enddo; enddo; enddo            !!
              
              do i=1,im; do j=1,in             !!
                kc=kc+1; ad(i,j)=par(kc)     !!
              enddo; enddo                   !! dispersal terms
        !*************************
        kc=0
        
        do iti=1,im
         do itj=1,in

          kc=kc+1
          
          iflag=mod(kc,in); if(iflag==0) iflag=in; cptsum=0.0           !!
           do itj1=1,in                                                 !! competition 
            if(itj1/=iflag) then                                       !! terms 
           ! print*,kc,itj1+(iti-1)*n                                  !! 
            cptsum=cptsum+acpt(iti,itj,itj1)*Y(itj1+(iti-1)*in)         !!
            endif
           enddo                                                       !!   

           cdsum=0.0                                                  !! 
           do iti1=1,im                                                 !! dispersal
            if(iti1/=iti) then                                         !! terms
           !  print*,iti1,iflag,kc,(iti1-1)*n+iflag                    !!
             cdsum=cdsum+ad(iti1,iflag)*Y((iti1-1)*in+iflag)            !!
            endif                                                      !!  
           enddo
           
          YDOT(kc)=ar(iti,itj)*Y(kc)*(1.0-(Y(kc)+cptsum)/akc(iti,itj))+(-ad(iti,itj)*Y(kc)+cdsum/float(im-1))

        enddo
       enddo 
       !pause
       
      RETURN
      END SUBROUTINE DERIVS

    END MODULE TESTER

!******************************************************************

    PROGRAM DEMOWANG

      USE TESTER
      USE random
!     Type declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: mlyap=0, NEQ=(mlyap+1)*imn
      real  :: h, T, TOUT, xrho1, xrho2, xvar1, xvar2, kppa
      integer:: IOUT, i, j, ntrans, niter
      real,allocatable:: Y(:), YDER(:), noi_spec(:)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! variables for the random number generator
      INTEGER, ALLOCATABLE  :: ix(:), seed(:) 
      REAL, ALLOCATABLE     :: xep(:,:), xptch(:), xmeanep(:,:), xmeanptch(:), covep(:,:), covptch(:), chol_fep(:,:), &
                               chol_fptch(:), xseed(:), xcov(:,:), xtemp(:), xmeantemp(:), xcovtemp(:), xchol(:)
      INTEGER               :: nep, nptch, ndf, n_binom, k, kep, kptch, ier, pos1, pos2, which, iindx
      REAL                  :: average, stdvn, shape, a, b, mu, p, pop_mean, pop_sd,   &
                               one = 1.0, zero = 0.0, sk, two = 2.0, pcntile(9), &
                               middle, xmax, xmin, xtseed
      LOGICAL               :: first
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Set the random number seed.

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k),xseed(k))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FIXING THE SEED
do i=1,k; call cpu_time(xtseed); xseed(i)=xtseed; enddo
!call random_number(xseed)
seed=int(5000000*xseed)
CALL RANDOM_SEED(put=seed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

xrho2=-0.4  !!print*, "Between patch correlation"; 
xrho1=-0.05 !!print*, "Between species correlation"; 
xvar1=0.04  !! species specific variance
xvar2=0.04  !! patch specific variance
kppa=8.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! allocation for the random number for species specific environmental response
      kep=in
      nep = kep*(kep+1)/2
       ALLOCATE( xep(im,kep), xmeanep(im,kep), covep(im,nep), chol_fep(im,nep), xtemp(kep), xmeantemp(kep), & 
                 xcovtemp(nep), xchol(nep) )
!!!!!! MEAN VECTOR
      xmeanep=0.0; xmeantemp=0.0
! ! 
!!!!!! COVARIANCE MATRIX FOR SPECIES SPECIFIC ENVIRONMENTAL RESPONSE
      DO j=1,im
         pos2 = 0
         DO i = 1, kep
          pos1 = pos2 + 1
          pos2 = pos2 + i
          covep(j,pos1:pos2)=xrho1*xvar1; covep(j,pos2)=xvar1
         END DO
      ENDDO
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ntrans=500000; niter=10000; allocate(Y(NEQ),YDER(NEQ),noi_spec(imn)); par=0.0

!     Set  the problem parameters:
      call random_number(par(1:imn)); par(1:imn)=0.5+par(1:imn)*1.0  !! M*N MAXIMAL GROWTH RATES
      call random_number(par(imn+1:2*imn)); par(imn+1:2*imn)=0.5+par(imn+1:2*imn)*1.0  !! M*N CARRYING CAPACITIES
      call random_number(par(2*imn+1:2*imn+in*in*im)); par(2*imn+1:2*imn+in*in*im)=0.25+par(2*imn+1:2*imn+in*in*im)*0.50 !! M*N*N COMPETITION COEFFICIENTS FOR M PATCHES
      call random_number(par(2*imn+in*in*im+1:2*imn+in*in*im+imn))                                     !! DISPERSAL RATES 
      par(2*imn+in*in*im+1:2*imn+in*in*im+imn)=0.050+par(2*imn+in*in*im+1:2*imn+in*in*im+imn)*0.10 !! WITHOUT SWITCHING OFF THE VALUES
      par(2*imn+in*in*im+1:2*imn+in*in*im+imn:2)=0.0                    !! EVERY ALTERNATE TERM SWITCHED OFF GIVING HALF GENERALISTS (SPECIALISTS) IN EACH PATCH 
do i=1,npar; print*,par(i); enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Set the initial conditions:
      call random_number(Y); YDER=0.0
!     Output points:
      T = 0.0
      h=0.01
      
!     Perform the integration:
       DO IOUT = 1, ntrans
!        print*,Y; pause
         call DERIVS(NEQ,T,Y,YDER)
          Y=Y+h*YDER; T=T+h
       ENDDO
       do i=1,imn; print*, Y(i); enddo
      first = .true.
      DO IOUT = 1, niter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! between species random number patch
        do i=1,im; xtemp=0.0
          !xmeantemp = xmeanep(i,:); xcovtemp = covep(i,:); xchol = chol_fep(i,:)!; print*,xmeantemp; pause; print*,xcovtemp; pause; print*,xchol; pause
          !CALL random_mvnorm(kep, xmeantemp, xcovtemp, xchol, first, xtemp, ier)
          CALL random_mvnorm(kep, xmeanep(i,:), covep(i,:), chol_fep(i,:), first, xep(i,:), ier)
           IF (ier .NE. 0) THEN
            WRITE(*, *) '** Covariance matrix1 is not +ve definite **'
            PAUSE; !EXIT replaced since the do loop was removed
           END IF
          !xep(i,:)=xtemp
        enddo
        !print*,IOUT,xep; pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! between patch random number
       ! print*, xmeanptch; pause; print*,covptch; pause; print*,chol_fptch; pause
        CALL random_mvnorm(kptch, xmeanptch, covptch, chol_fptch, first, xptch, ier)
        IF (ier .NE. 0) THEN
          WRITE(*, *) '** Covariance matrix2 is not +ve definite **'
          PAUSE; !EXIT replaced since the do loop was removed
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i=1,im
         pos1=(i-1)*in+1; pos2=i*in
         noi_spec(pos1:pos2)=xep(i,:)
        enddo
!!        noi_spec=abs(noi_spec); xptch=abs(xptch)   !!! DO WE NEED THIS STEP?
        
        call DERIVS(NEQ,T,Y,YDER)
        
        Y=Y+h*YDER+sqrt(h)*noi_spec !! between species noise added

        do i=1,im
         pos1=(i-1)*in+1; pos2=i*in
         Y(pos1:pos2)=Y(pos1:pos2)+Y(pos1:pos2)*sqrt(h)*xptch(i) !! between patch noise added
         Y(pos1:pos2)=Y(pos1:pos2)+Y(pos1:pos2)*h*kppa*xptch(i)*noi_spec(pos1:pos2) !! multiplicative response added
        enddo         
         T=T+h
        write(11,24) T,Y(1:4)
        !write(12,24) T, noi_spec,xptch
        write(12,24) T,Y(6:10)
        write(13,24) T,Y(11:15)
        write(14,24) T,Y(16:20)

        first = .false.
      END DO ! IOUT
      !first = .false.
      
!!! OUTPUT FILE FORMAT
      24  Format(40F20.10)

END PROGRAM DEMOWANG

!!!! CODE CURRENTLY CRASHES DUE TO IMPROPER RESTRICTIONS ON POPULATION VALUES
