MODULE WANG   
   implicit none
   integer, parameter:: iid=2,id=iid*2 !! n is the number of species in each patch, m is the number of patches in the metacommunity
   integer, parameter:: npar=iid*(iid+1)+3
   real(kind=4):: par(npar)

    CONTAINS

      SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        REAL T, Y, YDOT
        DIMENSION Y(NEQ), YDOT(NEQ)
        integer(kind=8):: itr,it1,itr2,it2,ib,itrr,j,ii,jj,ticount
        real(kind=8):: csum,cdsum
        real(kind=8)::aa,cc,bb(iid),tadj(iid,iid),ep !!!parameters
        aa=par(1); cc=par(2)           !!aa,cc
        !*************************
        do ii=1,iid                   !!bb
         bb(ii)=par(2+ii)
        enddo
        !*************************
        ticount=0                      !!adj
        do ii=1,iid
         do jj=1,iid
          ticount=ticount+1
          tadj(ii,jj)=par(2+iid+ticount)
         enddo
       enddo
       !*************************
       ep=par(iid*(iid+1)+3)        !! epsilon
      !print*,par; pause
      do itr=1,id,2
       csum=0.d0
       it1=(itr+1)/2
       do itr2=1,id,2
         it2=(itr2+1)/2
         csum=csum+tadj(it1,it2)*(Y(2*it2-1)-Y(2*it1-1))
       enddo
       ib=(itr+1)/2
       YDOT(itr)=Y(itr)*(aa-Y(itr))*(Y(itr)-1.d0)-Y(itr+1)+ep*csum
       YDOT(itr+1)=bb(ib)*Y(itr)-cc*Y(itr+1)
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
      INTEGER, PARAMETER :: mlyap=0, NEQ=(mlyap+1)*id
      real  :: h,t,TOUT,conn_sum,xe1,xrho,xvar,xsd
      integer:: IOUT,icount,j,looparm,i,ntrans,niter,kc
      real,allocatable:: Y(:),YDER(:),adj(:,:),noi_ftzn(:)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! variables for the random number generator
      INTEGER, ALLOCATABLE  :: ix(:), seed(:) 
      REAL, ALLOCATABLE     :: x(:), xmean(:), cov(:), chol_f(:), xseed(:)
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

xrho=0.0; xvar=1E-8; xsd=1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! allocation for the random number for species specific environmental response
      kep=iid
      nep = kep*(kep+1)/2
       ALLOCATE( x(kep), xmean(kep), cov(nep), chol_f(nep))
!!!!!! MEAN VECTOR
      xmean=0.0
!!!!!! COVARIANCE MATRIX FOR SPECIES SPECIFIC ENVIRONMENTAL RESPONSE
         pos2 = 0
         DO i = 1, kep
          pos1 = pos2 + 1
          pos2 = pos2 + i
          cov(pos1:pos2)=xrho*xvar; cov(pos2)=xvar
         END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ntrans=10000000; niter=10000000; allocate(Y(NEQ),YDER(NEQ),adj(iid,iid),noi_ftzn(NEQ))
   par=0.0; adj=0.0; noi_ftzn=0.0
   
      adj(1,2)=1.d0; adj(2,1)=1.d0
!     Set  the problem parameters:
      par(1)=-0.025794d0 !!aa
      par(2)=0.02d0  !!cc
      par(3)=0.0065d0 !!bb1
      par(4)=0.0135d0 !!bb2     
     icount=0
       do i=1,iid
        conn_sum=0.d0
        do j=1,iid
          conn_sum=conn_sum+adj(i,j)
        enddo
        do j=1,iid
         icount=icount+1
         par(2+iid+icount)=adj(i,j)!/conn_sum  (already taken care of in this case)
        enddo
       enddo
      
      !do i=1,51; print*,i
      xe1=0.126 !+float(i-1)/50.0*(0.129-0.126)
!     Set the initial conditions:
      call random_number(Y); YDER=0.0
!     Output points:
      T = 0.0
      h=0.005
      !!!!!!!!!! loop for parameters
      par(npar)=xe1
!     Perform the integration:
       DO IOUT = 1, ntrans
         call DERIVS(NEQ,T,Y,YDER)
          Y=Y+h*YDER
          T=T+h
        !  write(11,24)T,Y
       ENDDO
!      
       first = .true.
      kc=0
      DO IOUT = 1, niter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! random number generation
          CALL random_mvnorm(kep, xmean, cov, chol_f, first, x, ier)
           IF (ier .NE. 0) THEN
            WRITE(*, *) '** Covariance matrix1 is not +ve definite **'
            PAUSE; !EXIT replaced since the do loop was removed
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call DERIVS(NEQ,T,Y,YDER)
        Y=Y+h*YDER; Y(1)=Y(1)+sqrt(h)*sqrt(xsd)*x(1); Y(3)=Y(3)+sqrt(h)*sqrt(xsd)*x(2)
          T=T+h; if((Y(1)+Y(3)/2)>0.4) kc=kc+1
        if (mod(IOUT,3)==0) write(11,24) T,Y
        if (mod(IOUT,3)==0) write(12,24) T,x
        first = .false.
      END DO ! IOUT
      print*,xe1,kc
!      write(10,24) xe1,float(kc)
!     enddo 
      
!!! OUTPUT FILE FORMAT
      24  Format(40F20.10)

END PROGRAM DEMOWANG
