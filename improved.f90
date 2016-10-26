MODULE test   
   implicit none
   integer, parameter:: in=10,im=2,igen=5 !! n is the number of species in each patch, m is the number of patches in the metacommunity
   integer, parameter:: imn=im*in, npar=imn*(imn+3)
   real(kind=4):: par(npar)

    CONTAINS

      SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        REAL T, Y, YDOT
        DIMENSION Y(NEQ), YDOT(NEQ)
        integer(kind=4)::i,j,iti,itj,kc,icptch,iflag
        real(kind=4):: cptsum,cdsum
        real(kind=4):: ar(imn),akc(imn),acpt(imn,imn),ad(imn) !! parameters

        !*************************
        kc=0; do i=1,imn                   !!       
                kc=kc+1; ar(i)=par(kc)     !!
              enddo                        !! maximal growth rates
           !print*,ar; pause   
              do i=1,imn                    !!
                kc=kc+1; akc(i)=par(kc)     !!
              enddo                         !! carrying capacities
           !print*,akc; pause   
              acpt=0.0                                   !! competition initialization
              do i=1,imn; do j=1,imn                     !!
                kc=kc+1; acpt(i,j)=par(kc)               !! 
                if(j==i) acpt(i,j)=0.0                   !! competition terms
                if((i-10)*(j-10)<0) acpt(i,j)=0.0        !!
                if((i-10)>0.and.(j-10)==0) acpt(i,j)=0.0 !!
              enddo; enddo                         !!
             ! do i=1,imn; do j=1,imn; print*,i,j,acpt(i,j); enddo; enddo; pause
              do i=1,imn                     !!
                kc=kc+1; ad(i)=par(kc)       !!
              enddo                          !! dispersal terms
              !print*,ad; pause
        !*************************

        do iti=1,imn
          
           cptsum=0.0                                                   !!
           do itj=1,imn                                                 !! competition
            cptsum=cptsum+acpt(iti,itj)*Y(itj)                          !!
           enddo                                                        !!   

           cdsum=0.0                                                           !!
           icptch=int(float(iti-1)/float(in))+1                                !!
           iflag=mod(iti,in); if(iflag==0) iflag=in                            !!  
           do itj=1,im                                                         !! dispersal
            if(itj/=icptch) cdsum=cdsum+ad(iti)*Y(iflag+(itj-1)*in)            !!
           ! if(itj/=icptch) print*,iti,itj,icptch,iflag,iflag+(itj-1)*in
           enddo

            YDOT(iti)=ar(iti)*Y(iti)*(1.0-(Y(iti)+cptsum)/akc(iti))+(-ad(iti)*Y(iti)+cdsum/float(im-1))
            
       enddo
       
      RETURN
      END SUBROUTINE DERIVS

    END MODULE WANG
!******************************************************************

    PROGRAM improve

      USE test
      USE random
!     Type declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: mlyap=0, NEQ=(mlyap+1)*imn
      real  :: h, T, TOUT, xrho1, xrho2, xvar1, xvar2, kppa
      real  :: spsum, alpd
      integer:: IOUT, i, j, ntrans, niter,kc
      real,allocatable:: Y(:),Yold(:), YDER(:), noi_spec(:),noi_ptch(:)
      real,allocatable:: loc_simp(:)
      
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

xrho2=0.4  !!print*, "Between patch correlation"; 
xrho1=-0.05 !!print*, "Between species correlation"; 
xvar1=0.04  !! species specific variance
xvar2=0.04  !! patch specific variance
kppa=8.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! allocation for the random number for species specific environmental response
      kep=imn-igen
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ntrans=40000; niter=10000; allocate(Y(NEQ),YDER(NEQ),noi_spec(imn),noi_ptch(imn),loc_simp(im)); par=0.0; noi_spec=0.0; noi_ptch=0.0

!     Set  the problem parameters:
      call random_number(par(1:imn)); 
      par(1:imn)=0.5+par(1:imn)*1.0                 !! M*N MAXIMAL GROWTH RATES
      do i=1,(im-1)
        pos1=i*in-igen+1; pos2=(i+1)*in-igen+1
        par(pos1:i*in)=par(pos2:(i+1)*in)
      enddo                                    !!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!
      
      call random_number(par(imn+1:2*imn))
      par(imn+1:2*imn)=0.5+par(imn+1:2*imn)*1.0  !! M*N CARRYING CAPACITIES
!       do i=1,(im-1)
!         pos1=(imn+1)+i*in-igen; pos2=(imn+1)+(i+1)*in-igen
!         par(pos1:imn+i*in)=par(pos2:imn+(i+1)*in)
!       enddo                                    !!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      call random_number(par(2*imn+1:2*imn+imn**2)); 
      par(2*imn+1:2*imn+imn**2)=0.25+par(2*imn+1:2*imn+imn**2)*0.5 !! M*N*N COMPETITION COEFFICIENTS FOR M PATCHES (NEED TO THINK ............!!!) **
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      call random_number(par(2*imn+imn**2+1+in-igen:2*imn+imn**2+in)) !! DISPERSAL RATES
      par(2*imn+imn**2+1+in-igen:2*imn+imn**2+in)=0.05+par(2*imn+imn**2+1+in-igen:2*imn+imn**2+in)*0.1
      do i=1,(im-1)
        pos1=(2*imn+imn**2+1)+i*in-igen; pos2=(2*imn+imn**2+1)+(i+1)*in-igen
        par(pos2:(2*imn+imn**2)+(i+1)*in)=par(pos1:(2*imn+imn**2)+i*in)
      enddo !; print*,par(2*imn+imn**2+1:2*imn+imn**2+imn)      ; pause
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Set the initial conditions:
      call random_number(Y); YDER=0.0
!     Output points:
      T = 0.0
      h=0.1
      
!     Perform the integration:
       DO IOUT = 1, ntrans
!        print*,Y; pause
         call DERIVS(NEQ,T,Y,YDER)
          Y=Y+h*YDER; T=T+h
       ENDDO
     !  do i=1,imn; print*, Y(i); enddo
       Y=Y+0.1
      first = .true.
      DO IOUT = 1, niter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! between species random number patch
          CALL random_mvnorm(kep, xmeanep, covep, chol_fep, first, xep, ier)
           IF (ier .NE. 0) THEN
            WRITE(*, *) '** Covariance matrix1 is not +ve definite **'
            PAUSE; !EXIT replaced since the do loop was removed
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! between patch random number
       ! print*, xmeanptch; pause; print*,covptch; pause; print*,chol_fptch; pause
        CALL random_mvnorm(kptch, xmeanptch, covptch, chol_fptch, first, xptch, ier)
        IF (ier .NE. 0) THEN
          WRITE(*, *) '** Covariance matrix2 is not +ve definite **'
          PAUSE; !EXIT replaced since the do loop was removed
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do i=1,imn-igen; noi_spec(i)=xep(i); enddo;         !!!!!!!!!!!!!noise arrangement for integration
         do i=1,(im-1)  !! making up the dummy noise terms
           noi_spec((i*in+1)+(in-igen):(i*in+1)+in)=noi_spec(((i-1)*in+1)+(in-igen):((i-1)*in+1)+in)
         enddo
         kc=0
         do i=1,im
          do j=1,in
           kc=kc+1
           noi_ptch(kc)=xptch(i)
          enddo
         enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
        
        call DERIVS(NEQ,T,Y,YDER)
        Y=Y+h*YDER+sqrt(h)*(noi_spec+noi_ptch)+kppa*h*noi_ptch*noi_spec !!  noise added
         do i=1,imn
          if (Y(i)<0.0) Y(i)=0.001
         enddo
!!!! local simpson indices         
         do i=1,im
         spsum=0.0
          for j=1,in
           spsum=spsum+(Y((i-1)*i+j)/sum(Y((i-1)*im+1:(i-1)*im+in)))**2
          enddo
          loc_simp(i)=1.0/spsum
         enddo 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        alpd=0.0
        do i=1,im
         alpd=alpd+sum(Y((i-1)*im+1:(i-1)*im+in))*loc_simp(i)
        enddo 
         T=T+h
!         write(11,24) T,Y(1:5)
!         !write(12,24) T, noi_spec,xptch
!         write(12,24) T,Y(6:10)
!         write(20,24) T,(Y(6:10)+Y(16:20))/2
!         write(13,24) T,Y(11:15)
!         write(14,24) T,Y(16:20)

        first = .false.
      END DO ! IOUT
      
!!! OUTPUT FILE FORMAT
      24  Format(40F20.10)

END PROGRAM improve
