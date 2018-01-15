   MODULE WANG    
   implicit none
   integer, parameter:: in=10,im=2,igen=5            !! 'in' is the number of species in each patch, 'im' is the number of patches in the metacommunity
   integer, parameter:: imn=im*in, npar=imn*(imn+3)  !! npar is the number of parameters
   INTEGER, PARAMETER :: mlyap=0,NEQ=(mlyap+1)*imn,niter=5000,ntrans=50000  !! ntrans is transient time, niter is iteration beyond transients
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
              !24  Format(10F10.7)
      RETURN
      END SUBROUTINE DERIVS

!      SUBROUTINE JACD(NEQ,T,Y,ML,MU,PD,NROWPD)
!     Subroutine to define the exact Jacobian for this problem
!        IMPLICIT NONE
!        INTEGER NEQ, ML, MU, NROWPD
!        DOUBLE PRECISION T, Y, PD
!        DIMENSION Y(NEQ), PD(NROWPD,NEQ)
!        PD=0.d0
        !N = Y(1)
        !PHI = Y(2)
        !PD(1,1) = -(ALPHA*PHI+BETA)
        !PD(2,1) = PHI*RHO + TAU
        !PD(1,2) = -N*ALPHA
        !PD(2,2) = RHO*N - SIGMA
!        RETURN
!      END SUBROUTINE JACD

    END MODULE WANG

!******************************************************************

    PROGRAM DEMOWANG

      USE WANG
      USE DVODE_F90_M

!     Type declarations:
      IMPLICIT NONE
      real(kind=8):: h,t,TOUT
      integer(kind=8):: i,j,k,pos1,pos2,kcout
      real(kind=8),allocatable:: Y(:),ctmp(:,:)

!!!!!!!!!!!!!! INTEGRATOR ENTRIES; USUALLY NOT TO BE TOUCHED !!!!!!     
      INTEGER ITASK, ISTATE, ISTATS, IOUT, SWITCH
      DOUBLE PRECISION RSTATS, RTOL, ATOL
      DIMENSION RSTATS(22), ISTATS(31)
      TYPE (VODE_OPTS) :: OPTIONS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(Y(NEQ),ctmp(imn,imn)); par=0.d0

!     Set  the problem parameters:
      call random_number(par(1:imn)); 
      par(1:imn)=0.5+par(1:imn)*1.0              !! M*N MAXIMAL GROWTH RATES
!      do i=1,(im-1)
!        pos1=i*in-igen+1; pos2=(i+1)*in-igen+1   !!!!! NOT SURE IF NEEDED
!        par(pos1:i*in)=par(pos2:(i+1)*in)
!      enddo                                    !!!!!!!!!!!!!!!!!!
      
      call random_number(par(imn+1:2*imn))
      par(imn+1:2*imn)=0.5+par(imn+1:2*imn)*1.0  !! M*N CARRYING CAPACITIES
!      do i=1,(im-1)
!        pos1=(imn+1)+i*in-igen; pos2=(imn+1)+(i+1)*in-igen   !!!!! NOT SURE IF NEEDED
!        par(pos1:imn+i*in)=par(pos2:imn+(i+1)*in)
!      enddo                                    !!!!!!!!!!!!!!!!!!

      ctmp=0.d0
          do i=1,im
           call random_number(ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in))
           ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in)=0.25d0+ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in)*0.5d0
           do j=(i-1)*in+1,i*in-1; do k=j+1,i*in  !!
            ctmp(k,j)=ctmp(j,k)                   !! REMOVE IF SYMMETRY IN COMPETITION PARAMETERS IS NOT REQUIRED
           enddo; enddo                           !!
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
        par(2*imn+imn**2+1+i*in-igen:2*imn+imn**2+i*in)=0.05+par(2*imn+imn**2+1+i*in:2*imn+imn**2+i*in)*0.1  !!
      enddo
!       do i=1,(im-1)
!         pos1=(2*imn+imn**2+1)+i*in-igen; pos2=(2*imn+imn**2+1)+(i+1)*in-igen    !! NOT SURE IF NEEDED
!         par(pos2:(2*imn+imn**2)+(i+1)*in)=par(pos1:(2*imn+imn**2)+i*in)
!       enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!! Set the integration parameters:
      RTOL = 1.0D-6
      ATOL = 1.0D-9
      IOUT = 1
      ITASK = 1
      ISTATE = 1
     ! WRITE (6,*) RTOL, ATOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Set the initial conditions:
      call random_number(Y); Y=2.d0*Y
!     Output points:
      t = 0.0D0
      TOUT=0.d0
      h=0.1d0
!     Start with nonstiff method:
      OPTIONS = SET_OPTS(RELERR=RTOL,ABSERR=ATOL)

!     Perform the integration:
      DO IOUT = 1, ntrans
        TOUT = TOUT + h
        CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)

!!!!! ERROR ANALYSIS: DO NOT TOUCH    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Gather and write the integration statistics for this problem:
        CALL GET_STATS(RSTATS,ISTATS)
!        Stop the integration if an error occurred:
        IF (ISTATE<0) THEN
          WRITE (6,24) ISTATE
          STOP
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Print the solution and write the plot data file:
      END DO ! IOUT
!       
      DO IOUT = 1, niter
        TOUT = TOUT + h
        CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)

        write(11,24) t,Y(1:10)
        write(12,24) t,Y(11:20)
       ! write(13,24) t,Y(11:15)
       ! write(14,24) t,Y(16:20)
        
!!!!! ERROR ANALYSIS: DO NOT TOUCH    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Gather and write the integration statistics for this problem:
        CALL GET_STATS(RSTATS,ISTATS)
!        Stop the integration if an error occurred:
        IF (ISTATE<0) THEN
          WRITE (6,24) ISTATE
          STOP
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Print the solution and write the plot data file: ......................
      END DO ! IOUT
       24  Format(11F15.7)
!     Write the integration final integration statistics:
      WRITE (6,*) ISTATS(11), ISTATS(12), ISTATS(13)
!     Format statements for this problem:

END PROGRAM DEMOWANG
