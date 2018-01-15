   MODULE WANG    
   implicit none
   integer, parameter:: in=10,im=2,igen=5            !! 'in' is the number of species in each patch, 'im' is the number of patches in the metacommunity
   integer, parameter:: imn=im*in, npar=imn*(imn+3)  !! npar is the number of parameters
   INTEGER, PARAMETER :: mlyap=0,NEQ=(mlyap+1)*imn,niter=5000,ntrans=1  !! ntrans is transient time, niter is iteration beyond transients
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
!     Type declarations:
      IMPLICIT NONE
      real(kind=8):: h,t,TOUT
      integer(kind=8):: IOUT,i,j,k,pos1,pos2,kcout
      real(kind=8),allocatable:: Y(:),YDER(:),ctmp(:,:)

      allocate(Y(NEQ),YDER(NEQ),ctmp(imn,imn)); par=0.d0

!     Set  the problem parameters:
      call random_number(par(1:imn)); 
      par(1:imn)=0.5+par(1:imn)*1.0              !! M*N MAXIMAL GROWTH RATES
      
      call random_number(par(imn+1:2*imn))
      par(imn+1:2*imn)=0.5+par(imn+1:2*imn)*1.0  !! M*N CARRYING CAPACITIES

      ctmp=0.d0
          do i=1,im
           call random_number(ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in))
           ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in)=0.25d0+ctmp((i-1)*in+1:i*in,(i-1)*in+1:i*in)*0.5d0
           do j=(i-1)*in+1,i*in-1; do k=j+1,i*in  !!!!!!!
            ctmp(k,j)=ctmp(j,k)                   !!!!! REMOVE IF SYMMETRY IN COMPETITION PARAMETERS IS NOT REQUIRED
           enddo; enddo                           !!!!!!!
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
       T = 0.0D0
       h=0.010
!       
! !     Perform the integration:
       DO IOUT = 1, ntrans
         call DERIVS(NEQ,T,Y,YDER)
         Y=Y+h*YDER; T=T+h
       END DO ! IOUT
!       
!       do i=1,imn; print*, Y(i); enddo
       DO IOUT = 1, niter
         call DERIVS(NEQ,T,Y,YDER)
         Y=Y+h*YDER; T=T+h
         write(11,24) T,Y(1:10)
         write(12,24) T,Y(11:20)
       END DO ! IOUT
!!! OUTPUT FILE FORMAT
      24  Format(11F15.7)

END PROGRAM DEMOWANG
