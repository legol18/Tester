 MODULE information
   integer:: nptch
   integer:: NEQ   !! Number of patches
   integer, parameter:: niter=20000,ntrans=50000, net_Tscl=1000+ntrans, id=2,ninic=50
   real(kind=8), parameter:: d_max=14000.d0, a_max=1000.d0, h=0.1d0 
   real(kind=8):: par(20000), xL
   
   CONTAINS
   
   SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER, intent(in):: NEQ
        real(kind=8), intent(in):: T, Y(NEQ) 
        real(kind=8), intent(out):: YDOT(NEQ)
        INTEGER(kind=8):: i,j,xk,idx1,idx2,idx3
        real(kind=8):: xsum,xsum_out,xsum_in,xeff,xcin,xcout,xout,xlmbda,xm
        real(kind=8):: xrm(nptch),Kcar(nptch),xci(nptch),xdi(nptch),xbi(nptch),xhi(nptch),xadj(nptch,nptch),xweight(nptch)
        
        !print*,'sub'
        xrm=par(1:nptch)               ! max growth rate
        Kcar=par(nptch+1:2*nptch)      ! carrying capacity
        xci=par(2*nptch+1:3*nptch)     ! C parameter
        xdi=par(3*nptch+1:4*nptch)     ! death parameter
        xbi=par(4*nptch+1:5*nptch)     ! conversion efficiency parameter
        xhi=par(5*nptch+1:6*nptch)     ! handling parameter
        xlmbda=par(6*nptch+1)          ! lambda 
        xm=par(6*nptch+2)              ! movement parameter
        
        xk=0
        do i=1,nptch
         do j=1,nptch
          xk=xk+1
          xadj(i,j)=par(6*nptch+2+xk)     ! adjacency matrix elements
         enddo
        enddo
        
        
        !print*,xadj; pause
        do i=1,nptch
!        xsum=0.d0; xcin=0.d0; xcout=0.d0
        !print*,xcin,sum(xadj(:,i)),xcout,sum(xadj(i,:)); pause
        !print*, sum(xadj(:,i)), sum(xadj(i,:)), dfloat(NEQ); pause
!         if (sum(xadj(:,i))>0.00001d0) then !! OUTGOING
!                 xsum=xsum/sum(xadj(:,i))
!         else
!                 xsum=0.d0
!         endif
         do j=1,nptch  !! calculation of weights for all patches
         idx2=2*j; idx1=idx2-1
          xweight(j)=xci(j)*xbi(j)*Y(idx1)/(1.d0+xci(j)*xhi(j)*Y(idx1))-xdi(j)
         enddo
         
         xsum_out=0.d0 !!! outgoing loss from source patch
         do j=1,nptch
          xsum_out=xsum_out+xadj(j,i)*dexp(xlmbda*(xweight(j)-xweight(i)))
         enddo 

         xsum_in=0.d0  !!! incoming contribution from other patches
         do j=1,nptch
          idx2=2*j; idx1=idx2-1
          xsum_in=xsum_in+Y(idx2)*xadj(i,j)*dexp(xlmbda*(xweight(i)-xweight(j)))
         enddo 
         
         if (sum(xadj(:,i))>0.00001d0) then !! AVERAGE OUTGOING
          xsum_out=xsum_out/sum(xadj(:,i))
         else
                 xsum_out=0.d0
         endif
         
         if (sum(xadj(i,:))>0.00001d0) then !! AVERAGE INCOMING
          xsum_in=xsum_in/sum(xadj(i,:))
         else
                 xsum_in=0.d0
         endif
         
         idx2=2*i; idx1=(idx2-1)
         YDOT(idx1)=Y(idx1)*(xrm(i)-Kcar(i)*Y(idx1))-xci(i)*Y(idx2)*Y(idx1)/(1.d0+xci(i)*xhi(i)*Y(idx1))
         YDOT(idx2)=Y(idx2)*xweight(i)-xm*Y(idx2)*xsum_out+xm*xsum_in
         
        enddo 
         
   RETURN
   END SUBROUTINE DERIVS
   
 END MODULE information

    PROGRAM KH
    use information
    use DVODE_F90_M
    
      IMPLICIT NONE
      integer:: i,ii,j,kij,ipar1,ipar2,inic,iter,ieff,iprob,indx1,indx2,kcount,irepchk,kimp
      real(kind=8):: p_waxmn, p_mc, T, TOUT, amn,amx,sumy,xprob,xsd
      real(kind=8):: xlbd, mov_par, xnormsq(id),xmoransum(id),xgearysum(id),mmoran(id),xmoran(id,ninic),&
      & fmmoran(id),mgeary(id),xgeary(id,ninic),fmgeary(id),m2moran(id),sd2moran(id),m2geary(id),sd2geary(id)
      real(kind=8),allocatable:: Y(:), ysum(:), ytemp1(:), ytemp(:,:), tempmat(:,:)      !! dynamical variables
      real(kind=8),allocatable:: adj(:,:), beta(:), dist_mat(:,:), r_max(:), car_cap(:), Pcar_cap(:), C_parm(:),&
      & death_parm(:), ceff_parm(:), h_parm(:), ymean(:), y2mean(:), ysd(:), yinicon(:), y2inicon(:), yiniconsd(:) 

!!!!!!!!!!!!!! INTEGRATOR ENTRIES; USUALLY NOT TO BE TOUCHED !!!!!!     
      INTEGER ITASK, ISTATE, ISTATS, IOUT, SWITCH!
      DOUBLE PRECISION RSTATS, RTOL, ATOL
      DIMENSION RSTATS(22), ISTATS(31)
      TYPE (VODE_OPTS) :: OPTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call initialize()
      ALLOCATE(Y(NEQ), adj(nptch,nptch), beta(nptch), dist_mat(nptch,nptch), r_max(nptch), car_cap(nptch), Pcar_cap(nptch),&
      & C_parm(nptch), death_parm(nptch), ceff_parm(nptch), h_parm(nptch), ymean(id), y2mean(id), ysd(id),&
      & yinicon(id), y2inicon(id), yiniconsd(id))

!      write(1,*) "digraph{"
       do i=1,nptch
        beta(i)=1.d0      ! AREA PROXIES for Waxman
       enddo

!!!!!!!!!!!!!!!!!! SYSTEM PARAMETERS
!!! XXX to discuss

      r_max=1.d0               !! intrinsic max growth rates for the species (can it depend on local conditions?) XXX
      par(1:nptch)=r_max(1:nptch)             !! parameter allocation
      
      Pcar_cap=1.d0
      car_cap=Pcar_cap         !! carrying capacities for each KH XXX
      par(nptch+1:2*nptch)=car_cap(1:nptch) !! parameter allocation
      
      !call random_number(C_parm)        !! C parameter XXX
      C_parm=(/(0.25d0+dfloat(ii)/(nptch+1), ii=1,nptch) /)
      par(2*nptch+1:3*nptch)=C_parm(1:nptch); print*, C_parm   !! parameter allocation      
      
      death_parm=0.02d0        !! death parameter XXX
      par(3*nptch+1:4*nptch)=death_parm(1:nptch)   !! parameter allocation

      ceff_parm=0.25d0         !! conversion efficiency parameter XXX
      par(4*nptch+1:5*nptch)=ceff_parm(1:nptch)   !! parameter allocation
      
      h_parm=3.d0           !! handling parameter XXX
      par(5*nptch+1:6*nptch)=h_parm(1:nptch)   !! parameter allocation
      
      xlbd=200.d0              !! lambda
      par(6*nptch+1)=xlbd    !! parameter allocation XXX
      
      mov_par=0.005d0         !! movement parameter XXX
      par(6*nptch+2)=mov_par    !! parameter allocation
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    dist_mat=0.d0
      do i=1,nptch
       do j=i+1,nptch
        call random_number(dist_mat(j,i)); dist_mat(i,j)=dist_mat(j,i)  !! dist_mat is the distance matrix for the KHs
        enddo
      enddo
      
        irepchk=0
11      kcount=0
        kimp=0
    yinicon=0.d0; y2inicon=0.d0; yiniconsd=0.d0
    
   do inic=1,ninic      
   
      ymean=0.d0; y2mean=0.d0

      adj=0.d0; kij=0
        do i=1,nptch                                             ! 
          do j=1,nptch
          kij=kij+1                                            !
            p_waxmn=(beta(i)*beta(j))*dexp(-dist_mat(i,j)/xL)  !(beta(i)*beta(j))*dexp(-dist_mat(i,j)/xL)  !
            !call random_number(p_mc)                           ! Adjacency matrix
            !if(p_mc<p_waxmn) adj(i,j)=1.d0
            adj(i,j)=1.d0
            if(i==j) adj(i,j)=0.d0                             !
             par(6*nptch+2+kij)=adj(i,j)
          enddo                                                ! 
         enddo
!!!!!!!! Set the integration parameters:
      RTOL = 1.0D-6
      ATOL = 1.0D-9
      IOUT = 1
      ITASK = 1
      ISTATE = 1
      !WRITE (6,*) RTOL, ATOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
!     Output points:
      T = 0.0D0
      TOUT=0.d0
      
      call random_number(Y) !! initial conditions
!     Start with nonstiff method:
      OPTIONS = SET_OPTS(RELERR=RTOL,ABSERR=ATOL)

      do iter=1,ntrans  !! iteration loop
       TOUT=TOUT+h
        CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)

!!!!! ERROR ANALYSIS: DO NOT TOUCH    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Gather and write the integration statistics for this problem:
        CALL GET_STATS(RSTATS,ISTATS)
!        Stop the integration if an error occurred:
        IF (ISTATE<0) THEN
          WRITE (6,*) ISTATE
          exit
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mod(iter,net_Tscl)==0) then  !! NETWORK RESET
       adj=0.d0; kij=0 !! INITIALIZE ADJACENCY MATRIX
         do i=1,nptch                                                             ! 
          do j=1,nptch
            kij=kij+1                                                             !
            p_waxmn=(beta(i)*beta(j))*dexp(-dist_mat(i,j)/xL)                   !
            call random_number(p_mc)                            ! Adjacency matrix
            if(p_mc<p_waxmn) adj(i,j)=1.d0                      !
            if(i==j) adj(i,j)=0.d0                              !
            par(6*nptch+2+kij)=adj(i,j)     
          enddo                                                 ! 
         enddo
     endif
      enddo !!! transient loop finishes

      
      do iter=1,niter  !! iteration loop

       IF (ISTATE<0) THEN
          kcount=kcount+1
          WRITE (6,*) ISTATE
          exit
        END IF
        
       TOUT=TOUT+h
        CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
!!!!! ERROR ANALYSIS: DO NOT TOUCH    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Gather and write the integration statistics for this problem:
        CALL GET_STATS(RSTATS,ISTATS)
!        Stop the integration if an error occurred:
        IF (ISTATE<0) THEN
          WRITE (6,*) ISTATE
          exit
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mod(iter,net_Tscl)==0) then  !! NETWORK RESET
       adj=0.d0; kij=0 !! INITIALIZE ADJACENCY MATRIX
         do i=1,nptch                                                             ! 
          do j=1,nptch
            kij=kij+1                                                             !
            p_waxmn=(beta(i)*beta(j))*dexp(-dist_mat(i,j)/xL)                   !
            call random_number(p_mc)                         ! Adjacency matrix
            if(p_mc<p_waxmn) adj(i,j)=1.d0                   !
            if(i==j) adj(i,j)=0.d0                           !
            par(6*nptch+2+kij)=adj(i,j)
          enddo                                              ! 
         enddo
       endif

    do i=1,nptch
     indx2=2*i; indx1=indx2-1

     ymean(1)=ymean(1)+Y(indx1)
     ymean(2)=ymean(2)+Y(indx2)
     
     y2mean(1)=y2mean(1)+Y(indx1)**2
     y2mean(2)=y2mean(2)+Y(indx2)**2

    enddo
    
    enddo !! end evolution loop

    ymean=ymean/dfloat(niter)
    y2mean=y2mean/dfloat(niter)
    Ysd=dsqrt(dabs(y2mean-ymean))

    if(ymean(1)>0.00000001d0.and.ymean(2)>0.00000001d0) then
       do i=1,id       
        yinicon(i)=yinicon(i)+Ysd(i)/ymean(i)
        y2inicon(i)=y2inicon(i)+(Ysd(i)/ymean(i))**2
       enddo      
    ! kimp=kimp+1
    endif 
    write(11,24) dfloat(nptch),ymean,y2mean,Ysd/ymean

   if(kcount>30.and.irepchk<5) then
     irepchk=irepchk+1
     goto 11
   endif
    
    enddo !! inic
    
        
    yinicon=yinicon/dfloat(ninic-kcount); y2inicon=y2inicon/dfloat(ninic-kcount)
    yiniconsd=dsqrt(dabs(y2inicon-yinicon))/sqrt(dfloat(ninic-kcount))

    write(12,24) dfloat(nptch), yinicon, yiniconsd, dfloat(kcount)
    write(*,24) dfloat(nptch), yinicon, yiniconsd, dfloat(kcount)
    
      24  Format(52F17.10)
!     Write the integration final integration statistics:
      WRITE (6,*) ISTATS(11), ISTATS(12), ISTATS(13)
!     Format statements for this problem:
      
    END PROGRAM KH  
      
subroutine initialize()
      use information
      implicit none
      read*,nptch, xL
      NEQ=id*nptch
end subroutine initialize
