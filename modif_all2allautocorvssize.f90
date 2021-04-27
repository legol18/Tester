 MODULE information   
   implicit none
   integer, parameter:: nptch=128, id=2, NEQ=id*nptch, ninic=50   !! Number of patches
   integer, parameter:: niter=20000,npar=nptch*(nptch+6)+2,ntrans=50000, net_Tscl=1000+ntrans
   real(kind=8), parameter:: d_max=14000.d0, a_max=1000.d0, xL=1.d0, h=0.1d0 
   real(kind=8):: par(npar)
   
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
      integer:: i,ii,j,kij,ipar1,ipar2,inic,iter,ieff,iprob,indx1,indx2,kcount
      real(kind=8):: p_waxmn, p_mc, T, TOUT, amn,amx,sumy,xprob,xsd
      real(kind=8):: xlbd, mov_par, xnormsq(id),xmoransum(id),xgearysum(id),mmoran(id),xmoran(id,ninic),&
      & fmmoran(id),mgeary(id),xgeary(id,ninic),fmgeary(id),m2moran(id),sd2moran(id),m2geary(id),sd2geary(id)
      real(kind=8),allocatable:: Y(:), ysum(:), ytemp1(:), ytemp(:,:), tempmat(:,:)      !! dynamical variables
      real(kind=8),allocatable:: adj(:,:), beta(:), dist_mat(:,:), r_max(:), car_cap(:), Pcar_cap(:), C_parm(:),&
      & death_parm(:), ceff_parm(:), h_parm(:), avin(:), avout(:), xwt_in(:,:), xwt_out(:,:) 

!!!!!!!!!!!!!! INTEGRATOR ENTRIES; USUALLY NOT TO BE TOUCHED !!!!!!     
      INTEGER ITASK, ISTATE, ISTATS, IOUT, SWITCH!
      DOUBLE PRECISION RSTATS, RTOL, ATOL
      DIMENSION RSTATS(22), ISTATS(31)
      TYPE (VODE_OPTS) :: OPTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE(Y(NEQ), adj(nptch,nptch), beta(nptch), dist_mat(nptch,nptch), r_max(nptch), car_cap(nptch), Pcar_cap(nptch),&
      & C_parm(nptch), death_parm(nptch), ceff_parm(nptch), h_parm(nptch), avin(nptch), avout(nptch), ysum(ninic),&
      & xwt_in(nptch,nptch), xwt_out(nptch,nptch), ytemp(id,nptch), tempmat(id,nptch))

!      write(1,*) "digraph{"
       do i=1,nptch
        beta(i)=1.d0!0.5d0!0.01d0+(1.d0-0.01d0)*dfloat(i-1)/dfloat(NEQ)  ! AREA PROXIES for Waxman
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
      kcount=0
   do inic=1,ninic      
      
       adj=0.d0; kij=0; xwt_out=0.d0; xwt_in=0.d0
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
        xwt_in=adj; xwt_out=adj
!!!! normalizing the connection/weight matrix as per total incoming/outgoing connections per node
     do i=1,nptch
      if (sum(xwt_in(i,:))>0.00001d0) xwt_in(i,:)=xwt_in(i,:)/sum(xwt_in(i,:))
      if (sum(xwt_out(:,i))>0.00001d0) xwt_out(:,i)=xwt_out(:,i)/sum(xwt_out(:,i))
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

      avin=0.d0; avout=0.d0 !! average in and out degree variables
      mmoran=0.d0; m2moran=0.d0; mgeary=0.d0; m2geary=0.d0 
      
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
       adj=0.d0; kij=0; xwt_in=0.d0; xwt_out=0.d0 !! INITIALIZE ADJACENCY MATRIX
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
        xwt_in=adj; xwt_out=adj
!!!! normalizing the connection/weight matrix as per total incoming/outgoing connections per node
      do i=1,nptch
       if (sum(xwt_in(i,:))>0.00001d0) xwt_in(i,:)=xwt_in(i,:)/sum(xwt_in(i,:))
       if (sum(xwt_out(:,i))>0.00001d0) xwt_out(:,i)=xwt_out(:,i)/sum(xwt_out(:,i))
      enddo        
     endif
!!!!!! Autocorrelation calculation     

     do i=1,nptch  !! creating vectors
      indx2=2*i; indx1=(indx2-1)
      ytemp(1,i)=Y(indx1)
      ytemp(2,i)=Y(indx2)
     enddo
     
     xgearysum=0.d0
     do i=1,nptch   !for Geary
      do j=1,nptch
      xgearysum(1)=xgearysum(1)+xwt_out(i,j)*(ytemp(1,j)-ytemp(1,i))**2
      xgearysum(2)=xgearysum(2)+xwt_out(i,j)*(ytemp(2,j)-ytemp(2,i))**2
      enddo
     enddo
     
     ytemp(1,:)=ytemp(1,:)-sum(ytemp(1,:))/dfloat(nptch) !! detrending x vector
     ytemp(2,:)=ytemp(2,:)-sum(ytemp(2,:))/dfloat(nptch) !! detrending y vector

     tempmat=0.d0
     do i=1,nptch
      do j=1,nptch
      tempmat(1,i)=tempmat(1,i)+xwt_out(i,j)*ytemp(1,j)
      tempmat(2,i)=tempmat(2,i)+xwt_out(i,j)*ytemp(2,j)
      enddo
     enddo

     xmoransum=0.d0; xnormsq=0.d0
      do j=1,nptch
      xmoransum(1)=xmoransum(1)+tempmat(1,j)*ytemp(1,j)
      xmoransum(2)=xmoransum(2)+tempmat(2,j)*ytemp(2,j)
      xnormsq(1)=xnormsq(1)+ytemp(1,j)**2
      xnormsq(2)=xnormsq(2)+ytemp(2,j)**2
      enddo

    xmoransum(1)=xmoransum(1)/xnormsq(1); xmoransum(2)=xmoransum(2)/xnormsq(2)
    xgearysum(1)=dfloat(nptch-1)*xgearysum(1)/xnormsq(1)/2.d0/dfloat(nptch)
    xgearysum(2)=dfloat(nptch-1)*xgearysum(2)/xnormsq(2)/2.d0/dfloat(nptch)

    mmoran=mmoran+xmoransum
    mgeary=mgeary+xgearysum
    
    enddo !! end evolution loop

    xmoran(:,inic)=mmoran/dfloat(niter); xgeary(:,inic)=mgeary/dfloat(niter)
   enddo !! inic
  fmmoran(1)=sum(xmoran(1,:))/dfloat(ninic-kcount); fmmoran(2)=sum(xmoran(2,:))/dfloat(ninic-kcount)
  fmgeary(1)=sum(xgeary(1,:))/dfloat(ninic-kcount); fmgeary(2)=sum(xgeary(2,:))/dfloat(ninic-kcount)
  
  do i=1,ninic
   m2moran=m2moran+xmoran(:,i)**2
   m2geary=m2geary+xgeary(:,i)**2
  enddo
  m2moran=m2moran/dfloat(ninic-kcount); m2geary=m2geary/dfloat(ninic-kcount)
  
  sd2moran=dsqrt(m2moran-fmmoran**2)/dsqrt(dfloat(ninic-kcount)) 
  sd2geary=dsqrt(m2geary-fmgeary**2)/dsqrt(dfloat(ninic-kcount))
    write(11,24)dfloat(nptch),fmmoran,sd2moran,fmgeary,sd2geary
  print*,kcount
      24  Format(52F17.10)
!     Write the integration final integration statistics:
      WRITE (6,*) ISTATS(11), ISTATS(12), ISTATS(13)
!     Format statements for this problem:
      
    END PROGRAM KH  
      
