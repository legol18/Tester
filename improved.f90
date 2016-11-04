MODULE WANG   
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
            
            do i=1,imn                    !!
                kc=kc+1; akc(i)=par(kc)     !!
            enddo                         !! carrying capacities
            
            acpt=0.0                                   !! competition initialization
            do i=1,imn; do j=1,imn                     !!
                kc=kc+1; acpt(i,j)=par(kc)               !! 
                if(j==i) acpt(i,j)=0.0                   !! competition terms
                if((i-10)*(j-10)<0) acpt(i,j)=0.0        !!
                if((i-10)>0.and.(j-10)==0) acpt(i,j)=0.0 !!
            enddo; enddo                         !!
            
            do i=1,imn                     !!
                kc=kc+1; ad(i)=par(kc)       !!
            enddo                          !! dispersal terms
            
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
            if(itj/=icptch) then
            cdsum=cdsum+ad(iti)*Y(iflag+(itj-1)*in)
            !print*,iti,icptch,iflag+(itj-1)*in,ad(iti); pause!!
            endif
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
    INTEGER, PARAMETER :: mlyap=0, NEQ=(mlyap+1)*imn
    real  :: h, T, TOUT, xrho1, xrho2, xvar1, xvar2, kppa
    real  :: spsum, svsum, ald, sald, btd, sbtd, gmd, sgmd, alcv, btcv, gmcv
    integer:: IOUT, i, ii, j, ntrans, niter, kc, inic, ninic, ncalc, ndisc
    real,allocatable:: Y(:), Yold(:), YDER(:), noi_spec(:), noi_ptch(:)
    real,allocatable:: lo_simp(:), omg(:), ymean(:), ysd(:), ycovar(:), ypopl(:,:), yprel(:,:), ysave(:,:,:)
    
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

xrho2=0.0  !!print*, "Between patch correlation"; 
xrho1=0.0   !!print*, "Between species correlation"; 
xvar1=0.04  !! species specific variance
xvar2=0.04  !! patch specific variance
kppa=8.0    !! kappa
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
ntrans=50000; niter=50000; ndisc=40000; ninic=2000; ncalc=niter-ndisc
allocate(Y(NEQ), YDER(NEQ), noi_spec(imn), noi_ptch(imn), lo_simp(im), omg(im), ypopl(im,in), yprel(im,in), ymean(im), &
ysd(im), ycovar(im), ysave(im,in,ncalc))
    
do inic=1,ninic; print*,inic
    par=0.0
!     Set  the problem parameters:
    call random_number(par(1:imn)); 
    par(1:imn)=0.5+par(1:imn)*1.0                 !! M*N MAXIMAL GROWTH RATES
!       do i=1,(im-1)
!         pos1=i*in-igen+1; pos2=(i+1)*in-igen+1
!         par(pos1:i*in)=par(pos2:(i+1)*in)
!       enddo                                    !!!!!!!!!!!!!!!!!!
!       !!!!!!!!!!!!!!!!!!!!!!!!!
    
    call random_number(par(imn+1:2*imn))
    par(imn+1:2*imn)=0.5+par(imn+1:2*imn)*1.0  !! M*N CARRYING CAPACITIES
                                    !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    par(2*imn+1:2*imn+imn**2)=0.75
    call random_number(par(2*imn+1:2*imn+imn**2)); 
    par(2*imn+1:2*imn+imn**2)=0.25+par(2*imn+1:2*imn+imn**2)*0.5 !! M*N*N COMPETITION COEFFICIENTS FOR M PATCHES (NEED TO THINK ............!!!) **
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call random_number(par(2*imn+imn**2+1+in-igen:2*imn+imn**2+in)) !! DISPERSAL RATES
    par(2*imn+imn**2+1+in-igen:2*imn+imn**2+in)=0.05+par(2*imn+imn**2+1+in-igen:2*imn+imn**2+in)*0.1
    do i=1,(im-1)
        pos1=(2*imn+imn**2+1)+i*in-igen; pos2=(2*imn+imn**2+1)+(i+1)*in-igen 
        !print*,pos1,(2*imn+imn**2)+i*in,pos2,(2*imn+imn**2)+(i+1)*in,npar; pause
        par(pos2:(2*imn+imn**2)+(i+1)*in)=par(pos1:(2*imn+imn**2)+i*in)
    enddo !; print*,par(2*imn+imn**2+1:2*imn+imn**2+imn)      ; pause
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*,par(2*imn+imn**2+1:2*imn+imn**2+imn); pause
!     Set the initial conditions:
    call random_number(Y); YDER=0.0; ; noi_spec=0.0; noi_ptch=0.0
!     Output points:
    T = 0.0
    h=0.1
    
!     Perform the integration:
    DO IOUT = 1, ntrans
        call DERIVS(NEQ,T,Y,YDER)
        Y=Y+h*YDER; T=T+h
    ENDDO
    
    Y=Y+0.1
    first = .true.
    
    DO IOUT = 1, niter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EVOLUTION BEGINS            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! between species random number patch
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
        do i=1,imn-igen; noi_spec(i)=xep(i); enddo;         !!!!!!!!!!!!!noise arrangement for integration
        do i=1,(im-1)  !! making up the dummy noise terms
        !print*,(i*in+1)+(in-igen),(i*in)+in,(i-1)*in+1+(in-igen),((i-1)*in)+in; pause
        noi_spec((i*in+1)+(in-igen):(i*in)+in)=noi_spec(((i-1)*in+1)+(in-igen):((i-1)*in)+in) !! ERROR IN INDICES CORRECTED
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
        Y=Y+h*YDER+sqrt(h)*Y*(noi_spec+noi_ptch+kppa*noi_ptch*noi_spec) !!  noise added

        do i=1,imn
        if (Y(i)<1E-10) Y(i)=0.0
        enddo
        kc=0
        do i=1,im
        do j=1,in
        kc=kc+1
        ypopl(i,j)=Y(kc)
        enddo
        enddo        
    if (IOUT>ndisc) then
        kc=0
        do i=1,im
        do ii=1,in
        kc=kc+1
        ysave(i,ii,(IOUT-ndisc))=Y(kc)
        enddo
        enddo
    endif
    T=T+h
    first = .false.
    END DO ! IOUT
    
    sald=0.0; sbtd=0.0; sgmd=0.0      
do IOUT=1,ncalc       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EVOLUTION DONE        
!!!! COMPUTATION OF ENS BASED MEASURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!local simpson indices
do i=1,im
    yprel(i,:)=ysave(i,:,IOUT)/sum(ysave(i,:,IOUT))
enddo 

        do i=1,im
        lo_simp(i)=sum(yprel(i,:)**2)
        enddo 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,im                                    !!!! weight of each patch
        omg(i)=sum(ysave(i,:,IOUT))/sum(ysave(:,:,IOUT))
        enddo 
        
        ald=1.0/sum(omg*lo_simp)             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! alpha_D (QC done)
        sald=sald+ald
!!!! global simpson index
        spsum=0.0
        do i=1,im
            do j=1,in; !print*,i, j, (i-1)*in+j                             QC done
            spsum=spsum+(yprel(i,j)*omg(i))**2          !(Yrelptch((i-1)*in+j)*omg(i))**2            
            enddo 
        enddo 
        gmd=1.0/spsum                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! gamma_D
        btd=gmd/ald                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! beta_D
        sbtd=sbtd+btd; sgmd=sgmd+gmd
enddo  !!IOUT 
sald=sald/float(ncalc); sbtd=sbtd/float(ncalc); sgmd=sgmd/float(ncalc) 
!         
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! !!!!! VARIABILITY BASED MEASURES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,im        !!!! temporal mean for each patch
        ymean(i)=sum(ysave(i,:,:))/float(ncalc)
        enddo          !                                       QC done
!         
        do i=1,im                                    !!!! standard dev. for each patch
        spsum=0.0
        do IOUT=1,ncalc                     !                                  QC done
            spsum=spsum+(sum(ysave(i,:,IOUT))-ymean(i))**2
        enddo       !                                                    QC done 
        ysd(i)=sqrt(spsum/float(ncalc))
        enddo
        alcv=(sum(ysd)/sum(ymean))**2
! 
        do i=1,im
        svsum=0.0
        do ii=i,im
            spsum=0.0
            do IOUT=1,ncalc
            spsum=spsum+(sum(ysave(i,:,IOUT))-ymean(i))*(sum(ysave(ii,:,IOUT))-ymean(ii))
            enddo
            if(spsum<0.0) then 
            spsum=-sqrt(abs(spsum/float(ncalc)))
            else
            spsum=sqrt(spsum/float(ncalc))
            endif 
            svsum=svsum+spsum
        enddo
        enddo 
        gmcv=(sqrt(svsum)/sum(ymean))**2
        btcv=gmcv/alcv
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(11,24) sald, sbtd, sgmd, alcv, btcv, gmcv 
!       endif   
!       enddo !! IOUT   
    enddo !!!  INIC
            
!!! OUTPUT FILE FORMAT
    24  Format(40F20.10)

END PROGRAM DEMOWANG
