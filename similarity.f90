     SUBROUTINE similarity(xx,yy,simfx)!! similarity function calculation. Can be used to estimate periodicity of the signal
                                       !! by providing the same signal as the xx,yy input. 
     implicit none
     real(kind=8), intent(in):: xx(niter),yy(niter)
     real(kind=8), intent(out):: simfx(imaxlag)
     real(kind=8):: simsum, sqm1,sqm2
     integer:: i,j
     
     simfx=0.d0

     do i=1,imaxlag
       simsum=0.d0; sqm1=0.d0; sqm2=0.d0
       do j=1,niter-i
        simsum=simsum+(xx(j)-yy(j+i))**2
        sqm1=sqm1+xx(j)**2/dfloat(niter-i)
        sqm2=sqm2+yy(j)**2/dfloat(niter-i)
       enddo
      simfx(i)=(simsum/dfloat(niter-i))/sqrt(sqm1*sqm2)
     enddo
     
     RETURN
     end SUBROUTINE similarity
