       subroutine hilb1   (wx,wp,s,dw,
      d              nw,nwp,
      o              x)

c 94.08.25
c hilbert transform X(w) = I[wmin,wmax] dw' S(w')/(w-w')
c wx = w, wp = w'
c NOTE: wp must be on uniform mesh but w can be anything

       implicit real*8 (a-h,o-z)
       dimension wx(nw),wp(nwp),s(nwp),
      o          x(nw)
       data tol/1.d-6/

*      write(6,*) 'in hilb'
*      do   iw  =  1,nwp
*      write(6,6000) wp(iw),s(iw)
*      enddo
* 6000 format (1x,6f14.6)
*      do   iw  =  1,nw
*      write(6,6000) wx(iw)
*      enddo


c constants
       odw        = 1.d0/dw
       hdw        = 0.5d0*dw

c loop over w
       do 111  iw = 1,nw
       w          = wx(iw)
       sum        = 0.d0

c w is outside wmin,wmax
       if (w .lt. wp(2)-hdw .or. w .gt. wp(nwp-1)+hdw) then
       do     iwp = 3,nwp-2
       sum        = sum + s(iwp)/(w-wp(iwp))
       enddo
       sum        = sum + 0.5d0*( s(2)  /(w-wp(2))
      .                          +s(nwp-1)/(w-wp(nwp-1)) )
       x(iw)      = dw*sum
       goto 111
       endif

c locate the point just before singularity
       w0         = w - wp(1)
       iw0        = idnint(odw*w0)

c integrate over w' to the point just before singularity
       do     iwp = 2,iw0-1
       sum        = sum + s(iwp)/(w-wp(iwp))
       enddo
       if (iw0 .ge. 2)
      .sum        = sum + 0.5d0*( s(1)/(w-wp(1)) + s(iw0)/(w-wp(iw0)) )

c take into account singularity
       sum        = sum -(s(iw0+2)-s(iw0))*odw

c integrate from the point just after the singularity to wmax
       do     iwp = iw0+3,nwp-1
       wwp        = 1.d0/(w - wp(iwp))
       sum        = sum + wwp*s(iwp)
       enddo
       if (iw0 .le. nwp-3)
      .sum        = sum + 0.5d0
      .           * ( s(iw0+2)/(w-wp(iw0+2)) + s(nwp)/(w-wp(nwp)) )

       x(iw)      = dw*sum
   111 continue

       return
       end
