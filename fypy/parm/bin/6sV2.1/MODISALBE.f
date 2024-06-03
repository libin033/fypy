      subroutine modisalbe(p1,p2,p3,
     s           brdfalb)
      real p1,p2,p3,brdfalb
      brdfalb=p1+p2*0.23062499-p3*1.6798152
      if (brdfalb.lt.0) then
      Write(6,*) "warning white sky albedo <0  , reset to 0."
      brdfalb=0.
      endif
      return
      end
