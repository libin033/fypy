      subroutine modisbrdffos(p1,p2,p3,mu,rm,
     s           rosur,wfisur,fisur)
      include "paramdef.inc"
      real p1,p2,p3,xmu,view
      real dts,dtv,dfs,dfv,dfi
      real rts,rtv,rfs,rfv,rfi,rpha
      real cts,ctv,cfi,cpha,ct0
      real sts,stv,sfi
      real tanti,tantv
      real cost,sint,tvar
      real rossthick,rosselt,lispars
      real angdist,angtemp,angover
      integer mu,np,k,j
      real rm(-mu:mu)
      real rosur(0:mu_p,mu_p,83),fisur(83),wfisur(83),psip
      real mu1,mu2
      
      
       pisp=acos(-1.)
       pi=atan(1.)*4.
       fac=180./pi
c      mu=mu_p
      call gauss(0,pisp,fisur,wfisur,83)
      do i=0,mu
      do j=1,mu
      do k=1,83
      if (i.eq.0) then
      mu1=rm(0)
      else
      mu1=rm(i)
      endif
      mu2=rm(j)
      rfi=fisur(k)+pisp
      rts=acos(mu1)
      rtv=acos(mu2)
      if ((rts*180./pi).gt.75.) rts=75.*pi/180.
      if ((rtv*180./pi).gt.65.) rtv=65.*pi/180.      
        rfi=abs(rfi)
	cts=cos(rts)
	ctv=cos(rtv)
	sts=sin(rts)
	stv=sin(rtv)
	cfi=cos(rfi)
	sfi=sin(rfi)
	cpha=cts*ctv+sts*stv*cfi
	if (cpha.gt.1.0)  cpha=1.0
	if (cpha.lt.-1.0) cpha=-1.0
	rpha=acos(cpha)
		
	rosselt=(pi/2-rpha)*cpha+sin(rpha)
	rossthick=(rosselt/(cts+ctv))-pi/4.
	
	tanti=tan(rts)
	tantv=tan(rtv)
	
	angdist=tanti*tanti+tantv*tantv-2.*tanti*tantv*cfi
	angdist=sqrt(angdist)
	
	angtemp=1./cts+1./ctv
	cost=2.*sqrt(angdist*angdist+tanti*tanti*tantv
     &	*tantv*sfi*sfi)
	cost=cost/angtemp
	if (cost.ge.1.) cost=1.
	if (cost.le.-1.) cost=-1.
	tvar=acos(cost)
	sint=sqrt(1.-cost*cost)
	angover=(tvar-sint*cost)*angtemp/pi
	lispars=angover-angtemp+0.5*(1.+cpha)/cts/ctv
      
      rosur(i,j,k)=p1+p2*rossthick+p3*lispars
c      write(6,*) "modisbrffos ",acos(mu1)*fac,acos(mu2)*fac,rfi*fac,
c     &  rosur(i,j,k)
      if (rosur(i,j,k).lt.0.) then 
          rosur(i,j,k)=0.
          write(6,*) "Warning negative reflectance from the BRDF model->set to 0."
	  endif
      enddo
      enddo
      enddo
      return
      end
