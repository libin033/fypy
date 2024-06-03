      Subroutine RAHMBRDFFOS (rho0,af,xk,mu,rm,rosur,wfisur,fisur)
C***********************************************************************
C
C           A semi-empirical bidirectional reflectance model
C
C  Purpose:
C
C  To generate a single bidirectional reflectance factor value for a
C  semi-infinite medium, given the illumination and viewing geometry,
C  the optical and the structural properties of the scatterers.
C
C  Definitions:
c   geometrical conditions
C     mu1          : Illumination zenith angle, in radians
C     mu2          : Observation zenith angle, in radians
C     fi           : Relative azimuth angle, in radians
C   optical characteristics of the scatterers:
C     Rho0         : Intensity of the reflectance of the surface cover.
C                    N/D value greater or equal to 0.0 (Rho_0)
C     af           : Phase function parameter:
C                    Asymmetry factor, N/D value between -1.0 and 1.0
C     xk           : Structural parameter of the medium
C
C  References:
C
C        [1] R. Rahman, Verstraete M., and Pinty B., 1993, `A coupled 
C            surface-atmosphere reflectance (CSAR) model. Part 1: Model
C            Description and Inversion on Synthetic Data', submitted to
C            JGR.
C        [2] R. Rahman, Pinty B., and Verstraete M., 1993, `A coupled 
C            surface-atmosphere reflectance (CSAR) model. Part 2: a 
C            semi-empirical model usable with NOAA/AVHRR data', 
C            submitted to JGR.
C
C***********************************************************************
C 
      include "paramdef.inc"
      integer mu,np,j,k
      real xk,af,rho0
      real rm(-mu:mu)
      real coef1,coef2,cospha,geofac
      real fi,mu1,mu2,phaang,phafun,pi,tante1,tante2
      real rosur(0:mu_p,mu_p,83),fisur(83),wfisur(83),psip
C
       pisp=acos(-1.)
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
      fi=fisur(k)+pisp
C Compute various trigonometric expressions:
      cospha=mu1*mu2+sqrt(1.-mu1*mu1)*sqrt(1.-mu2*mu2)*cos(fi)
      if (cospha.le.-1.0) cospha=-1.0
      if (cospha.ge.1.0) cospha=1.0
      tante1=sqrt(1.-mu1*mu1)/mu1
      tante2=sqrt(1.-mu2*mu2)/mu2
      geofac=sqrt(tante1*tante1+tante2*tante2-2.0*tante1*tante2*cos(fi))
C Compute the first term
      coef1=(mu1**(xk-1.))*(mu2**(xk-1.))/((mu1+mu2)**(1.-xk))
C Compute the phase function:
      phafun=(1.0-af*af)/((1.0+af*af-2.0*af*(-cospha))**1.5)

C Compute the opposition (hot spot) function:
      coef2=1.+(1.-rho0)/(1.+geofac)
C Compute the bidirectional reflectance factor:
      rosur(i,j,k)=rho0*coef1*phafun*coef2
c      write(6,*) "rosur", mu1,mu2,fi,rosur(i,j,k)
      enddo
      enddo
      enddo
C
      return
      end

