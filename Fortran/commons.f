      parameter (nx=2000)
c
      dimension fm(nx), dm(nx), dm12(nx),
     &          r (nx), dr(nx), dr12(nx), rold(nx),
     &          v (nx), a (nx), at12(nx), ak12(nx),
     &          u (nx), p (nx), rho (nx), eps(nx),
     &          w (nx), wt12(nx),
     &          aux(nx)
c
      logical cartesian
c
      common fm,dm,dm12,r,dr,dr12,rold,v,a,at12,
     &       ak12,u,p,rho,eps,w,wt12,
     &       aux,pi,gamma,t,dt,dt12,tmax,q,fmratio,
     &       cartesian,cflfactor,dfactor,efac,il,istep

