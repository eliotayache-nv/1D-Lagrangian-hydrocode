      program lh1
c
c     This Lagrangian code follows the adiabatic expansion 
c     of a hot spherical bubble in a uniform ambient medium.
c     The present setup is ajusted to spherical geometry.
c     The results may be plotted in Supermomgo environment with the help
c     of commands contained in files pltmod and plte.
c     The code can also work in Cartesian geometry. To switch from 
c     spherical to Cartesian geometry go to subroutine inicond and
c     set cartesian = .true.
c
      include 'commons.f'
c
c                                  open the input-data file
c
      open(15,file='lh1.dat',form='formatted')
c
c                                  open the file in which 
c                                  the final model will be stored 
c
      open(16,file='lh1.out',form='formatted')
c
c                                  open the file in which values
c                                  of thermal, kinetic and total
c                                  energy will be stored at every
c                                  10th time-step
c
      open(17,file='lh2.out',form='formatted')
c
      read(15,*)
      read(15,*)nsteps
      read(15,*)tmax   
      read(15,*)efac   
      read(15,*)il
      read(15,*)q
      read(15,*)fmratio
c
c                                  set initial conditions
c
      call inicond
c
c
c                                  begin loop over time steps
c
c
      do istep = 1, nsteps
c
c                                  determine time step
c
        call tmsc
c
c                                  update velocities
c
        do ix=2,nx
          u(ix) = u(ix)-a(ix)*(p(ix)-p(ix-1))*dt/dm(ix)
     &           -0.5*(w(ix  )*(3.*ak12(ix  )-a(ix))
     &                -w(ix-1)*(3.*ak12(ix-1)-a(ix)))
     &               *dt/dm(ix)
        end do
        if (cartesian) then
          u(1) = u(2)
        else
          u(1) = 0.
        end if
c
c                                 update radii, surfaces and volumes
c
        do ix=1,nx
          rold(ix) = r(ix)
        end do
        do ix=1,nx
          r(ix)    = rold(ix)+u(ix)*dt12
        end do
        do ix=1,nx-1
          dr12(ix) = r(ix+1)-r(ix)
        end do
        if (cartesian) then
          do ix=1,nx
            v   (ix) = r(ix)
          end do
        else
          do ix=1,nx
            at12(ix) = 4.*pi*(0.5*(r(ix)+rold(ix)))**2
            a   (ix) = 4.*pi*r(ix)**2 
            v   (ix) = 4./3.*pi*r(ix)**3 
          end do
          do ix=1,nx-1
            ak12(ix) = 0.5*(at12(ix+1)+at12(ix))
          end do
        end if
c
c                                 update densities
c
        do ix=1,nx-1
          rho(ix)  = dm12(ix)/(v(ix+1)-v(ix))
        end do
          rho(nx)  = rho(nx-1)
c
c                                 artificial viscosity
c
        do ix=1,nx-1
          w  (ix)  =-q**2*rho(ix)*abs(u(ix+1)-u(ix))
     &              *(u(ix+1)*(1.-at12(ix+1)/3./ak12(ix))
     &               -u(ix  )*(1.-at12(ix  )/3./ak12(ix)))
        end do
        do ix=1,nx-1
          if (u(ix+1).gt.u(ix)) w(ix)=0.
        end do
c
c                                 update internal energies and pressures
c
        do ix=1,nx-1
          aux(ix)  = eps(ix)-p(ix)
     &              *( at12(ix+1)*u(ix+1)
     &                -at12(ix  )*u(ix  ))*dt12/dm12(ix)
        end do
        do ix=1,nx-1
          p(ix)    = 0.5*(p(ix)+(gamma-1.)*rho(ix)*aux(ix))
        end do
        do ix=1,nx-1
          eps(ix)  = eps(ix)-p(ix)
     &              *( at12(ix+1)*u(ix+1)
     &                -at12(ix  )*u(ix  ))*dt12/dm12(ix)
        end do
c
c contribution from artificial viscosity
c
        do ix=1,nx-1
          eps(ix)  = eps(ix)-0.5*w(ix)*dt12/dm12(ix)
     &              *(u(ix+1)*(3.*ak12(ix)-at12(ix+1))
     &               -u(ix  )*(3.*ak12(ix)-at12(ix  )))
        end do 
c
        do ix=1,nx-1
          p(ix)    = (gamma-1.)*rho(ix)*eps(ix)
        end do
        p  (nx) = p  (nx-1)
        eps(nx) = eps(nx-1)
c
c                                 end time step
c
        t = t+dt12  
c
c                                 check energy conservstion
c
        ethe = 0.
        ekin = 0.
        do ix=2,nx-1
          ethe = ethe + (eps(ix))*dm(ix)
          ekin = ekin + 0.5*(0.5*(u(ix+1)+u(ix)))**2*dm(ix)
        end do
        etot = ethe + ekin
        if (istep.eq.1) etot0 = etot
        etot = etot/etot0
        ethe = ethe/etot0
        ekin = ekin/etot0
        if (mod(istep,10).eq.0)  write(17,200)istep,t,etot,
     &                   ethe, ekin
        if (mod(istep,100).eq.0) write(* ,100)istep,t,etot,
     &                   ethe, ekin
c
c
c                                      end loop over time steps

        if (t.eq.tmax) goto 10
c
      end do
c
c
c                                      final printout
c
   10 continue

      umax=0.
      rhomax=0.
      pmax=0.
      emax=0.
      wmax=0.

      do ix=1,nx
         if (u(ix).gt.0) umax  = max(u(ix),umax)
         rhomax = max(rho(ix),rhomax)
         pmax   = max(p(ix),pmax)
         emax   = max(eps(ix),emax)
         wmax   = max(w(ix),wmax)
      end do

      epsi=1.e-20

      do ix=1,nx
         u(ix)   = u(ix)/(epsi+umax)
         rho(ix) = rho(ix)/(epsi+rhomax)
         p(ix)   = p(ix)/(epsi+pmax)
         eps(ix) = eps(ix)/(epsi+emax)
         w(ix)   = w(ix)/(epsi+wmax)
      end do

      write(16,101)(ix,fm(ix),r(ix)/r(nx),u(ix),
     &      rho(ix),p(ix),eps(ix),w(ix),ix=1,nx)
  100 format(1x,i5,'; t:',1pe10.2,';  etot:',1pe10.2,
     &          ';  eth:',1pe10.2,';  ekin:',1pe10.2)
  101 format(1x,i5,1p7e12.4) 
  200 format(1x,i5,1p4e11.3)
c
      stop
      end
c
c
c
c
      subroutine inicond
      include 'commons.f'
c
      cartesian = .false.
c
c     Courant factor, must be smaller than 1
      cflfactor = 0.1
c
c     diffusion factor, must be greater that 1
      dfactor   = 1.5
      dfactor   = dfactor**2
c
      gamma     = 5.0/3.0
c
      pi        = asin(1.0)*2.
c
c                                         high pressure ejecta
c
      rhoej  =  ((float(nx)/float(il))**3-1.)*fmratio
      dmej   =  4./3.*pi*(float(il)/float(nx))**3*rhoej/float(il)
      do ix = 1,il
        dm12(ix) = dmej
        rho (ix) = rhoej
        eps (ix) = efac
        p   (ix) = (gamma-1.)*rho(ix)*eps(ix)
        u   (ix) = 0.
      end do
c     
c                                         low pressure ambient medium
c
      rhoamb = 1.
      dmamb  = 4./3.*pi*(1.-(float(il)/float(nx))**3)/float(nx-il)
     &         *rhoamb
      do ix=il+1,nx
        dm12(ix) = dmamb
        rho (ix) = rhoamb
        eps (ix) = 1.
        p   (ix) = (gamma-1.)*rho(ix)*eps(ix)
        u   (ix) = 0. 
      end do
c
      fm(1) = dm12(1) 
      do ix=2,nx
        fm (ix) = fm(ix-1)+dm12(ix)
      end do
      do ix=2,nx
        dm(ix) = 0.5*(dm12(ix)+dm12(ix-1))
      end do
c
c
      r   (1) = 0.
      v   (1) = 0.
      do ix=2,nx
        v(ix) = v(ix-1)+dm12(ix-1)/rho(ix-1)
        if (cartesian) then
          r(ix)    = v(ix)
          a(ix)    = 1.
          at12(ix) = 1.
        else
          r(ix) = (v(ix)/(4./3.*pi))**(1./3.)
          a(ix) = 4.*pi*r(ix)**2
        end if
      end do
      do ix=1,nx-1
        dr12(ix) = r(ix+1)-r(ix)
      end do
c
      t    = 0.
c
      return
      end
c
c
c
c
      subroutine tmsc
      include 'commons.f'
c
      dtc = 1.e30
      do ix=1,nx-1
        dtc = min(dtc,dr12(ix)
     &         /(abs(u(ix))+sqrt(gamma*eps(ix))))
      end do
c
      dtc  = cflfactor*dtc
      if (t+dtc.gt.tmax) dtc=tmax-t 
c
c     diffusion limit
c
      dtd = 1.e-30
      do ix=1,nx-1
        dtd = max(dtd,abs(at12(ix+1)*u(ix+1)-at12(ix)*u(ix))
     &                  /(   v(ix+1)        -v(ix)         )
     &           )
      end do
      dtd = 0.5/dtd/dfactor
c
      dtc = min(dtc,dtd)
c
      dt   = 0.5*(dt12+dtc)
      dt12 = dtc
c
      return
      end
