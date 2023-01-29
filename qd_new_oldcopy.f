      program qdot  
      use omp_lib
      implicit none

      include 'sizes.f'

      character*3 prefix

      complex*16 one,ione,zero
      parameter (one= (1.d0,0.d0))
      parameter (ione=(0.d0,1.d0))
      parameter (zero=(0.d0,0.d0))

      real*8     pi,twopi
      parameter (pi=3.14159265358979334d0)
      parameter (twopi=6.28318530717958668d0)

c_____declarations for bath common block
      complex*16 aa(nb),gg(nb)

c_____declarations for syst common block
      complex*16 xx(n,n),hh(n,n)

c_____declarations for dmat common block
      complex*16 zz(n,n,nb)  ! for non-Markovian propagation
      complex*16 zs(n,n)     ! for Markov propagation

c_____declarations for laser common block       
      real*8     pulse_par(6)
      real*8     lambda,delta1,delta2

      complex*16 corr    ! for bath symmetry testing
      complex*16 ccc,sss

c_____temperature and time stuff
       integer   ntime
       integer   i,j,k

       integer   kk,nruns


       parameter (nruns=48)
c      parameter (nruns=8*96)
c      parameter (nruns=1)

       real*8     signal(nruns,3),efield(nruns),scanvar(nruns)

       complex*16 w,c1,dw,jjc,jjcxx,x,nbeta

       real*8     kT,kT_K,beta,wc,wc_eV
       real*8     maxtime,time,tend,timestep
       real*8     phisec,tau0,taup,alpha,w0,e0,ep
       real*8     delta1_au,delta2_au,w0_au
       real*8     fr,fi
       real*8     dd,om,oo,nx,nz
       complex*16 sin_cos



       parameter (timestep=4.0d0/0.024189d0)

       real*8 xold,h

       common/coup/xx
       common/bath/beta,aa,gg   
       COMMON/CONDO5/XOLD,H
       !$OMP THREADPRIVATE(/coup/,/bath/,/CONDO5/)


c____________________________________________________
      read(5,'(a)') prefix
      read(5,*) w0,e0,tau0,phisec           ! pump pulse in fs
      read(5,*) lambda                      ! coupling strength 
      read(5,*) delta1                       ! excitation energy 
      read(5,*) delta2                       ! excitation energy 
      read(5,*) kT_K                        ! temperature       
      read(5,*) wc_eV                       ! cut-off energy    

c______writing parameter file
      open (44,file=prefix//'.param')
      write(44,'(a)') prefix
      write(44,*) 'central frequency w0',w0
      write(44,*) 'field amplitude e0'  ,e0
      write(44,*) 'pulse length FT pulse tau0' ,tau0
      write(44,*) 'chirp paramater phisec', phisec
      write(44,*) 'coupling strength paramater lambda',  lambda
      write(44,*) 'delta1    ',  delta1
      write(44,*) 'delta2    ',  delta2
      write(44,*) 'temperature ',  kT_K   
      write(44,*) 'cut-off    ',  wc_eV
      close(44)

c      open(67,file=prefix//'.dat')


c     call omp_set_num_threads( 6 )  
c     call omp_set_num_threads( 24 )  
      call omp_set_num_threads( 48 )  
c     call omp_set_num_threads( 1 )  

      !$omp parallel default (private) 
     . shared (wc_eV,kT_K,w0,e0,tau0,taup,
     .         phisec,lambda,delta1,delta2,efield,scanvar,signal)


c______conversion into a.u.
       wc=wc_eV/27.211d0
c      beta=1.d0/(kT_K*3.16686d-6)
       beta=1.d0/(kT_K/11604.d0/27.211d0)

c_____conversion from cm-1 to a.u.
       w0_au   =w0/219474.d0
       delta1_au=delta1/219474.d0
       delta2_au=delta2/219474.d0

c_____here: everything in fs !!, expressions from J. Degert, these, eq. I.11
       taup=tau0*dsqrt(1.d0+(2.d0*phisec/tau0**2)**2)  ! taup in fs
       alpha=2.d0*phisec/(tau0**4+(2.d0*phisec)**2)    ! alpha in fs^-2
       ep=e0*dsqrt(tau0/taup)

c_____everything into atomic units
       taup=taup/0.024189d0
       alpha=alpha*0.024189d0**2

c      ep=e0*pi/dsqrt(pi*tau0*taup)  ! here: e0 is in multipleof Pi

       pulse_par(1)=ep
       pulse_par(2)=w0_au
       pulse_par(3)=taup
       pulse_par(4)=alpha
       pulse_par(5)=delta1_au
       pulse_par(6)=delta2_au

       call init_coup(n,xx)
       call init_bath(nb,wc,beta,aa,gg)




c______check bath correlation function
c      do i=1,100
c      time=dble(i-1)/dble(100)*200000.d0
c      call check_bath(33,nb,wc,beta,time,aa,gg)
c      enddo
c      stop

c      time=2000.d0

c________________________
      do i=1,nb
      aa(i)=lambda*aa(i)
      enddo

c     time=00.d0
c     call check_bath(33,nb,wc,beta,lambda,time,aa,gg)

      
c      oo=0.d0
c      write(*,*) 'sss ',sss(nb,oo,aa,gg)
c      write(*,*) 'ccc ',ccc(nb,oo,aa,gg)

c      stop


       oo=wc/3.d0
       write(*,*) '-----'
c      write(*,*) 'sss ',46487200.d0*sss(nb,oo,aa,gg)
       write(*,*) 'sss ',sss(nb,oo,aa,gg)
       c1=dimag(sss(nb,oo,aa,gg))
       write(*,*) 'c1',c1
c      write(*,*) 'corr ',time,corr(nb,aa,gg,time)
c      write(*,*) pi*46487200.d0*oo**3*dexp(-oo**3/wc**3)
c      write(*,*) pi*46487200.d0*oo**3*dexp(-oo**3/wc**3)
c    .        *(dexp(0.5d0*beta*oo)+dexp(-0.5d0*beta*oo))
c    .        /(dexp(0.5d0*beta*oo)-dexp(-0.5d0*beta*oo))
       write(*,*) 'jjc',pi/2*46487200.d0*jjc(dcmplx(oo,0.d0),wc)
c      write(*,*) pi/2*46487200.d0*oo**3*dexp(-oo**2/wc**2)

       write(*,*) '-----'
c      write(*,*) 'ccc ',wc,beta,46487200.d0*ccc(nb,oo,aa,gg)
       write(*,*) 'ccc ',wc,beta,ccc(nb,oo,aa,gg)
       write(*,*) 'jjc coth',pi/2*46487200.d0*jjc(dcmplx(oo,0.d0),wc)
     .        *(dexp(0.5d0*beta*oo)+dexp(-0.5d0*beta*oo))
     .        /(dexp(0.5d0*beta*oo)-dexp(-0.5d0*beta*oo))
c      write(*,*) pi/2*46487200.d0*oo**3*dexp(-oo**2/wc**2)
c    .        *(dexp(0.5d0*beta*oo)+dexp(-0.5d0*beta*oo))
c    .        /(dexp(0.5d0*beta*oo)-dexp(-0.5d0*beta*oo))

c      stop



c________________________

      !$omp do
       do kk=1,nruns


c______here we scan !!!!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       efield(kk)=ep*dble(kk)/dble(nruns)  ! for multiple runs
       pulse_par(1)=efield(kk)

       write(*,*) 'runs',kk,omp_get_thread_num(),efield(kk)

c______for pulse lengtyh scaw w/o touching anything else !!
c      pulse_par(3)=taup*1.5d0*dble(kk)/dble(nruns)
c      scanvar(kk)=pulse_par(3)



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c______initial state
       do i=1,n
       do j=1,n
       zs(i,j)=dcmplx(0.d0,0.d0)
       do k=1,nb
       zz(i,j,k)=dcmplx(0.d0,0.d0)
       enddo
       enddo
       enddo
       zz(1,1,1)=dcmplx(1.d0,0.d0)    ! initial state is ground state
       zs(1,1)  =dcmplx(1.d0,0.d0)                 ! initial state is ground state

c______control: add seeding
c      zz(1,1,1)=0.999d0    ! initial state is ground state
c      zz(2,1,1)=0.001d0    ! initial state is ground state

c      ntime=6.d0*taup/timestep+1
       ntime=6.d0*pulse_par(3)/timestep+1





c______main time loop
       do i=1,ntime

c      time=-3.d0*taup+dble(i-1)*timestep
       time=-3.d0*pulse_par(3)+dble(i-1)*timestep

       tend=time+timestep

       call PROP_CH(time,tend,zz,pulse_par)   ! uses rk_ch   

c       if (kk.eq.nruns) write(67,88) time,dreal(zz(1,1,1)),
c    .dreal(zz(2,2,1)),dreal(zz(3,3,1)),pulse_par(4)

       enddo  ! end of time loop

c       close(67)

       signal(kk,1)=dreal(zz(1,1,1))
       signal(kk,2)=dreal(zz(2,2,1))
       signal(kk,3)=dreal(zz(4,4,1))

       enddo    ! end of kk loop
       !$omp end do
       !$omp end parallel




c______output 
      open(45,file=prefix//'.result')
      do kk=1,nruns
c     write(45,*) efield(kk)*dsqrt(taup),signal(kk,1),signal(kk,2)
      write(45,89) efield(kk)*dsqrt(taup),signal(kk,1),signal(kk,2),
     .signal(kk,3)
c     write(45,*) 0.024189*scanvar(kk),signal(kk,1),signal(kk,2)
89     format(11F28.14)

      enddo
      close(45)

88     format(11F28.14)

998    continue

       end


       function sss(nb,oo,aa,gg)
       implicit none

       integer    nb,i
       real*8     oo
       complex*16 sum,sss,ig,aa(nb),gg(nb)

       sum=(0.d0,0.d0)
       do i=2,nb
       ig=-(0.d0,1.d0)*gg(i)
       sum=sum+aa(i)*oo/(oo**2+ig**2)
       enddo

       sss=sum

       return
       end


       function ccc(nb,oo,aa,gg)
       implicit none

       integer    nb,i
       real*8     oo
       complex*16 sum,ccc,ig,aa(nb),gg(nb)

       sum=(0.d0,0.d0)
       do i=2,nb
       ig=-(0.d0,1.d0)*gg(i)
       sum=sum+aa(i)*ig/(oo**2+ig**2)
       enddo

       ccc=sum

       return
       end








       function corr(nb,aa,gg,time)
       implicit none

       integer nb,i

       real*8     time,beta

       complex*16 sum,corr,aa(nb),gg(nb)

       complex*16 one,ione,zero
       parameter (one= (1.d0,0.d0))
       parameter (ione=(0.d0,1.d0))
       parameter (zero=(0.d0,0.d0))

       sum=(0.0,0.0)
       do i=2,nb
       sum=sum+aa(i)*cdexp(ione*gg(i)*time)
       enddo

       corr=sum

       return
       end










      subroutine init_coup(n,xx)
      implicit none

      integer    i,j,n
      complex*16 xx(n,n)

c_____coupling
      do i=1,n
      do j=1,n
      xx(i,j)=dcmplx(0.d0,0.d0)
      enddo
      enddo
      xx(2,2)=dcmplx(1.d0,0.d0)
      xx(4,4)=dcmplx(2.d0,0.d0)   ! four


      return
      end









      subroutine init_bath(nb,wc,beta,aa,gg)
      implicit none

      real*8     pi,twopi
c     parameter (pi=3.14159265358979334d0)
c     parameter (twopi=6.28318530717958668d0)

      complex*16 one,ione,zero
      parameter (one= (1.d0,0.d0))
      parameter (ione=(0.d0,1.d0))
      parameter (zero=(0.d0,0.d0))

      integer    i,nb,firstmatsu
      complex*16 aa(nb),gg(nb)
      complex*16 zz1,zc1,zz2,zc2,jjc,jjcxx,nbeta,xx,c1
      real*8     beta,wc,nu_k

      include 'poles.f'

      pi=dacos(-1.d0)
      twopi=2.d0*pi

      aa(1)=dcmplx(0.d0,0.d0)  ! for system density matrix
      gg(1)=dcmplx(0.d0,0.d0)  ! for system density matrix

      c1=(0.d0,0.d0)

      do i=1,nfit

      zz1=z1(i)
      zc1=dconjg(z1(i))
      zz2=z2(i)
      zc2=dconjg(z2(i))

c     write(*,*) i,zz1,zz2

      gg((i-1)*nfit+2)= wc*zz1
      gg((i-1)*nfit+3)=-wc*zc1
      gg((i-1)*nfit+4)= wc*zz2
      gg((i-1)*nfit+5)=-wc*zc2

      aa((i-1)*nfit+2)=c(i)*nbeta(beta,wc*zz1) *ione*pi*wc**4*zz1**2
     .                /(zz1**2-zz2**2)/(zz1**2-zc1**2)/(zz1**2-zc2**2)

      aa((i-1)*nfit+3)=c(i)*nbeta(beta,-wc*zc1)*ione*pi*wc**4*zc1**2
     .                /(zc1**2-zz1**2)/(zc1**2-zz2**2)/(zc1**2-zc2**2)

      aa((i-1)*nfit+4)=c(i)*nbeta(beta,wc*zz2) *ione*pi*wc**4*zz2**2
     .                /(zz2**2-zz1**2)/(zz2**2-zc1**2)/(zz2**2-zc2**2)

      aa((i-1)*nfit+5)=c(i)*nbeta(beta,-wc*zc2)*ione*pi*wc**4*zc2**2
     .                /(zc2**2-zz1**2)/(zc2**2-zz2**2)/(zc2**2-zc1**2)

      enddo

      firstmatsu=4*nfit+2

      do i=firstmatsu,nb
        nu_k=twopi/beta*dble(i-firstmatsu+1)
        xx=dcmplx(0.d0,nu_k/wc)
        gg(i)=dcmplx(0.d0,nu_k)
        aa(i)=dcmplx(0.d0,twopi)*wc**3/beta*jjcxx(xx)
      enddo

      return
      end



       subroutine check_bath(iunit,nb,wc,beta,lambda,time,aa,gg)
       implicit none

       integer nb,i,iunit
       real*8  wc,beta,lambda,dw,time
       complex*16 c0,c1,w,x,corr

c_____declarations for bath common block
      complex*16 aa(nb),gg(nb)
      complex*16 nbeta,jjc,jjcxx

      write(*,*) ' checking bath....',beta,wc,time 

       c0=dcmplx(0.d0,0.d0)
       c1=dcmplx(0.d0,0.d0)

       do i=-4400000,4400000
       if (i.eq.0) goto 899
       w=wc*dcmplx(dble(i)*0.00001d0,0.d0)
       dw=wc*0.00001d0

c      if (iunit.ne.0) 
c    . write(iunit,*) dreal(w),dreal(nbeta(beta,w)*jjc(w,wc))
c    . write(iunit,*) dreal(w),dreal(jjc(w,wc))

       x=w/wc

       c0=c0+nbeta(beta,w)*lambda*w**3*cdexp(-w**2/wc**2)
     .    *cdexp(w*dcmplx(0.d0,time))*dw
       c1=c1+nbeta(beta,w)*lambda*wc**3*jjcxx(x)
     .    *cdexp(w*dcmplx(0.d0,time))*dw

899    continue
       enddo

       write(*,*) time,dreal(c0),dimag(c0)
       write(*,*) time,dreal(c1),dimag(c1)
       write(*,*) time,dreal(corr(nb,aa,gg,time))
     .                 ,dimag(corr(nb,aa,gg,time))

       return
       end



      function jjcxx(xx)
      implicit  none

      integer    k
      complex*16 xx,jjcxx

      include 'poles.f' 

      jjcxx=dcmplx(0.d0,0.d0)
      do k=1,nfit
      jjcxx=jjcxx+c(k)*xx**3
     . /(xx-z1(k))/(xx+z1(k))/(xx-dconjg(z1(k)))/(xx+dconjg(z1(k)))
     . /(xx-z2(k))/(xx+z2(k))/(xx-dconjg(z2(k)))/(xx+dconjg(z2(k)))
      enddo

      return
      end


      function jjc(w,wc)
      implicit  none

      integer    k
      real*8     wc
      complex*16 w,xx,jjc

      include 'poles.f' 

      jjc=dcmplx(0.d0,0.d0)
      do k=1,nfit
      xx=w/wc
      jjc=jjc+c(k)*wc**3*xx**3
     . /(xx-z1(k))/(xx+z1(k))/(xx-dconjg(z1(k)))/(xx+dconjg(z1(k)))
     . /(xx-z2(k))/(xx+z2(k))/(xx-dconjg(z2(k)))/(xx+dconjg(z2(k)))
      enddo

      return
      end


      function nbeta(beta,w)
      implicit none
      real*8     beta
      complex*16 nbeta,w

c     if (dreal(w).lt.0.d0) then
c     nbeta=1.d0/(cdexp(beta*w)-1.d0)
c     else
      nbeta=cdexp(-beta*w)/(1.d0-cdexp(-beta*w))
c     endif


      return
      end
      


      subroutine get_pulse(time,fr,fi,pulse_par)
      implicit none

      real*8     time,e,omega,tau,alpha,delta1,delta2,pulse_par(6)
      real*8     fr,fi
      complex*16 f,comm  ! for later use: control

      complex*16 one,ione,zero
      parameter (one= (1.d0,0.d0))
      parameter (ione=(0.d0,1.d0))
      parameter (zero=(0.d0,0.d0))

      e=    pulse_par(1)
      omega=pulse_par(2)
      tau=  pulse_par(3)
      alpha=pulse_par(4)
      delta1=pulse_par(5)
      delta2=pulse_par(6)

c_____here is the chirped pulse
      if ((time.gt.-3.d0*tau).and.(time.lt.3.d0*tau)) then
         f=e*dexp(-(time/tau)**2)*cdexp(ione*(omega*time+alpha*time**2))  ! chirp
      else
         f=dcmplx(0.d0,0.d0)
      endif

      fr=dreal(f)
      fi=dimag(f)

c_____here we can add the control pulse

      return
      end


      subroutine get_ham(time,hh,pulse_par)
      implicit none

      include 'sizes.f'

      real*8     time,e,omega,tau,alpha,delta1,delta2,pulse_par(6)
      real*8     fr,fi
      complex*16 hh(n,n),f,comm  ! for later use: control

      complex*16 one,ione,zero
      parameter (one= (1.d0,0.d0))
      parameter (ione=(0.d0,1.d0))
      parameter (zero=(0.d0,0.d0))

      e=    pulse_par(1)
      omega=pulse_par(2)
      tau=  pulse_par(3)
      alpha=pulse_par(4)
      delta1=pulse_par(5)
      delta2=pulse_par(6)

      hh=dcmplx(0.d0,0.d0)

c________RWA_______________
      hh(1,1)= (omega+2.d0*alpha*time-delta1)
      hh(2,2)= (0.d0,0.d0)     
      hh(3,3)= (0.d0,0.d0)     
      hh(4,4)=-(omega+2.d0*alpha*time-delta2)

      hh(1,2)= e/2.d0*dexp(-(time/tau)**2)
      hh(2,1)= e/2.d0*dexp(-(time/tau)**2)
      hh(2,4)= e/2.d0*dexp(-(time/tau)**2)
      hh(4,2)= e/2.d0*dexp(-(time/tau)**2)

c_____no RWA_______________
c     hh(1,1)= dcmplx(-delta,0.d0)     
c     hh(1,2)= dcmplx(e*dexp(-(time/tau)**2)
c    .  *dcos(omega*time+alpha*time**2),0.d0)
c     hh(2,1)= dcmplx(e*dexp(-(time/tau)**2)
c    .  *dcos(omega*time+alpha*time**2),0.d0)
c     hh(2,2)= dcmplx(0.d0,0.d0)     

c_____control______________
c     use pulse_par(1) as field scaling, and pulse_par(4) as field
c     hh(1,1)= (0.d0,0.d0)     
c     hh(1,2)= e*dexp(-(time/tau)**2)*alpha
c     hh(2,1)= e*dexp(-(time/tau)**2)*alpha
c     hh(2,2)= delta           

      return
      end


      subroutine PROP_MARKOV(t,tend,y,RPAR)
c     use omp_lib
      implicit none

      include 'sizes.f'

      integer     n2
      parameter (n2=2*n*n)

c_____prop stuff
      integer     I,LWORK,LIWORK,IDID
      parameter (LWORK=8*n2+21,LIWORK=21)
      integer     ITOL,IOUT,IWORK(LIWORK)
      real*8      RTOL,ATOL,WORK(LWORK)
      real*8      t,tend,y(n2)

c_____parameter passing
      integer     IPAR(5)
      real*8      RPAR(6)

      external    G,SOLOUT

      do i=1,5
      ipar(i)=0
      enddo

C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
c       RTOL=1.0D-13
        RTOL=1.0D-8
        ATOL=RTOL
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.D0
C --- OUTPUT
        IOUT=0  ! no call to SOLOUT


C --- CALL OF THE SUBROUTINE DOPRI5
        CALL DOPRI5(n2,G,t,y,tend,
     &              RTOL,ATOL,ITOL,
     &              SOLOUT,IOUT,
     &              WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

c       CALL DOPRI5(N,FAREN,X,Y,XEND,
c    &              RTOL,ATOL,ITOL,
c    &              SOLOUT,IOUT,
c    &              WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)


      return
      end



      subroutine PROP_CH(t,tend,y,RPAR)
      use omp_lib
      implicit none

      include 'sizes.f'

c_____prop stuff
      integer     I,LWORK,LIWORK,IDID
      parameter (LWORK=8*bign+21,LIWORK=21)
      integer     ITOL,IOUT,IWORK(LIWORK)
      real*8      RTOL,ATOL,WORK(LWORK)
      real*8      t,tend,y(bign)

c_____parameter passing
      integer     IPAR(5)
      real*8      RPAR(6)

      external    FCH,SOLOUT    

      do i=1,5
      ipar(i)=0
      enddo

C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
c       RTOL=1.0D-12
        RTOL=1.0D-10
        ATOL=RTOL
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.D0
C --- OUTPUT
        IOUT=0  ! no call to SOLOUT


       
C --- CALL OF THE SUBROUTINE DOPRI5   
        CALL DOPRI5(bign,FCH,t,y,tend,
     &              RTOL,ATOL,ITOL,
     &              SOLOUT,IOUT,
     &              WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

c       CALL DOPRI5(N,FAREN,X,Y,XEND,
c    &              RTOL,ATOL,ITOL,
c    &              SOLOUT,IOUT,
c    &              WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)


      return
      end


      subroutine FCH(ndummy,time,zz,zp,rpar,ipar)  ! 
      use omp_lib
      implicit none

      include 'sizes.f'

       integer    ndummy
       integer    ipar(5)
       real*8     rpar(6)

c______declarations for bath common block
       complex*16 aa(nb),gg(nb)

       complex*16 xx(n,n),hh(n,n)

c______declarations for dmat common block
       complex*16 zz(n,n,nb),zp(n,n,nb)
       complex*16 zh(n,n)

       integer    o,m,k,p
       real*8     time,beta
       real*8     fr,fi

       complex*16 one,ione,zero,carg,cc
       parameter (one= (1.d0,0.d0))
       parameter (ione=(0.d0,1.d0))
       parameter (zero=(0.d0,0.d0))

       common/coup/xx
       common/bath/beta,aa,gg
       !$OMP THREADPRIVATE(/coup/,/bath/)
   

c_____control
c     rpar(4)=dimag(zz(1,2,1))

      call get_ham(time,hh,rpar)


c      write(*,*) 'xx',xx

c============== aux density matrices

       do k=2,nb

       do o=1,n
       do m=1,n

       cc=ione*gg(k)*zz(o,m,k)
       
       do p=1,n
       cc=cc
     . -ione*(hh(o,p)*zz(p,m,k)-zz(o,p,k)*hh(p,m))
     . +ione*aa(k)*xx(o,p)*zz(p,m,1)
c    . -ione*aa(k)*zz(o,p,1)*xx(p,m)*cdexp(gg(k)*beta)     ! varaint (46)
        enddo
       
        zp(o,m,k)=cc

        enddo
        enddo

       enddo   ! end of k-loop for aux density matrices


c====================== system part


c_____initializing  to zero
       do o=1,n
       do m=1,n
       zp(o,m,1)=dcmplx(0.d0,0.d0)
       enddo
       enddo

c______Hamiltonian part_________________________________
       do o=1,n
       do m=1,n
       do p=1,n

       zp(o,m,1)=zp(o,m,1)-ione*(hh(o,p)*zz(p,m,1)-zz(o,p,1)*hh(p,m))

       enddo
       enddo
       enddo

c______dissipative part: summing up aux for system 
       do k=2,nb

       do o=1,n
       do m=1,n
       zh(o,m)=zz(o,m,k)+dconjg(zz(m,o,k))    ! variant (41)
c      zh(o,m)=zz(o,m,k)                      ! variant (46)
       enddo
       enddo

       do o=1,n
       do m=1,n
       do p=1,n
       zp(o,m,1)=zp(o,m,1)+ione*(xx(o,p)*zh(p,m)-zh(o,p)*xx(p,m))
       enddo
       enddo
       enddo

       enddo   ! end of k-loop


       return
       end




      subroutine G(ndummy,time,zz,zp,pulse_par,ipar)
      use omp_lib
      implicit none

      include 'sizes.f'

      integer ndummy
      integer ipar(5)
c     real*8  rpar(5)

      real*8     pi,twopi
      parameter (pi=3.14159265358979334d0)
      parameter (twopi=6.28318530717958668d0)

c______declarations for bath common block
       complex*16 aa(nb),gg(nb)

c______declarations for syst common block
       complex*16 xx(n,n),dd(n,n),dc(n,n),hh(n,n)

c______aux variables for dissipative operator
       real*8     e,omega,tau,alpha,delta
       real*8     dddd,om,oo
       complex*16 ccc,sss,cc0,cc1,ss1         

c______declarations for dmat common block
       complex*16 zz(n,n),zp(n,n),zh(n,n)
       complex*16 coco,brr,bcc

c_____declarations for laser common block
       real*8     pulse_par(5)

       integer    o,m,k,p,pp
       real*8     time,beta
       real*8     fr,fi


c______variables for ramsay comparison
       real*8     ram_a,ram_b
       real*8     lambda,wc
       complex*16 jjc


       complex*16 one,ione,zero,carg,cc
       parameter (one= (1.d0,0.d0))
       parameter (ione=(0.d0,1.d0))
       parameter (zero=(0.d0,0.d0))

       common/coup/xx
       common/bath/beta,aa,gg
       !$OMP THREADPRIVATE(/coup/,/bath/)

c_____create dissipative operator dd

      e=    pulse_par(1)
      omega=pulse_par(2)
      tau=  pulse_par(3)
      alpha=pulse_par(4)
      delta=pulse_par(5)

c_____no RWA
c     hh(2,1)= e*dexp(-(time/tau)**2)*dcos(omega*time+alpha*time**2)
c     hh(2,2)= delta

      dddd=-delta+omega-2.d0*alpha*time
c     dddd=0.d0

      om=  e*dexp(-(time/tau)**2)
      oo=  dsqrt(om**2+dddd**2)



c_____RWA
      hh(1,1)= (0.d0,0.d0)
      hh(1,2)= om/2.d0
      hh(2,1)= om/2.d0
      hh(2,2)= dddd 

c_____no RWA
c     hh(1,1)= (0.d0,0.d0)
c     hh(2,1)= e*dexp(-(time/tau)**2)*dcos(omega*time+alpha*time**2)
c     hh(1,2)= e*dexp(-(time/tau)**2)*dcos(omega*time+alpha*time**2)
c     hh(2,2)= delta

      cc0=dcmplx(dreal(ccc(nb,0.d0,aa,gg)),0.d0)
      cc1=dcmplx(dreal(ccc(nb,oo,aa,gg)),0.d0)
      ss1=dcmplx(0.d0,dimag(sss(nb,oo,aa,gg)))

c     cc0=1.d0*(ccc(nb,0.d0,aa,gg))
c     cc1=1.d0*(ccc(nb,oo,aa,gg))
c     ss1=1.d0*(sss(nb,oo,aa,gg))

      dd(1,1)=  0.5d0*om**2/oo**2*(cc0-cc1)

      dd(1,2)=  0.5d0*(om*dddd/oo**2*(cc0-cc1)-ione*om/oo*ss1)
      dd(2,1)=  0.5d0*(om*dddd/oo**2*(cc0-cc1)+ione*om/oo*ss1)

      dd(2,2)= -0.5d0*om**2/oo**2*(cc0-cc1)

      dc(1,1)=dconjg(dd(1,1))
      dc(2,1)=dconjg(dd(1,2))
      dc(1,2)=dconjg(dd(2,1))
      dc(2,2)=dconjg(dd(2,2))

c--------------
       do o=1,n
       do m=1,n

       cc=(0.d0,0.d0)

       do p=1,n
       cc=cc-ione*(hh(o,p)*zz(p,m)-zz(o,p)*hh(p,m))  
       enddo

       do pp=1,n
       do p=1,n
       cc=cc
     . -xx(o,p)*dd(p,pp)*zz(pp,m)  ! dissi
     . +xx(o,p)*zz(p,pp)*dc(pp,m)  ! dissi
     . +dd(o,p)*zz(p,pp)*xx(pp,m)  ! dissi
     . -zz(o,p)*dc(p,pp)*xx(pp,m)  ! dissi
       enddo
       enddo


       zp(o,m)=cc

       enddo
       enddo


c______comparison with ramsay

c     lambda=46487200.d0
c     wc=0.00144d0/27.211d0

c     ram_a=pi/8.d0*46487200.d0*jjc(dcmplx(oo,0.d0),wc)
c     ram_b=pi/8.d0*46487200.d0*jjc(dcmplx(oo,0.d0),wc)
c    .        *(dexp(0.5d0*beta*oo)+dexp(-0.5d0*beta*oo))
c    .        /(dexp(0.5d0*beta*oo)-dexp(-0.5d0*beta*oo))


c      zp(1,2)=zp(1,2)-ram_b*4.d0*zz(1,2)
c      zp(2,1)=zp(2,1)-ram_b*4.d0*zz(2,1)

c      zp(1,2)=zp(1,2)-ram_b*4.d0*zz(1,2)-2.d0*ram_a*(zz(1,1)+zz(2,2))
c      zp(2,1)=zp(2,1)-ram_b*4.d0*zz(2,1)-2.d0*ram_a*(zz(1,1)+zz(2,2))

       return
       end

