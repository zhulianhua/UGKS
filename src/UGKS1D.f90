!-----------------------------------------------------------------------------
!UGKS1D - shock structure calculation using Unified Gas-Kinetic Scheme
!-----------------------------------------------------------------------------
!Copyright (C) 2012 Ruijie Wang <lainme993@gmail.com>
!
!This program is free software: you can redistribute it and/or modify it under
!the terms of the GNU General Public License as published by the Free Software
!Foundation, either version 3 of the License, or (at your option) any later
!version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT
!ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
!FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
!details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

!--------------------------------------------------
!>store the global variables
!--------------------------------------------------
module global_data
    !--------------------------------------------------
    !kind selection
    !--------------------------------------------------
    integer,parameter :: RKD = 8

    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(kind=RKD),parameter :: PI = 4.0*atan(1.0) !Pi
    real(kind=RKD),parameter :: SMV = tiny(real(1.0,8)) !small value to avoid 0/0
    real(kind=RKD),parameter :: UP = 1.0 !used in sign() function
    real(kind=RKD) :: cfl !global CFL number
    real(kind=RKD) :: dt !global time step
    real(kind=RKD) :: sim_time !current simulation time
    real(kind=RKD) :: max_time !maximum simulation time
    character(len=9),parameter :: HSTFILENAME = "shock.hst" !history file name
    character(len=9),parameter :: RSTFILENAME = "shock.dat" !result file name
    integer :: iter !iteration
    integer :: method_output !output the solution with normalized value or not

    !--------------------------------------------------
    !gas properties
    !--------------------------------------------------
    real(kind=RKD) :: gam !ratio of specific heat
    real(kind=RKD) :: omega !temperature dependence index in HS/VHS/VSS model
    real(kind=RKD) :: pr !Prandtl number
    real(kind=RKD) :: mu_ref !viscosity coefficient in reference state
    integer :: ck !internal degree of freedom

    !--------------------------------------------------
    !macros for a readable code
    !--------------------------------------------------
    !I/O
    integer,parameter :: HSTFILE = 20 !history file ID
    integer,parameter :: RSTFILE = 21 !result file ID
    !output method
    integer,parameter :: ORIGINAL = 0 !do not normalize the solution
    integer,parameter :: NORMALIZE = 1 !normalize the solution

    !--------------------------------------------------
    !basic derived type
    !--------------------------------------------------
    !cell center
    type :: cell_center
        !geometry
        real(kind=RKD) :: x !cell center coordinates
        real(kind=RKD) :: length !length
        !flow field
        real(kind=RKD) :: w(3) !density, x-momentum,total energy
        real(kind=RKD) :: primary(2) !density,temperature
        real(kind=RKD),allocatable,dimension(:) :: h,b !distribution function
        real(kind=RKD),allocatable,dimension(:) :: sh,sb !slope of distribution function
    end type cell_center

    !cell interface
    type :: cell_interface
        real(kind=RKD) :: flux(3) !mass flux, x momentum flux, energy flux
        real(kind=RKD),allocatable,dimension(:) :: flux_h,flux_b !flux of distribution function
    end type cell_interface

    !--------------------------------------------------
    !flow field
    !--------------------------------------------------
    !index method
    !     ----------------
    !  (i)|      (i)     |(i+1)
    !     ----------------
    integer :: ixmin,ixmax !index range
    type(cell_center),allocatable,dimension(:) :: ctr !cell centers
    type(cell_interface),allocatable,dimension(:) :: vface !vertical interfaces

    !--------------------------------------------------
    !discrete velocity space
    !--------------------------------------------------
    integer :: unum !number of velocity points
    real(kind=RKD) :: umax !maximum micro velocity
    real(kind=RKD),allocatable,dimension(:) :: uspace !discrete velocity space
    real(kind=RKD),allocatable,dimension(:) :: weight !weight at velocity u_k
end module global_data

!--------------------------------------------------
!>define some commonly used functions/subroutines
!--------------------------------------------------
module tools
    use global_data
    implicit none
    contains
        !--------------------------------------------------
        !>obtain discretized Maxwellian distribution
        !>@param[out] h,b  :distribution function
        !>@param[in]  prim :primary variables
        !--------------------------------------------------
        subroutine discrete_maxwell(h,b,prim)
            real(kind=RKD),dimension(:),intent(out) :: h,b
            real(kind=RKD),intent(in) :: prim(3)

            h = prim(1)*(prim(3)/PI)**(1.0/2.0)*exp(-prim(3)*(uspace-prim(2))**2)
            b = h*ck/(2.0*prim(3))
        end subroutine discrete_maxwell

        !--------------------------------------------------
        !>calculate the Shakhov part H^+, B^+
        !>@param[in]  H,B           :Maxwellian distribution function
        !>@param[in]  qf            :heat flux
        !>@param[in]  prim          :primary variables
        !>@param[out] H_plus,B_plus :Shakhov part
        !--------------------------------------------------
        subroutine shakhov_part(H,B,qf,prim,H_plus,B_plus)
            real(kind=RKD),dimension(:),intent(in) :: H,B
            real(kind=RKD),intent(in) :: qf
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD),dimension(:),intent(out) :: H_plus,B_plus

            H_plus = 0.8*(1-pr)*prim(3)**2/prim(1)*&
                     (uspace-prim(2))*qf*(2*prim(3)*(uspace-prim(2))**2+ck-5)*H
            B_plus = 0.8*(1-pr)*prim(3)**2/prim(1)*&
                     (uspace-prim(2))*qf*(2*prim(3)*(uspace-prim(2))**2+ck-3)*B
        end subroutine shakhov_part

        !--------------------------------------------------
        !>convert primary variables to conservative variables
        !>@param[in] prim          :primary variables
        !>@return    get_conserved :conservative variables
        !--------------------------------------------------
        function get_conserved(prim)
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD) :: get_conserved(3)

            get_conserved(1) = prim(1)
            get_conserved(2) = prim(1)*prim(2)
            get_conserved(3) = 0.5*prim(1)/prim(3)/(gam-1.0)+0.5*prim(1)*prim(2)**2
        end function get_conserved

        !--------------------------------------------------
        !>convert conservative variables to primary variables
        !>@param[in] w           :conservative variables
        !>@return    get_primary :conservative variables
        !--------------------------------------------------
        function get_primary(w)
            real(kind=RKD),intent(in) :: w(3)
            real(kind=RKD) :: get_primary(3) !primary variables

            get_primary(1) = w(1)
            get_primary(2) = w(2)/w(1)
            get_primary(3) = 0.5*w(1)/(gam-1.0)/(w(3)-0.5*w(2)**2/w(1))
        end function get_primary

        !--------------------------------------------------
        !>obtain ratio of specific heat
        !>@param[in] ck        :internal degree of freedom
        !>@return    get_gamma :ratio of specific heat
        !--------------------------------------------------
        function get_gamma(ck)
            integer,intent(in) :: ck
            real(kind=RKD) :: get_gamma

            get_gamma = float(ck+3)/float(ck+1)
        end function get_gamma

        !--------------------------------------------------
        !>obtain speed of sound
        !>@param[in] prim    :primary variables
        !>@return    get_sos :speed of sound
        !--------------------------------------------------
        function get_sos(prim)
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD) :: get_sos !speed of sound

            get_sos = sqrt(0.5*gam/prim(3))
        end function get_sos
        
        !--------------------------------------------------
        !>calculate collision time
        !>@param[in] prim    :primary variables
        !>@return    get_tau :collision time
        !--------------------------------------------------
        function get_tau(prim)
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD) :: get_tau

            get_tau = mu_ref*2*prim(3)**(1-omega)/prim(1)
        end function get_tau

        !--------------------------------------------------
        !>get heat flux
        !>@param[in] h,b           :distribution function
        !>@param[in] prim          :primary variables
        !>@return    get_heat_flux :heat flux in normal and tangential direction
        !--------------------------------------------------
        function get_heat_flux(h,b,prim)
            real(kind=RKD),dimension(:),intent(in) :: h,b
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD) :: get_heat_flux !heat flux in normal and tangential direction

            get_heat_flux = 0.5*(sum(weight*(uspace-prim(2))*(uspace-prim(2))**2*h)+sum(weight*(uspace-prim(2))*b)) 
        end function get_heat_flux

        !--------------------------------------------------
        !>get temperature
        !>@param[in] h,b             :distribution function
        !>@param[in] prim            :primary variables
        !>@return    get_temperature :temperature
        !--------------------------------------------------
        function get_temperature(h,b,prim)
            real(kind=RKD),dimension(:),intent(in) :: h,b
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD) :: get_temperature !pressure

            get_temperature = 2.0*(sum(weight*(uspace-prim(2))**2*h)+sum(weight*b))/(ck+1)/prim(1)
        end function get_temperature

        !--------------------------------------------------
        !>get the nondimensionalized viscosity coefficient
        !>@param[in] kn          :Knudsen number
        !>@param[in] alpha,omega :index related to HS/VHS/VSS model
        !>@return    get_mu      :nondimensionalized viscosity coefficient
        !--------------------------------------------------
        function get_mu(kn,alpha,omega)
            real(kind=RKD),intent(in) :: kn,alpha,omega
            real(kind=RKD) :: get_mu

            get_mu = 5*(alpha+1)*(alpha+2)*sqrt(PI)/(4*alpha*(5-2*omega)*(7-2*omega))*kn
        end function get_mu
end module tools

!--------------------------------------------------
!>flux calculation
!--------------------------------------------------
module flux
    use global_data
    use tools
    implicit none

    integer,parameter :: MNUM = 6 !number of normal velocity moments
    integer,parameter :: MTUM = 4 !number of tangential velocity moments

    contains
        !--------------------------------------------------
        !>calculate flux of inner interface
        !>@param[in]    cell_L :cell left to the target interface
        !>@param[inout] face   :the target interface
        !>@param[in]    cell_R :cell right to the target interface
        !--------------------------------------------------
        subroutine calc_flux(cell_L,face,cell_R)
            type(cell_center),intent(in) :: cell_L,cell_R
            type(cell_interface),intent(inout) :: face
            real(kind=RKD),allocatable,dimension(:) :: h,b !distribution function at the interface
            real(kind=RKD),allocatable,dimension(:) :: H0,B0 !Maxwellian distribution function
            real(kind=RKD),allocatable,dimension(:) :: H_plus,B_plus !Shakhov part of the equilibrium distribution
            real(kind=RKD),allocatable,dimension(:) :: sh,sb !slope of distribution function at the interface
            integer,allocatable,dimension(:) :: delta !Heaviside step function
            real(kind=RKD) :: w(3),prim(3) !conservative and primary variables at the interface
            real(kind=RKD) :: qf !heat flux in normal and tangential direction
            real(kind=RKD) :: sw(3) !slope of W
            real(kind=RKD) :: aL(3),aR(3),aT(3) !micro slope of Maxwellian distribution, left,right and time.
            real(kind=RKD) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<\xi^l>
            real(kind=RKD) :: Mau_0(3),Mau_L(3),Mau_R(3),Mau_T(3) !<u\psi>,<aL*u^n*\psi>,<aR*u^n*\psi>,<A*u*\psi>
            real(kind=RKD) :: tau !collision time
            real(kind=RKD) :: Mt(5) !some time integration terms
            integer :: i,j

            !--------------------------------------------------
            !prepare
            !--------------------------------------------------
            !allocate array
            allocate(delta(unum))
            allocate(h(unum))
            allocate(b(unum))
            allocate(sh(unum))
            allocate(sb(unum))
            allocate(H0(unum))
            allocate(B0(unum))
            allocate(H_plus(unum))
            allocate(B_plus(unum))

            !Heaviside step function
            delta = (sign(UP,uspace)+1)/2

            !--------------------------------------------------
            !reconstruct initial distribution
            !--------------------------------------------------
            h = (cell_L%h+0.5*cell_L%length*cell_L%sh)*delta+&
                (cell_R%h-0.5*cell_R%length*cell_R%sh)*(1-delta)
            b = (cell_L%b+0.5*cell_L%length*cell_L%sb)*delta+&
                (cell_R%b-0.5*cell_R%length*cell_R%sb)*(1-delta)
            sh = cell_L%sh*delta+cell_R%sh*(1-delta)
            sb = cell_L%sb*delta+cell_R%sb*(1-delta)

            !--------------------------------------------------
            !obtain macroscopic variables
            !--------------------------------------------------
            !conservative variables w_0
            w(1) = sum(weight*h)
            w(2) = sum(weight*uspace*h)
            w(3) = 0.5*(sum(weight*uspace**2*h)+sum(weight*b))

            !convert to primary variables
            prim = get_primary(w)

            !heat flux
            qf = get_heat_flux(h,b,prim) 

            !--------------------------------------------------
            !calculate a^L,a^R
            !--------------------------------------------------
            sw = (w-cell_L%w)/(0.5*cell_L%length) !left slope of W
            aL = micro_slope(prim,sw) !calculate a^L

            sw = (cell_R%w-w)/(0.5*cell_R%length) !right slope of W
            aR = micro_slope(prim,sw) !calculate a^R

            !--------------------------------------------------
            !calculate time slope of W and A
            !--------------------------------------------------
            !<u^n>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
            call calc_moment_u(prim,Mu,Mxi,Mu_L,Mu_R) 

            Mau_L = moment_au(aL,Mu_L,Mxi,1) !<aL*u*\psi>_{>0}
            Mau_R = moment_au(aR,Mu_R,Mxi,1) !<aR*u*\psi>_{<0}

            sw = -prim(1)*(Mau_L+Mau_R) !time slope of W
            aT = micro_slope(prim,sw) !calculate A

            !--------------------------------------------------
            !calculate collision time and some time integration terms
            !--------------------------------------------------
            tau = get_tau(prim)

            Mt(4) = tau*(1.0-exp(-dt/tau))
            Mt(5) = -tau*dt*exp(-dt/tau)+tau*Mt(4)
            Mt(1) = dt-Mt(4)
            Mt(2) = -tau*Mt(1)+Mt(5) 
            Mt(3) = dt**2/2.0-tau*Mt(1)

            !--------------------------------------------------
            !calculate the flux of conservative variables related to g0
            !--------------------------------------------------
            Mau_0 = moment_uv(Mu,Mxi,1,0) !<u*\psi>
            Mau_L = moment_au(aL,Mu_L,Mxi,2) !<aL*u^2*\psi>_{>0}
            Mau_R = moment_au(aR,Mu_R,Mxi,2) !<aR*u^2*\psi>_{<0}
            Mau_T = moment_au(aT,Mu,Mxi,1) !<A*u*\psi>

            face%flux = Mt(1)*prim(1)*Mau_0+Mt(2)*prim(1)*(Mau_L+Mau_R)+Mt(3)*prim(1)*Mau_T

            !--------------------------------------------------
            !calculate the flux of conservative variables related to g+ and f0
            !--------------------------------------------------
            !Maxwellian distribution H0 and B0
            call discrete_maxwell(H0,B0,prim)

            !Shakhov part H+ and B+
            call shakhov_part(H0,B0,qf,prim,H_plus,B_plus)

            !macro flux related to g+ and f0
            face%flux(1) = face%flux(1)+Mt(1)*sum(weight*uspace*H_plus)+Mt(4)*sum(weight*uspace*h)-Mt(5)*sum(weight*uspace**2*sh)
            face%flux(2) = face%flux(2)+Mt(1)*sum(weight*uspace**2*H_plus)+Mt(4)*sum(weight*uspace**2*h)-Mt(5)*sum(weight*uspace**3*sh)
            face%flux(3) = face%flux(3)+&
                           Mt(1)*0.5*(sum(weight*uspace*uspace**2*H_plus)+sum(weight*uspace*B_plus))+&
                           Mt(4)*0.5*(sum(weight*uspace*uspace**2*h)+sum(weight*uspace*b))-&
                           Mt(5)*0.5*(sum(weight*uspace**2*uspace**2*sh)+sum(weight*uspace**2*sb))

            !--------------------------------------------------
            !calculate flux of distribution function
            !--------------------------------------------------
            face%flux_h = Mt(1)*uspace*(H0+H_plus)+&
                          Mt(2)*uspace**2*(aL(1)*H0+aL(2)*uspace*H0+0.5*aL(3)*(uspace**2*H0+B0))*delta+&
                          Mt(2)*uspace**2*(aR(1)*H0+aR(2)*uspace*H0+0.5*aR(3)*(uspace**2*H0+B0))*(1-delta)+&
                          Mt(3)*uspace*(aT(1)*H0+aT(2)*uspace*H0+0.5*aT(3)*(uspace**2*H0+B0))+&
                          Mt(4)*uspace*h-Mt(5)*uspace**2*sh

            face%flux_b = Mt(1)*uspace*(B0+B_plus)+&
                          Mt(2)*uspace**2*(aL(1)*B0+aL(2)*uspace*B0+0.5*aL(3)*(uspace**2*B0+Mxi(2)*H0))*delta+&
                          Mt(2)*uspace**2*(aR(1)*B0+aR(2)*uspace*B0+0.5*aR(3)*(uspace**2*B0+Mxi(2)*H0))*(1-delta)+&
                          Mt(3)*uspace*(aT(1)*B0+aT(2)*uspace*B0+0.5*aT(3)*(uspace**2*B0+Mxi(2)*H0))+&
                          Mt(4)*uspace*b-Mt(5)*uspace**2*sb
        end subroutine calc_flux
        
        !--------------------------------------------------
        !>calculate micro slope of Maxwellian distribution
        !>@param[in] prim        :primary variables
        !>@param[in] sw          :slope of W
        !>@return    micro_slope :slope of Maxwellian distribution
        !--------------------------------------------------
        function micro_slope(prim,sw)
            real(kind=RKD),intent(in) :: prim(3),sw(3)
            real(kind=RKD) :: micro_slope(3)

            micro_slope(3) = 4.0*prim(3)**2/(ck+1)/prim(1)*(2.0*sw(3)-2.0*prim(2)*sw(2)+sw(1)*(prim(2)**2-0.5*(ck+1)/prim(3)))

            micro_slope(2) = 2.0*prim(3)/prim(1)*(sw(2)-prim(2)*sw(1))-prim(2)*micro_slope(3)
            micro_slope(1) = sw(1)/prim(1)-prim(2)*micro_slope(2)-0.5*(prim(2)**2+0.5*(ck+1)/prim(3))*micro_slope(3)
        end function micro_slope

        !--------------------------------------------------
        !>calculate moments of velocity
        !>@param[in] prim :primary variables
        !>@param[out] Mu        :<u^n>
        !>@param[out] Mxi       :<\xi^l>
        !>@param[out] Mu_L,Mu_R :<u^n>_{>0},<u^n>_{<0}
        !--------------------------------------------------
        subroutine calc_moment_u(prim,Mu,Mxi,Mu_L,Mu_R)
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD),intent(out) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM)
            real(kind=RKD),intent(out) :: Mxi(0:2)
            integer :: i

            !moments of normal velocity
            Mu_L(0) = 0.5*erfc(-sqrt(prim(3))*prim(2))
            Mu_L(1) = prim(2)*Mu_L(0)+0.5*exp(-prim(3)*prim(2)**2)/sqrt(PI*prim(3))
            Mu_R(0) = 0.5*erfc(sqrt(prim(3))*prim(2))
            Mu_R(1) = prim(2)*Mu_R(0)-0.5*exp(-prim(3)*prim(2)**2)/sqrt(PI*prim(3))

            do i=2,MNUM
                Mu_L(i) = prim(2)*Mu_L(i-1)+0.5*(i-1)*Mu_L(i-2)/prim(3)
                Mu_R(i) = prim(2)*Mu_R(i-1)+0.5*(i-1)*Mu_R(i-2)/prim(3)
            end do

            Mu = Mu_L+Mu_R

            !moments of \xi
            Mxi(0) = 1.0 !<\xi^0>
            Mxi(1) = 0.5*ck/prim(3) !<\xi^2>
            Mxi(2) = (ck**2+2.0*ck)/(4.0*prim(3)**2) !<\xi^4>
        end subroutine calc_moment_u

        !--------------------------------------------------
        !>calculate <u^\alpha*\xi^\delta*\psi>
        !>@param[in] Mu        :<u^\alpha>
        !>@param[in] Mxi       :<\xi^l>
        !>@param[in] alpha     :exponential index of u
        !>@param[in] delta     :exponential index of \xi
        !>@return    moment_uv :moment of <u^\alpha*\xi^\delta*\psi>
        !--------------------------------------------------
        function moment_uv(Mu,Mxi,alpha,delta)
            real(kind=RKD),intent(in) :: Mu(0:MNUM),Mxi(0:2)
            integer,intent(in) :: alpha,delta
            real(kind=RKD) :: moment_uv(3)

            moment_uv(1) = Mu(alpha)*Mxi(delta/2)
            moment_uv(2) = Mu(alpha+1)*Mxi(delta/2)
            moment_uv(3) = 0.5*(Mu(alpha+2)*Mxi(delta/2)+Mu(alpha)*Mxi((delta+2)/2))
        end function moment_uv

        !--------------------------------------------------
        !>calculate <a*u^\alpha*\psi>
        !>@param[in] a         :micro slope of Maxwellian
        !>@param[in] Mu        :<u^\alpha>
        !>@param[in] Mxi       :<\xi^l>
        !>@param[in] alpha     :exponential index of u
        !>@return    moment_au :moment of <a*u^\alpha*\psi>
        !--------------------------------------------------
        function moment_au(a,Mu,Mxi,alpha)
            real(kind=RKD),intent(in) :: a(3)
            real(kind=RKD),intent(in) :: Mu(0:MNUM),Mxi(0:2)
            integer,intent(in) :: alpha
            real(kind=RKD) :: moment_au(3)

            moment_au = a(1)*moment_uv(Mu,Mxi,alpha+0,0)+&
                        a(2)*moment_uv(Mu,Mxi,alpha+1,0)+&
                        0.5*a(3)*moment_uv(Mu,Mxi,alpha+2,0)+&
                        0.5*a(3)*moment_uv(Mu,Mxi,alpha+0,2)
        end function moment_au
end module flux

!--------------------------------------------------
!>UGKS solver
!--------------------------------------------------
module solver
    use global_data
    use tools
    use flux
    implicit none
    contains
        !--------------------------------------------------
        !>calculate time step
        !--------------------------------------------------
        subroutine timestep()
            real(kind=RKD) :: tmax !max 1/dt allowed
            real(kind=RKD) :: sos !speed of sound
            real(kind=RKD) :: prim(3) !primary variables
            integer :: i

            !set initial value
            tmax = 0.0

            !$omp parallel 
            !$omp do private(i,sos,prim) reduction(max:tmax)
            do i=ixmin,ixmax
                !convert conservative variables to primary variables
                prim = get_primary(ctr(i)%w)

                !sound speed
                sos = get_sos(prim)

                !maximum velocity
                prim(2) = max(umax,abs(prim(2)))+sos

                !maximum 1/dt allowed
                tmax = max(tmax,prim(2)/ctr(i)%length)
            end do 
            !$omp end do
            !$omp end parallel

            !time step
            dt = cfl/tmax
        end subroutine timestep

        !--------------------------------------------------
        !>calculate the slope of distribution function
        !--------------------------------------------------
        subroutine interpolation()
            integer :: i

            call interp_boundary(ctr(ixmin),ctr(ixmin),ctr(ixmin+1))
            call interp_boundary(ctr(ixmax),ctr(ixmax-1),ctr(ixmax))

            !$omp parallel
            !$omp do
            do i=ixmin+1,ixmax-1
                call interp_inner(ctr(i-1),ctr(i),ctr(i+1))
            end do
            !$omp end do nowait
            !$omp end parallel
        end subroutine interpolation

        !--------------------------------------------------
        !>calculate the flux across the interfaces
        !--------------------------------------------------
        subroutine evolution()
            integer :: i

            !$omp parallel
            !$omp do
            do i=ixmin,ixmax+1 !with ghost cell
                call calc_flux(ctr(i-1),vface(i),ctr(i))
            end do
            !$omp end do nowait
            !$omp end parallel 
        end subroutine evolution

        !--------------------------------------------------
        !>update cell averaged values
        !--------------------------------------------------
        subroutine update()
            real(kind=RKD),allocatable,dimension(:) :: H_old,B_old !equilibrium distribution at t=t^n
            real(kind=RKD),allocatable,dimension(:) :: H,B !equilibrium distribution at t=t^{n+1}
            real(kind=RKD),allocatable,dimension(:) :: H_plus,B_plus !Shakhov part
            real(kind=RKD) :: w_old(3) !conservative variables at t^n
            real(kind=RKD) :: prim_old(3),prim(3) !primary variables at t^n and t^{n+1}
            real(kind=RKD) :: tau_old,tau !collision time and t^n and t^{n+1}
            real(kind=RKD) :: qf
            integer :: i

            !allocate arrays
            allocate(H_old(unum))
            allocate(B_old(unum))
            allocate(H(unum))
            allocate(B(unum))
            allocate(H_plus(unum))
            allocate(B_plus(unum))

            !set initial value
            do i=ixmin,ixmax
                !--------------------------------------------------
                !store W^n and calculate H^n,B^n,\tau^n
                !--------------------------------------------------
                w_old = ctr(i)%w !store W^n
                
                prim_old = get_primary(w_old) !convert to primary variables
                call discrete_maxwell(H_old,B_old,prim_old) !calculate Maxwellian
                tau_old = get_tau(prim_old) !calculate collision time \tau^n

                !--------------------------------------------------
                !update W^{n+1} and calculate H^{n+1},B^{n+1},\tau^{n+1}
                !--------------------------------------------------
                ctr(i)%w = ctr(i)%w+(vface(i)%flux-vface(i+1)%flux)/ctr(i)%length !update W^{n+1}

                prim = get_primary(ctr(i)%w)
                call discrete_maxwell(H,B,prim)
                tau = get_tau(prim)

                !--------------------------------------------------
                !Shakhov part
                !--------------------------------------------------
                !heat flux at t=t^n
                qf = get_heat_flux(ctr(i)%h,ctr(i)%b,prim_old) 

                !h^+ = H+H^+ at t=t^n
                call shakhov_part(H_old,B_old,qf,prim_old,H_plus,B_plus) !H^+ and B^+
                H_old = H_old+H_plus !h^+
                B_old = B_old+B_plus !b^+

                !h^+ = H+H^+ at t=t^{n+1}
                call shakhov_part(H,B,qf,prim,H_plus,B_plus)
                H = H+H_plus
                B = B+B_plus

                !--------------------------------------------------
                !update distribution function
                !--------------------------------------------------
                ctr(i)%h = (ctr(i)%h+(vface(i)%flux_h-vface(i+1)%flux_h)/ctr(i)%length+&
                                  0.5*dt*(H/tau+(H_old-ctr(i)%h)/tau_old))/(1.0+0.5*dt/tau)
                ctr(i)%b = (ctr(i)%b+(vface(i)%flux_b-vface(i+1)%flux_b)/ctr(i)%length+&
                                  0.5*dt*(B/tau+(B_old-ctr(i)%b)/tau_old))/(1.0+0.5*dt/tau)
            end do
        end subroutine update

        !--------------------------------------------------
        !>one-sided interpolation of the boundary cell
        !>@param[inout] cell_N :the target boundary cell
        !>@param[inout] cell_L :the left cell
        !>@param[inout] cell_R :the right cell
        !--------------------------------------------------
        subroutine interp_boundary(cell_N,cell_L,cell_R)
            type(cell_center),intent(inout) :: cell_N
            type(cell_center),intent(inout) :: cell_L,cell_R

            cell_N%sh = (cell_R%h-cell_L%h)/(0.5*cell_R%length+0.5*cell_L%length)
            cell_N%sb = (cell_R%b-cell_L%b)/(0.5*cell_R%length+0.5*cell_L%length)
        end subroutine interp_boundary

        !--------------------------------------------------
        !>interpolation of the inner cells
        !>@param[in]    cell_L :the left cell
        !>@param[inout] cell_N :the target cell
        !>@param[in]    cell_R :the right cell
        !--------------------------------------------------
        subroutine interp_inner(cell_L,cell_N,cell_R)
            type(cell_center),intent(in) :: cell_L,cell_R
            type(cell_center),intent(inout) :: cell_N
            real(kind=RKD),allocatable,dimension(:) :: sL,sR

            !allocate array
            allocate(sL(unum))
            allocate(sR(unum))

            sL = (cell_N%h-cell_L%h)/(0.5*cell_N%length+0.5*cell_L%length)
            sR = (cell_R%h-cell_N%h)/(0.5*cell_R%length+0.5*cell_N%length)
            cell_N%sh = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

            sL = (cell_N%b-cell_L%b)/(0.5*cell_N%length+0.5*cell_L%length)
            sR = (cell_R%b-cell_N%b)/(0.5*cell_R%length+0.5*cell_N%length)
            cell_N%sb = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)
        end subroutine interp_inner
end module solver

!--------------------------------------------------
!>input and output
!--------------------------------------------------
module io
    use global_data
    use tools
    implicit none
    contains
        !--------------------------------------------------
        !>main initialization subroutine
        !--------------------------------------------------
        subroutine init()
            real(kind=RKD) :: init_gas(3) !initial condition
            real(kind=RKD) :: alpha_ref,omega_ref !molecule model coefficient in referece state
            real(kind=RKD) :: kn !Knudsen number in reference state
            real(kind=RKD) :: Ma !Mach number in front of shock
            real(kind=RKD) :: xlength !length of computational domain
            real(kind=RKD) :: umin !smallest and largest discrete velocity (for newton-cotes)
            real(kind=RKD) :: xscale !cell size/mean free path

            !control
            cfl = 0.95 !CFL number
            max_time = 250 !output time interval
            method_output = NORMALIZE !normalize the variables by (V-V1)/(V2-V1), where V1,V2 are upstream and downstream values

            !gas
            ck = 2 !internal degree of freedom
            gam = get_gamma(ck) !ratio of specific heat
            pr = 2.0/3.0 !Prandtl number
            omega = 0.72 !temperature dependence index in VHS model
            kn = 1.0 !Knudsen number in reference state
            alpha_ref = 1.0 !coefficient in HS model
            omega_ref = 0.5 !coefficient in HS model
            mu_ref = get_mu(kn,alpha_ref,omega_ref) !reference viscosity coefficient

            !velocity space
            unum = 100
            umin = -15.0
            umax = +15.0

            !geometry
            xlength = 50.0
            xscale = 0.5

            !Mach number
            Ma = 8.0

            call init_geometry(xlength,xscale) !initialize the geometry
            call init_velocity_newton(unum,umin,umax) !initialize discrete velocity space
            call init_allocation() !allocate arrays
            call init_flow_field(Ma) !set the initial condition
        end subroutine init

        !--------------------------------------------------
        !>initialize the geometry
        !>@param[inout] xlength :domain length
        !>@param[in]    xscale  :cell size/mean free path
        !--------------------------------------------------
        subroutine init_geometry(xlength,xscale)
            real(kind=RKD),intent(inout) :: xlength
            real(kind=RKD),intent(in) :: xscale
            real(kind=RKD) :: xnum !number of cells
            real(kind=RKD) :: dx !cell length
            integer :: i

            !adjust values
            xnum = xlength/xscale
            xlength = xnum*xscale

            !cell index range
            ixmin = 1
            ixmax = xnum

            !allocation
            allocate(ctr(ixmin-1:ixmax+1)) !cell center (with ghost cell)
            allocate(vface(ixmin:ixmax+1)) !vertical and horizontal cell interface

            !cell length and area
            dx = xlength/(ixmax-ixmin+1)

            !cell center (with ghost cell)
            forall(i=ixmin-1:ixmax+1) 
                ctr(i)%x = (i-0.5)*dx
                ctr(i)%length = dx
            end forall
        end subroutine init_geometry

        !--------------------------------------------------
        !>set discrete velocity space using Newton–Cotes formulas
        !>@param[inout] num_u :number of velocity points
        !>@param[in]    min_u :smallest discrete velocity
        !>@param[in]    max_u :largest discrete velocity
        !--------------------------------------------------
        subroutine init_velocity_newton(num_u,min_u,max_u)
            integer,intent(inout) :: num_u
            real(kind=RKD),intent(in) :: min_u,max_u
            real(kind=RKD) :: du !spacing in u velocity space
            integer :: i

            !modify unum if not appropriate
            unum = (num_u/4)*4+1
    
            !allocate array
            allocate(uspace(unum))
            allocate(weight(unum))

            !spacing in u and v velocity space
            du = (max_u-min_u)/(unum-1)

            !velocity space
            forall(i=1:unum)
                uspace(i) = min_u+(i-1)*du
                weight(i) = (newton_coeff(i,unum)*du)
            end forall

            !maximum micro velocity
            umax = max_u

            contains
                !--------------------------------------------------
                !>calculate the coefficient for newton-cotes formula
                !>@param[in] idx          :index in velocity space
                !>@param[in] num          :total number in velocity space
                !>@return    newton_coeff :coefficient for newton-cotes formula
                !--------------------------------------------------
                pure function newton_coeff(idx,num)
                    integer,intent(in) :: idx,num
                    real(kind=RKD) :: newton_coeff

                    if (idx==1 .or. idx==num) then 
                        newton_coeff = 14.0/45.0
                    else if (mod(idx-5,4)==0) then
                        newton_coeff = 28.0/45.0
                    else if (mod(idx-3,4)==0) then
                        newton_coeff = 24.0/45.0
                    else
                        newton_coeff = 64.0/45.0
                    end if
                end function newton_coeff
        end subroutine init_velocity_newton

        !--------------------------------------------------
        !>allocate arrays
        !--------------------------------------------------
        subroutine init_allocation()
            integer :: i

            !cell center (with ghost cell)
            do i=ixmin-1,ixmax+1
                allocate(ctr(i)%h(unum))
                allocate(ctr(i)%b(unum))
                allocate(ctr(i)%sh(unum))
                allocate(ctr(i)%sb(unum))
            end do

            !cell interface
            do i=ixmin,ixmax+1
                allocate(vface(i)%flux_h(unum))
                allocate(vface(i)%flux_b(unum))
            end do
        end subroutine init_allocation

        !--------------------------------------------------
        !>set the initial condition
        !>@param[in] Ma_L :Mach number in front of shock
        !--------------------------------------------------
        subroutine init_flow_field(Ma_L)
            real(kind=RKD),intent(in) :: Ma_L
            real(kind=RKD) :: Ma_R !Mach number after shock
            real(kind=RKD) :: ratio_T !T2/T1
            real(kind=RKD) :: prim_L(3), prim_R(3) !primary variables before and after shock
            real(kind=RKD) :: w_L(3), w_R(3) !conservative variables before and after shock
            real(kind=RKD),allocatable,dimension(:) :: H_L,B_L !distribution function before shock
            real(kind=RKD),allocatable,dimension(:) :: H_R,B_R !distribution function before shock
            integer :: i

            !allocation
            allocate(H_L(unum))
            allocate(B_L(unum))
            allocate(H_R(unum))
            allocate(B_R(unum))

            !upstream condition (before shock)
            prim_L(1) = 1.0 !density
            prim_L(2) = Ma_L*sqrt(gam/2.0) !velocity
            prim_L(3) = 1.0 !lambda=1/temperature

            !downstream condition (after shock)
            Ma_R = sqrt((Ma_L**2*(gam-1)+2)/(2*gam*Ma_L**2-(gam-1)))
            ratio_T = (1+(gam-1)/2*Ma_L**2)*(2*gam/(gam-1)*Ma_L**2-1)/(Ma_L**2*(2*gam/(gam-1)+(gam-1)/2))

            prim_R(1) = prim_L(1)*(gam+1)*Ma_L**2/((gam-1)*Ma_L**2+2)
            prim_R(2) = Ma_R*sqrt(gam/2.0)*sqrt(ratio_T)
            prim_R(3) = prim_L(3)/ratio_T

            !conservative variables and distribution function
            w_L = get_conserved(prim_L)
            w_R = get_conserved(prim_R)
            call discrete_maxwell(H_L,B_L,prim_L)
            call discrete_maxwell(H_R,B_R,prim_R)

            !initialize field (with ghost cell)
            forall(i=ixmin-1:(ixmin+ixmax)/2)
                ctr(i)%w = w_L
                ctr(i)%h = H_L
                ctr(i)%b = B_L
            end forall
            
            forall(i=(ixmin+ixmax)/2+1:ixmax+1)
                ctr(i)%w = w_R
                ctr(i)%h = H_R
                ctr(i)%b = B_R
            end forall 

            !initialize slope of distribution function at ghost cell
            ctr(ixmin-1)%sh = 0.0
            ctr(ixmin-1)%sb = 0.0
            ctr(ixmax+1)%sh = 0.0
            ctr(ixmax+1)%sb = 0.0
        end subroutine init_flow_field

        !--------------------------------------------------
        !>write result
        !--------------------------------------------------
        subroutine output()
            real(kind=RKD) :: rho_avg !average density
            real(kind=RKD) :: xmid !location of average density
            integer :: i

            !--------------------------------------------------
            !preparation
            !--------------------------------------------------
            !solution (with ghost cell)
            do i=ixmin-1,ixmax+1
                ctr(i)%primary(1) = ctr(i)%w(1) !density
                ctr(i)%primary(2) = get_temperature(ctr(i)%h,ctr(i)%b,get_primary(ctr(i)%w)) !temperature
            end do

            !find middle location - the location of average density
            rho_avg = 0.5*(ctr(ixmin-1)%w(1)+ctr(ixmax+1)%w(1))

            do i=ixmin,ixmax
                if ((ctr(i)%w(1)-rho_avg)*(ctr(i+1)%w(1)-rho_avg)<=0) then
                    xmid = ctr(i)%x+(ctr(i+1)%x-ctr(i)%x)/(ctr(i+1)%w(1)-ctr(i)%w(1))*(rho_avg-ctr(i)%w(1))
                end if
            end do

            !normalization
            if (method_output==NORMALIZE) then
                ctr%primary(1) = (ctr%primary(1)-ctr(ixmin-1)%primary(1))/(ctr(ixmax+1)%primary(1)-ctr(ixmin-1)%primary(1))
                ctr%primary(2) = (ctr%primary(2)-ctr(ixmin-1)%primary(2))/(ctr(ixmax+1)%primary(2)-ctr(ixmin-1)%primary(2))
            end if

            !--------------------------------------------------
            !write to file
            !--------------------------------------------------
            !open result file and write header
            open(unit=RSTFILE,file=RSTFILENAME,status="replace",action="write")
            write(RSTFILE,*) "VARIABLES = X, RHO, T"
            write(RSTFILE,*) 'ZONE  T="Time: ',sim_time,'", I = ',ixmax-ixmin+1,', DATAPACKING=BLOCK'

            !write geometry (cell-centered)
            write(RSTFILE,"(6(ES23.16,2X))") ctr(ixmin:ixmax)%x-xmid

            !write solution (cell-centered)
            do i=1,2
                write(RSTFILE,"(6(ES23.16,2X))") ctr(ixmin:ixmax)%primary(i)
            end do
    
            !close file
            close(RSTFILE)
        end subroutine output
end module io

!--------------------------------------------------
!>main program
!--------------------------------------------------
program main
    use global_data
    use tools
    use solver
    use flux
    use io
    implicit none

    !initialization
    call init() 

    !set initial value
    iter = 1 !number of iteration
    sim_time = 0.0 !simulation time

    !open file and write header
    open(unit=HSTFILE,file=HSTFILENAME,status="replace",action="write") !open history file
    write(HSTFILE,*) "VARIABLES = iter, sim_time, dt" !write header

    !iteration
    do while(.true.)
        call timestep() !calculate time step
        call interpolation() !calculate the slope of distribution function
        call evolution() !calculate flux across the interfaces
        call update() !update cell averaged value

        !check if output
        if (sim_time>=max_time) exit

        !write iteration situation every 10 iterations
        if (mod(iter,10)==0) then
            write(*,"(A18,I15,2E15.7)") "iter,sim_time,dt:",iter,sim_time,dt
            write(HSTFILE,"(I15,2E15.7)") iter,sim_time,dt
        end if

        iter = iter+1
        sim_time = sim_time+dt
    end do

    !close history file
    close(HSTFILE)

    !output solution
    call output()
end program main

! vim: set ft=fortran tw=0: 
