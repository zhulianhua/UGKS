!-----------------------------------------------------------------------------
!UGKS2D - lid-driven cavity calculation using Unified Gas-Kinetic Scheme
!-----------------------------------------------------------------------------
!Copyright (C) 2012 Wang Ruijie <lainme993@gmail.com>
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
    real(kind=RKD) :: res(4) !residual
    real(kind=RKD) :: eps !convergence criteria
    real(kind=RKD) :: sim_time !current simulation time
    character(len=255),parameter :: HSTFILENAME = "cavity.hst" !history file name
    character(len=255),parameter :: RSTFILENAME = "cavity.rst" !result file name
    integer :: iter !iteration

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
    !direction
    integer,parameter :: IDIRC = 1 !i direction
    integer,parameter :: JDIRC = 2 !j direction
    !rotation
    integer,parameter :: RN = 1 !no frame rotation
    integer,parameter :: RY = -1 !with frame rotation
    !I/O
    integer,parameter :: HSTFILE = 20 !history file ID
    integer,parameter :: RSTFILE = 21 !result file ID

    !--------------------------------------------------
    !basic derived type
    !--------------------------------------------------
    !cell center
    type :: cell_center
        !geometry
        real(kind=RKD) :: x,y !x and y coordinates of cell center
        real(kind=RKD) :: area !cell area
        real(kind=RKD) :: length(2) !length in i and j direction
        !flow field
        real(kind=RKD) :: w(4) !density, x-momentum, y-momentum, total energy
        real(kind=RKD),allocatable,dimension(:,:) :: h,b !distribution function
        real(kind=RKD),allocatable,dimension(:,:,:) :: sh,sb !slope of distribution function in i and j direction
    end type cell_center

    !cell interface
    type :: cell_interface
        !geometry
        real(kind=RKD) :: length !length of cell interface
        real(kind=RKD) :: cosx,cosy !directional cosine
        !flow flux
        real(kind=RKD) :: flux(4) !mass flux, x and y momentum flux, energy flux
        real(kind=RKD),allocatable,dimension(:,:) :: flux_h,flux_b !flux of distribution function
    end type cell_interface

    !--------------------------------------------------
    !flow field
    !--------------------------------------------------
    !index method
    !          (i,j+1)
    !     ----------------
    !     |              |
    !     |              |
    !     |              |
    !(i,j)|     (i,j)    |(i+1,j)
    !     |              |
    !     |              |
    !     |              |
    !     ----------------
    !           (i,j)
    integer :: ixmin,ixmax,iymin,iymax !index range in i and j direction
    integer :: ngrid !number of total grids
    type(cell_center),allocatable,dimension(:,:) :: ctr !cell centers
    type(cell_interface),allocatable,dimension(:,:) :: vface,hface !vertical and horizontal interfaces
    real(kind=RKD) :: bc_W(4),bc_E(4),bc_S(4),bc_N(4) !boundary conditions at west,east,south and north boundary

    !--------------------------------------------------
    !discrete velocity space. Inefficient implementation but clear code
    !--------------------------------------------------
    integer :: unum,vnum !number of velocity points for u and v
    real(kind=RKD) :: umax,vmax !maximum micro velocity
    real(kind=RKD),allocatable,dimension(:,:) :: uspace,vspace !u and v discrete velocity space
    real(kind=RKD),allocatable,dimension(:,:) :: weight !weight at velocity u_k and v_l
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
        !>@param[out] h,b   :distribution function
        !>@param[in]  vn,vt :normal and tangential velocity
        !>@param[in]  prim  :primary variables
        !--------------------------------------------------
        subroutine discrete_maxwell(h,b,vn,vt,prim)
            real(kind=RKD),dimension(:,:),intent(out) :: h,b
            real(kind=RKD),dimension(:,:),intent(in) :: vn,vt
            real(kind=RKD),intent(in) :: prim(4)

            h = prim(1)*(prim(4)/PI)*exp(-prim(4)*((vn-prim(2))**2+(vt-prim(3))**2))
            b = h*ck/(2.0*prim(4))
        end subroutine discrete_maxwell

        !--------------------------------------------------
        !>calculate the Shakhov part H^+, B^+
        !>@param[in]  H,B           :Maxwellian distribution function
        !>@param[in]  vn,vt         :normal and tangential velocity
        !>@param[in]  qf            :heat flux
        !>@param[in]  prim          :primary variables
        !>@param[out] H_plus,B_plus :Shakhov part
        !--------------------------------------------------
        subroutine shakhov_part(H,B,vn,vt,qf,prim,H_plus,B_plus)
            real(kind=RKD),dimension(:,:),intent(in) :: H,B
            real(kind=RKD),dimension(:,:),intent(in) :: vn,vt
            real(kind=RKD) :: qf(2)
            real(kind=RKD) :: prim(4)
            real(kind=RKD),dimension(:,:),intent(out) :: H_plus,B_plus

            H_plus = 0.8*(1-pr)*prim(4)**2/prim(1)*&
                     ((vn-prim(2))*qf(1)+(vt-prim(3))*qf(2))*(2*prim(4)*((vn-prim(2))**2+(vt-prim(3))**2)+ck-5)*H
            B_plus = 0.8*(1-pr)*prim(4)**2/prim(1)*&
                     ((vn-prim(2))*qf(1)+(vt-prim(3))*qf(2))*(2*prim(4)*((vn-prim(2))**2+(vt-prim(3))**2)+ck-3)*B
        end subroutine shakhov_part

        !--------------------------------------------------
        !>convert primary variables to conservative variables
        !>@param[in] prim          :primary variables
        !>@return    get_conserved :conservative variables
        !--------------------------------------------------
        function get_conserved(prim)
            real(kind=RKD),intent(in) :: prim(4)
            real(kind=RKD) :: get_conserved(4)

            get_conserved(1) = prim(1)
            get_conserved(2) = prim(1)*prim(2)
            get_conserved(3) = prim(1)*prim(3)
            get_conserved(4) = 0.5*prim(1)/prim(4)/(gam-1.0)+0.5*prim(1)*(prim(2)**2+prim(3)**2)
        end function get_conserved

        !--------------------------------------------------
        !>convert conservative variables to primary variables
        !>@param[in] w           :conservative variables
        !>@return    get_primary :conservative variables
        !--------------------------------------------------
        function get_primary(w)
            real(kind=RKD),intent(in) :: w(4)
            real(kind=RKD) :: get_primary(4) !primary variables

            get_primary(1) = w(1)
            get_primary(2) = w(2)/w(1)
            get_primary(3) = w(3)/w(1)
            get_primary(4) = 0.5*w(1)/(gam-1.0)/(w(4)-0.5*(w(2)**2+w(3)**2)/w(1))
        end function get_primary

        !--------------------------------------------------
        !>convert macro variables from local frame to global
        !>@param[in] w            :macro variables in local frame
        !>@param[in] cosx,cosy    :directional cosine
        !>@return    global_frame :macro variables in global frame
        !--------------------------------------------------
        function global_frame(w,cosx,cosy)
            real(kind=RKD),intent(in) :: w(4)
            real(kind=RKD),intent(in) :: cosx,cosy
            real(kind=RKD) :: global_frame(4)

            global_frame(1) = w(1)
            global_frame(2) = w(2)*cosx-w(3)*cosy
            global_frame(3) = w(2)*cosy+w(3)*cosx
            global_frame(4) = w(4)
        end function global_frame

        !--------------------------------------------------
        !>convert macro variables from global frame to local
        !>@param[in] w            :macro variables in global frame
        !>@param[in] cosx,cosy    :directional cosine
        !>@return    local_frame  :macro variables in local frame
        !--------------------------------------------------
        function local_frame(w,cosx,cosy)
            real(kind=RKD),intent(in) :: w(4)
            real(kind=RKD),intent(in) :: cosx,cosy
            real(kind=RKD) :: local_frame(4)

            local_frame(1) = w(1)
            local_frame(2) = w(2)*cosx+w(3)*cosy
            local_frame(3) = w(3)*cosx-w(2)*cosy
            local_frame(4) = w(4)
        end function local_frame

        !--------------------------------------------------
        !>obtain ratio of specific heat
        !>@param[in] ck        :internal degree of freedom
        !>@return    get_gamma :ratio of specific heat
        !--------------------------------------------------
        function get_gamma(ck)
            integer,intent(in) :: ck
            real(kind=RKD) :: get_gamma

            get_gamma = float(ck+4)/float(ck+2)
        end function get_gamma

        !--------------------------------------------------
        !>obtain speed of sound
        !>@param[in] prim    :primary variables
        !>@return    get_sos :speed of sound
        !--------------------------------------------------
        function get_sos(prim)
            real(kind=RKD),intent(in) :: prim(4)
            real(kind=RKD) :: get_sos !speed of sound

            get_sos = sqrt(0.5*gam/prim(4))
        end function get_sos
        
        !--------------------------------------------------
        !>calculate collision time
        !>@param[in] prim    :primary variables
        !>return     get_tau :collision time
        !--------------------------------------------------
        function get_tau(prim)
            real(kind=RKD),intent(in) :: prim(4)
            real(kind=RKD) :: get_tau

            get_tau = mu_ref*2*prim(4)**(1-omega)/prim(1)
        end function get_tau

        !--------------------------------------------------
        !>get heat flux
        !>@param[in] h,b   :distribution function
        !>@param[in] vn,vt :normal and tangential velocity
        !>@param[in] prim  :primary variables
        !--------------------------------------------------
        function get_heat_flux(h,b,vn,vt,prim)
            real(kind=RKD),dimension(:,:),intent(in) :: h,b
            real(kind=RKD),dimension(:,:),intent(in) :: vn,vt
            real(kind=RKD),intent(in) :: prim(4)
            real(kind=RKD) :: get_heat_flux(2) !heat flux in normal and tangential direction

            get_heat_flux(1) = 0.5*(sum(weight*(vn-prim(2))*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*(vn-prim(2))*b)) 
            get_heat_flux(2) = 0.5*(sum(weight*(vt-prim(3))*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*(vt-prim(3))*b)) 
        end function get_heat_flux

        !--------------------------------------------------
        !>get pressure
        !>@param[in] h,b   :distribution function
        !>@param[in] vn,vt :normal and tangential velocity
        !>@param[in] prim  :primary variables
        !--------------------------------------------------
        function get_pressure(h,b,vn,vt,prim)
            real(kind=RKD),dimension(:,:),intent(in) :: h,b
            real(kind=RKD),dimension(:,:),intent(in) :: vn,vt
            real(kind=RKD),intent(in) :: prim(4)
            real(kind=RKD) :: get_pressure !pressure

            get_pressure = (sum(weight*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*b))/(ck+2)
        end function get_pressure

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
        !>@param[in]    idx    :index indicating i or j direction
        !--------------------------------------------------
        subroutine calc_flux(cell_L,face,cell_R,idx)
            type(cell_center),intent(in) :: cell_L,cell_R
            type(cell_interface),intent(inout) :: face
            integer,intent(in) :: idx
            real(kind=RKD),allocatable,dimension(:,:) :: vn,vt !normal and tangential micro velocity
            real(kind=RKD),allocatable,dimension(:,:) :: h,b !distribution function at the interface
            real(kind=RKD),allocatable,dimension(:,:) :: H0,B0 !Maxwellian distribution function
            real(kind=RKD),allocatable,dimension(:,:) :: H_plus,B_plus !Shakhov part of the equilibrium distribution
            real(kind=RKD),allocatable,dimension(:,:) :: sh,sb !slope of distribution function at the interface
            integer,allocatable,dimension(:,:) :: delta !Heaviside step function
            real(kind=RKD) :: w(4),prim(4) !conservative and primary variables at the interface
            real(kind=RKD) :: qf(2) !heat flux in normal and tangential direction
            real(kind=RKD) :: sw(4) !slope of W
            real(kind=RKD) :: aL(4),aR(4),aT(4) !micro slope of Maxwellian distribution, left,right and time.
            real(kind=RKD) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM),Mv(0:MTUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<v^m>,<\xi^l>
            real(kind=RKD) :: Mau_0(4),Mau_L(4),Mau_R(4),Mau_T(4) !<u\psi>,<aL*u^n*\psi>,<aR*u^n*\psi>,<A*u*\psi>
            real(kind=RKD) :: tau !collision time
            real(kind=RKD) :: Mt(5) !some time integration terms
            integer :: i,j

            !--------------------------------------------------
            !prepare
            !--------------------------------------------------
            !allocate array
            allocate(vn(unum,vnum))
            allocate(vt(unum,vnum))
            allocate(delta(unum,vnum))
            allocate(h(unum,vnum))
            allocate(b(unum,vnum))
            allocate(sh(unum,vnum))
            allocate(sb(unum,vnum))
            allocate(H0(unum,vnum))
            allocate(B0(unum,vnum))
            allocate(H_plus(unum,vnum))
            allocate(B_plus(unum,vnum))

            !convert the micro velocity to local frame
            vn = uspace*face%cosx+vspace*face%cosy
            vt = vspace*face%cosx-uspace*face%cosy

            !Heaviside step function
            delta = (vn+abs(vn)+1)/(2*abs(vn)+1)

            !--------------------------------------------------
            !reconstruct initial distribution
            !--------------------------------------------------
            h = (cell_L%h+0.5*cell_L%length(idx)*cell_L%sh(:,:,idx))*delta+&
                (cell_R%h-0.5*cell_R%length(idx)*cell_R%sh(:,:,idx))*(1-delta)
            b = (cell_L%b+0.5*cell_L%length(idx)*cell_L%sb(:,:,idx))*delta+&
                (cell_R%b-0.5*cell_R%length(idx)*cell_R%sb(:,:,idx))*(1-delta)
            sh = cell_L%sh(:,:,idx)*delta+cell_R%sh(:,:,idx)*(1-delta)
            sb = cell_L%sb(:,:,idx)*delta+cell_R%sb(:,:,idx)*(1-delta)

            !--------------------------------------------------
            !obtain macroscopic variables (local frame)
            !--------------------------------------------------
            !conservative variables w_0
            w(1) = sum(weight*h)
            w(2) = sum(weight*vn*h)
            w(3) = sum(weight*vt*h)
            w(4) = 0.5*(sum(weight*(vn**2+vt**2)*h)+sum(weight*b))

            !convert to primary variables
            prim = get_primary(w)

            !heat flux
            qf = get_heat_flux(h,b,vn,vt,prim) 

            !--------------------------------------------------
            !calculate a^L,a^R
            !--------------------------------------------------
            sw = (w-local_frame(cell_L%w,face%cosx,face%cosy))/(0.5*cell_L%length(idx)) !left slope of W
            aL = micro_slope(prim,sw) !calculate a^L

            sw = (local_frame(cell_R%w,face%cosx,face%cosy)-w)/(0.5*cell_R%length(idx)) !right slope of W
            aR = micro_slope(prim,sw) !calculate a^R

            !--------------------------------------------------
            !calculate time slope of W and A
            !--------------------------------------------------
            !<u^n>,<v^m>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
            call calc_moment_u(prim,Mu,Mv,Mxi,Mu_L,Mu_R) 

            Mau_L = moment_au(aL,Mu_L,Mv,Mxi,1,0) !<aL*u*\psi>_{>0}
            Mau_R = moment_au(aR,Mu_R,Mv,Mxi,1,0) !<aR*u*\psi>_{<0}

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
            Mau_0 = moment_uv(Mu,Mv,Mxi,1,0,0) !<u*\psi>
            Mau_L = moment_au(aL,Mu_L,Mv,Mxi,2,0) !<aL*u^2*\psi>_{>0}
            Mau_R = moment_au(aR,Mu_R,Mv,Mxi,2,0) !<aR*u^2*\psi>_{<0}
            Mau_T = moment_au(aT,Mu,Mv,Mxi,1,0) !<A*u*\psi>

            face%flux = Mt(1)*prim(1)*Mau_0+Mt(2)*prim(1)*(Mau_L+Mau_R)+Mt(3)*prim(1)*Mau_T

            !--------------------------------------------------
            !calculate the flux of conservative variables related to g+ and f0
            !--------------------------------------------------
            !Maxwellian distribution H0 and B0
            call discrete_maxwell(H0,B0,vn,vt,prim)

            !Shakhov part H+ and B+
            call shakhov_part(H0,B0,vn,vt,qf,prim,H_plus,B_plus)

            !macro flux related to g+ and f0
            face%flux(1) = face%flux(1)+Mt(1)*sum(weight*vn*H_plus)+Mt(4)*sum(weight*vn*h)-Mt(5)*sum(weight*vn**2*sh)
            face%flux(2) = face%flux(2)+Mt(1)*sum(weight*vn*vn*H_plus)+Mt(4)*sum(weight*vn*vn*h)-Mt(5)*sum(weight*vn*vn**2*sh)
            face%flux(3) = face%flux(3)+Mt(1)*sum(weight*vt*vn*H_plus)+Mt(4)*sum(weight*vt*vn*h)-Mt(5)*sum(weight*vt*vn**2*sh)
            face%flux(4) = face%flux(4)+&
                           Mt(1)*0.5*(sum(weight*vn*(vn**2+vt**2)*H_plus)+sum(weight*vn*B_plus))+&
                           Mt(4)*0.5*(sum(weight*vn*(vn**2+vt**2)*h)+sum(weight*vn*b))-&
                           Mt(5)*0.5*(sum(weight*vn**2*(vn**2+vt**2)*sh)+sum(weight*vn**2*sb))

            !--------------------------------------------------
            !calculate flux of distribution function
            !--------------------------------------------------
            face%flux_h = Mt(1)*vn*(H0+H_plus)+&
                          Mt(2)*vn**2*(aL(1)*H0+aL(2)*vn*H0+aL(3)*vt*H0+0.5*aL(4)*((vn**2+vt**2)*H0+B0))*delta+&
                          Mt(2)*vn**2*(aR(1)*H0+aR(2)*vn*H0+aR(3)*vt*H0+0.5*aR(4)*((vn**2+vt**2)*H0+B0))*(1-delta)+&
                          Mt(3)*vn*(aT(1)*H0+aT(2)*vn*H0+aT(3)*vt*H0+0.5*aT(4)*((vn**2+vt**2)*H0+B0))+&
                          Mt(4)*vn*h-Mt(5)*vn**2*sh

            face%flux_b = Mt(1)*vn*(B0+B_plus)+&
                          Mt(2)*vn**2*(aL(1)*B0+aL(2)*vn*B0+aL(3)*vt*B0+0.5*aL(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*delta+&
                          Mt(2)*vn**2*(aR(1)*B0+aR(2)*vn*B0+aR(3)*vt*B0+0.5*aR(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*(1-delta)+&
                          Mt(3)*vn*(aT(1)*B0+aT(2)*vn*B0+aT(3)*vt*B0+0.5*aT(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
                          Mt(4)*vn*b-Mt(5)*vn**2*sb

            !--------------------------------------------------
            !final flux
            !--------------------------------------------------
            !convert to global frame
            face%flux = global_frame(face%flux,face%cosx,face%cosy) 
            !total flux
            face%flux = face%length*face%flux
            face%flux_h = face%length*face%flux_h
            face%flux_b = face%length*face%flux_b
        end subroutine calc_flux
        
        !--------------------------------------------------
        !>calculate flux of boundary interface, assuming left wall
        !>@param[in]    bc   :boundary condition
        !>@param[inout] face :the boundary interface
        !>@param[in]    cell :cell next to the boundary interface
        !>@param[in]    idx  :index indicating i or j direction
        !>@param[in]    rot  :indicating rotation
        !--------------------------------------------------
        subroutine calc_flux_boundary(bc,face,cell,idx,rot) 
            real(kind=RKD),intent(in) :: bc(4)
            type(cell_interface),intent(inout) :: face
            type(cell_center),intent(in) :: cell
            integer,intent(in) :: idx,rot
            real(kind=RKD),allocatable,dimension(:,:) :: vn,vt !normal and tangential micro velocity
            real(kind=RKD),allocatable,dimension(:,:) :: h,b !distribution function
            real(kind=RKD),allocatable,dimension(:,:) :: H0,B0 !Maxwellian distribution function at the wall
            integer,allocatable,dimension(:,:) :: delta !Heaviside step function
            real(kind=RKD) :: prim(4) !boundary condition in local frame
            real(kind=RKD) :: SF,SG

            !--------------------------------------------------
            !prepare
            !--------------------------------------------------
            !allocate array
            allocate(vn(unum,vnum))
            allocate(vt(unum,vnum))
            allocate(delta(unum,vnum))
            allocate(h(unum,vnum))
            allocate(b(unum,vnum))
            allocate(H0(unum,vnum))
            allocate(B0(unum,vnum))

            !convert the micro velocity to local frame
            vn = uspace*face%cosx+vspace*face%cosy
            vt = vspace*face%cosx-uspace*face%cosy

            !Heaviside step function. The rotation accounts for the right wall
            delta = (vn*rot+abs(vn)+1)/(2*abs(vn)+1)

            !boundary condition in local frame
            prim = local_frame(bc,face%cosx,face%cosy)

            !--------------------------------------------------
            !obtain h^{in} and b^{in}, rotation accounts for the right wall
            !--------------------------------------------------
            h = cell%h-rot*0.5*cell%length(idx)*cell%sh(:,:,idx)
            b = cell%b-rot*0.5*cell%length(idx)*cell%sb(:,:,idx)

            !--------------------------------------------------
            !calculate wall density and Maxwellian distribution
            !--------------------------------------------------
            SF = sum(weight*vn*h*(1-delta))
            SG = (prim(4)/PI)*sum(weight*vn*exp(-prim(4)*((vn-prim(2))**2+(vt-prim(3))**2))*delta)

            prim(1) = -SF/SG

            call discrete_maxwell(H0,B0,vn,vt,prim)

            !--------------------------------------------------
            !distribution function at the boundary interface
            !--------------------------------------------------
            h = H0*delta+h*(1-delta)
            b = B0*delta+b*(1-delta)

            !--------------------------------------------------
            !calculate flux
            !--------------------------------------------------
            face%flux(1) = sum(weight*vn*h)
            face%flux(2) = sum(weight*vn*vn*h)
            face%flux(3) = sum(weight*vn*vt*h)
            face%flux(4) = 0.5*sum(weight*vn*((vn**2+vt**2)*h+b))

            face%flux_h = vn*h
            face%flux_b = vn*b

            !--------------------------------------------------
            !final flux
            !--------------------------------------------------
            !convert to global frame
            face%flux = global_frame(face%flux,face%cosx,face%cosy) 
            !total flux
            face%flux = dt*face%length*face%flux
            face%flux_h = dt*face%length*face%flux_h
            face%flux_b = dt*face%length*face%flux_b
        end subroutine calc_flux_boundary

        !--------------------------------------------------
        !>calculate micro slope of Maxwellian distribution
        !>@param[in] prim        :primary variables
        !>@param[in] sw          :slope of W
        !>@return    micro_slope :slope of Maxwellian distribution
        !--------------------------------------------------
        function micro_slope(prim,sw)
            real(kind=RKD),intent(in) :: prim(4),sw(4)
            real(kind=RKD) :: micro_slope(4)

            micro_slope(4) = 4.0*prim(4)**2/(ck+2)/prim(1)*(2.0*sw(4)-2.0*prim(2)*sw(2)-2.0*prim(3)*sw(3)+sw(1)*(prim(2)**2+prim(3)**2-0.5*(ck+2)/prim(4)))

            micro_slope(3) = 2.0*prim(4)/prim(1)*(sw(3)-prim(3)*sw(1))-prim(3)*micro_slope(4)
            micro_slope(2) = 2.0*prim(4)/prim(1)*(sw(2)-prim(2)*sw(1))-prim(2)*micro_slope(4)
            micro_slope(1) = sw(1)/prim(1)-prim(2)*micro_slope(2)-prim(3)*micro_slope(3)-0.5*(prim(2)**2+prim(3)**2+0.5*(ck+2)/prim(4))*micro_slope(4)
        end function micro_slope

        !--------------------------------------------------
        !>calculate moments of velocity
        !>@param[in] prim :primary variables
        !>@param[out] Mu,Mv     :<u^n>,<v^m>
        !>@param[out] Mxi       :<\xi^l>
        !>@param[out] Mu_L,Mu_R :<u^n>_{>0},<u^n>_{<0}
        !--------------------------------------------------
        subroutine calc_moment_u(prim,Mu,Mv,Mxi,Mu_L,Mu_R)
            real(kind=RKD),intent(in) :: prim(4)
            real(kind=RKD),intent(out) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM)
            real(kind=RKD),intent(out) :: Mv(0:MTUM)
            real(kind=RKD),intent(out) :: Mxi(0:2)
            integer :: i

            !moments of normal velocity
            Mu_L(0) = 0.5*erfc(-sqrt(prim(4))*prim(2))
            Mu_L(1) = prim(2)*Mu_L(0)+0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))
            Mu_R(0) = 0.5*erfc(sqrt(prim(4))*prim(2))
            Mu_R(1) = prim(2)*Mu_R(0)-0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))

            do i=2,MNUM
                Mu_L(i) = prim(2)*Mu_L(i-1)+0.5*(i-1)*Mu_L(i-2)/prim(4)
                Mu_R(i) = prim(2)*Mu_R(i-1)+0.5*(i-1)*Mu_R(i-2)/prim(4)
            end do

            Mu = Mu_L+Mu_R

            !moments of tangential velocity
            Mv(0) = 1.0
            Mv(1) = prim(3)

            do i=2,MTUM
                Mv(i) = prim(3)*Mv(i-1)+0.5*(i-1)*Mv(i-2)/prim(4)
            end do

            !moments of \xi
            Mxi(0) = 1.0 !<\xi^0>
            Mxi(1) = 0.5*ck/prim(4) !<\xi^2>
            Mxi(2) = (ck**2+2.0*ck)/(4.0*prim(4)**2) !<\xi^4>
        end subroutine calc_moment_u

        !--------------------------------------------------
        !>calculate <u^\alpha*v^\beta*\xi^\delta*\psi>
        !>@param[in] Mu,Mv      :<u^\alpha>,<v^\beta>
        !>@param[in] Mxi        :<\xi^l>
        !>@param[in] alpha,beta :exponential index of u and v
        !>@param[in] delta      :exponential index of \xi
        !--------------------------------------------------
        function moment_uv(Mu,Mv,Mxi,alpha,beta,delta)
            real(kind=RKD),intent(in) :: Mu(0:MNUM),Mv(0:MTUM),Mxi(0:2)
            integer,intent(in) :: alpha,beta 
            integer,intent(in) :: delta
            real(kind=RKD) :: moment_uv(4)

            moment_uv(1) = Mu(alpha)*Mv(beta)*Mxi(delta/2)
            moment_uv(2) = Mu(alpha+1)*Mv(beta)*Mxi(delta/2)
            moment_uv(3) = Mu(alpha)*Mv(beta+1)*Mxi(delta/2)
            moment_uv(4) = 0.5*(Mu(alpha+2)*Mv(beta)*Mxi(delta/2)+Mu(alpha)*Mv(beta+2)*Mxi(delta/2)+Mu(alpha)*Mv(beta)*Mxi((delta+2)/2))
        end function moment_uv

        !--------------------------------------------------
        !>calculate <a*u^\alpha*v^\beta*\psi>
        !>@param[in] a          :micro slope of Maxwellian
        !>@param[in] Mu,Mv      :<u^\alpha>,<v^\beta>
        !>@param[in] Mxi        :<\xi^l>
        !>@param[in] alpha,beta :exponential index of u and v
        !--------------------------------------------------
        function moment_au(a,Mu,Mv,Mxi,alpha,beta)
            real(kind=RKD),intent(in) :: a(4)
            real(kind=RKD),intent(in) :: Mu(0:MNUM),Mv(0:MTUM),Mxi(0:2)
            integer,intent(in) :: alpha,beta
            real(kind=RKD) :: moment_au(4)

            moment_au = a(1)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+0,0)+&
                        a(2)*moment_uv(Mu,Mv,Mxi,alpha+1,beta+0,0)+&
                        a(3)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+1,0)+&
                        0.5*a(4)*moment_uv(Mu,Mv,Mxi,alpha+2,beta+0,0)+&
                        0.5*a(4)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+2,0)+&
                        0.5*a(4)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+0,2)
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
            real(kind=RKD) :: prim(4) !primary variables
            integer :: i,j

            !set initial value
            tmax = 0.0

            !$omp parallel 
            !$omp do private(i,j,sos,prim) reduction(max:tmax)
            do j=iymin,iymax
                do i=ixmin,ixmax
                    !convert conservative variables to primary variables
                    prim = get_primary(ctr(i,j)%w)

                    !sound speed
                    sos = get_sos(prim)

                    !maximum velocity
                    prim(2) = max(umax,abs(prim(2)))+sos
                    prim(3) = max(vmax,abs(prim(3)))+sos

                    !maximum 1/dt allowed
                    tmax = max(tmax,(ctr(i,j)%length(2)*prim(2)+ctr(i,j)%length(1)*prim(3))/ctr(i,j)%area)
                end do
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
            integer :: i,j
                !$omp parallel
                !--------------------------------------------------
                !i direction
                !--------------------------------------------------
                !$omp do
                do j=iymin,iymax
                    call interp_boundary(ctr(ixmin,j),ctr(ixmin,j),ctr(ixmin+1,j),IDIRC) !the last argument indicating i direction
                    call interp_boundary(ctr(ixmax,j),ctr(ixmax-1,j),ctr(ixmax,j),IDIRC)
                end do
                !$omp end do nowait

                !$omp do
                do j=iymin,iymax
                    do i=ixmin+1,ixmax-1
                        call interp_inner(ctr(i-1,j),ctr(i,j),ctr(i+1,j),IDIRC)
                    end do
                end do
                !$omp end do nowait

                !--------------------------------------------------
                !j direction
                !--------------------------------------------------
                !$omp do
                do i=ixmin,ixmax
                    call interp_boundary(ctr(i,iymin),ctr(i,iymin),ctr(i,iymin+1),JDIRC)
                    call interp_boundary(ctr(i,iymax),ctr(i,iymax-1),ctr(i,iymax),JDIRC)
                end do
                !$omp end do nowait

                !$omp do
                do j=iymin+1,iymax-1
                    do i=ixmin,ixmax
                        call interp_inner(ctr(i,j-1),ctr(i,j),ctr(i,j+1),JDIRC)
                    end do
                end do
                !$omp end do nowait
                !$omp end parallel
        end subroutine interpolation

        !--------------------------------------------------
        !>calculate the flux across the interfaces
        !--------------------------------------------------
        subroutine evolution()
            integer :: i,j

            !$omp parallel
            !$omp do
            do j=iymin,iymax
                call calc_flux_boundary(bc_W,vface(ixmin,j),ctr(ixmin,j),IDIRC,RN) !RN means no frame rotation
                call calc_flux_boundary(bc_E,vface(ixmax+1,j),ctr(ixmax,j),IDIRC,RY) !RY means with frame rotation
            end do
            !$omp end do nowait

            !$omp do
            do j=iymin,iymax
                do i=ixmin+1,ixmax
                    call calc_flux(ctr(i-1,j),vface(i,j),ctr(i,j),IDIRC)
                end do
            end do
            !$omp end do nowait

            !$omp do
            do i=ixmin,ixmax
                call calc_flux_boundary(bc_S,hface(i,iymin),ctr(i,iymin),JDIRC,RN)
                call calc_flux_boundary(bc_N,hface(i,iymax+1),ctr(i,iymax),JDIRC,RY)
            end do
            !$omp end do nowait

            !$omp do
            do j=iymin+1,iymax
                do i=ixmin,ixmax
                    call calc_flux(ctr(i,j-1),hface(i,j),ctr(i,j),JDIRC)
                end do
            end do
            !$omp end do nowait
            !$omp end parallel 
        end subroutine evolution

        !--------------------------------------------------
        !>update cell averaged values
        !--------------------------------------------------
        subroutine update()
            real(kind=RKD),allocatable,dimension(:,:) :: H_old,B_old !equilibrium distribution at t=t^n
            real(kind=RKD),allocatable,dimension(:,:) :: H,B !equilibrium distribution at t=t^{n+1}
            real(kind=RKD),allocatable,dimension(:,:) :: H_plus,B_plus !Shakhov part
            real(kind=RKD) :: w_old(4) !conservative variables at t^n
            real(kind=RKD) :: prim_old(4),prim(4) !primary variables at t^n and t^{n+1}
            real(kind=RKD) :: tau_old,tau !collision time and t^n and t^{n+1}
            real(kind=RKD) :: sum_res(4),sum_avg(4)
            real(kind=RKD) :: qf(2)
            integer :: i,j

            !allocate arrays
            allocate(H_old(unum,vnum))
            allocate(B_old(unum,vnum))
            allocate(H(unum,vnum))
            allocate(B(unum,vnum))
            allocate(H_plus(unum,vnum))
            allocate(B_plus(unum,vnum))

            !set initial value
            res = 0.0
            sum_res = 0.0
            sum_avg = 0.0

            do j=iymin,iymax
                do i=ixmin,ixmax
                    !--------------------------------------------------
                    !store W^n and calculate H^n,B^n,\tau^n
                    !--------------------------------------------------
                    w_old = ctr(i,j)%w !store W^n
                    
                    prim_old = get_primary(w_old) !convert to primary variables
                    call discrete_maxwell(H_old,B_old,uspace,vspace,prim_old) !calculate Maxwellian
                    tau_old = get_tau(prim_old) !calculate collision time \tau^n

                    !--------------------------------------------------
                    !update W^{n+1} and calculate H^{n+1},B^{n+1},\tau^{n+1}
                    !--------------------------------------------------
                    ctr(i,j)%w = ctr(i,j)%w+(vface(i,j)%flux-vface(i+1,j)%flux+hface(i,j)%flux-hface(i,j+1)%flux)/ctr(i,j)%area !update W^{n+1}

                    prim = get_primary(ctr(i,j)%w)
                    call discrete_maxwell(H,B,uspace,vspace,prim)
                    tau = get_tau(prim)

                    !--------------------------------------------------
                    !record residual
                    !--------------------------------------------------
                    sum_res = sum_res+(w_old-ctr(i,j)%w)**2
                    sum_avg = sum_avg+abs(ctr(i,j)%w)

                    !--------------------------------------------------
                    !Shakhov part
                    !--------------------------------------------------
                    !heat flux at t=t^n
                    qf = get_heat_flux(ctr(i,j)%h,ctr(i,j)%b,uspace,vspace,prim_old) 

                    !h^+ = H+H^+ at t=t^n
                    call shakhov_part(H_old,B_old,uspace,vspace,qf,prim_old,H_plus,B_plus) !H^+ and B^+
                    H_old = H_old+H_plus !h^+
                    B_old = B_old+B_plus !b^+

                    !h^+ = H+H^+ at t=t^{n+1}
                    call shakhov_part(H,B,uspace,vspace,qf,prim,H_plus,B_plus)
                    H = H+H_plus
                    B = B+B_plus

                    !--------------------------------------------------
                    !update distribution function
                    !--------------------------------------------------
                    ctr(i,j)%h = (ctr(i,j)%h+(vface(i,j)%flux_h-vface(i+1,j)%flux_h+hface(i,j)%flux_h-hface(i,j+1)%flux_h)/ctr(i,j)%area+&
                                        0.5*dt*(H/tau+(H_old-ctr(i,j)%h)/tau_old))/(1.0+0.5*dt/tau)
                    ctr(i,j)%b = (ctr(i,j)%b+(vface(i,j)%flux_b-vface(i+1,j)%flux_b+hface(i,j)%flux_b-hface(i,j+1)%flux_b)/ctr(i,j)%area+&
                                        0.5*dt*(B/tau+(B_old-ctr(i,j)%b)/tau_old))/(1.0+0.5*dt/tau)
                end do
            end do

            !final residual
            res = sqrt(ngrid*sum_res)/(sum_avg+SMV)
        end subroutine update

        !--------------------------------------------------
        !>one-sided interpolation of the boundary cell
        !>@param[inout] cell_N :the target boundary cell
        !>@param[inout] cell_L :the left cell
        !>@param[inout] cell_R :the right cell
        !>@param[in]    idx    :the index indicating i or j direction
        !--------------------------------------------------
        subroutine interp_boundary(cell_N,cell_L,cell_R,idx)
            type(cell_center),intent(inout) :: cell_N
            type(cell_center),intent(inout) :: cell_L,cell_R
            integer,intent(in) :: idx

            cell_N%sh(:,:,idx) = (cell_R%h-cell_L%h)/(0.5*cell_R%length(idx)+0.5*cell_L%length(idx))
            cell_N%sb(:,:,idx) = (cell_R%b-cell_L%b)/(0.5*cell_R%length(idx)+0.5*cell_L%length(idx))
        end subroutine interp_boundary

        !--------------------------------------------------
        !>interpolation of the inner cells
        !>@param[in]    cell_L :the left cell
        !>@param[inout] cell_N :the target cell
        !>@param[in]    cell_R :the right cell
        !>@param[in]    idx    :the index indicating i or j direction
        !--------------------------------------------------
        subroutine interp_inner(cell_L,cell_N,cell_R,idx)
            type(cell_center),intent(in) :: cell_L,cell_R
            type(cell_center),intent(inout) :: cell_N
            integer,intent(in) :: idx
            real(kind=RKD),allocatable,dimension(:,:) :: sL,sR

            !allocate array
            allocate(sL(unum,vnum))
            allocate(sR(unum,vnum))

            sL = (cell_N%h-cell_L%h)/(0.5*cell_N%length(idx)+0.5*cell_L%length(idx))
            sR = (cell_R%h-cell_N%h)/(0.5*cell_R%length(idx)+0.5*cell_N%length(idx))
            cell_N%sh(:,:,idx) = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

            sL = (cell_N%b-cell_L%b)/(0.5*cell_N%length(idx)+0.5*cell_L%length(idx))
            sR = (cell_R%b-cell_N%b)/(0.5*cell_R%length(idx)+0.5*cell_N%length(idx))
            cell_N%sb(:,:,idx) = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)
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
            real(kind=RKD) :: init_gas(4) !initial condition
            real(kind=RKD) :: alpha !another index in molecule model
            real(kind=RKD) :: kn !Knudsen number in reference state
            real(kind=RKD) :: xlength,ylength !length of computational domain in x and y direction
            integer :: xnum,ynum !number of cells in x and y direction

            !control
            cfl = 0.8 !CFL number
            eps = 1.0E-5 !convergence criteria

            !gas
            ck = 1 !internal degree of freedom
            gam = get_gamma(ck) !ratio of specific heat
            pr = 2.0/3.0 !Prandtl number
            omega = 0.81 !temperature dependence index in VHS model

            !reference viscosity coefficient
            alpha = 1.0 !another index in VHS model
            kn = 1.0 !Knudsen number in reference state
            !mu_ref = get_mu(kn,alpha,omega) !reference viscosity coefficient
            mu_ref = get_mu(kn,real(1.0,RKD),real(0.5,RKD)) !reference viscosity coefficient

            !geometry
            xlength = 1.0
            ylength = 1.0
            xnum = 61
            ynum = 61

            !initial condition (density,u-velocity,v-velocity,lambda=1/temperature)
            init_gas = [1.0, 0.0, 0.0, 1.0]

            !set boundary condition (density,u-velocity,v-velocity,lambda)
            bc_W = [1.0, 0.0, 0.0, 1.0] !west
            bc_E = [1.0, 0.0, 0.0, 1.0] !east
            bc_S = [1.0, 0.0, 0.0, 1.0] !south
            bc_N = [1.0, 0.15, 0.0, 1.0] !north

            call init_geometry(xlength,ylength,xnum,ynum) !initialize the geometry
            call init_velocity() !initialize discrete velocity space using Gaussian-Hermite quadrature
            call init_allocation() !allocate arrays
            call init_flow_field(init_gas) !set the initial condition
        end subroutine init

        !--------------------------------------------------
        !>initialize the geometry
        !>@param[in] xnum,ynum       :number of cells in x and y direction
        !>@param[in] xlength,ylength :domain length in x and y direction
        !--------------------------------------------------
        subroutine init_geometry(xlength,ylength,xnum,ynum)
            real(kind=RKD),intent(in) :: xlength,ylength
            integer,intent(in) :: xnum,ynum
            real(kind=RKD) :: dx,dy !cell length in x and y direction
            real(kind=RKD) :: area !cell area
            integer :: i,j

            !cell index range
            ixmin = 1
            iymin = 1
            ixmax = xnum
            iymax = ynum

            !total number of cell
            ngrid = (ixmax-ixmin+1)*(iymax-iymin+1)

            !allocation
            allocate(ctr(ixmin:ixmax,iymin:iymax)) !cell center
            allocate(vface(ixmin:ixmax+1,iymin:iymax),hface(ixmin:ixmax,iymin:iymax+1)) !vertical and horizontal cell interface

            !cell length and area
            dx = xlength/(ixmax-ixmin+1)
            dy = ylength/(iymax-iymin+1)
            area = dx*dy

            forall(i=ixmin:ixmax,j=iymin:iymax) !cell center
                ctr(i,j)%x = (i-0.5)*dx
                ctr(i,j)%y = (j-0.5)*dy
                ctr(i,j)%length(1) = dx
                ctr(i,j)%length(2) = dy
                ctr(i,j)%area = area
            end forall

            forall(i=ixmin:ixmax+1,j=iymin:iymax) !vertical interface
                vface(i,j)%length = dy
                vface(i,j)%cosx = 1.0
                vface(i,j)%cosy = 0.0
            end forall

            forall(i=ixmin:ixmax,j=iymin:iymax+1) !horizontal interface
                hface(i,j)%length = dx
                hface(i,j)%cosx = 0.0
                hface(i,j)%cosy = 1.0
            end forall
        end subroutine init_geometry

        !--------------------------------------------------
        !>set discrete velocity space using Gaussian-Hermite quadrature
        !--------------------------------------------------
        subroutine init_velocity()
            real(kind=RKD) :: vcoords(28), weights(28) !velocity points and weight for 28 points (symmetry)
            integer :: i,j

            !set velocity points and weight
            vcoords = [ -6.59160566329956050,-5.85701465606689450,-5.24328517913818360,-4.69075632095336910,&
                        -4.17663669586181640,-3.68913412094116210,-3.22111201286315920,-2.76779532432556150,&
                        -2.32574987411499020,-1.89236044883728030,-1.46553730964660640,-1.04353523254394530,&
                        -0.62483674287796021,-0.20806738734245300,+0.20806738734245300,+0.62483674287796021,&
                        +1.04353523254394530,+1.46553730964660640,+1.89236044883728030,+2.32574987411499020,&
                        +2.76779532432556150,+3.22111201286315920,+3.68913412094116210,+4.17663669586181640,&
                        +4.69075632095336910,+5.24328517913818360,+5.85701465606689450,+6.59160566329956050]

            weights = [ +0.84476017951965332,+0.65798896551132202,+0.57779419422149658,+0.53077453374862671,&
                        +0.49934458732604980,+0.47681638598442078,+0.46000820398330688,+0.44718948006629944,&
                        +0.43733271956443787,+0.42979142069816589,+0.42414394021034241,+0.42011165618896484,&
                        +0.41751345992088318,+0.41624009609222412,+0.41624009609222412,+0.41751345992088318,&
                        +0.42011165618896484,+0.42414394021034241,+0.42979142069816589,+0.43733271956443787,&
                        +0.44718948006629944,+0.46000820398330688,+0.47681638598442078,+0.49934458732604980,&
                        +0.53077453374862671,+0.57779419422149658,+0.65798896551132202,+0.84476017951965332]

            !set grid number for u-velocity and v-velocity
            unum = 28
            vnum = 28

            !allocate discrete velocity space
            allocate(uspace(unum,vnum)) !x direction
            allocate(vspace(unum,vnum)) !y direction
            allocate(weight(unum,vnum)) !weight at u_k and v_l

            !set velocity space and weight
            forall(i=1:unum,j=1:vnum)
                uspace(i,j) = vcoords(i)
                vspace(i,j) = vcoords(j)
                weight(i,j) = weights(i)*weights(j)
            end forall

            !store the maximum micro velocity
            umax = maxval(abs(uspace(:,1)))
            vmax = maxval(abs(vspace(1,:)))
        end subroutine init_velocity

        !--------------------------------------------------
        !>allocate arrays
        !--------------------------------------------------
        subroutine init_allocation()
            integer :: i,j

            !cell center
            do j=iymin,iymax
                do i=ixmin,ixmax
                    allocate(ctr(i,j)%h(unum,vnum))
                    allocate(ctr(i,j)%b(unum,vnum))
                    allocate(ctr(i,j)%sh(unum,vnum,2))
                    allocate(ctr(i,j)%sb(unum,vnum,2))
                end do
            end do

            !cell interface
            do j=iymin,iymax
                do i=ixmin,ixmax+1
                    allocate(vface(i,j)%flux_h(unum,vnum))
                    allocate(vface(i,j)%flux_b(unum,vnum))
                end do
            end do

            do j=iymin,iymax+1
                do i=ixmin,ixmax
                    allocate(hface(i,j)%flux_h(unum,vnum))
                    allocate(hface(i,j)%flux_b(unum,vnum))
                end do
            end do
        end subroutine init_allocation

        !--------------------------------------------------
        !>set the initial condition
        !>@param[in] init_gas :initial condition
        !--------------------------------------------------
        subroutine init_flow_field(init_gas)
            real(kind=RKD),intent(in) :: init_gas(4)
            real(kind=RKD),allocatable,dimension(:,:) :: H,B !reduced Maxwellian distribution functions
            real(kind=RKD) :: w(4) !conservative variables
            integer :: i,j

            !allocate array
            allocate(H(unum,vnum))
            allocate(B(unum,vnum))

            !convert primary variables to conservative variables
            w = get_conserved(init_gas)

            !obtain discretized Maxwellian distribution H and B
            call discrete_maxwell(H,B,uspace,vspace,init_gas)

            !initial condition
            forall(i=ixmin:ixmax,j=iymin:iymax)
                ctr(i,j)%w = w
                ctr(i,j)%h = H
                ctr(i,j)%b = B
            end forall
        end subroutine init_flow_field

        !--------------------------------------------------
        !>write result
        !--------------------------------------------------
        subroutine output()
            real(kind=RKD) :: prim(4) !primary variables
            real(kind=RKD) :: pressure,temperature,qf(2) !pressure,temperature,heat flux
            integer :: i,j

            !open result file
            open(unit=RSTFILE,file=RSTFILENAME,status="replace",action="write")

            !write header
            write(RSTFILE,*) "VARIABLES = X, Y, RHO, U, V, T, P, QX, QY"
            write(RSTFILE,*) "ZONE  I = ",ixmax-ixmin+1,", J = ",iymax-iymin+1

            !write solution
            do j=iymin,iymax
                do i=ixmin,ixmax
                    prim = get_primary(ctr(i,j)%w)
                    pressure = get_pressure(ctr(i,j)%h,ctr(i,j)%b,uspace,vspace,prim)
                    temperature = 2.0*pressure/prim(1)
                    qf = get_heat_flux(ctr(i,j)%h,ctr(i,j)%b,uspace,vspace,prim)

                    write(RSTFILE,*) ctr(i,j)%x,ctr(i,j)%y,prim(1),prim(2),prim(3),temperature,pressure,qf(1),qf(2)
                end do
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
    write(HSTFILE,*) "VARIABLES = iter, sim_time, dt, res_rho, res_ru, res_rv, res_re" !write header

    !iteration
    do while(.true.)
        call timestep() !calculate time step
        call interpolation() !calculate the slope of distribution function
        call evolution() !calculate flux across the interfaces
        call update() !update cell averaged value

        !check if exit
        if (all(res<eps)) exit

        !write iteration situation every 10 iterations
        !if (mod(iter,10)==0) then
            write(*,"(A18,I15,2E15.7)") "iter,sim_time,dt:",iter,sim_time,dt
            write(*,"(A18,4E15.7)") "res:",res
            write(HSTFILE,"(I15,6E15.7)") iter,sim_time,dt,res
        !end if

        iter = iter+1
        sim_time = sim_time+dt
    end do

    !close history file
    close(HSTFILE)

    !output solution
    call output()
end program main

! vim: set ft=fortran tw=0: 
