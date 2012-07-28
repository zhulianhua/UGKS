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
    integer :: iter !iteration

    !--------------------------------------------------
    !gas properties
    !--------------------------------------------------
    real(kind=RKD) :: kn !Knudsen number
    real(kind=RKD) :: gam !ratio of specific heat
    real(kind=RKD) :: omega !temperature dependence index
    real(kind=RKD) :: pr !Prandtl number
    integer :: ck !internal degree of freedom

    !--------------------------------------------------
    !macros for a readable code
    !--------------------------------------------------
    !direction
    integer,parameter :: WEST = 1
    integer,parameter :: EAST = 2 
    integer,parameter :: SOUTH = 3
    integer,parameter :: NORTH = 4
    integer,parameter :: IDIRC = 1 !i direction
    integer,parameter :: JDIRC = 2 !j direction

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

    !discrete velocity space
    type :: quad
        real(kind=RKD) :: s !speed
        real(kind=RKD) :: w !weight
    end type quad
    
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
    !discrete velocity space
    !--------------------------------------------------
    integer :: unum,vnum !number of velocity points for u and v
    real(kind=RKD) :: umax,vmax !maximum micro velocity
    type(quad),allocatable,dimension(:) :: uspace,vspace !u and v discrete velocity space
end module global_data

!--------------------------------------------------
!>define some commonly used functions/subroutines
!--------------------------------------------------
module tools
    use global_data
    implicit none
    contains
        !--------------------------------------------------
        !>convert primary variables to conservative variables
        !>@param[in] prim primary variables
        !>@param[out] w conservative variables
        !--------------------------------------------------
        subroutine convert_prim_w(prim,w)
            real(kind=RKD),intent(in) :: prim(4)
            real(kind=RKD),intent(out) :: w(4)

            w(1) = prim(1)
            w(2) = prim(1)*prim(2)
            w(3) = prim(1)*prim(3)
            w(4) = 0.5*prim(1)/prim(4)/(gam-1.0)+0.5*prim(1)*(prim(2)**2+prim(3)**2)
        end subroutine convert_prim_w

        !--------------------------------------------------
        !>obtain discretized Maxwellian distribution
        !>@param[in] prim primary variables in global frame
        !>@param[out] h,b reduced Maxwellian distribution
        !--------------------------------------------------
        subroutine discrete_maxwell(prim,h,b)
            real(kind=8),intent(in) :: prim(4)
            real(kind=8),dimension(:,:),intent(out) :: h,b
            real(kind=8) :: temp
            integer :: i,j

            temp = prim(1)*(prim(4)/PI)
            forall(i=1:unum,j=1:vnum)
                h(i,j) = temp*exp(-prim(4)*((uspace(i)%s-prim(2))**2+(vspace(j)%s-prim(3))**2))
                b(i,j) = h(i,j)*ck/(2.0*prim(4))
            end forall
        end subroutine discrete_maxwell

        !--------------------------------------------------
        !>convert macro variables from global frame to local
        !>@param[in] w_local macro variables in local frame
        !>@param[out] w_global macro variables in global frame
        !>@param[in] cosx,cosy directional cosine
        !--------------------------------------------------
        subroutine global_frame(w_local,w_global,cosx,cosy)
            real(kind=8),intent(in) :: w_local(2) !variables in local frame
            real(kind=8),intent(out) :: w_global(2) !variables in global frame
            real(kind=8),intent(in) :: cosx,cosy !directional cosine
            real(kind=8) :: wx,wy !variables in x and y direction

            w_global(1) = w_local(1)*cosx-w_local(2)*cosy
            w_global(2) = w_local(1)*cosy+w_local(2)*cosx
        end subroutine global_frame

        !--------------------------------------------------
        !>obtain ratio of specific heat
        !>@param[in] ck internal degree of freedom
        !>@return get_gamma ratio of specific heat
        !--------------------------------------------------
        function get_gamma(ck)
            integer,intent(in) :: ck
            real(kind=RKD) :: get_gamma

            get_gamma = float(ck+4)/float(ck+2)
        end function get_gamma
end module tools

!--------------------------------------------------
!>UGKS solver
!--------------------------------------------------
module solver
    use global_data
    use tools
    implicit none
    contains
        !--------------------------------------------------
        !>calculate time step
        !--------------------------------------------------
        subroutine timestep()
            real(kind=8) :: tmax !max 1/dt allowed
            real(kind=8) :: sos !speed of sound
            real(kind=8) :: prim(4) !primary variables
            integer :: i,j

            !set initial value
            tmax = 0.0

            !$omp parallel 
            !$omp do private(i,j,sos,prim) reduction(max:tmax)
            do j=iymin,iymax
                do i=ixmin,ixmax
                    !convert conservative variables to primary variables
                    call convert_w_prim(ctr(i,j)%w,prim)

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
                    call interp_boundary(ctr(ixmin,j),ctr(ixmin,j),ctr(ixmin+1,j),1)
                    call interp_boundary(ctr(ixmax,j),ctr(ixmax-1,j),ctr(ixmax,j),1)
                end do
                !$omp end do nowait

                !$omp do
                do j=ixmin,ixmax
                    do i=ixmin+1,ixmax-1
                        call interp_inner(ctr(i-1,j),ctr(i,j),ctr(i+1,j),1)
                    end do
                end do
                !$omp end do nowait

                !--------------------------------------------------
                !j direction
                !--------------------------------------------------
                !$omp do
                do i=ixmin,ixmax
                    call interp_boundary(ctr(i,iymin),ctr(i,iymin),ctr(i,iymin+1),2)
                    call interp_boundary(ctr(i,iymax),ctr(i,iymax-1),ctr(i,iymax),2)
                end do
                !$omp end do nowait

                !$omp do
                do j=iymin+1,iymax-1
                    do i=ixmin,ixmax
                        call interp_inner(ctr(i,j-1),ctr(i,j),ctr(i,j+1),2)
                    end do
                end do
                !$omp end do nowait
                !$omp end parallel
        end subroutine interpolation

        !--------------------------------------------------
        !>one-sided interpolation of the boundary cell
        !>@param[in] cell_N the target boundary cell
        !>@param[in] cell_L the left cell
        !>@param[in] cell_R the right cell
        !>@param[in] idx  the index indicating i or j direction
        !--------------------------------------------------
        subroutine interp_boundary(cell_N,cell_L,cell_R,idx)
            type(cell_center),intent(in) :: cell_N,cell_L,cell_R
            integer,intent(in) :: idx

            cell_N%sh(:,:,idx) = (cell_R%h-cell_L%h)/(0.5*cell_R%length(idx)+0.5*cell_L%length(idx))
            cell_N%sb(:,:,idx) = (cell_R%b-cell_L%b)/(0.5*cell_R%length(idx)+0.5*cell_L%length(idx))
        end subroutine interp_boundary

        !--------------------------------------------------
        !>interpolation of the inner cells
        !>@param[in] cell_L the left cell
        !>@param[in] cell_N the target cell
        !>@param[in] cell_R the right cell
        !>@param[in] idx  the index indicating i or j direction
        !--------------------------------------------------
        subroutine interp_inner(cell_L,cell_N,cell_R,idx)
            type(cell_center),intent(in) :: cell_N,cell_L,cell_R
            integer,intent(in) :: idx
            real(kind=8),allocatable,dimension(:,:) :: sL,sR

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
!>flux calculation
!--------------------------------------------------
module flux
end module flux

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
            real(kind=RKD) :: xlength,ylength !length of computational domain in x and y direction
            integer :: xnum,ynum !number of cells in x and y direction

            !control
            cfl = 0.8 !CFL number
            eps = 1.0E-5 !convergence criteria
            !gas
            ck = 1 !internal degree of freedom
            kn = 1.0 !Knudsen number
            omega = 0.81 !temperature dependence index in HS/VHS model
            pr = 2.0/3.0 !Prandtl number
            gam = get_gamma(ck) !ratio of specific heat
            !geometry
            xlength = 1.0
            ylength = 1.0
            xnum = 61
            ynum = 61
            !initial condition (density,u-velocity,v-velocity,lambda=1/temperature)
            init_gas = [1.0, 0.0, 0.0, 1.0]
            !set boundary condition (density,u-velocity,v-velocity,lambda)
            bc_W = [1.0, 0.00, 0.0, 1.0] !west
            bc_E = [1.0, 0.00, 0.0, 1.0] !east
            bc_S = [1.0, 0.00, 0.0, 1.0] !south
            bc_N = [1.0, 0.15, 0.0, 1.0] !north

            call init_geometry(xlength,ylength,xnum,ynum) !initialize the geometry
            call init_velocity() !initialize discrete velocity space using Gaussian-Hermite quadrature
            call init_allocation() !allocate arrays
            call init_flow_field(init_gas) !set the initial condition
        end subroutine init

        !--------------------------------------------------
        !>initialize the geometry
        !>@param[in] xnum,ynum number of cells in x and y direction
        !>@param[in] xlength,ylength domain length in x and y direction
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
            real(kind=RKD) :: coords(14), weight(14) !half space velocity points and weight for 28 points (symmetry)
            integer :: i

            !set velocity points and weight
            coords = [0.208067382691,0.624836719505,1.043535273750,1.465537263460,1.892360496840,2.325749842660,2.76779535291,3.221112076560,3.689134238460,4.176636742130,4.690756523940,5.243285373200,5.857014641380,6.591605442370]
            weight = [0.416240109246,0.417513453286,0.420111646094,0.424143944203,0.429791424953,0.437332725764,0.44718947117,0.460008206624,0.476816375219,0.499344597245,0.530774551498,0.577794173639,0.657988990298,0.844760204047]

            !set grid number for u-velocity and v-velocity
            unum = 28
            vnum = 28

            !allocate discrete velocity space
            allocate(uspace(unum)) !x direction
            allocate(vspace(vnum)) !y direction

            !set grid points and weight for u-velocity
            do i=unum/2+1,unum
                uspace(i)%s = coords(i-unum/2)
                uspace(i)%w = weight(i-unum/1)
                uspace(unum+1-i)%s = -coords(i-unum/2)
                uspace(unum+1-i)%w = weight(i-unum/2)
            end do

            !set grid points and weight for v-velocity (the same)
            vspace = uspace

            !store the maximum micro velocity
            umax = maxval(abs(uspace%s))
            vmax = maxval(abs(vspace%s))
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
        !>@param[in] init_gas initial condition
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
            call convert_prim_w(init_gas,w)

            !obtain discretized Maxwellian distribution H and B
            call discrete_maxwell(init_gas,H,B)

            !initial condition
            forall(i=ixmin:ixmax,j=iymin:iymax)
                ctr(i,j)%w = w
                ctr(i,j)%h = H
                ctr(i,j)%b = B
            end forall
        end subroutine init_flow_field
end module io

!--------------------------------------------------
!>main program
!--------------------------------------------------
program main
    use global_data
    use tools
    use solver
    use io
    implicit none

    !initialization
    call init() 

    !set initial value
    iter = 1 !number of iteration
    sim_time = 0.0 !simulation time

    !iteration
    do while(.true.)
        call timestep() !calculate time step
        call interpolation() !calculate the slope of distribution function
        call evolution() !calculate flux across the interfaces
        call update() !update cell averaged value

        !check if exit
        if (all(res<eps)) exit

        !write iteration situation every 10 iterations
        if (mod(iter,10)==0) then
            write(*,"(A18,I15,2E15.7)") "iter,sim_time,dt:",iter,sim_time,dt
            write(*,"(A18,4E15.7)") "res:",res
        end if

        iter = iter+1
        sim_time = sim_time+dt
    end do
end program main
