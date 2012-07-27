!-----------------------------------------------------------------------------
!UGKS1D - shock structure calculation at M=8 using Unified Gas-Kinetic Scheme
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
!store the global variables
!--------------------------------------------------
module global_data
end module global_data

!--------------------------------------------------
!define some commonly used functions/subroutines
!--------------------------------------------------
module tools
end module tools

!--------------------------------------------------
!UGKS solver
!--------------------------------------------------
module solver
end module solver

!--------------------------------------------------
!input and output
!--------------------------------------------------
module io
    use global_data
    use tools
    implicit none
    contains
        !--------------------------------------------------
        !initialization
        !--------------------------------------------------
        subroutine init()
            !main initialization subroutine
        end subroutine init

end module io

!--------------------------------------------------
!main program
!--------------------------------------------------
program main
    use global_data
    use tools
    use solver
    use io
    implicit none

    !///initialization>
    call init() 
end program main

! vim: set ft=fortran tw=78: 
