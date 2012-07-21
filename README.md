Description
===============
This program provides two codes for solving gas flow problem based on the
Unified Gas-Kinetic Scheme. For detail, please refer to

1. K. Xu and J.C. Huang, "A unified gas-kinetic scheme for continuum and
rarefied flows", J. Comput. Physics, vol. 229 (2010), pp. 7747-7764 (October).
2. J.C. Huang, K. Xu, and P.B. Yu, "A Unified Gas-kinetic Scheme for Continuum
and Rarefied Flows II: Multi-dimensional Cases", Communications in
Computational Physics, vol. 3, No. 3, pp. 662-690, September (2012).

For the source file

1. UGKS1D.f90 is for 1D shock structure calculation at Ma=8
2. UGKS2D.f90 is for 2D lid-driven cavity flow calculation

Usage
===============
The Makefile is provided for those who have gnu make installed. If you are
using any IDE (e.g. Visual Studio), use the compiling function provided by the
IDE

Note: openmp and intel fortran compiler is used by default.

1. make both 1D and 2D program

    make all

2. make 1D only

    make 1D

3. make 2D only

    make 2D

4. make both 1D and 2D WITHOUT openmp, and WITH gfortran

    make all OMP=no FC=gfortran

5. clean
    
    make clean

LICENSE
===============
Copyright (C) 2012 Wang Ruijie <lainme993@gmail.com>

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

/* vim: set ft=md tw=78: */
