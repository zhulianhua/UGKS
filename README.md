Description
================
This program provides two codes for solving gas flow problem based on the
Unified Gas-Kinetic Scheme. For detail, please refer to the manual in doc
directory.

For the source file

1. UGKS1D.f90 is for 1D shock structure calculation at Ma=8
2. UGKS2D.f90 is for 2D lid-driven cavity flow calculation

Pre-requirements
================
1. Fortran compiler: `ifort` or `gfortran` supporting Fortran 2003
2. Latex installation: only for re-compilation of the manual, requiring
`hyperref`, `parskip`, `amsmath`, `amssymb`, `fullpage` and `appendix`
packages. Also requires `bibtex`, `dvips` and `ps2pdf`

Usage
================
The Makefile is provided for those who have gnu make installed. If you are
using any IDE (e.g. Visual Studio), use the compiling function provided by the
IDE

Note: openmp and Intel Fortran compiler is used by default. **DO NOT** type the
prompt symbol $

1. Make both 1D and 2D program

        $ make

2. Make 1D only

        $ make 1D

3. Make 2D only

        $ make 2D

4. Make both 1D and 2D **WITHOUT** openmp, and **WITH** gfortran

        $ make OMP=no FC=gfortran

5. Re-generate manual

        $ make manual

6. Clean
 
        $ make clean

License
================
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
