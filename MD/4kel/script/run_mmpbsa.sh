#!/bin/bash
mpirun -np 12 MMPBSA.py.MPI -O -i mmpbsa.in -o MMPBSA1.dat -sp com_sol.prmtop -cp com.prmtop -rp rec.prmtop -lp lig.prmtop -y md1.nc
mpirun -np 12 MMPBSA.py.MPI -O -i mmpbsa.in -o MMPBSA2.dat -sp com_sol.prmtop -cp com.prmtop -rp rec.prmtop -lp lig.prmtop -y md2.nc
mpirun -np 12 MMPBSA.py.MPI -O -i mmpbsa.in -o MMPBSA3.dat -sp com_sol.prmtop -cp com.prmtop -rp rec.prmtop -lp lig.prmtop -y md3.nc
