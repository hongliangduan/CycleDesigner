#!/bin/bash
bash AmberMdPrep.sh -p com_sol.prmtop -c com_sol.inpcrd --temp 300
pmemd.cuda -O -i md1.in -p com_sol.prmtop -c final.1.ncrst -o md1.out -r md1.rst7 -x md1.nc
bash AmberMdPrep.sh -p com_sol.prmtop -c com_sol.inpcrd --temp 300
pmemd.cuda -O -i md1.in -p com_sol.prmtop -c final.1.ncrst -o md2.out -r md2.rst7 -x md1.nc
bash AmberMdPrep.sh -p com_sol.prmtop -c com_sol.inpcrd --temp 300
pmemd.cuda -O -i md1.in -p com_sol.prmtop -c final.1.ncrst -o md3.out -r md3.rst7 -x md3.nc
