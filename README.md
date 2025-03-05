# CycleDesigner
![image](https://github.com/user-attachments/assets/68a7e619-7e37-4f83-b780-7b80ebe14438)
CycleDesigner is a novol workflow for de novo cyclic peptide design. The modified powerful RFdiffusion model(CycRFdiffusion) allows the generation of cyclic peptide backbone and integrated it with ProteinMPNN and HighFold to design binders for specific targets. 

## Instructions
### Installation
You can install the [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) at first, and then copy the source codes of CycleDesigner into the installed RFdifusion preject.
### Run inference
To run the generation, type
```
./rfd_zch_inference.sh /YOUR/PDB/PATH/pdb.pdb
```
The inference adapt for the targets with binders only, if there is no binders in pdb file, please type follows to run the generation.
```
./scripts/run_inference.py 'contigmap.contigs=[A1-n/0 a-a]' 'ppi.hotspot_res=[A1,A2]' inference.deterministic=Ture inference.num_designs=5 inference.model_directory_path=/YOUR/MODEL/PATH inference.input_pdb=/YOUR/PDB/PATH/pdb.pdb inference.output_prefix=/YOUR/OUTPUT/CYCPEP/NAME > /PATH/NAME.log 2>&1
```
### MD traj files
The output trajectory files of molecular dynamics simulation prosess have uploaded at [Zenodo](https://zenodo.org/records/14955049).
