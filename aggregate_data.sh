#!/bin/bash
rm directory_key_smiles.dat
for dir in $(ls -d ../l1* | cut -d"/" -f2 | sort -V) $(ls -d ../l2* | cut -d"/" -f2 | sort -V) $(ls -d ../l3* | cut -d"/" -f2 | sort -V) 
do
	grep "smiles" ../$dir/nematic_dft.py | head -1 | cut -d"\"" -f2 | awk -v dir=$dir '{print dir,$1}' >> directory_key_smiles.dat
	#smiles=$(grep "smiles" ../$dir/nematic_dft.py | head -1 | cut -d"\"" -f2)
	#python polymerize_dp.py "$smiles" 2 | awk -v dir=$dir '{print dir,$1}' >> directory_key_smiles_dimerized.dat
	#mkdir $dir
	#cp ../$dir/minim*xyz $dir/
done
