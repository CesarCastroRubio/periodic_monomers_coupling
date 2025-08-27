import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolops

# adjust sys.path to include tools directory
script_dir = os.path.dirname(__file__)
sys.path.append(os.path.normpath(os.path.join(script_dir, '..', 'tools')))

class Polymerizer:
    def __init__(self, p_smiles, degree=2):
        self.p_smiles = p_smiles
        self.degree = degree

    def preprocess_monomer(self, mol):
        rw_mol = Chem.RWMol(mol)
        star_atoms = [atom for atom in rw_mol.GetAtoms() if atom.GetSymbol() == '*']
        if len(star_atoms) < 2:
            raise ValueError("Monomer must have at least two connection points denoted by [*]")
        star_atoms[0].SetAtomicNum(29)
        star_atoms[1].SetAtomicNum(29)
        return rw_mol

    def polymerize(self):
        mol = Chem.MolFromSmiles(self.p_smiles)
        if mol is None:
            raise ValueError("Invalid p-SMILES string: %s" % self.p_smiles)

        monomer = self.preprocess_monomer(mol)
        degree = self.degree 

        copies = []
        connection_info = []

        for _ in range(degree):
            copy_mol = Chem.Mol(monomer)
            dummy_atoms = [atom.GetIdx() for atom in copy_mol.GetAtoms() if atom.GetAtomicNum() == 29]
            if len(dummy_atoms) != 2:
                raise ValueError("Monomer must have exactly two dummy atoms (atomic num 29)")
            head_dummy, tail_dummy = sorted(dummy_atoms)

            head_neighbors = copy_mol.GetAtomWithIdx(head_dummy).GetNeighbors()
            tail_neighbors = copy_mol.GetAtomWithIdx(tail_dummy).GetNeighbors()
            if len(head_neighbors) != 1 or len(tail_neighbors) != 1:
                raise ValueError("Each dummy atom must be connected to exactly one atom")

            head_neighbor = head_neighbors[0].GetIdx()
            tail_neighbor = tail_neighbors[0].GetIdx()
            connection_info.append((head_dummy, head_neighbor, tail_dummy, tail_neighbor))
            copies.append(copy_mol)

        combined = copies[0]
        offsets = [0]
        for mol_copy in copies[1:]:
            offsets.append(combined.GetNumAtoms())
            combined = Chem.CombineMols(combined, mol_copy)

        emol = Chem.EditableMol(combined)
        for i in range(degree - 1):
            tail_idx = connection_info[i][3] + offsets[i]
            head_idx = connection_info[i+1][1] + offsets[i+1]
            emol.AddBond(tail_idx, head_idx, Chem.rdchem.BondType.SINGLE)

        polymer_mol = emol.GetMol()
        dummy_idxs = sorted([a.GetIdx() for a in polymer_mol.GetAtoms() if a.GetAtomicNum() == 29], reverse=True)
        rw = Chem.RWMol(polymer_mol)
        for idx in dummy_idxs[1:-1]:
            rw.RemoveAtom(idx)

        final = rw.GetMol()
        Chem.SanitizeMol(final)
        p_smiles_e = Chem.MolToSmiles(final).replace('Cu', '*')
        #print(p_smiles_e)
        return p_smiles_e

if __name__ == '__main__':
    import sys
    p_smiles = sys.argv[1]
    degree = int(sys.argv[2])
    polymerizer = Polymerizer(p_smiles, degree=degree)
    try:
        print(polymerizer.polymerize())
    except ValueError as e:
        print(f"Error: {e}")
