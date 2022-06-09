from rdkit import Chem
from rdkit.Chem import AllChem

def generate_monomer_sdf():
    naa_symbols = ["G", "A", "R", "N", "D", "C", "Q", "E", "H", "I",\
                 "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    nnaa_symbols = ["Me"+char for char in naa_symbols if char != "P"]
    print("Number of natural amino acids ", len(naa_symbols))
    print("Number of non-natural amino acids ", len(nnaa_symbols))

    naa_mols = []
    for symbol in naa_symbols:
        mol = Chem.MolFromSequence(symbol)
        if mol:
            naa_mols.append(mol)

    # N-methylation
    rxn = AllChem.ReactionFromSmarts("OC([C:1][NX3;H2:2])=O>>OC([C:1][NX3;H1:2]C)=O")
    nnaa_mols = []
    for mol in naa_mols:
        try:
            ps = rxn.RunReactants([mol])
            nnaa = ps[0][0]
            if nnaa:
                nnaa_mols.append(nnaa)
        except:
            # Proline fails but no problem
            pass

    test_mols = naa_mols + nnaa_mols
    test_symbols = naa_symbols + nnaa_symbols
    w = Chem.SDWriter("test_monomers.sdf")
    for mol, symbol in zip(test_mols, test_symbols):
        mol.SetProp("symbol", symbol)
        w.write(mol)
    w.close()


if __name__ == "__main__":
    generate_monomer_sdf()    