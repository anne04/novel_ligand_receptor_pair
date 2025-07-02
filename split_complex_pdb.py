from Bio.PDB import PDBParser, PDBIO, Select

class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.id == self.chain_id

# Load the unrelaxed complex
ligand_name = 'TGFB2' #'RETN'  #'CDH1'
rec_name = 'TGFBR2' #'CD44' #'ITGB6'
pdb_path =  'TGFB2_TGFBR2_unrelaxed_model_1_multimer_v3_pred_0.pdb' #"RETN_CD44_unrelaxed_model_1_multimer_v3_pred_0.pdb" #"/cluster/projects/schwartzgroup/fatema/LRbind/ParallelFold-main/output/old_lrbind/lrbind_CDH1_ITGB6/unrelaxed_model_1_multimer_v3_pred_0.pdb"
parser = PDBParser(QUIET=True)
structure = parser.get_structure("complex", pdb_path)

for model in structure:
    for chain in model:
        print("Chain", chain.id)
        for res in chain:
            print(res.resname, end=" ")
        print("\n---")

# Save chain A (example: ligand)
io = PDBIO()
io.set_structure(structure)
io.save(ligand_name+"_holo.pdb", select=ChainSelect("B"))

# Save chain B (example: receptor)
io.save(rec_name+"_holo.pdb", select=ChainSelect("C"))
