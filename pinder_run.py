from pinder.core import get_pinder_location
get_pinder_location()
from pinder.core import get_index, get_pinder_location
from pinder.core.index.utils import IndexEntry

index = get_index()
# Define UniProt IDs
prot1 = "P78556"  # CCL19
prot2 = "P32248"  # CCR7

prot1 = "P61812" #TGFB2
#prot2 = "Q03167" #TGFBR3
prot2 = "P37173" #TGFBR2

hits = index[((index["uniprot_L"] == prot1) & (index["uniprot_R"] == prot2))]
hits_L = index[(index["uniprot_L"] == prot1)]
hits_R = index[(index["uniprot_R"] == prot2)]


hits['holo_R_pdb']
hits['holo_L_pdb']

hits['apo_R_pdb']
hits['apo_L_pdb']



hits_L["predicted_L_pdb"]
hits_R["predicted_R_pdb"]

from pinder.scoring import score_holo
result = score_holo("unrelaxed_model_1_multimer_v3_pred_0.pdb")
print(result)
