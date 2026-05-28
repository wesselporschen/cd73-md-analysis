import os
import numpy as np
from pathlib import Path

import MDAnalysis as mda 
from MDAnalysis.analysis import rms, align, diffusionmap 

from trajectory_hull_class import Trajectory


BINDING_SITE_RESIDUES = [
    (366, "A"), (456, "B"), (543, "A"), (543, "B"), (484, "A"), (484, "B"), (547, "B"),     # Major H-bond ints in top hits, Fig 6 Rahimova 2018 and text  
    (473, "A"), (473, "B"), (545, "A"), (545, "B"), (471, "A"),                             # Charged residues with non bonding ints, described in text
    (455, "A"), (542, "A"),                                                                 # H-bond ints in Fig 6 with RR6, RR4
    (453, "A"), (453, "B"), (454, "A"), (454, "B"), (544, "A"), (544, "B"),                 # Glycines hydrophobic described in text, Fig 7
    (366, "A"), (531, "A"),                                                                 # D366 makes H-bond, Y531 Pi-stack (Grosjean 2025)    
]

CWD = Path(os.getcwd()).parent

paths = {
        '4h2i': CWD / '4h2i' / 'unbiased',
        '6twa': CWD / '6twa' / 'unbiased',
}

allosteric_bindingsite_trajectories = {
        '4h2i': Trajectory(paths['4h2i'], paths['4h2i'] / 'nosol.pdb', BINDING_SITE_RESIDUES),
        '6twa': Trajectory(paths['6twa'], paths['6twa'] / 'nosol.pdb', BINDING_SITE_RESIDUES),
}


#for stepsize in [1000, 200, 100, 50, 10]:
#    print(stepsize)
#    matrix_4h2i_bindingsite = diffusionmap.DistanceMatrix(allosteric_bindingsite_trajectories['4h2i'].residue_atoms)
#    matrix_4h2i_bindingsite.run(start=1, stop=50000, step=stepsize)
#    np.save(f"matrix_4h2i_bindingsite_{stepsize}", matrix_4h2i_bindingsite.results.dist_matrix)
#
    
for stepsize in [1000, 200, 100, 50, 10]:
    print(stepsize)
    matrix_6twa_bindingsite = diffusionmap.DistanceMatrix(allosteric_bindingsite_trajectories['6twa'].residue_atoms)
    matrix_6twa_bindingsite.run(start=1, stop=50000, step=stepsize)
    np.save(f"matrix_6twa_bindingsite_{stepsize}", matrix_6twa_bindingsite.results.dist_matrix)

    
