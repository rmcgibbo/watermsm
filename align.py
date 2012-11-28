"""
Load a TRR, align all the frames to the first frame, save the results
as lh5
"""

import numpy as np
from msmbuilder import Trajectory
from  lprmsd import LPRMSD

t = Trajectory.load_from_trr('/home/jweber1/CNT_for_Robert/runCNT_300K.trr', PDBFilename='/home/jweber1/CNT_for_Robert/CNT300.pdb')

C_indices = np.where(t['AtomNames'] == 'C')[0]
all_indices = np.arange(t['XYZList'].shape[1])

lp = LPRMSD(C_indices, altindices=all_indices)
pt = lp.prepare_trajectory(t)
distances, aligned = lp.one_to_all_aligned(pt, pt, 0)

t['XYZList']  = aligned
t.save('runCNT_300K.lh5')
