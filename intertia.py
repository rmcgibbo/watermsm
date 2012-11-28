import numpy as np

def align_to_moments(xyz, thresh=1e-3):
    """
    Align a molecule to its moments of inertia
    
    Parameters
    ----------
    xyz : np.ndarray, shape=(n_atoms, 3)
        xyz coordinates of each atom in a single frame
    thresh : float, default=1e-3
        Trigger an error if the eigenvalues are too small
        
    Note
    ----
    This is a refactored version of code from Lee-Ping's LPRMSD
    
    Returns
    -------
    rotated_xyz : np.ndarray, shape=(n_atoms, 3)
        A rotated version of xyz
    rotation_matrix : np.ndarray, shape=(3,3)
        The rotation matrix
    """
    I = np.zeros((3,3))
    for i, xi in enumerate(xyz):
        I += np.dot(xi, xi) * np.eye(3) - np.outer(xi, xi)
    
    A, B = np.linalg.eig(I)
    # Sort eigenvectors by eigenvalue
    BB = B[:, np.argsort(A)]
    determ = np.linalg.det(BB)

    if np.abs(determ - 1.0) > thresh:
        if np.abs(determ + 1.0) > thresh:
            raise RuntimeError("determinant is % .3f", determ)
        BB[:,2] *= -1

    rotated_xyz = np.dot(BB.T, xyz.T).T.copy()

    return rotated_xyz, BB


def center(xyzlist):
    """
    Center a collection of frames INPLACE
    
    Parameters
    ----------
    xyzlist : np.ndarray, shape=(n_frames, n_atoms, 3)
    
    Notes
    -----
    This acts in place on xyzlist
    
    Returns
    -------
    centers : np.ndarray, shape=(n_frames, 3)
        The mean position in each frame that was subtracted from each
        atom
    """
    centers = np.zeros((xyzlist.shape[0], xyzlist.shape[2]))
    
    for i in xrange(xyzlist.shape[0]):
        X = xyzlist[i].astype(np.float64)
        centers[i] = X.mean(0)
        X -= centers[i]
        xyzlist[i] = X

    return centers
    
if __name__ == '__main__':
    from msmbuilder import Trajectory
    t = Trajectory.load_trajectory_file('short_traj.lh5')
    xyz = t['XYZList'][:, :320, :]
    
    centers = center(xyz)
    rotations = np.zeros((len(xyz), 3, 3))
    
    for i in range(len(xyz)):
        frame, B = align_to_moments(xyz[i])
        xyz[i] = frame
        rotations[i] = B
    
    print 
    print rotations