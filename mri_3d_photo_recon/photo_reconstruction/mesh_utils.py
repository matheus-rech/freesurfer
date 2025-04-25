import os
import numpy as np
from nibabel.freesurfer import read_geometry


# Function that reads, reorients, and centers mesh
def read_and_reorient_mesh(filename, vertices_string, fsprefix, output_directory, hemi):
    print('Trying to read mesh with nibabel')
    try:
        P, T, meta = read_geometry(filename, read_metadata=True)
        print('Success!')
    except:
        print('Nibabel could not read surface; let us convert to freesurfer format first')
        a = os.system(fsprefix + ' mris_convert ' + filename + ' ' + output_directory + '/temp.surf >/dev/null')
        if a > 0:
            raise Exception('mris_convert failed; exiting')
        P, T, meta = read_geometry(output_directory + '/temp.surf', read_metadata=True)
        os.system('rm -rf ' + output_directory + '/temp.surf >/dev/null')

    # A bit dumb if the mesh has a proper header but I just do not trust it...
    meta["valid"] = "1  # volume info valid"
    meta["filename"] = ""
    meta["volume"] = np.array([256, 256, 256]).astype(int)
    meta["voxelsize"] = np.array([1.0, 1.0, 1.0])
    meta["xras"] = np.array([-1.0, 0.0, 0.0])
    meta["yras"] = np.array([0.0, 0.0, -1.0])
    meta["zras"] = np.array([0.0, 1.0, 0.0])
    meta["cras"] = np.array([0.0, 0.0, 0.0])

    # Reorienting mesh with provided vertices")
    idx = np.zeros(3).astype(int)
    aux = vertices_string.split(",")
    for i in range(len(idx)):
        idx[i] = int(aux[i])
    K = P[idx, :]
    K = K - np.mean(K, axis=0)
    Kref = np.array([[0, 85, -20], [0, -80, -25], [0, -5, 45]]).astype(float) # rough RAS aligment, already demeaned!
    H = np.transpose(Kref) @ K
    U, S, Vt = np.linalg.svd(H)
    if np.linalg.det(np.transpose(Vt) @ U) > 0:
        R = np.transpose(Vt) @ np.transpose(U)
    else:
        E = np.eye(3)
        E[2, 2] = -1
        R = np.transpose(Vt) @ (E @ np.transpose(U))
    P -= np.mean(P, axis=0)
    P = P @ R
    # nib.freesurfer.write_geometry(output_directory + '/temp.surf', P, T, volume_info=meta)

    # Let's check if the mesh is likely to be flipped
    xmid = np.mean(P[idx,0])
    score_left = np.sum(P[:,0] < (xmid - 5)) / P.shape[0]
    score_right = np.sum(P[:, 0] > (xmid + 5)) / P.shape[0]
    if ((hemi=='left') and (score_right > score_left)) or ((hemi=='right') and (score_right < score_left)):
        print('  ***IMPORTANT***  Given your selection of hemisphere, the mesh seems to be flipped')
        print('  ***IMPORTANT***  Please review output carefully')
        P[:, 0] = -P[:,0]

    return P, T, meta