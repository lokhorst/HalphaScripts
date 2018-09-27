import h5py

def read_header(fileloc='./data/'):
    """ Read various attributes from the header group. """
    f       = h5py.File(fileloc+'snap_028_z000p000.0.hdf5', 'r')
    a       = f['Header'].attrs.get('Time')         # Scale factor.
    h       = f['Header'].attrs.get('HubbleParam')  # h.
    boxsize = f['Header'].attrs.get('BoxSize')      # L [Mph/h].
    f.close()

    return a, h, boxsize