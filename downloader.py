from pdb import set_trace
import os 
from os.path import exists, join, dirname, basename


def download(remote_path=None, download_dir=None, force=False):

    assert download_dir != None 

    os.system('mkdir -p %s' % download_dir)

    save_path = join(download_dir, basename(remote_path)) 

    if not exists(save_path) or force==True: 
        os.system('wget -O %s %s' % (save_path, remote_path))

    return save_path

    
def unzip(local_filename, target_dir=None):
    if target_dir == None: 
        target_dir = dirname(local_filename)

    os.system( "unzip -q -o %s -d %s" % (local_filename, target_dir))

