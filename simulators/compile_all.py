import os
import subprocess

fnames = ['lru', 'drrip', 'popt-8b', 'opt-ideal']
for fname in fnames:

    origDir = os.getcwd()
    os.chdir(origDir + '/' + fname)
    subprocess.call('make', shell=True)
    os.chdir(origDir)
    print('**************************************')
    print('~~~~~~ Compiled ' + fname + ' ~~~~~~~~')
    print('**************************************')

subprocess.call('unset CXX', shell=True)
