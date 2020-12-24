import subprocess
import os

fnames = ['baseline', 'popt', 'opt-ideal']

for fname in fnames:
    origDir = os.getcwd()
    os.chdir(origDir + '/' + fname)
    subprocess.call('make clean; make', shell=True)
    os.chdir(origDir)
    print('**************************************')
    print('~~~~~~ Compiled ' + fname + ' ~~~~~~~~')
    print('**************************************')
