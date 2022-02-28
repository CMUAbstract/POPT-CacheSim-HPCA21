import os, subprocess

subprocess.call('wget https://software.intel.com/sites/landingpage/pintool/downloads/pin-3.21-98484-ge7cd811fd-gcc-linux.tar.gz', shell=True)
subprocess.call('tar -xvzf pin-3.21-98484-ge7cd811fd-gcc-linux.tar.gz', shell=True)
subprocess.call('rm -rf ../pin-3.21; mv pin-3.21-98484-ge7cd811fd-gcc-linux ../pin-3.21/', shell=True)
subprocess.call('rm pin-3.21-98484-ge7cd811fd-gcc-linux.tar.gz', shell=True)

