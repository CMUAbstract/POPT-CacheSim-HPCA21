import os, subprocess

subprocess.call('wget http://software.intel.com/sites/landingpage/pintool/downloads/pin-2.14-71313-gcc.4.4.7-linux.tar.gz', shell=True)
subprocess.call('tar -xvzf pin-2.14-71313-gcc.4.4.7-linux.tar.gz', shell=True)
subprocess.call('rm -rf ../pin-2.14; mv pin-2.14-71313-gcc.4.4.7-linux ../pin-2.14/', shell=True)
subprocess.call('rm pin-2.14-71313-gcc.4.4.7-linux.tar.gz', shell=True)

