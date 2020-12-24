import os, subprocess

## The dbpedia dataset used to be available on konect.cc/networks but that website is currently unresponsive 
## Last snapshot (May 6, 2020) -- https://web.archive.org/web/20200506193607/http://konect.cc/networks/ 
## Apparently, Konect datasets stopped being freely available from Oct 2020 -- https://twitter.com/kunegis/status/1323231324036124672?s=20

links  = ['https://suitesparse-collection-website.herokuapp.com/MM/LAW/uk-2002.tar.gz', \
          'https://suitesparse-collection-website.herokuapp.com/MM/DIMACS10/hugebubbles-00020.tar.gz'] 
graphs = ['uk-2002', 'hugebubbles-00020'] 


## Download and extract graph files
subprocess.call('mkdir -p ../input-graphs/', shell=True)
origDir = os.getcwd()
os.chdir(origDir + '/../input-graphs') 

for link in links:
    subprocess.call('wget ' + link, shell=True) 

print('*************************')
print('[DOWNLOADED RAW FILES]')
print('*************************')

for graph in graphs:
    subprocess.call('tar -xvzf ' + graph + '.tar.gz', shell=True)

print('*************************')
print('[EXTRACTED EDGELISTS]')
print('*************************')


## Randomize vertex ordering to isolate any locality artifacts from reordering
origDir = os.getcwd()
os.chdir(origDir + '/../applications/baseline')
subprocess.call('make', shell=True)
os.chdir(origDir)

for graph in graphs:
    subprocess.call('../applications/baseline/randomizer -f ' + graph + '/' + graph + '.mtx -b ' + graph + '.sg', shell=True)

print('*************************')
print('[BUILT REAL GRAPHS]')
print('*************************')

## Cleanup
for graph in graphs:
    subprocess.call('rm ' + graph + '.tar.gz', shell=True)
    subprocess.call('rm -rf ' + graph, shell=True)

## Now build synthetic graphs
graphs = ['kron25-d4', 'urand25-d4']
for graph in graphs:
    if 'kron25-d4' in graph:
        subprocess.call('../applications/baseline/randomizer -g 25 -k 2 -b ' + graph + '.sg', shell=True)
    elif 'urand25-d4' in graph:
        subprocess.call('../applications/baseline/randomizer -u 25 -k 2 -b ' + graph + '.sg', shell=True)

print('*************************')
print('[BUILT SYNTHETIC GRAPHS]')
print('*************************')

