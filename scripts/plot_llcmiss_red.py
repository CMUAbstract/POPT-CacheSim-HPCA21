import matplotlib.pyplot as plt
import numpy as np
import compiler

app         = 'pr'
appNick     = 'PAGERANK'
policies    = ['lru', 'drrip', 'popt-8b', 'opt-ideal']
policyNicks = ['LRU', 'DRRIP', 'P-OPT', 'T-OPT']
versions    = ['baseline', 'baseline', 'popt', 'opt-ideal']
graphs      = ['uk-2002', 'kron25-d4', 'urand25-d4', 'hugebubbles-00020']
graphNicks  = ['UK-02', 'KRON', 'URND', 'HBBL']

## Collect data
data = {}
for graph in graphs:
    data[graph] = {}
    for p in range(len(policies)):
        policy  = policies[p]
        version = versions[p]
        data[graph][policy] = compiler.getLLCMisses(app, graph, policy, version)

## Normalize results
for graph in graphs:
    base = data[graph][policies[0]]
    for policy in policies:
        data[graph][policy] = base / data[graph][policy]

## Plot results
fig = plt.figure(figsize = (5, 3))
bwidth = 1 / (1.25 + len(policies))
maxVal = 0
for p in range(len(policies)):
    policy = policies[p]
    vals   = []
    for graph in graphs:
        vals.append(data[graph][policy])
    maxVal = max(maxVal, max(vals))
    ind = np.arange(len(vals))
    plt.bar(ind + p * bwidth, vals, width = bwidth, label = policyNicks[p])

ax = plt.gca()
ax.set_title('App - ' + appNick, fontweight = 'bold')
ax.set_xlabel('Input Graphs', fontweight = 'bold')
ax.set_ylabel('LLC Miss Reduction', fontweight = 'bold')
ax.set_xticks(np.arange(len(graphs)) + 1.5 * bwidth)
ax.set_xticklabels(graphNicks)
ax.set_yticks(np.arange(0, 1.1 * maxVal, 0.2))
ax.grid(alpha = 0.64, axis = 'y', linestyle = '--')
ax.legend(loc = 'lower right', framealpha = 1.0, ncol = 2)
ax.axhline(y = 1, linestyle = '--', color = 'k')

plt.tight_layout()
plt.savefig('llcmiss-red.pdf')

print('*********************************')
print('[OUTPUT SAVED TO llcmiss-red.pdf]')
print('*********************************')
