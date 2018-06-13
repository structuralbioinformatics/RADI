import sys
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['svg.fonttype'] = 'none'

prefix = sys.argv[1]
try:
    todo = sys.argv[2]
except IndexError:
    todo = 'show'
cc = pd.read_csv('{}cmap.dat'.format(prefix),    sep='\t')
mi = pd.read_csv('{}MI_top.dat'.format(prefix),  sep='\t')
di = pd.read_csv('{}DI_top.dat'.format(prefix),  sep='\t')
oc = pd.read_csv('{}cmapOff.dat'.format(prefix), sep='\t')

rc = max(max(cc.max()['#pos1'], cc.max()['pos2']),
         max(mi.max()['#pos1'], mi.max()['pos2']),
         max(di.max()['#pos1'], di.max()['pos2']),
         max(oc.max()['#pos1'], oc.max()['pos2']))

fig = plt.figure(figsize=(12, 12))
ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)
ax.scatter(cc['#pos1'], cc['pos2'],
           label='close contacts', color='lightgrey')
ax.scatter(oc['#pos1'], oc['pos2'],
           label='off-side contacts', color='lightblue')
ax.scatter(mi['#pos1'], mi['pos2'],
           label='MI', color='red', marker='X')
ax.scatter(di['#pos1'], di['pos2'],
           label='DI', color='blue', marker='*')
ax.set_ylim(1, rc)
ax.set_xlim(1, rc)
ax.set_xlabel('sequence')
ax.set_ylabel('sequence')
ax.set_title('contact map')
ax.legend(ncol=2)
ax.grid(color='black', linestyle=':', linewidth=2)
ax.set_axisbelow(True)
plt.tight_layout()

if todo == 'png' or todo == 'svg':
    plt.savefig(".".join([prefix, todo]))
else:
    plt.show()
