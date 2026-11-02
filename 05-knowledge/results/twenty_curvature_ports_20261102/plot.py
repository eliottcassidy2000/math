"""Scientific comparison from the exact central-pair formula."""
from pathlib import Path
from math import acosh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                     'svg.hashsalt':'twenty-curvature-ports'})
h=list(range(1,21))
fig,ax=plt.subplots(figsize=(9,5.2))
fig.subplots_adjust(left=.10,right=.98,bottom=.19,top=.90)
ax.plot(h,[2*k for k in h],color='#777777',linestyle='--',label='Binary tree distance 2h')
for base,color in ((2,'#d17829'),(3,'#276aa7'),(4,'#387d50')):
    distances=[acosh(1+((base-2)*base**(k-1)+1)**2/2) for k in h]
    ax.plot(h,distances,color=color,linewidth=2.4,label=f'Hyperbolic distance, scale {base}')
ax.set(xlabel='Depth h of opposite boundary leaves',ylabel='Distance',
       title='A positive boundary gap repairs the loss of tree distance',xlim=(1,20),ylim=(0,60))
ax.set_xticks([1,5,10,15,20])
ax.grid(alpha=.18)
ax.legend(frameon=False,loc='upper left')
fig.text(.50,.025,'Each curve uses its stated metric; no Euclidean unit-distance equivalence is asserted.',
         ha='center',fontsize=9,color='#555555')
fig.savefig(HERE/'gap_distance.png',dpi=170,metadata={'Software':'Matplotlib'})
fig.savefig(HERE/'gap_distance.svg',metadata={'Date':None,'Creator':'plot.py'})
print('Wrote gap_distance.png and gap_distance.svg')
