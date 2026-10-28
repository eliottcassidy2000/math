"""Draw the exact arithmetic geometry recorded in normal.json."""
from pathlib import Path
import json
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

p = Path(__file__).parent
d = json.loads((p/'normal.json').read_text(encoding='utf-8'))
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':12})
fig = plt.figure(figsize=(15,10), facecolor='#f4f6f9')
gs = fig.add_gridspec(2,2, height_ratios=[1.15,1], hspace=.32, wspace=.2)
blue, teal, purple = '#214a80', '#007d80', '#6b46a2'

ax = fig.add_subplot(gs[0,0]); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('7 primes: the arithmetic Fano plane', loc='left', weight='bold', color=blue, pad=18)
pos = {1:(0,1), 2:(-math.sqrt(3)/2,-.5), 4:(math.sqrt(3)/2,-.5),
       3:(-math.sqrt(3)/4,.25), 5:(math.sqrt(3)/4,.25), 6:(0,-.5), 7:(0,0)}
for a,b in [(1,2),(1,4),(2,4),(1,6),(2,5),(4,3)]:
    ax.plot([pos[a][0],pos[b][0]],[pos[a][1],pos[b][1]],color='#a5b7cc',lw=2,zorder=1)
ax.add_patch(Circle((0,0),.5, fill=False,edgecolor=teal,lw=2))
for v,(x,y) in pos.items():
    ax.add_patch(Circle((x,y),.105, facecolor='white',edgecolor=blue,lw=1.7,zorder=2))
    ax.text(x,y,str(d['prime_bank7'][v-1]['prime']),ha='center',va='center',weight='bold',color=blue,zorder=3)
ax.set_xlim(-1.2,1.2); ax.set_ylim(-.8,1.3)
ax.text(0,-.75,'7 three-point lines; 14 circuits in total',ha='center',color=blue)

ax = fig.add_subplot(gs[0,1]); ax.axis('off')
ax.set_title('15 primes: one pentad of 5 square-total triples',loc='left',weight='bold',color=teal,pad=18)
for i,L in enumerate(d['geometry']['six_pentads'][0]):
    y=.92-i*.165
    ax.add_patch(FancyBboxPatch((.04,y-.10),.92,.125,boxstyle='round,pad=.012',
                             facecolor='white',edgecolor='#b8d6d5',transform=ax.transAxes))
    ax.text(.5,y-.04,'{ ' + ',  '.join(str(n) for n in L) + ' }',ha='center',va='center',fontsize=18,
            color=teal,transform=ax.transAxes)
ax.text(.5,.02,'Each row is a three-prime circuit.\nFor each row: sigma(product of its primes) is a square.',
        ha='center',va='bottom',fontsize=11,transform=ax.transAxes,color='#405568')

ax = fig.add_subplot(gs[1,0]); ax.axis('off')
ax.set_title('Relations at one scale become points at the next',loc='left',weight='bold',color=blue,pad=18)
entries=[('7 original primes','16 square-total products'),
         ('Remove the identity product 1','15 nonzero codewords'),
         ('Choose a code basis + new prime labels','15 points; 35 three-term relations')]
for i,(label,value) in enumerate(entries):
    y=.94-i*.285
    ax.text(.04,y,label,color='#405568',transform=ax.transAxes,va='top')
    ax.text(.04,y-.085,value,color=blue,fontsize=18,weight='bold',transform=ax.transAxes,va='top')
ax.text(.04,.045,'Group law:  x * y = xy / gcd(x,y)^2',fontsize=14,color=teal,transform=ax.transAxes)

ax = fig.add_subplot(gs[1,1]); ax.axis('off')
ax.set_title('Parity hides two different integer mechanisms',loc='left',weight='bold',color=purple,pad=18)
ax.text(.03,.93,'{2, 5, 7}',fontsize=20,weight='bold',color=purple,transform=ax.transAxes)
ax.text(.03,.81,'Integer relation; Smith factors (1, 1)',transform=ax.transAxes)
ax.text(.03,.63,'{7, 19, 89}',fontsize=20,weight='bold',color=purple,transform=ax.transAxes)
ax.text(.03,.51,'Torsion relation; Smith factors (1, 1, 6)',transform=ax.transAxes)
ax.text(.03,.27,'Fano symmetries: 168  →  2',fontsize=20,weight='bold',color=purple,transform=ax.transAxes)
ax.text(.03,.14,'after retaining every subset’s integer module type',fontsize=11,transform=ax.transAxes)

fig.suptitle('Arithmetic dependencies beyond ordinary graph cycles',fontsize=23,weight='bold',color='#183246',y=.98)
fig.text(.5,.025,'Explicit finite witnesses and scoped proofs • relation_gluing_7_15_20261028 • markings are part of the maps',
         ha='center',fontsize=11,color='#506879')
fig.subplots_adjust(left=.06,right=.96,top=.90,bottom=.09)
fig.savefig(p/'overview.png',dpi=150,facecolor=fig.get_facecolor())
plt.close(fig)
