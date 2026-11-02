"""Scientific figure generated from the exact receipt; no layout claims about edge lengths."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'verification.json').read_text(encoding='utf-8'))
fig,axes=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1.3,1]})
ax=axes[0];r=list(range(27))
ax.plot(r,d['inverse_ball_counts'],color='#153e75',lw=2.5,label='Integers reaching 5')
ax.plot(r,[1+3*x*(x+1)//2 for x in r],color='#a84b00',lw=2,label='Honeycomb capacity')
ax.plot(r,[1+3*x*(x+1) for x in r],color='#387b50',lw=2,label='Triangular capacity')
ax.scatter([21,26],[704,2263],color='#153e75',zorder=4)
ax.annotate('Radius 21: 704 > 694',(21,704),(9,1150),arrowprops={'arrowstyle':'->','color':'#555'},fontsize=10)
ax.annotate('Radius 26: 2263 > 2107',(26,2263),(8,2300),arrowprops={'arrowstyle':'->','color':'#555'},fontsize=10)
ax.set(xlabel='Ordinary Collatz steps / lattice radius',ylabel='Distinct vertices',xlim=(0,27),ylim=(0,2550))
ax.set_title('A finite obstruction inside a convergent basin',fontsize=12,pad=12)
ax.legend(loc='upper left',frameon=False,fontsize=9);ax.grid(alpha=.18)
ax.spines[['top','right']].set_visible(False)
ax=axes[1]
pos={13:(-.85,1),53:(1.3,3),160:(1.3,2),80:(.7,1),40:(0,0),20:(0,-1),10:(0,-2),5:(0,-3)}
edges=[(13,40),(53,160),(160,80),(80,40),(40,20),(20,10),(10,5)]
for u,v in edges:
    ax.annotate('',xy=pos[v],xytext=pos[u],arrowprops={'arrowstyle':'-|>','color':'#777','shrinkA':16,'shrinkB':16,'lw':1.7})
for n,(x,y) in pos.items():
    col='#153e75' if n in(13,53) else ('#a84b00' if n==40 else '#387b50' if n==5 else '#eee')
    ax.scatter([x],[y],s=700,color=col,zorder=3,edgecolors='white')
    ax.text(x,y,str(n),ha='center',va='center',color='white' if col!='#eee' else '#222',fontsize=11,zorder=4)
ax.text(.5,0,'First native meeting\nclocks 1 and 3',va='center',fontsize=10)
ax.text(.5,-3,'Odd ancestor',va='center',fontsize=10)
ax.text(-1.25,-3.9,'Arrow diagram; lengths are not lattice distances.',fontsize=9,color='#555')
ax.set(xlim=(-1.5,2.65),ylim=(-4.2,3.7));ax.axis('off')
ax.set_title('Odd-tree edges share their even tails',fontsize=12,pad=12)
fig.tight_layout(pad=2)
fig.savefig(HERE/'growth_and_ports.png',dpi=180)
fig.savefig(HERE/'growth_and_ports.svg')
