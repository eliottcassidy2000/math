#!/usr/bin/env python3
"""Illustration only. Exact geometry is verified in the Penrose/audit checkers."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/session71-mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/session71-cache')
os.makedirs(os.environ['XDG_CACHE_HOME'],exist_ok=True)
from pathlib import Path
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
matplotlib.rcParams['svg.hashsalt']='session71-penrose-guards'
phi=(1+math.sqrt(5))/2;q=2-phi
D=[(0,0),(1,0),(q,q),(0,1)]
K=[(1,1),(1,0),(q,q),(0,1)]
good=[(x,y-1) for x,y in K]
bad=[((1-phi)*x-y+phi,x-1) for x,y in K]
def xy(p):
 x,y=p;return ((x+y)*math.cos(math.pi/5),(x-y)*math.sin(math.pi/5))
def tile(ax,ps,colors,face):
 pts=[xy(p) for p in ps]
 ax.add_patch(Polygon(pts,closed=True,facecolor=face,edgecolor='#23354d',lw=1.8))
 cx=sum(x for x,y in pts)/4;cy=sum(y for x,y in pts)/4
 for (x,y),color in zip(pts,colors):
  # Distinct dots slightly inside each incident tile preserve coincident markings.
  dx,dy=cx-x,cy-y;l=math.hypot(dx,dy)
  ax.plot(x+.065*dx/l,y+.065*dy/l,'o',ms=8,mfc=color,mec='#23354d',mew=1.1,zorder=6)
 return pts
fig,axes=plt.subplots(1,2,figsize=(10,5.2))
for ax,ps,title,color in zip(axes,[good,bad],['Accepted local attachment','Rejected local attachment'],['#19704a','#b33c38']):
 tile(ax,D,['black','white','black','white'],'#e0edf8')
 tile(ax,ps,['white','black','white','black'],'#f7e8bc')
 a,b=xy((0,0)),xy((1,0));ax.plot([a[0],b[0]],[a[1],b[1]],color='#7c4c9b',lw=3,zorder=5)
 ax.text(a[0]-.15,a[1]-.06,'0',fontsize=12,color='#7c4c9b')
 ax.text(b[0]+.05,b[1]+.02,'U',fontsize=12,color='#7c4c9b')
 ax.set_title(title,color=color,fontsize=15,pad=12)
 ax.set_aspect('equal');ax.axis('off');ax.set_xlim(-1.05,1.05);ax.set_ylim(-.72,1.55)
fig.suptitle('Same kite + dart, unit seam, and adjacency graph',fontsize=17,y=.98)
fig.text(.5,.055,'Endpoint marks agree on the left and disagree on the right. Neither placement overlaps.\nFinite local matching test only; no claim about extension to an infinite tiling.',ha='center',fontsize=10,color='#334155')
fig.subplots_adjust(top=.85,bottom=.16,wspace=.18)
out=Path(__file__).resolve().parent
fig.savefig(out/'penrose-guards.svg',metadata={'Date':None,'Creator':'Session71 exact-coordinate illustration'})
fig.savefig(out/'penrose-guards.png',dpi=160,metadata={'Software':'Session71 exact-coordinate illustration'})
print('Wrote penrose-guards.svg and penrose-guards.png; exact predicates are checked separately.')
