import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from math import sqrt
from pathlib import Path
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--png",action="store_true",help="also render local PNG for visual review")
args=parser.parse_args()
OUT=Path(__file__).resolve().parent
plt.rcParams['svg.fonttype']='none'
plt.rcParams['svg.hashsalt']='session61-hn-map'
fig=plt.figure(figsize=(13,6),facecolor='#f6f4ed')
ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,13),ylim=(0,6));ax.axis('off')
ax.text(.45,5.55,'FROM SIX NEIGHBORS TO FOUR SURFACE CODES',fontsize=19,weight='bold',color='#173b40')
ax.text(.45,5.13,'The bridge requires a quotient, triangular faces and a marking. Counts alone do not suffice.',fontsize=11,color='#526365')
# seven branch-label patch, with physical120degree basis
points=[(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]
for a,b in points:
 x,y=1.75+.75*(a-b/2),3.45+.75*sqrt(3)*b/2
 label=(a+2*b)%7
 ax.scatter(x,y,s=440,c='#137e82' if label else '#d58b35',zorder=3,edgecolors='white',linewidths=2)
 ax.text(x,y,str(label),ha='center',va='center',color='white',fontsize=12,weight='bold',zorder=4)
ax.text(1.75,2.35,'Seven dual branch labels',ha='center',fontsize=11,weight='bold',color='#173b40')
ax.text(1.75,1.98,'center + six Eisenstein units',ha='center',fontsize=10,color='#526365')

def box(x,y,w,h,title,body,col='#137e82'):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.13,rounding_size=.12',edgecolor=col,facecolor='white',linewidth=1.6))
 ax.text(x+w/2,y+h-.35,title,ha='center',va='center',fontsize=12,weight='bold',color=col)
 ax.text(x+w/2,y+h/2-.18,body,ha='center',va='center',fontsize=10,color='#354d50',linespacing=1.6)
def arrow(x1,y1,x2,y2,label=None):
 ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=15,color='#49676a',linewidth=1.5))
 if label:ax.text((x1+x2)/2,(y1+y2)/2+.15,label,fontsize=9,ha='center',color='#49676a')
box(3.6,2.65,2.3,1.65,'Choose an ideal','Two chiral quotients\na + 2b or a + 4b (mod 7)')
box(6.6,2.65,2.35,1.65,'Keep the triangles','Torus: 7 vertices, 14 faces\nTwo Fano face types')
box(9.65,2.65,2.65,1.65,'Transport the marking','4 free Campedelli codes\nOne class; Aut = C₂ × A₄')
arrow(2.65,3.45,3.43,3.45);arrow(6.05,3.45,6.43,3.45);arrow(9.1,3.45,9.48,3.45)
box(.65,.48,5.3,.85,'Source graph: 3 colors','Its K₇ quotient needs 7. A lower bound does not lift backwards.','#b24e38')
box(6.55,.48,5.7,.85,'Classical plane upper bound: 7 colors','Retain Voronoi diameter and intercell separation to restore the metric.','#137e82')
fig.savefig(OUT/'six-seven-hn-map.svg',bbox_inches='tight',pad_inches=.08,metadata={'Date':None})
if args.png:
    fig.savefig(OUT/'six-seven-hn-map.png',dpi=140,bbox_inches='tight',pad_inches=.08)
