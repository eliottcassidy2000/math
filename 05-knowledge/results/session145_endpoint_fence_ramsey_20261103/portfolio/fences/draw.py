"""Illustration only; exact certificates are in verify.py and PROOF.md."""
import os, tempfile
os.environ.setdefault('MPLCONFIGDIR', os.path.join(tempfile.gettempdir(), 'math-geometry-plot-cache'))
os.environ.setdefault('XDG_CACHE_HOME', os.path.join(tempfile.gettempdir(), 'math-geometry-plot-cache'))
from pathlib import Path
from math import sqrt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
r=sqrt(13); d=sqrt((r-1)/6); a=(r+1)/4; b=3*d/2; h=d*(r+7)/4
xy=[(-1,0),(-a,b),(-.5,h),(.5,h),(a,b),(1,0)]
fig=plt.figure(figsize=(10,6),facecolor='#fcfaf5')
ax=fig.add_axes([.06,.13,.43,.76]); ax.set_aspect('equal'); ax.axis('off')
ax.add_patch(Polygon(xy,facecolor='#d7e8e7',edgecolor='#253d4c',linewidth=2.2))
ax.add_patch(Polygon([(x,-y) for x,y in xy],facecolor='#e9dfef',edgecolor='#253d4c',linewidth=2.2))
ax.scatter([0],[0],s=22,c='#253d4c',zorder=3)
ax.set_xlim(-1.45,1.45); ax.set_ylim(-2,2)
for x in (-.5,.5): ax.text(x,.12,'1',ha='center',fontsize=13,color='#253d4c')
ax.text(0,.85,'area a',ha='center',fontsize=15,color='#253d4c')
ax.text(0,-.9,'area a',ha='center',fontsize=15,color='#253d4c')
fig.text(.06,.94,'At 12 edges, sharing a path beats every cactus',fontsize=19,color='#253d4c')
fig.text(.55,.77,'12 unit edges\nTheta path lengths: 5, 5, 2',fontsize=15,linespacing=1.8,color='#253d4c')
fig.text(.55,.55,'Known product\n(89 + 26√13)/16 ≈ 11.42152082',fontsize=15,linespacing=1.8,color='#253d4c')
fig.text(.55,.32,'Exact cactus maximum\n6 + 3√3 ≈ 11.19615242',fontsize=15,linespacing=1.8,color='#253d4c')
fig.text(.06,.045,"Reconstruction of Friedman's known example. No global optimality claim.",fontsize=11,color='#253d4c')
p=Path(__file__).parent
fig.savefig(p/'theta12.svg',facecolor=fig.get_facecolor())
fig.savefig(p/'theta12.png',dpi=160,facecolor=fig.get_facecolor())
