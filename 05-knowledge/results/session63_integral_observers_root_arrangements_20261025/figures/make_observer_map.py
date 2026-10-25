#!/usr/bin/env python3
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
p=argparse.ArgumentParser();p.add_argument('--png',action='store_true');args=p.parse_args();OUT=Path(__file__).resolve().parent
plt.rcParams.update({'svg.fonttype':'none','svg.hashsalt':'session63-integral-observers'})
fig=plt.figure(figsize=(13,6),facecolor='#f6f3eb');ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,13),ylim=(0,6));ax.axis('off')
ax.text(.45,5.52,'KEEP THE CONSTRUCTION; DISTINGUISH ITS OBSERVERS',fontsize=18,weight='bold',color='#173e46')
ax.text(.45,5.08,'Proposed framework: integral carrier + operation bank + realizability cone',fontsize=12,color='#586c6d')
items=[('Invariant lifting','Odd group; exterior powers\nof a permutation lattice','Commutes with every\ncoefficient ring','No new invariant classes', '#187d84'),('Multiplication','C₃ exterior invariants','x ∧ y = 3z','Product degenerates mod 3','#a36632'),('Scalar jets','25 gaps; 50 vertices','Depths: 5 over Q\n21 mod 5; 25 mod 2','One invariant-gap toggle','#5968a5'),('Effectivity','Hourglass open base (good char.)','T² − 4T + 7\nT = q − 1','Not a split-torus sum','#9c4e50')]
for i,(title,subtitle,formula,guard,col) in enumerate(items):
 x=.48+i*3.16;w=2.72;y=2.05;hh=2.45
 ax.add_patch(FancyBboxPatch((x,y),w,hh,boxstyle='round,pad=.1,rounding_size=.12',facecolor='white',edgecolor=col,linewidth=1.7))
 ax.text(x+w/2,y+hh-.35,title,ha='center',fontsize=13,weight='bold',color=col)
 ax.text(x+w/2,y+hh-.78,subtitle,ha='center',fontsize=9.7,color='#52666a')
 ax.text(x+w/2,y+1.02,formula,ha='center',va='center',fontsize=12,color='#243f45',linespacing=1.6)
 ax.text(x+w/2,y+.3,guard,ha='center',fontsize=9.5,color=col)
ax.add_patch(FancyBboxPatch((.6,.45),11.8,.95,boxstyle='round,pad=.1,rounding_size=.1',facecolor='#e7eee9',edgecolor='#54736a'))
ax.text(6.5,1.08,'Open root cover = multiplicative power map pulled back along the additive branch map',ha='center',fontsize=12,weight='bold',color='#173e46')
ax.text(6.5,.69,'Fresh twist means, fixed-family extensions and ordinary motives are different operations.',ha='center',fontsize=10.5,color='#52666a')
fig.savefig(OUT/'integral-observer-map.svg',bbox_inches='tight',pad_inches=.08,metadata={'Date':None})
if args.png:fig.savefig(OUT/'integral-observer-map.png',dpi=140,bbox_inches='tight',pad_inches=.08)
