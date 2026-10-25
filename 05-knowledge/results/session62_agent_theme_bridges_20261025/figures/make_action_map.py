#!/usr/bin/env python3
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
p=argparse.ArgumentParser();p.add_argument('--png',action='store_true');args=p.parse_args()
OUT=Path(__file__).resolve().parent
plt.rcParams.update({'svg.fonttype':'none','svg.hashsalt':'session62-fixed-space-map'})
fig=plt.figure(figsize=(12.5,6),facecolor='#f5f3ec');ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,12.5),ylim=(0,6));ax.axis('off')
ax.text(.45,5.5,'ONE ACTION, TWO DIFFERENT FIXED-SPACE QUESTIONS',fontsize=18,weight='bold',color='#183e45')
ax.text(.45,5.08,'Keep the whole permutation action. A minimum over elements is not a common invariant space.',fontsize=10.6,color='#52686a')
def box(x,y,w,h,title,body,color='#187f86'):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.12,rounding_size=.12',facecolor='white',edgecolor=color,linewidth=1.7))
 ax.text(x+w/2,y+h-(.23 if h<1 else .32),title,ha='center',va='center',fontsize=12,weight='bold',color=color)
 ax.text(x+w/2,y+(.22 if h<1 else .46),body,ha='center',va='center',fontsize=11,color='#304c50',linespacing=1.55)
def arrow(a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=15,linewidth=1.6,color='#557275'))
box(.65,2.85,3.0,1.35,'Permutation module','Odd-order A acts on n gaps')
box(4.55,3.45,3.15,1.05,'One element at a time','d = min dim Fix(g)')
box(4.55,1.8,3.15,1.05,'All elements together','r = dim Fix(A)')
box(8.55,3.45,3.15,1.05,'Clone cancellation','vanishing order at −1 = d')
box(8.55,1.8,3.15,1.05,'Connected, holonomy A','dim H¹ = n(β − 1) + r')
arrow((3.8,3.65),(4.4,3.95));arrow((3.8,3.2),(4.4,2.3));arrow((7.85,3.95),(8.4,3.95));arrow((7.85,2.3),(8.4,2.3))
box(.65,.42,11.05,.85,'Regular C₃² on nine gaps; connected gluing over F₂','r = 1      d = 3      β = 2      dim H¹(V) = 10      dim H¹(ker parity) = 8',color='#a55e32')
ax.text(.7,2.1,'The gluing graph is additional data.\nNeither dimension is a plane-coloring bound.',fontsize=10,color='#52686a',linespacing=1.6)
fig.savefig(OUT/'common-action-map.svg',bbox_inches='tight',pad_inches=.08,metadata={'Date':None})
if args.png:fig.savefig(OUT/'common-action-map.png',dpi=140,bbox_inches='tight',pad_inches=.08)
