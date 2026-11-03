#!/usr/bin/env python3
"""Illustrate the proved coordinate maps; matplotlib is only a rendering dependency."""
from pathlib import Path
from math import isqrt,cos,sin,pi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def prime(n): return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
base=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,3,figsize=(16,6),gridspec_kw={'width_ratios':[1.35,1,1]})
fig.patch.set_facecolor('#f8fafc')
for ax in axes:
 ax.set_facecolor('#f8fafc')
 ax.spines[['top','right']].set_visible(False)
# Prime-pair coordinate chart.
ax=axes[0]; xs=[]; ys=[]; tx=[]; ty=[]; rx=[]; ry=[]
for h in range(3,161):
 for d in range(h-2):
  if prime(h-d) and prime(h+d) and (h-d)%2 and (h+d)%2:
   xs.append(d);ys.append(h)
   if d==1:tx.append(d);ty.append(h)
   if h==145:rx.append(d);ry.append(h)
ax.scatter(xs,ys,s=7,color='#94a3b8',alpha=.7)
ax.scatter(tx,ty,s=17,color='#dc2626',zorder=4,label='Twin primes: d = 1')
ax.axhline(145,color='#0369a1',lw=1,alpha=.6)
ax.scatter(rx,ry,s=30,color='#0369a1',zorder=5,label='Sum 290: 10 prime pairs')
ax.plot([0,157],[3,160],color='#cbd5e1',lw=1)
ax.set(xlim=(-4,161),ylim=(0,164),xlabel='Half-difference d = (q − p)/2',ylabel='Half-sum h = (p + q)/2')
ax.set_title('Prime-pair triangle',loc='left',fontsize=15,pad=15,fontweight='bold')
ax.legend(loc='lower right',frameon=True,fontsize=9)
ax.text(.02,1.01,'p = h − d and q = h + d must both be prime',transform=ax.transAxes,fontsize=9,color='#475569')
# Divisor trial coordinate chart. Same odd pair, different predicate.
ax=axes[1]; xs=[];ys=[];dx=[];dy=[]
for p in range(3,26,2):
 for q in range(3,26,2):
  if p+q<=28:
   xs.append(p);ys.append(q)
   if (q-1)%(p-1)==0:dx.append(p);dy.append(q)
ax.scatter(xs,ys,s=17,color='#cbd5e1')
ax.scatter(dx,dy,s=33,color='#15803d',label='Trial divisor succeeds')
route=[(2*d+1,2*(12-d)+1) for d in range(11,5,-1)]
for (x,y),(xx,yy) in zip(route,route[1:]):
 ax.annotate('',(xx,yy),(x,y),arrowprops={'arrowstyle':'->','color':'#7c3aed','lw':1.6,'shrinkA':4,'shrinkB':4})
ax.scatter(*zip(*route),s=40,facecolors='none',edgecolors='#7c3aed',zorder=5)
ax.annotate('n = 12: trials 11 → 6\nfirst divisor at (13, 13)',xy=(13,13),xytext=(10,22),fontsize=9,color='#6d28d9',arrowprops={'arrowstyle':'-','color':'#7c3aed'})
ax.set(xlim=(0,27),ylim=(0,27),xlabel='Ordered odd port p = 2d + 1',ylabel='Ordered odd port q = 2(n − d) + 1',xticks=[3,7,11,15,19,23],yticks=[3,7,11,15,19,23])
ax.set_title('Conway divisor triangle',loc='left',fontsize=15,pad=15,fontweight='bold')
ax.text(.02,1.01,'Branch guard: (p − 1) divides (q − 1)',transform=ax.transAxes,fontsize=9,color='#475569')
ax.legend(loc='upper right',bbox_to_anchor=(1,-.14),frameon=False,fontsize=9)
# One signed cycle on the units.
ax=axes[2]; ax.axis('off');ax.set_aspect('equal')
coords=[]
for j in range(28):
 angle=pi/2-2*pi*j/28;coords.append((cos(angle),sin(angle)))
for j in range(28):
 x,y=coords[j];xx,yy=coords[(j+1)%28]
 ax.annotate('',(xx,yy),(x,y),arrowprops={'arrowstyle':'->','color':'#94a3b8','lw':.8,'shrinkA':4,'shrinkB':4})
 r=pow(2,j,145);col='#dc2626' if j in(0,14) else '#334155'
 ax.scatter([x],[y],s=28 if j in(0,14) else 9,color=col)
 ax.text(1.17*x,1.17*y,str(r),ha='center',va='center',fontsize=9,color=col,fontweight='bold' if j in(0,14) else 'normal')
ax.text(0,.16,'2¹⁴ = −1',ha='center',fontsize=17,color='#dc2626')
ax.text(0,-.12,'2²⁸ = 1',ha='center',fontsize=17,color='#334155')
ax.text(0,-.42,'Folding ± gives period 14',ha='center',fontsize=10,color='#475569')
ax.set(xlim=(-1.48,1.48),ylim=(-1.48,1.48))
ax.set_title('One of four unit cycles',loc='left',fontsize=15,pad=15,fontweight='bold')
ax.text(.02,1.01,'Multiply by 2 modulo 145',transform=ax.transAxes,fontsize=9,color='#475569')
fig.suptitle('Three exact maps — three different guards and clocks',x=.045,ha='left',fontsize=20,fontweight='bold',y=.98,color='#0f172a')
fig.text(.045,.035,'The same row sum 290 means Goldbach center h = 145, but Conway candidate n = 144.  The unit-cycle clock is neither one.',fontsize=11,color='#475569')
fig.subplots_adjust(left=.055,right=.98,top=.80,bottom=.22,wspace=.35)
fig.savefig(base/'triangle-and-wheel.png',dpi=170,facecolor=fig.get_facecolor())
fig.savefig(base/'triangle-and-wheel.svg',facecolor=fig.get_facecolor(),metadata={'Date':None})
