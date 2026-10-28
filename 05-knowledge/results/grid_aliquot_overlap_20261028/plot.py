"""Render exact coordinate diagrams; requires matplotlib. No numerical inference."""
from pathlib import Path
from math import sqrt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

ROOT=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                     'axes.spines.right':False,'svg.fonttype':'none'})
blue='#21618C'; orange='#C65D22'; teal='#187B6B'; grey='#89939C'
fig,axs=plt.subplots(2,2,figsize=(12,9.5),layout='constrained')
ax=axs[0,0]
first=[(x,sqrt(3)*y) for x in range(-2,3) for y in range(-2,3)]
second=[(x+.5,sqrt(3)*(y+.5)) for x in range(-3,3) for y in range(-2,2)]
for y in range(-2,3): ax.axhline(sqrt(3)*y,color='#e5e9ed',lw=.7,zorder=0)
for x in range(-2,3): ax.axvline(x,color='#e5e9ed',lw=.7,zorder=0)
ax.scatter(*zip(*first),s=36,color=blue,label=r'$L=\mathbb{Z}[\sqrt{-3}]$')
ax.scatter(*zip(*second),s=36,color=orange,label=r'$h+L$')
ax.add_patch(Circle((0,0),1,fill=False,ls='--',color=grey,lw=1))
ax.annotate(r'$h=(1/2,\sqrt{3}/2)$',(.5,sqrt(3)/2),(1,1.25),arrowprops={'arrowstyle':'-','color':orange},fontsize=10)
ax.set(xlim=(-2.5,2.5),ylim=(-2.05,2.05),xlabel='Real coordinate',ylabel='Imaginary coordinate',
       title='The half-shift completes the Eisenstein lattice')
ax.set_aspect('equal'); ax.legend(loc='lower left',fontsize=9,framealpha=.95)
ax=axs[0,1]
for b in range(5):
    for a in range(5-b):
        n=16*b+a; inside=n<=50
        ax.scatter(b,a,s=280,facecolors=teal if inside else 'white',edgecolors=teal if inside else grey,zorder=3)
        ax.annotate(str(n),(b,a),ha='center',va='center',color='white' if inside else grey,fontsize=10)
ax.plot([0,4],[4,0],color=grey,ls='--',lw=1,zorder=1)
ax.set(xlim=(-.45,4.45),ylim=(-.4,4.45),xticks=range(5),yticks=range(5),
       xlabel='b: summands equal to 16',ylabel='a: summands equal to 1',
       title='Four fourth powers: 14 values through 50')
ax.grid(alpha=.16); ax.text(2.2,3.45,r'$a+b\leq4$'+'\n13 positive values\nplus zero',ha='center',fontsize=11)
for ax,q in zip(axs[1],[3,5]):
    n=4*q
    for j in (0,1): ax.plot([0,2],[j,j],color='#cbd3d9',zorder=0)
    for i in range(3): ax.plot([i,i],[0,1],color='#cbd3d9',zorder=0)
    for i in range(3):
        for j in range(2):
            v=2**i*q**j; top=v==n
            ax.scatter(i,j,s=800,facecolor='white' if top else blue,edgecolor=grey if top else blue,zorder=2)
            ax.text(i,j,str(v),color=grey if top else 'white',ha='center',va='center',fontsize=14)
            if top: ax.plot([i-.18,i+.18],[j-.14,j+.14],color=grey,zorder=3)
    total=sum(2**i*q**j for i in range(3) for j in range(2))-n
    ax.set(xlim=(-.45,2.55),ylim=(-.45,1.55),xticks=range(3),yticks=[0,1],
           xlabel='Exponent of 2',ylabel='Exponent of '+str(q),
           title=f'{n} = 2² × {q}: the same 3 × 2 divisor box')
    ax.text(1,-.31,f's({n}) = {total}; excess = {total-n}',ha='center',fontsize=12,color=teal)
fig.suptitle('Coordinate geometry, summand budgets, and divisor weights',fontsize=16)
fig.savefig(ROOT/'grids-and-divisors.png',dpi=170)
plt.close(fig)
