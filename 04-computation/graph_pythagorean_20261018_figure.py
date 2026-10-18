"""Draw the proved five-pole correspondence, without numerical inference."""
from pathlib import Path
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, RED, INK = "#2166ac", "#bd4859", "#263445"
fig, axes = plt.subplots(1, 3, figsize=(14, 5.2))
fig.patch.set_facecolor("white")
fig.suptitle("Petersen's isosceles triangle: equal, disjoint boundary supports",
             fontsize=17, color=INK, y=.97)


def ring(radius):
    return {i: (radius * math.cos(math.pi/2 + 2*math.pi*i/5),
                radius * math.sin(math.pi/2 + 2*math.pi*i/5)) for i in range(5)}


def edge(ax, pos, u, v, color, width=2.1):
    ax.plot([pos[u][0],pos[v][0]], [pos[u][1],pos[v][1]],
            color=color, lw=width, zorder=1)


def nodes(ax, pos, labels=True):
    for i, (x,y) in pos.items():
        ax.scatter([x],[y],s=250,facecolor="white",edgecolor=INK,linewidth=1.4,zorder=3)
        if labels:
            ax.text(x,y,str(i),ha="center",va="center",fontsize=10,color=INK,zorder=4)


for ax in axes:
    ax.set_aspect("equal")
    ax.axis("off")

p = ring(1)
for i in range(5):
    edge(axes[0], p, i, (i+1)%5, BLUE, 2.7)
    edge(axes[0], p, i, (i+2)%5, RED, 2.2)
nodes(axes[0],p)
axes[0].set_title(r"Boundary-pair space: $K_5$",fontsize=13,color=INK,pad=15)
axes[0].text(0,-1.38,"Two complementary five-cycles\n5 blue pairs + 5 red pairs",ha="center",fontsize=11,color=INK)
axes[0].set_xlim(-1.4,1.4)
axes[0].set_ylim(-1.65,1.2)

outer, inner = ring(1.1), ring(.5)
for i in range(5):
    edge(axes[1],outer,i,(i+1)%5,BLUE,2.7)
    edge(axes[1],inner,i,(i+2)%5,RED,2.2)
    axes[1].plot([outer[i][0],inner[i][0]],[outer[i][1],inner[i][1]],
                 color="#89929e",lw=1.8,ls="--",zorder=0)
nodes(axes[1],outer)
nodes(axes[1],inner,False)
axes[1].set_title("Petersen: cut the five spokes",fontsize=13,color=INK,pad=15)
axes[1].text(0,-1.38,"Each boundary pair has 6 colour states\n30 states on each side; no common state",ha="center",fontsize=11,color=INK)
axes[1].set_xlim(-1.4,1.4)
axes[1].set_ylim(-1.65,1.2)

ax=axes[2]
ax.plot([0,1],[0,0],lw=3,color=BLUE)
ax.plot([0,0],[0,1],lw=3,color=RED)
ax.plot([0,1],[1,0],lw=2,color=INK)
ax.plot([.09,.09,0],[0,.09,.09],lw=1,color=INK)
ax.text(.5,-.13,r"$f:\ \sqrt{30}$",ha="center",fontsize=13,color=BLUE)
ax.text(-.14,.5,r"$g:\ \sqrt{30}$",ha="center",va="center",rotation=90,fontsize=13,color=RED)
ax.text(.62,.63,r"$\sqrt{60}$",ha="center",fontsize=13,color=INK)
ax.text(.5,-.40,r"$\langle f,g\rangle=0$"+"\n"+r"$\|f\|^2=\|g\|^2=30$",ha="center",fontsize=13,color=INK)
ax.set_title("The boundary vectors are perpendicular",fontsize=13,color=INK,pad=15)
ax.set_xlim(-.35,1.3)
ax.set_ylim(-.55,1.15)

fig.text(.5,.035,"The vectors count extensions of labelled boundary colours. Identifying the cyclic orders instead gives the colourable pentagonal prism.",
         ha="center",fontsize=10,color=INK)
fig.subplots_adjust(left=.035,right=.975,top=.84,bottom=.17,wspace=.34)
dest=Path(__file__).resolve().parents[1]/"05-knowledge"/"results"/"graph_pythagorean_20261018_boundary.png"
fig.savefig(dest,dpi=160,facecolor="white")
print(dest)
