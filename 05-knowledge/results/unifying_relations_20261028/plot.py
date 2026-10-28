"""Render the two exact dependence models and the distinct tower layers."""
from pathlib import Path
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle

HERE=Path(__file__).resolve().parent
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11})
fig=plt.figure(figsize=(13,9),facecolor="#f7f9fc")
grid=fig.add_gridspec(2,2,height_ratios=[1.4,1],hspace=.20,wspace=.15)
ink="#18314d"; blue="#276aa5"; amber="#a56716"
fig.suptitle("Arithmetic dependence survives; special properties can change",
             x=.06,ha="left",fontsize=19,color=ink,fontweight="bold",y=.97)

ax=fig.add_subplot(grid[0,0]); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Six primes give an ordinary multigraph",loc="left",color=ink,pad=18)
pts={0:(0,0),2:(2,0),3:(1,1.5),7:(3.5,0)}
def edge(u,v,label,rad=0,offset=(0,0)):
    a,b=pts[u],pts[v]
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle="-",connectionstyle=f"arc3,rad={rad}",
                                linewidth=2,color=blue))
    x,y=(a[0]+b[0])/2+offset[0],(a[1]+b[1])/2+offset[1]
    ax.text(x,y,label,color=blue,ha="center",va="center",
            bbox={"facecolor":"#f7f9fc","edgecolor":"none","pad":2})
edge(0,3,"p=2",.22,(-.35,0))
edge(0,3,"p=11",-.22,(.35,0))
edge(2,3,"p=5",0,(.15,.08))
edge(0,2,"p=7",0,(0,-.18))
edge(2,7,"p=13",0,(0,-.18))
ax.add_patch(Circle((-.30,-.17),.28,fill=False,color=blue,lw=2))
ax.text(-.65,-.55,"p=3",color=blue)
for v,(x,y) in pts.items():
    ax.scatter([x],[y],s=290,c="white",edgecolors=ink,zorder=5,linewidths=1.6)
    ax.text(x,y,str(v),ha="center",va="center",color=ink,zorder=6)
ax.text(-.65,-1.03,"66 = loop + parallel pair\n70 = triangle\n13 lies on a bridge, so no zero-sum subset uses it.",
        color=ink,fontsize=11,linespacing=1.6,va="top")
ax.set_xlim(-.9,4); ax.set_ylim(-1.85,2.05)

ax=fig.add_subplot(grid[0,1]); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Seven primes give the nongraphic Fano matroid",loc="left",color=ink,pad=18)
rt=math.sqrt(3)
q={1:(0,1),2:(-rt/2,-.5),4:(rt/2,-.5),
   3:(-rt/4,.25),5:(rt/4,.25),6:(0,-.5),7:(0,0)}
for a,b in [(1,2),(1,4),(2,4),(1,6),(2,5),(4,3)]:
    ax.plot([q[a][0],q[b][0]],[q[a][1],q[b][1]],color=blue,lw=1.7)
ax.add_patch(Circle((0,0),.5,fill=False,color=blue,lw=1.7))
primes={1:7,2:2,4:19,3:5,5:89,6:59,7:29}
offsets={1:(0,.14),2:(-.10,-.16),4:(.11,-.16),3:(-.23,.06),
         5:(.25,.06),6:(0,-.17),7:(.05,-.15)}
for k,(x,y) in q.items():
    ax.scatter([x],[y],s=80,c=ink,zorder=5)
    dx,dy=offsets[k]
    ax.text(x+dx,y+dy,f"p={primes[k]}",color=ink,ha="center",va="center",fontsize=11,
            bbox={"facecolor":"#f7f9fc","edgecolor":"none","pad":1})
ax.text(0,-.90,"Columns are sf(p+1) = 2, 3, 5, 6, 10, 15, 30.\n16 zero-sum subsets:  1 + 7z³ + 7z⁴ + z⁷\nThe [7,4,3] Hamming code.",ha="center",va="top",color=ink,linespacing=1.6)
ax.set_xlim(-1.5,1.5); ax.set_ylim(-1.65,1.4)

ax=fig.add_subplot(grid[1,:]);ax.axis("off")
ax.set_title("Equal divisor sums persist along both doubling towers",loc="left",color=ink,pad=14,fontsize=14)
rows=[
["0","33, 35","48","(2, 2, 2)","Neither special property"],
["1","66, 70","144 = 12²","(6, 3, 6)","Only square-total layer"],
["2","132, 140","336","(10, 3, 7)","Only F = U + S layer"],
["3","264, 280","720","(14, 3, 7)","Neither special property"]]
tab=ax.table(cellText=rows,colLabels=["a","2ᵃ × {33,35}","Common σ","(F, U, S)","Property"],
             colWidths=[.06,.19,.17,.19,.39],cellLoc="left",loc="upper left",bbox=[0,.14,1,.76])
tab.auto_set_font_size(False);tab.set_fontsize(12)
for (r,c),cell in tab.get_celld().items():
    cell.set_edgecolor("white");cell.set_linewidth(2)
    cell.set_facecolor(ink if r==0 else "#edf2f8")
    cell.set_text_props(color="white" if r==0 else ink,weight="bold" if r==0 else "normal")
    if r==2:cell.set_facecolor("#dcecf9")
    if r==3:cell.set_facecolor("#f6e8ce")
ax.text(0,-.01,"For every a ≥ 0:  σ(2ᵃ·33) = σ(2ᵃ·35) = 48(2ᵃ⁺¹ − 1).  Both exceptional layers are proved unique.",
        color=ink,fontsize=11,transform=ax.transAxes)
fig.subplots_adjust(left=.06,right=.97,top=.88,bottom=.06)
fig.savefig(HERE/"relations.png",dpi=160,facecolor=fig.get_facecolor())
