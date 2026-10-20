#!/usr/bin/env python3
"""Static research figures drawn directly from the exact n=5 matrices."""
import os
os.environ.setdefault("MPLCONFIGDIR","/tmp/session14-mpl")
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

ROOT=Path(__file__).resolve().parent
D=json.loads((ROOT/"census_n5.json").read_text())
N=D["nodes"]
BLUE="#2968b2"
INK="#23324a"
GOLD="#fff0c9"
PALE="#eaf2fc"


def style(ax):
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-4.4,4.4)
    ax.set_ylim(-3.3,3.5)


def node(ax,u,xy,mode="blue"):
    r=N[u]
    title=f"u{u}"
    if u==0:
        title+=" · transitive"
    if u==7:
        title+=" · regular C₅"
    if mode=="blue":
        text=title+f"\nc₃={r['c3']}   M={r['mass']}   B={r['blue']}"
    else:
        text=title+f" · {r['type']}\nc₃={r['c3']}   M={r['mass']}   B={r['blue']}\ndegree in Q = {r['reflection_pairs']}"
    fc=GOLD if r["type"]=="PB" else PALE if r["type"]=="MX" else "#edf0f3"
    ax.text(*xy,text,ha="center",va="center",fontsize=11 if mode=="blue" else 10,
            color=INK,zorder=3,linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.55",facecolor=fc,
                      edgecolor="#b77a05" if u==7 else INK,lw=2 if u==7 else 1.2))


def blue():
    pos={7:(0,1.15),3:(-1.65,0),5:(0,-1.15),6:(1.65,0),
         9:(0,2.8),1:(-3.35,0),0:(0,-2.8),8:(3.35,0)}
    assert set(pos)=={r["node"]for r in N if r["blue"]}
    deg={u:0 for u in pos}
    fig,ax=plt.subplots(figsize=(13.5,9.4),facecolor="white")
    style(ax)
    for u,v,w in D["blue_edges"]:
        assert w==1 and u!=v
        deg[u]+=w;deg[v]+=w
        ax.plot([pos[u][0],pos[v][0]],[pos[u][1],pos[v][1]],
                color=BLUE,lw=2.8,zorder=1)
    assert all(deg[u]==N[u]["blue"]for u in pos)
    for u,p in pos.items():
        node(ax,u,p)
    fig.suptitle("Pure-blue does not mean pendant",fontsize=23,fontweight="bold",color=INK,y=.97)
    fig.text(.5,.917,"Exact n=5 blue subgraph · 8 vertices · 8 distinct lines · no loops or parallel lines",
             ha="center",fontsize=12,color=INK)
    fig.text(.5,.085,"The regular C₅ vertex has three blue tilings and three different neighbors.",
             ha="center",fontsize=14,fontweight="bold",color=INK)
    fig.text(.5,.052,"Gold: pure-blue    Pale blue: mixed    M: merged tiling count    B: blue degree    c₃: cyclic triangles",
             ha="center",fontsize=11,color=INK)
    fig.text(.5,.025,"All three neighbors of regular C₅ have c₃=4 and score sequence (1,2,2,2,3), yet are nonisomorphic.",
             ha="center",fontsize=11,color=INK)
    fig.subplots_adjust(top=.89,bottom=.12,left=.03,right=.97)
    for suffix in ("svg","png"):
        fig.savefig(ROOT/f"n5_blue_counterexample.{suffix}",dpi=180)
    plt.close(fig)


def paired_black():
    pos={1:(-1.6,1.15),4:(-3.4,2.45),2:(0,0),
         5:(1.6,1.15),9:(3.4,0),6:(-1.6,-1.15),3:(1.6,-1.15)}
    assert set(pos)=={r["node"]for r in N if r["reflection_pairs"]}
    deg={u:0 for u in pos}
    fig,ax=plt.subplots(figsize=(13.5,9.4),facecolor="white")
    style(ax)
    for u,v,w in D["paired_black_edges"]:
        assert w==1
        deg[u]+=w;deg[v]+=w
        if u==v:
            x,y=pos[u]
            ax.add_patch(Arc((x,y-.58),1.15,1.5,theta1=115,theta2=425,
                            color=INK,lw=2.3,zorder=1))
        else:
            bridge=[u,v]in D["graphs"]["paired_black"]["bridges"]
            ax.plot([pos[u][0],pos[v][0]],[pos[u][1],pos[v][1]],
                    color="#bf6b09"if bridge else INK,lw=3.8 if bridge else 2,zorder=1)
    assert all(deg[u]==N[u]["reflection_pairs"]for u in pos)
    for u,p in pos.items():
        node(ax,u,p,"Q")
    fig.suptitle("Remove the forced black-edge doubling",fontsize=23,fontweight="bold",color=INK,y=.97)
    fig.text(.5,.917,"Exact n=5 paired-black quotient Q · 7 active vertices · 12 edges including 2 loops · cycle rank 6",
             ha="center",fontsize=12,color=INK)
    fig.text(.5,.085,"Every drawn edge represents two reflected black lines in the original metagraph.",
             ha="center",fontsize=13,fontweight="bold",color=INK)
    fig.text(.5,.052,"Orange: the unique bridge    Pale blue: mixed    Gray: non-self-converse (NS)    Loops count twice toward degree",
             ha="center",fontsize=10.5,color=INK)
    fig.text(.5,.025,"The three pure-blue vertices are isolated in Q and omitted. Each NS vertex has odd degree in Q.",
             ha="center",fontsize=11,color=INK)
    fig.subplots_adjust(top=.89,bottom=.12,left=.03,right=.97)
    for suffix in ("svg","png"):
        fig.savefig(ROOT/f"n5_paired_black_quotient.{suffix}",dpi=180)
    plt.close(fig)


if __name__=="__main__":
    blue()
    paired_black()
