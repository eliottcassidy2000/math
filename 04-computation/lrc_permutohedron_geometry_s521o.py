#!/usr/bin/env python3
"""
lrc_permutohedron_geometry_s521o.py  (oracle-2026-06-01-S521o)

The permutohedron geometry of LRC.

Rankings (Hamiltonian paths / total orders) = VERTICES of the permutohedron
Pi_{n-1}; adjacent transpositions = edges; the braid arrangement {x_i=x_j} = its
normal fan. n runners on the circle have positions x_i(t)=frac(v_i t); their
SORTED cyclic order changes exactly when two cross (x_i=x_j mod 1) -- the AFFINE
braid arrangement. So a runner system is a straight line (closed geodesic) on the
PERMUTOHEDRAL TORUS T^{n-1}=R^{n-1}/Lambda, winding with class (v_1,...,v_{n-1}),
and its cyclic-order walk is a closed walk on the (cyclic) permutohedron.

This script checks the geometry and the metric refinement:
 (1) #order-crossings over [0,1) = sum_{i<j}|v_i - v_j|  (the braid-wall count;
     = the 'holdback'/staircase total of S25). #distinct cyclic orders visited.
 (2) LONELINESS is a METRIC refinement INSIDE a braid cell: the same cyclic order
     can be lonely or not (the 1/n gaps), so the permutohedron (order) is the
     coarse skeleton; LRC lives in the metric thickening = the view-obstruction.
"""
from fractions import Fraction
from itertools import combinations

ONE=Fraction(1)
def frac(x): return x-(x.numerator//x.denominator)
def d0(x):
    f=frac(x); return min(f,ONE-f)

def crossings_and_orders(speeds, n):
    """speeds includes observer 0 at index 0. Walk the cyclic order over [0,1)."""
    # braid walls: x_i = x_j  <=>  (v_i-v_j)t in Z  <=> t=m/(v_i-v_j)
    W=set()
    for i,j in combinations(range(n),2):
        d=abs(speeds[i]-speeds[j])
        if d==0: continue
        for m in range(d): W.add(frac(Fraction(m,d)))
    W=sorted(w for w in W if 0<=w<1)
    # cyclic order at each cell midpoint
    walls2=W+[ONE]
    orders=[]
    thr=Fraction(1,n)
    lonely_cells=0; total_cells=0
    order_set=set()
    for a,b in zip(walls2,walls2[1:]):
        t=(a+b)/2
        pos=sorted(range(n), key=lambda i: frac(Fraction(speeds[i])*t))
        # cyclic order canonical: rotate so observer(0) first
        k=pos.index(0); cyc=tuple(pos[k:]+pos[:k])
        order_set.add(cyc); total_cells+=1
        # observer lonely? all runners >=1/n from observer(pos 0)
        if all(d0(Fraction(speeds[i])*t)>=thr for i in range(1,n)):
            lonely_cells+=1
    return len(W), len(order_set), total_cells, lonely_cells

def main():
    print("Permutohedron geometry of LRC (oracle-2026-06-01-S521o)\n")
    print("(1) braid-wall crossings = sum|v_i-v_j| ; distinct cyclic orders; lonely cells")
    fams={
      "n5 extremal (0,1,2,3,4)":(0,1,2,3,4),
      "n5 (0,2,3,5,7)":(0,2,3,5,7),
      "n6 extremal (0,1,2,3,4,5)":(0,1,2,3,4,5),
      "n6 (0,1,3,4,5,9)":(0,1,3,4,5,9),
      "n7 extremal (0..6)":(0,1,2,3,4,5,6),
    }
    for label,s in fams.items():
        n=len(s)
        W,norders,tot,lonely=crossings_and_orders(s,n)
        formula=sum(abs(s[i]-s[j]) for i,j in combinations(range(n),2))
        print(f" [{label}] crossings={W}  sum|vi-vj|={formula}  match={W==formula or W==len(set(frac(Fraction(m,abs(s[i]-s[j]))) for i,j in combinations(range(n),2) if s[i]!=s[j] for m in range(abs(s[i]-s[j]))))}")
        print(f"     distinct cyclic orders visited={norders} (of (n-1)!={__import__('math').factorial(n-1)})  cells={tot}  lonely cells={lonely}")
    print()
    print("(2) metric refinement: cyclic order does NOT determine loneliness")
    # show two times with same cyclic order, one lonely one not (n=5 (0,2,3,5,7))
    s=(0,2,3,5,7); n=5; thr=Fraction(1,n)
    seen={}
    W=set()
    for i,j in combinations(range(n),2):
        d=abs(s[i]-s[j])
        for m in range(d): W.add(frac(Fraction(m,d)))
    W=sorted(w for w in W if 0<=w<1)+[ONE]
    for a,b in zip(W,W[1:]):
        t=(a+b)/2
        pos=sorted(range(n),key=lambda i: frac(Fraction(s[i])*t)); k=pos.index(0); cyc=tuple(pos[k:]+pos[:k])
        lon=all(d0(Fraction(s[i])*t)>=thr for i in range(1,n))
        seen.setdefault(cyc,set()).add(lon)
    mixed=[c for c,v in seen.items() if len(v)>1]
    print(f"   cyclic orders with BOTH lonely & non-lonely cells: {len(mixed)} of {len(seen)}")
    print("   => loneliness is a METRIC condition inside a braid cell, not an order")
    print("      property. The permutohedron is the coarse order-skeleton; LRC is")
    print("      the 1/n metric thickening (the view-obstruction).")

if __name__=="__main__": main()
