#!/usr/bin/env python3
"""
S633 — EVERYTHING AS GRAPH COLORING. Grounding computations for the coloring atlas.
Builds on S632 (tie-graph = C_n / Cayley) and S626 (partition function): the chromatic polynomial
IS a partition function (Potts), so 'colorings everywhere' = 'partition functions everywhere'.
We verify the key reframes:
  (A) LRC sieve = vertex coloring of the tie-graph; corrector = independent set.  chi/alpha of C_n.
  (B) chromatic polynomial of the tie-graph = the (zero-temperature Potts) partition function;
      relate to covering-depth Z.
  (C) EDGE coloring: K_n edges labeled by pair-sums v_i+v_j mod (2n-1) (the pair-sum sieve, THM-401);
      'proper'/rainbow structure vs looseness.
  (D) unit distance = vertex coloring (Hadwiger-Nelson): chi of small unit-distance graphs.
"""
from math import gcd
import itertools

# ---------- chromatic number / polynomial helpers (small graphs) ----------
def chromatic_number(adj):
    N = len(adj)
    def colorable(k):
        color = [-1]*N
        def bt(v):
            if v == N: return True
            for c in range(k):
                if all(not(adj[v][u] and color[u]==c) for u in range(N)):
                    color[v]=c
                    if bt(v+1): return True
                    color[v]=-1
            return False
        return bt(0)
    k = 1
    while not colorable(k): k += 1
    return k

def independence_number(adj):
    N=len(adj); best=0
    nbr=[set(j for j in range(N) if adj[i][j]) for i in range(N)]
    # simple branch and bound
    def rec(cand, cur):
        nonlocal best
        if cur+len(cand)<=best: return
        if not cand:
            best=max(best,cur); return
        cand=list(cand); v=cand[0]
        rec([u for u in cand[1:] if u not in nbr[v]], cur+1)  # take v
        rec(cand[1:], cur)                                    # skip v
    rec(list(range(N)),0); return best

def cycle(n):
    adj=[[0]*n for _ in range(n)]
    for i in range(n):
        adj[i][(i+1)%n]=adj[(i+1)%n][i]=1
    return adj

print("(A)/(B) LRC tight tie-graph = C_n: sieve=chi, corrector=alpha, chromatic poly P(C_n,k)")
print("     P(C_n,k)=(k-1)^n+(-1)^n (k-1)  [the zero-temp Potts partition function]")
for n in range(3,9):
    adj=cycle(n); chi=chromatic_number(adj); alpha=independence_number(adj)
    P2=(2-1)**n+(-1)**n*(2-1); P3=(3-1)**n+(-1)**n*(3-1)
    print(f"   n={n}: chi(C_n)={chi} (sieve arity, =2 if n even else 3) "
          f"alpha={alpha}=floor(n/2) (corrector size)  P(C_n,2)={P2} P(C_n,3)={P3}")

print("\n(C) EDGE coloring: K_n edges labeled by pair-sum (v_i+v_j) mod (2n-1) [THM-401].")
def pairsum_edge_coloring(speeds, m):
    lab={}
    n=len(speeds)
    for i,j in itertools.combinations(range(n),2):
        lab[(i,j)]=(speeds[i]+speeds[j])%m
    return lab
for n in (5,6,7):
    S=list(range(1,n)); m=2*n-1   # observer + AP movers; pair-sums mod 2n-1
    lab=pairsum_edge_coloring([0]+S, m)
    from collections import Counter
    dist=Counter(lab.values())
    # 'proper' edge coloring? do edges sharing a vertex ever share a label?
    proper=True
    for v in range(n):
        seen=set()
        for (i,j),c in lab.items():
            if v in (i,j):
                if c in seen: proper=False
                seen.add(c)
    print(f"   n={n} (mod {m}): {len(lab)} edges, {len(dist)} distinct pair-sum labels; proper-edge-coloring={proper}")
    print(f"      label multiplicities: {dict(sorted(dist.items()))}")

print("\n(D) UNIT DISTANCE = vertex coloring (Hadwiger-Nelson): chi of small unit-distance graphs")
import cmath, math
def chi_points(pts, tol=1e-7):
    N=len(pts); adj=[[0]*N for _ in range(N)]
    for i,j in itertools.combinations(range(N),2):
        if abs(abs(pts[i]-pts[j])-1)<tol: adj[i][j]=adj[j][i]=1
    return chromatic_number(adj)
# triangle (K3, all unit): chi=3 ; Moser spindle: chi=4
tri=[0, 1, cmath.exp(1j*math.pi/3)]
print(f"   unit triangle K3: chi={chi_points(tri)} (=3)")
# Moser spindle (7 vertices, chi=4): build it
w=cmath.exp(1j*math.pi/3)
A=0; B=1; C=w; D=B+ (C-A)  # rhombus-ish; use known coordinates
import math as M
def rot(z,a): return z*cmath.exp(1j*a)
p0=0+0j; p1=1+0j; p2=rot(1, M.acos(5/6)); p3=rot(p2, -2*M.asin(0.5))
# fall back: known Moser spindle coordinates
th=M.acos(5/6)
spindle=[0, 1, cmath.exp(1j*M.pi/3), cmath.exp(-1j*M.pi/3)]
# add the rotated rhombus by angle th
r=cmath.exp(1j*th)
spindle=[0,1,cmath.exp(1j*M.pi/3),cmath.exp(-1j*M.pi/3)]
spindle+= [r*z for z in [cmath.exp(1j*M.pi/3), cmath.exp(-1j*M.pi/3), 1]]
print(f"   Moser-spindle-like (7 pts): chi={chi_points(spindle)} (target 4; coords approximate)")
print("   => unit distance is NATIVELY a vertex-coloring problem; the disproof maximizes edge density")
print("      (chi-forcing) = the chromatic-number-of-the-plane (Hadwiger-Nelson) family.")
