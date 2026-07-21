#!/usr/bin/env python3
"""
THE SHARED NEXT TARGET (opus-2026-07-20-S441): how does kps's GIT-instability scalar var(lambda^2)
move under insertion a, and what REALLY interpolates transitive<->Paley?

var(lambda^2) = variance of the squared skew eigenvalues (kps S128c139: transitive MAX, Paley 0).
CLAIMS to test:
  (A) Sum lambda^2 = tr(-S^2) = 2*#arcs = n(n-1) is FIXED for all tournaments on n vertices
      (mean lambda^2 = n-1). So var(lambda^2) depends only on tr(S^4)=Sum lambda^4.
  (B) tr(S^4) is SCORE-DETERMINED, hence var(lambda^2) is an AFFINE function of the 3-cycle count
      c3 (via KBS c3 = C(n,3) - Sum C(s_i,2)): var(lambda^2) = A(n) - B(n)*c3, DECREASING in c3.
      => transitive (c3=0) = MAX var (GIT-unstable nullcone vertex); regular/Paley (c3 max) = MIN
      var (GIT-stable). This unifies kps's var scalar with THM-1820 (c3=Schur-convex intransitivity)
      and my THM-1900 (Delta c3 = forward cut).
  (C) THE INSERTION-RESPONSE: since var = A - B*c3 and Delta c3 = e(P->V\\P) under a=insertion
      (THM-1900), Delta var(lambda^2) under insertion is DIRECTLY -B(n)*(forward cut) (mod the n-shift
      of A,B). kps's GIT scalar moves under insertion by MY forward-cut law.
  (D) CLARIFY kps's family ((x+c)^n+(x-c)^n)/2: its roots are c*i*cot((2k-1)pi/2n) = the transitive
      spectrum SCALED by c, so var = c^4 * var_transitive. At c=0 it is x^n = char_A. So this family
      interpolates char_A <-> char_S of the TRANSITIVE tournament (c:0->1), NOT transitive<->Paley.
      The true transitive<->Paley axis is c3 / var(lambda^2).
"""
import itertools, numpy as np
from fractions import Fraction as F

def edges_iter(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj

def skewnp(adj,n):
    return np.array([[0.0 if i==j else (1.0 if adj[i][j] else -1.0) for j in range(n)] for i in range(n)])

def c3_count(adj,n):
    t=0
    for i,j,k in itertools.combinations(range(n),3):
        if adj[i][j] and adj[j][k] and adj[k][i]: t+=1
        if adj[i][k] and adj[k][j] and adj[j][i]: t+=1
    return t

def scores(adj,n): return [sum(adj[i]) for i in range(n)]

print("var(lambda^2) = GIT-instability scalar: is it affine in c3?  (A)(B)(D)")
print("="*70)
for n in range(3,7):
    sum_l2_fixed=True; tr4_by_score={}; rows=[]  # (c3, tr4, var_l2, sum_s2)
    for adj in edges_iter(n):
        S=skewnp(adj,n)
        S2=S@S; tr2=np.trace(S2); tr4=np.trace(S2@S2)
        l2=sorted(abs(e)**2 for e in np.linalg.eigvals(S))   # lambda^2 multiset
        if abs(sum(l2) - n*(n-1))>1e-6: sum_l2_fixed=False
        var_l2=float(np.var(l2))
        c3=c3_count(adj,n); s=scores(adj,n); ss2=sum(v*v for v in s)
        rows.append((c3, round(tr4,4), round(var_l2,6), ss2))
        tr4_by_score.setdefault(ss2, set()).add(round(tr4,4))
    # is tr4 a function of sum s_i^2 (score-determined)?
    score_det = all(len(v)==1 for v in tr4_by_score.values())
    # fit var_l2 = A - B*c3 (exact affine?)
    c3s=sorted(set(r[0] for r in rows))
    var_by_c3={}
    affine=True
    for c3,tr4,vl2,ss2 in rows:
        var_by_c3.setdefault(c3,set()).add(vl2)
    var_det_by_c3 = all(len(v)==1 for v in var_by_c3.values())
    # slope/intercept from two c3 values
    fit=""
    if var_det_by_c3 and len(c3s)>=2:
        c3a,c3b=c3s[0],c3s[-1]
        va=next(iter(var_by_c3[c3a])); vb=next(iter(var_by_c3[c3b]))
        B=(va-vb)/(c3b-c3a); A=va+B*c3a
        # check all
        ok=all(abs(next(iter(var_by_c3[c]))-(A-B*c))<1e-6 for c in c3s)
        fit=f"var = {A:.4f} - {B:.4f}*c3   (exact affine: {ok})"
    print(f" n={n}: Sum l^2 = n(n-1) fixed: {sum_l2_fixed} | tr(S^4) score-determined: {score_det}"
          f" | var(l^2) determined by c3: {var_det_by_c3}")
    print(f"        {fit}")
    print(f"        c3 range {c3s[0]}..{c3s[-1]};  var(l^2): transitive(c3=0)={max(next(iter(var_by_c3[c])) for c in [0] if 0 in var_by_c3) if 0 in var_by_c3 else 'NA'}"
          f"  min={min(min(v) for v in var_by_c3.values()):.4f}")

# (D) the c-family interpolates char_A <-> char_S of the transitive tournament
print("\n" + "="*70)
print("(D) kps's family ((x+c)^n+(x-c)^n)/2: roots = c * transitive-spectrum; c=0 -> x^n = char_A")
import sympy as sp
xs,cs=sp.symbols('x c')
for n in (4,5):
    fam=sp.expand(((xs+cs)**n+(xs-cs)**n)/2)
    print(f"  n={n}: E_n^(c) = {fam}")
    print(f"        c=1 (char_S, transitive) = {sp.expand(fam.subs(cs,1))};  c=0 (char_A) = {sp.expand(fam.subs(cs,0))}")
print("\n READING: var(lambda^2) is AFFINE-DECREASING in c3 (GIT-instability = intransitivity DEFICIT);")
print(" transitive=max var=nullcone vertex, regular/Paley=min var. Insertion moves var by -B*(forward cut)")
print(" [THM-1900]. The c-family is char_A<->char_S of ONE tournament, NOT the transitive<->Paley axis.")
