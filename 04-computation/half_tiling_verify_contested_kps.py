#!/usr/bin/env python3
"""
half_tiling_verify_contested_kps.py
Independent re-verification (kind-pasteur 2026-06-20) of the THREE workflow claims
whose adversarial verifier died on an API error, before touching canon:
  (1) c3-parity dichotomy: phi-self-converse tournaments have c3 EVEN at even n,
      both parities at odd n  (HYP-2687 / candidate THM).
  (2) Fix_anti(phi)_full = 2^(half(n)+floor(n/2))  (HYP-2686 corrected exponent).
  (3) SC_n (self-converse iso classes) vs A000568(n-1) -- the claimed REFUTATION of
      'SC(n)=A000568(n-1)' in two-models-staircase-recursion.md.
All from first principles, no reuse of the workflow agents' scripts.
"""
import itertools
from math import comb, floor
from collections import Counter

def tiles(n):
    out=[]
    for y in range(1,n-1):
        for x in range(n,y+1,-1):
            out.append((x,y))
    return out

def rho(a,b,n): return (n+1-b, n+1-a)

def gridsym_tilings(n):
    """yield each phi-self-converse tiling as dict tile->bit (one free bit per orbit)."""
    tl=tiles(n); rmap={p:rho(*p,n) for p in tl}
    reps=[]; seen=set()
    for p in tl:
        if p in seen: continue
        q=rmap[p]; seen.add(p); seen.add(q); reps.append(p)
    for bits in itertools.product([0,1],repeat=len(reps)):
        t={}
        for i,r in enumerate(reps):
            t[r]=bits[i]; t[rmap[r]]=bits[i]
        yield t

def tiling_adj(t,n):
    A=[[0]*(n+1) for _ in range(n+1)]
    for k in range(n,1,-1): A[k][k-1]=1
    for (a,b),bit in t.items():
        if bit==0: A[a][b]=1
        else: A[b][a]=1
    return A

def c3_count(A,n):
    """directed 3-cycles = C(n,3) - sum_v C(outdeg(v),2)."""
    out=[sum(A[v][w] for w in range(1,n+1) if w!=v) for v in range(n+1)]
    return comb(n,3)-sum(comb(out[v],2) for v in range(1,n+1))

# ---- (1) c3-parity dichotomy ----
print("="*70); print("(1) c3-parity of phi-self-converse (grid-sym) tilings")
for n in range(3,8):
    par=set()
    for t in gridsym_tilings(n):
        A=tiling_adj(t,n)
        par.add(c3_count(A,n)%2)
    print(f"  n={n}: grid-sym c3 parities present = {sorted(par)}  "
          f"-> {'BOTH (odd n)' if par=={0,1} else 'EVEN only' if par=={0} else par}")
print("  predicted: EVEN-only for even n, BOTH for odd n.")

# ---- (2) Fix_anti(phi)_full = 2^(half+floor(n/2)) ----
print("="*70); print("(2) Fix_anti(phi) in full labeled space {0,1}^C(n,2)")
def half(n): return floor((n-1)**2/4)
for n in range(3,8):
    pairs=[(u,v) for u in range(1,n+1) for v in range(u+1,n+1)]
    phi=lambda i: n+1-i
    # count labeled tournaments with phi an anti-automorphism: T(phi u,phi v)=T(v,u).
    # orbit method on ordered constraints -> count = 2^(#free orbits) if consistent.
    # brute force for n<=6, orbit-count formula cross-check for all.
    # ORBIT FORMULA: #phi-orbits on unordered pairs = self-paired {u,n+1-u} + free 2-orbits.
    selfp=0; seen=set(); orb=0
    for (u,v) in pairs:
        key=frozenset([(u,v)])
        pu,pv=phi(u),phi(v)
        img=frozenset([min(pu,pv),max(pu,pv)])
        cur=frozenset([u,v])
        if cur in seen: continue
        seen.add(cur); seen.add(img); orb+=1
        if img==cur: selfp+=1
    pred_exp=orb
    formula_exp=half(n)+floor(n/2)
    line=f"  n={n}: #phi-orbits on pairs={orb} (self-paired={selfp}); half+floor(n/2)={formula_exp}"
    if n<=6:
        cnt=0
        for bits in itertools.product([0,1],repeat=len(pairs)):
            T={}
            for i,(u,v) in enumerate(pairs):
                T[(u,v)]=bits[i]; T[(v,u)]=1-bits[i]
            ok=all(T[(phi(u),phi(v))]==T[(v,u)] for (u,v) in pairs)
            if ok: cnt+=1
        line+=f"; BRUTE Fix_anti={cnt}=2^{cnt.bit_length()-1}"
        ok = (cnt==2**pred_exp==2**formula_exp)
        line+=f"  {'OK' if ok else 'MISMATCH'}"
    else:
        line+=f"  (predicted Fix_anti=2^{pred_exp})"
    print(line)

# ---- (3) SC_n self-converse iso classes vs A000568(n-1) ----
print("="*70); print("(3) SC_n (self-converse iso classes) vs A000568(n-1)")
A000568={0:1,1:1,2:1,3:2,4:4,5:12,6:56,7:456}  # tournament iso classes on k vertices
def canon(bits_pairs,n,pairs,perms):
    """canonical bit-int over sorted pairs, min over all relabelings."""
    arc=set()
    for i,(u,v) in enumerate(pairs):
        arc.add((u,v) if bits_pairs[i]==1 else (v,u))
    best=None
    pidx={p:i for i,p in enumerate(pairs)}
    for perm in perms:
        val=0
        for (a,b) in arc:
            pa,pb=perm[a-1],perm[b-1]
            uu,vv=(pa,pb) if pa<pb else (pb,pa)
            bit=1 if (pa,pb)==(uu,vv) else 0
            val|=bit<<pidx[(uu,vv)]
        if best is None or val<best: best=val
    return best
def converse_canon(canon_val,n,pairs,perms):
    """canonical form of the converse (reverse all arcs)."""
    arc=set()
    for i,(u,v) in enumerate(pairs):
        bit=(canon_val>>i)&1
        arc.add((u,v) if bit==1 else (v,u))
    rev={(b,a) for (a,b) in arc}
    pidx={p:i for i,p in enumerate(pairs)}
    best=None
    for perm in perms:
        val=0
        for (a,b) in rev:
            pa,pb=perm[a-1],perm[b-1]
            uu,vv=(pa,pb) if pa<pb else (pb,pa)
            bit=1 if (pa,pb)==(uu,vv) else 0
            val|=bit<<pidx[(uu,vv)]
        if best is None or val<best: best=val
    return best
for n in range(3,7):
    pairs=[(u,v) for u in range(1,n+1) for v in range(u+1,n+1)]
    perms=list(itertools.permutations(range(1,n+1)))
    classes=set()
    for bits in itertools.product([0,1],repeat=len(pairs)):
        classes.add(canon(bits,n,pairs,perms))
    sc=0
    for c in classes:
        if converse_canon(c,n,pairs,perms)==c: sc+=1
    print(f"  n={n}: total iso classes={len(classes)} (A000568={A000568[n]}, {'OK' if len(classes)==A000568[n] else 'CHECK'}); "
          f"SC_n={sc}; A000568(n-1)={A000568[n-1]}; "
          f"{'EQUAL' if sc==A000568[n-1] else 'DIFFERS -> refutes SC(n)=A000568(n-1)'}")
print("DONE.")
