#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_lambda_tanner_kpswf5.py   (kind-pasteur 2026-06-21, THREAD B)

The TANNER GRAPH of the LRC relation code Lambda(E), and whether its
EXPANSION (girth, spectral gap, low-weight-codeword profile) gives a usable
bound on the cover error corr(E) = Sum_{0!=n in Lambda(E)} K(n).

Lambda(E) = { n in Z^k : sum_i n_i e_i = 0 }  -- the OFFSET RELATION LATTICE,
the LRC twin of the K_n cycle space (mac-mini HYP-2719a). As a linear object
it is a rank-(k-1) lattice; we treat its LOW-WEIGHT generating relations as
the parity CHECKS of an LDPC code:

  VARIABLE nodes  = the k offset elements e_1..e_k         (k of them)
  CHECK nodes     = a chosen set of primitive relations r  (support>=2,|coef|<=B)
  EDGE (v,c)      <=> e_v appears (nonzero coef) in relation c

This is DISTINCT from codex's lrc14_delsarte_tanner_carrier (which is the
Tanner graph of the THM-534 DUAL certificate: rows=moment-checks,
variables=missed-depth atoms). Here the variables ARE the offsets and the
checks ARE the relations -- the code's own Tanner graph.

THREE QUESTIONS (the prompt):
 (1) build Tanner(Lambda(E)); girth, spectral gap, support/weight distribution.
     Does AP/consec give low-girth / poor-expansion (anti-MDS, hardest = HYP-2723)
     while Sidon/arc gives high-girth / good expansion?
 (2) does Tanner EXPANSION (large girth / few low-weight codewords) BOUND |corr|?
     Binding case = low expansion = AP.
 (3) circulant/QR_7 (Paley) relation code -- weakly-regular? small cases.

EXACT measS7/corr (reused from lrc_q108_relation_code_mds_kps). Tanner spectra
in float (eigenvalues of a real symmetric incidence -- fine for the gap).
Save output to 05-knowledge/results/lrc_lambda_tanner_kpswf5.out
"""
import itertools
from fractions import Fraction as Fr
from math import comb, gcd, factorial, sin, pi
from functools import reduce
import numpy as np

# ============================================================ measS7 / corr (EXACT)
def stirling2(n, k):
    if k == 0:
        return 1 if n == 0 else 0
    S = [[0]*(k+1) for _ in range(n+1)]
    S[0][0] = 1
    for i in range(1, n+1):
        for j in range(1, min(i, k)+1):
            S[i][j] = j*S[i-1][j] + S[i-1][j-1]
    return S[n][k]

def iid_k(k, p=7):
    return Fr(factorial(p) * stirling2(k, p), p**k)

def measS7(E, p=7):
    E = [int(e) for e in E]
    bps = set([Fr(0), Fr(1)])
    for e in E:
        if e == 0:
            continue
        for t in range(0, p*e):
            bps.add(Fr(t, p*e))
    pts = sorted(bps)
    total = Fr(0)
    for a, b in zip(pts, pts[1:]):
        mid = (a + b) / 2
        sectors = set()
        for e in E:
            y = (e*mid) % 1
            s = int(p*y)
            sectors.add(s)
            if len(sectors) == p:
                break
        if len(sectors) == p:
            total += (b - a)
    return total

def corr(E, p=7):
    return float(measS7(E, p) - iid_k(len(E), p))

# ===================================================== relation enumeration
def primitive_relations(E, B=2, max_support=4):
    """All PRIMITIVE nonzero integer relations sum n_i e_i = 0 with |n_i|<=B,
       support in [2,max_support]. Returns list of full-length-k coef tuples
       (canonical sign: first nonzero positive), de-duplicated by primitivity."""
    E = [int(e) for e in E]
    k = len(E)
    seen = set()
    rels = []
    for s in range(2, max_support+1):
        for combo in itertools.combinations(range(k), s):
            for coefs in itertools.product(range(-B, B+1), repeat=s):
                if any(c == 0 for c in coefs):
                    continue
                if sum(c*E[i] for c, i in zip(coefs, combo)) != 0:
                    continue
                g = reduce(gcd, [abs(c) for c in coefs])
                prim = tuple(c//g for c in coefs)
                if prim[0] < 0:
                    prim = tuple(-c for c in prim)
                full = [0]*k
                for c, i in zip(prim, combo):
                    full[i] = c
                key = tuple(full)
                if key in seen:
                    continue
                seen.add(key)
                rels.append(key)
    return rels

# ===================================================== Tanner graph + spectra
def tanner_incidence(rels, k):
    """0/1 variable x check incidence (presence of variable in relation)."""
    C = len(rels)
    H = np.zeros((k, C), dtype=float)
    for c, r in enumerate(rels):
        for v in range(k):
            if r[v] != 0:
                H[v, c] = 1.0
    return H

def bipartite_spectral_gap(H):
    """Bipartite graph on k variables + C checks, biadjacency H (0/1).
       Return (lambda2 of normalized adjacency, second singular value of
       normalized H). Spectral expansion of the bipartite Tanner graph.
       Uses the normalized bipartite Laplacian / singular values of
       D_v^{-1/2} H D_c^{-1/2}.  The bipartite expander quality is governed
       by the SECOND-largest singular value sigma_2 in (0,1); small=good."""
    k, C = H.shape
    dv = H.sum(axis=1)   # variable degrees
    dc = H.sum(axis=0)   # check degrees
    if (dv <= 0).any() or (dc <= 0).any():
        return None, None, dv, dc  # isolated node -> not connected
    Dv = np.diag(1.0/np.sqrt(dv))
    Dc = np.diag(1.0/np.sqrt(dc))
    M = Dv @ H @ Dc      # normalized biadjacency
    svals = np.linalg.svd(M, compute_uv=False)
    sigma1 = svals[0]              # =1 for connected bipartite
    sigma2 = svals[1] if len(svals) > 1 else 0.0
    return sigma1, sigma2, dv, dc

def tanner_girth(H, max_g=12):
    """Girth of the bipartite Tanner graph (shortest cycle length, even).
       BFS from each node; returns smallest even cycle, or inf."""
    k, C = H.shape
    # build adjacency: nodes 0..k-1 = variables, k..k+C-1 = checks
    N = k + C
    adj = [[] for _ in range(N)]
    for v in range(k):
        for c in range(C):
            if H[v, c]:
                adj[v].append(k+c)
                adj[k+c].append(v)
    best = float('inf')
    for src in range(N):
        # BFS tracking parent to find shortest cycle through src
        dist = [-1]*N
        par = [-1]*N
        dist[src] = 0
        from collections import deque
        q = deque([src])
        while q:
            u = q.popleft()
            if dist[u]*2 >= best:
                continue
            for w in adj[u]:
                if dist[w] == -1:
                    dist[w] = dist[u]+1
                    par[w] = u
                    q.append(w)
                elif par[u] != w:
                    cyc = dist[u]+dist[w]+1
                    if cyc < best:
                        best = cyc
        if best <= 4:
            break
    return best

# ===================================================== low-weight codeword profile
def lowweight_codeword_profile(E, B=2, max_support=6):
    """Weight (=support) distribution of low-coef relations -> A_s spectrum
       and dmin (minimum support). The 'code' weight enumerator (low part)."""
    E = [int(e) for e in E]
    k = len(E)
    counts = {s: 0 for s in range(2, max_support+1)}
    seen = set()
    for s in range(2, max_support+1):
        for combo in itertools.combinations(range(k), s):
            for coefs in itertools.product(range(-B, B+1), repeat=s):
                if any(c == 0 for c in coefs):
                    continue
                if sum(c*E[i] for c, i in zip(coefs, combo)) != 0:
                    continue
                g = reduce(gcd, [abs(c) for c in coefs])
                prim = tuple(c//g for c in coefs)
                if prim[0] < 0:
                    prim = tuple(-c for c in prim)
                key = (combo, prim)
                if key in seen:
                    continue
                seen.add(key)
                counts[s] += 1
    nz = [s for s in counts if counts[s] > 0]
    dmin = min(nz) if nz else None
    return counts, dmin

# ===================================================== weakly-regular test
def is_weakly_regular_bipartite(dv, dc):
    """Weakly regular (biregular): all variable degrees equal AND all check
       degrees equal (possibly two different constants)."""
    vu = set(np.round(dv).astype(int).tolist())
    cu = set(np.round(dc).astype(int).tolist())
    return (len(vu) == 1, len(cu) == 1, sorted(vu), sorted(cu))

# ===================================================== batteries
def is_sidon(E):
    sums = {}
    E = list(E)
    for i in range(len(E)):
        for j in range(i, len(E)):
            s = E[i]+E[j]
            if s in sums:
                return False
            sums[s] = 1
    return True

def main():
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a)
        out.append(s)
        print(s)

    p = 7
    P("="*86)
    P("THREAD B: the Tanner graph of the relation code Lambda(E) and its expansion vs corr(E)")
    P("="*86)
    P("VARIABLES = k offsets; CHECKS = primitive relations (support>=2, |coef|<=2, support<=4)")
    P("Tanner expansion = (girth, sigma_2 of normalized biadjacency); want corr small <=> good exp?")
    P("")

    # k=8 battery (0 = pinned observer)
    battery = {
        "consec/AP {0..7}":        [0,1,2,3,4,5,6,7],
        "AP step2":                [0,2,4,6,8,10,12,14],
        "AP step3":                [0,3,6,9,12,15,18,21],
        "dyadic {0,1,2,4,..64}":   [0,1,2,4,8,16,32,64],
        "Sidon (Mian-Chowla)":     [0,1,3,7,12,20,30,44],
        "Sidon b=2 perfectdiff":   [0,1,3,7,15,31,63,127],
        "wide near-consec":        [0,1,2,3,4,5,6,40],
        "two-block consec":        [0,1,2,3,40,41,42,43],
        "random wide":             [0,5,9,14,22,33,41,50],
        "QR7-lift {0,1,2,4}+block": [0,1,2,4,7,8,9,11],
    }
    P(f"{'set':<26}{'corr':>9}{'dmin':>6}{'#chk':>6}{'girth':>7}{'sig2':>9}{'A2':>5}{'A3':>5}{'A4':>5}{'biReg?':>8}")
    rows = []
    for name, E in battery.items():
        c = corr(E, p)
        rels = primitive_relations(E, B=2, max_support=4)
        H = tanner_incidence(rels, len(E))
        counts, dmin = lowweight_codeword_profile(E, B=2, max_support=4)
        if len(rels) == 0:
            girth = float('inf'); sig2 = None; bireg=(False,False,[],[])
            dv = H.sum(axis=1); dc = np.array([])
        else:
            sig1, sig2, dv, dc = bipartite_spectral_gap(H)
            girth = tanner_girth(H)
            bireg = is_weakly_regular_bipartite(dv, dc) if sig2 is not None else (False,False,[],[])
        sid = is_sidon([e for e in E if e != 0])
        gstr = ("inf" if girth==float('inf') else str(girth))
        s2str = ("-" if sig2 is None else f"{sig2:.4f}")
        rows.append(dict(name=name,E=E,corr=c,dmin=dmin,nchk=len(rels),
                         girth=girth,sig2=sig2,counts=counts,sidon=sid,bireg=bireg))
        P(f"{name:<26}{c:>9.4f}{str(dmin):>6}{len(rels):>6}{gstr:>7}{s2str:>9}"
          f"{counts.get(2,0):>5}{counts.get(3,0):>5}{counts.get(4,0):>5}"
          f"{('Y' if (bireg[0] and bireg[1]) else 'n'):>8}")

    P("")
    P("READING:")
    P("  HYP-2723/2724 predicts AP/consec = anti-MDS = small dmin, many low-weight checks,")
    P("  DENSE Tanner, low girth (=4), LARGE sigma_2 (poor expansion), LARGEST corr.")
    P("  Sidon/arc = MDS = large dmin, FEW checks, sparse Tanner, large girth, sigma_2 small")
    P("  (good expansion), corr ~ 0.")
    P("")

    # correlations
    import statistics
    def pearson(xs, ys):
        xs2 = [x for x,y in zip(xs,ys) if x is not None and y is not None]
        ys2 = [y for x,y in zip(xs,ys) if x is not None and y is not None]
        n = len(xs2)
        if n < 3 or len(set(xs2))<2 or len(set(ys2))<2: return None
        mx=sum(xs2)/n; my=sum(ys2)/n
        num=sum((a-mx)*(b-my) for a,b in zip(xs2,ys2))
        den=(sum((a-mx)**2 for a in xs2)*sum((b-my)**2 for b in ys2))**0.5
        return num/den
    cr   = [r['corr'] for r in rows]
    sig2 = [r['sig2'] for r in rows]
    nchk = [float(r['nchk']) for r in rows]
    A3   = [float(r['counts'].get(3,0)) for r in rows]
    dmin = [float(r['dmin']) if r['dmin'] else None for r in rows]
    P(f"  Pearson(corr, sigma_2)      = {fmtp(pearson(cr,sig2))}   (HYP: + ; poor expansion => large corr)")
    P(f"  Pearson(corr, #checks)      = {fmtp(pearson(cr,nchk))}   (HYP: + ; dense relation set => hard)")
    P(f"  Pearson(corr, A3)           = {fmtp(pearson(cr,A3))}   (canon HYP-2724: +0.93)")
    P(f"  Pearson(corr, dmin)         = {fmtp(pearson(cr,dmin))}   (HYP: - ; large min distance => easy)")
    return out, rows

def fmtp(x):
    return "  n/a" if x is None else f"{x:+.3f}"

if __name__ == "__main__":
    out, rows = main()
    import os
    os.makedirs("05-knowledge/results", exist_ok=True)
    with open("05-knowledge/results/lrc_lambda_tanner_kpswf5.out","w") as f:
        f.write("\n".join(out)+"\n")
