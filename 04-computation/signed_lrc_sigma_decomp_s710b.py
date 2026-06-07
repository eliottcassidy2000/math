"""Signed-LRC: NON-BRUTE deficiency via the sigma-decomposition of silent(D).  (monad-explorer-S710)

KEY IDEA.  For a move D (subset of magnitudes flipped), silent(D) = {eps : eps^D ~ eps}
splits by the A.B lemma (S708b): for all freq t, A_t=0 OR B_t=0, where
  B_t(eps) = sum_{i in D}   eps_i sin(2pi t i /C)   (depends on eps|_D =: sigma only)
  A_t(eps) = sum_{i in D^c} eps_i sin(2pi t i /C)   (depends on eps|_{D^c} only).
Partition by sigma = eps|_D (2^|D| patterns).  Given sigma, the move is silent iff
A_t = 0 for every t with B_t(sigma) != 0.  So
  |silent(D)| = sum_{sigma in {+-1}^D}  N(Req(sigma)),
  Req(sigma) = { t in 1..(C-1)/2 : B_t(sigma) != 0 },
  N(R) = #{ eps|_{D^c} in {+-1} : sum_{i in D^c} eps_i sin(2pi t i/C)=0  for all t in R }.

CONJECTURE (validated here): N(R) = 2^{ |D^c| - rank_R(Mat) }  (the +-1 solutions of these
SINE systems form an F2-coset of full real-nullity dimension; value-pairing structure THM-413).

This avoids the 2^{n-2} cut enumeration: each move needs 2^|D| small SVDs, |D| small for the
order-block moves.  Reaches C=63,75,81,99,105 which brute force cannot.

We VALIDATE silent(D) (sigma-decomp) vs BRUTE for every move D in V, all composite C<=27
(and 45 if feasible), THEN deploy the deficiency = Moebius-over-V framework.
"""
import sys
import numpy as np
from itertools import product
from collections import defaultdict, Counter

TOL = 1e-7

def factor_str(C):
    f, m, d = [], C, 2
    while d * d <= m:
        while m % d == 0:
            f.append(d); m //= d
        d += 1
    if m > 1:
        f.append(m)
    return "x".join(map(str, f)) if f else str(C)

def proper_divisors(C):
    return [d for d in range(2, C) if C % d == 0]

def order_d_halfsystem(C, d):
    K = set(((C // d) * j) % C for j in range(d))
    return frozenset(x for x in K if 1 <= x <= (C - 1) // 2)

def sine_row(C, mags, t):
    return np.array([np.sin(2*np.pi*t*i/C) for i in mags])

def nullity_pm1(C, Dc_mags, Req):
    """real-nullity of the |Req| x |Dc| sine matrix -> 2^nullity (conjectured +-1 solution count)."""
    k = len(Dc_mags)
    if not Req:
        return k                      # no constraints -> all 2^k free
    M = np.array([[np.sin(2*np.pi*t*i/C) for i in Dc_mags] for t in Req])
    if M.size == 0:
        return k
    s = np.linalg.svd(M, compute_uv=False)
    smax = s[0] if len(s) else 0.0
    rank = int(np.sum(s > TOL*max(smax,1.0)))
    return k - rank

def silent_count_sigma(C, n, Dmask, free_to_mag):
    """sigma-decomposition count of |silent(D)| where Dmask is over free coords (mags>=2)."""
    nfree = n - 2
    Dmags = [free_to_mag[k] for k in range(nfree) if (Dmask>>k)&1]
    Dcmags = [free_to_mag[k] for k in range(nfree) if not (Dmask>>k)&1]
    allT = list(range(1, (C-1)//2 + 1))
    # Precompute sin(t*i) for i in Dmags
    total = 0
    dlen = len(Dmags)
    # cache nullity by frozenset(Req)
    cache = {}
    for sbits in range(1 << dlen):
        sigma = [1 if (sbits>>j)&1 else -1 for j in range(dlen)]
        Req = []
        for t in allT:
            B = 0.0
            for j,i in enumerate(Dmags):
                B += sigma[j]*np.sin(2*np.pi*t*i/C)
            if abs(B) > TOL:
                Req.append(t)
        key = tuple(Req)
        if key not in cache:
            cache[key] = nullity_pm1(C, Dcmags, Req)
        total += (1 << cache[key])
    return total

# ---------- brute for validation ----------
def diff_multiset(S, C):
    cnt = [0]*C
    for a in S:
        for b in S:
            if a!=b: cnt[(a-b)%C]+=1
    return tuple(cnt)

def brute_sigs(n):
    C=2*n-1; nfree=n-2; runners=list(range(1,n))
    bits2sig={}
    for bits in range(1<<nfree):
        eps=[1]*(n-1)
        for b in range(nfree):
            if (bits>>b)&1: eps[b+1]=-1
        S=sorted((eps[i]*runners[i])%C for i in range(n-1))
        bits2sig[bits]=diff_multiset(S,C)
    return bits2sig

def V_elements(C, n):
    nfree=n-2; runners=list(range(1,n))
    mag_to_free={runners[k]:k-1 for k in range(1,n-1)}
    divs=proper_divisors(C)
    gens=[]
    for d in divs:
        Hd=order_d_halfsystem(C,d)
        if 1 in Hd: continue
        msk=0
        for x in Hd: msk|=1<<mag_to_free[x]
        gens.append((d,msk))
    basis=[]
    for _,msk in gens:
        cur=msk
        for b in basis: cur=min(cur,cur^b)
        if cur: basis.append(cur); basis.sort(reverse=True)
    elems=[0]
    for b in basis: elems=elems+[x^b for x in elems]
    return sorted(set(elems)), gens, len(basis)

def validate(n):
    C=2*n-1; nfree=n-2; runners=list(range(1,n))
    free_to_mag={k-1:runners[k] for k in range(1,n-1)}
    elems,gens,dimV=V_elements(C,n)
    bits2sig=brute_sigs(n)
    print(f"\n=== validate C={C}={factor_str(C)} n={n} dimV={dimV} ===")
    ok=True
    for D in elems:
        if D==0: continue
        brute=sum(1 for b in range(1<<nfree) if bits2sig[b]==bits2sig[b^D])
        sigma_c=silent_count_sigma(C,n,D,free_to_mag)
        mags=[free_to_mag[k] for k in range(nfree) if (D>>k)&1]
        tag=next((f"=H_{d}" for d,m in gens if m==D), "")
        match = (brute==sigma_c)
        ok = ok and match
        print(f"  D={mags}{tag}: brute={brute} sigma={sigma_c} {'OK' if match else '*** MISMATCH ***'}")
    print(f"  ALL MATCH: {ok}")
    return ok

if __name__=="__main__":
    allok=True
    for n in [5,8,11,13,14,17,18,20]:    # C=9,15,21,25,27,33,35,39
        allok = validate(n) and allok
    print(f"\n\n###### VALIDATION (sigma-decomp == brute) ALL C<=39: {allok} ######")
