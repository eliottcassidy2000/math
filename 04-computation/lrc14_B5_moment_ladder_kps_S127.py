# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont27: the Lean obligation hB5 IS the seven-sector moment ladder, one resolution finer.
# bandCount(v,q,p) = #{i: (v_i p mod q) NOT in [q/14, 13q/14]} = discrete empty-count N at scale 1/14 over ruler q.
# B5(v,q) = S0 - S1 + S2 - S3 + S4 - S5,  S_d = sum_p C(bandCount(p), d)  (alternating factorial-moment Bonferroni).
# THM-671: B5 <= liveCount(q); and B5 > 0 => a live multiplier p exists => t=p/q is a loneliness witness.
# hB5: EVERY residual covering family has some q with B5(v,q) > 0.  Confirmed empirically over 849 sets, not a-priori.
#
# THIS SESSION: (1) verify B5 = the moment ladder (same alternating-moment structure as THM-703's Phi-majorant);
# (2) find WHERE B5 is tight (which families barely clear 0) and whether depth-5 suffices; (3) the q-structure.
from math import comb

def bandCount(v, q, p):
    c = 0
    for vi in v:
        r = (vi * p) % q
        if not (q <= 14 * r <= 13 * q): c += 1
    return c
def momentS(v, q, d):
    return sum(comb(bandCount(v, q, p), d) for p in range(1, q))
def B5(v, q):
    return sum((-1)**d * momentS(v, q, d) for d in range(6))
def liveCount(v, q):
    return sum(1 for p in range(1, q) if bandCount(v, q, p) == 0)
def bestq(v, qmax=400):
    best = (-10**9, 0)
    for q in range(14, qmax+1):
        b = B5(v, q)
        if b > best[0]: best = (b, q)
    return best

def main():
    # canonical hard families + residual-type covering families (distinct, some |v|>=23, non-AP)
    fams = {
        'AP {1..13}':            list(range(1,14)),
        'GW (Goddyn-Wong)':      [1,2,3,4,5,6,7,8,9,10,11,12,24],
        'near-AP {1..12,26}':    list(range(1,13))+[26],
        'covering {1..11,13,26}':list(range(1,12))+[13,26],
        'covering {1..11,13,36}':list(range(1,12))+[13,36],
        'resid A':               [1,2,3,4,5,6,7,8,9,10,11,23,25],
        'resid B':               [1,2,3,4,5,7,9,11,13,15,17,23,29],
        'resid C':               [2,3,5,7,11,13,17,19,23,29,31,37,41],  # primes (dissociated)
    }
    print(" family                     bestB5  @q     liveCount@q   bandCount-hist@bestq (N -> #p)")
    tight=[]
    for name,v in fams.items():
        b,q = bestq(v)
        lc = liveCount(v,q)
        # histogram of bandCount over p in [1,q)
        hist={}
        for p in range(1,q):
            n=bandCount(v,q,p); hist[n]=hist.get(n,0)+1
        hs=' '.join(f'{k}:{hist[k]}' for k in sorted(hist)[:6])
        print(f"  {name:26s} {b:6d}  {q:4d}     {lc:6d}       {hs}")
        tight.append((b,name,q))
    print()
    tight.sort()
    print(f"  TIGHTEST family (min best-B5): {tight[0][1]} with B5={tight[0][0]} @q={tight[0][2]}")
    print(f"  => all clear B5>0 ? {'YES' if tight[0][0]>0 else 'NO -- '+tight[0][1]+' needs higher depth or larger q'}")
    # depth sufficiency: for the tightest, show S0-S1+...+(-1)^d S_d by depth d (does depth-5 = the right truncation?)
    b,name,q = tight[0]
    v = fams[name]
    print(f"\n  DEPTH LADDER for tightest ({name} @q={q}):  partial Bonferroni B_d = sum_{{i<=d}}(-1)^i S_i")
    part=0
    for d in range(0,7):
        part += (-1)**d * momentS(v,q,d)
        tag = '  <- B5 (the Lean depth)' if d==5 else ('  = liveCount (exact)' if d>=max(bandCount(v,q,p) for p in range(1,q)) else '')
        print(f"     B_{d} = {part}{tag}")
    print(f"     liveCount@q = {liveCount(v,q)} (the exact target; even-depth over-counts, odd-depth under-counts)")

if __name__ == '__main__':
    main()
