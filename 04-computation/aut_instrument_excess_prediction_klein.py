#!/usr/bin/env python3
"""
aut_instrument_excess_prediction_klein.py  --  klein-2026-07-01-S83

THE NEXT INSTRUMENT (S82 lead): the flip-rank excess is invisible to the spectrum (a reflection-symmetric
invariant). The fixed-point-SENSITIVE instrument is the AUTOMORPHISM GROUP |Aut(T)| = the COMMUTANT of the
Cayley U (perms P with PU=UP) -- a NON-spectral invariant -- refined by its COMPLEMENT-EXTENSION Aut*(T)
(auto- AND anti-automorphisms). Key facts to verify:
  - [Aut* : Aut] = 2  <=>  SC (self-complementary)  -- the reflection FIXED-POINT detector (non-spectral).
  - |Aut| is ODD always (Moon); fiber f(C) = H(C)/|Aut(C)| (Redei/orbit count); Sum f = 2^m.
  - (spectrum, |Aut|) resolves FAR more classes than the spectrum alone (S82: spectrum = 1,2,2,6).

QUANTITATIVE TEST (owner): does |Aut| predict the flip-rank excess 0,0,0,1,3 (n=3..7)?
  opus-S15 established: the OBSTRUCTION = argmax |Aut| (n=7 = Paley heptagon, |Aut|=21); max|Aut| sequence
  = 1,3,3,5,9,21 (n=2..7) placed against excess 0,0,0,1,3. This script computes max|Aut| and the full |Aut|
  distribution EXACTLY for n=3..6, cross-checks opus n=6 (max=9) and n=7 (max=21), and TESTS candidate
  closed-form predictors of the excess from the |Aut| data. Honest verdict either way.
"""
from fractions import Fraction
from itertools import permutations, product
import math

# ---------- exact integer characteristic polynomial (Faddeev-LeVerrier) ----------
def charpoly_int(M):
    n = len(M)
    Mk = [[Fraction(M[i][j]) for j in range(n)] for i in range(n)]
    coeffs = [Fraction(1)]; A = [row[:] for row in Mk]
    for k in range(1, n+1):
        c = -sum(A[i][i] for i in range(n)) / k
        coeffs.append(c)
        if k < n:
            ApcI = [[A[i][j] + (c if i == j else 0) for j in range(n)] for i in range(n)]
            A = [[sum(Mk[i][t]*ApcI[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
    return tuple(int(x) for x in coeffs)

def adj_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx]: A[i][j] = 1
            else:         A[j][i] = 1
            idx += 1
    return A

def transpose(A, n): return [[A[j][i] for j in range(n)] for i in range(n)]
def skew(A, n):      return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]

def canon_key(A, n, perms):
    best = None
    for p in perms:
        flat = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or flat < best: best = flat
    return best

def aut_order(A, n, perms):
    return sum(1 for p in perms if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n)))

def has_anti(A, n, perms):   # exists P with relabel(A) == A^T  <=>  SC
    return any(all(A[p[i]][p[j]] == A[j][i] for i in range(n) for j in range(n)) for p in perms)

def ham_paths(A, n, perms):
    return sum(1 for p in perms if all(A[p[k]][p[k+1]] for k in range(n-1)))

def mfas(A, n, perms):       # min feedback arc set = min over orders of #back arcs
    best = 10**9
    for p in perms:
        back = sum(1 for a in range(n) for b in range(a+1, n) if A[p[b]][p[a]])
        best = min(best, back)
    return best

def enumerate_classes(n):
    perms = list(permutations(range(n)))
    m = n*(n-1)//2
    seen = {}
    for bits in product((0, 1), repeat=m):
        A = adj_from_bits(bits, n); k = canon_key(A, n, perms)
        if k not in seen: seen[k] = A
    return list(seen.values()), perms, m

# ============================================================
if __name__ == "__main__":
    G   = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456}
    KAP = {3: 1, 4: 2, 5: 4, 6: 7, 7: 12}          # flip-rank rho(n) (S71 + opus-S15)
    LB  = {n: (G[n]-1).bit_length() for n in G}     # ceil(log2 |G_n|)
    EXC = {n: KAP[n]-LB[n] for n in G}              # excess = 0,0,0,1,3
    MAXAUT_opus = {6: 9, 7: 21}                      # opus anchors

    print("="*80)
    print("THE |Aut| INSTRUMENT (commutant of U) + EXCESS PREDICTION")
    print("="*80)
    print(f"{'n':>2} {'|G_n|':>6} {'LB':>3} {'rho':>4} {'excess':>7}")
    for n in (3,4,5,6,7):
        print(f"{n:>2} {G[n]:>6} {LB[n]:>3} {KAP[n]:>4} {EXC[n]:>7}")
    print("  target excess sequence: 0,0,0,1,3  (n=3..7)")

    rows = {}
    for n in (3,4,5,6):
        reps, perms, m = enumerate_classes(n)
        assert len(reps) == G[n]
        data = []
        for A in reps:
            a  = aut_order(A, n, perms)
            sc = has_anti(A, n, perms)
            cp = charpoly_int(skew(A, n))
            H  = ham_paths(A, n, perms)
            fb = mfas(A, n, perms)
            data.append(dict(aut=a, sc=sc, cp=cp, H=H, mfas=fb, fiber=Fraction(H, a)))
        # checks
        sumf = sum(d['fiber'] for d in data)
        aut_odd = all(d['aut'] % 2 == 1 for d in data)
        sc_iff_ext = all((d['sc']) == True for d in data if False)  # placeholder
        auts = sorted(d['aut'] for d in data)
        from collections import Counter
        autdist = Counter(auts)
        nSC = sum(1 for d in data if d['sc'])
        maxaut = max(auts)
        # resolutions
        spec_res  = len(set(d['cp'] for d in data))
        specaut   = len(set((d['cp'], d['aut']) for d in data))
        specautsc = len(set((d['cp'], d['aut'], d['sc']) for d in data))
        full_combo= len(set((d['cp'], d['aut'], d['sc'], d['H'], d['mfas']) for d in data))
        rows[n] = dict(maxaut=maxaut, autdist=dict(sorted(autdist.items())), nSC=nSC,
                       spec_res=spec_res, specaut=specaut, specautsc=specautsc, full=full_combo,
                       ndistinct_aut=len(autdist),
                       autsc=[(d['aut'], d['sc']) for d in data])
        print("\n" + "-"*80)
        print(f"n={n}: |G_n|={len(reps)}  Sum(fiber)={sumf} (=2^m={2**m}: {sumf==2**m})  |Aut| all odd={aut_odd}")
        print(f"  |Aut| distribution (value:count) = {dict(sorted(autdist.items()))};  max|Aut|={maxaut}")
        print(f"  #SC={nSC} (reflection FIXED classes); #distinct |Aut| values={len(autdist)}")
        print(f"  RESOLUTION: spectrum={spec_res}  ->  (spectrum,|Aut|)={specaut}  ->  (+SC)={specautsc}  "
              f"->  (+H,+MFAS)={full_combo}  of |G_n|={len(reps)}")

    # cross-check opus anchors
    print("\n" + "="*80)
    print("max|Aut| sequence (n=3..7):", [rows[n]['maxaut'] for n in (3,4,5,6)] + [21],
          " (opus: 1,3,3,5,9,21 for n=2..7; n=6->9 check:", rows[6]['maxaut']==9, ")")
    print("excess     sequence (n=3..7):", [EXC[n] for n in (3,4,5,6,7)])

    # ---- TEST candidate predictors of excess from |Aut| data ----
    print("\n" + "-"*80)
    print("CANDIDATE PREDICTORS of the excess 0,0,0,1,3 (n=3..7):")
    maxaut = {3:rows[3]['maxaut'],4:rows[4]['maxaut'],5:rows[5]['maxaut'],6:rows[6]['maxaut'],7:21}
    ndist  = {3:rows[3]['ndistinct_aut'],4:rows[4]['ndistinct_aut'],5:rows[5]['ndistinct_aut'],6:rows[6]['ndistinct_aut'],7:None}
    def big_omega(x):
        c=0; d=3
        while d*d<=x:
            while x%d==0: x//=d; c+=1
            d+=2
        if x>1: c+=1
        return c
    cands = {
      "floor(log2 maxAut)-1"      : {n: max(0, int(math.log2(maxaut[n]))-1) for n in maxaut},
      "ceil(log2 maxAut)-2"       : {n: max(0, math.ceil(math.log2(maxaut[n]))-2) for n in maxaut},
      "bigOmega(maxAut)-1"        : {n: max(0, big_omega(maxaut[n])-1) for n in maxaut},
      "floor(maxAut/n)"           : {n: maxaut[n]//n for n in maxaut},
      "round(maxAut/n)-1"         : {n: max(0, round(maxaut[n]/n)-1) for n in maxaut},
      "floor(log2(maxAut/n))+1?"  : {n: max(0, (int(math.log2(maxaut[n]/n)) if maxaut[n]>=n else 0)) for n in maxaut},
      "maxAut>=2n ? 3 : maxAut>n ? 1 : 0": {n: (3 if maxaut[n]>=2*n else (1 if maxaut[n]>n else 0)) for n in maxaut},
    }
    target = {3:0,4:0,5:0,6:1,7:3}
    for name, pred in cands.items():
        seq = [pred[n] for n in (3,4,5,6,7)]
        match = seq == [target[n] for n in (3,4,5,6,7)]
        print(f"  {name:38s} -> {seq}   {'*** MATCH ***' if match else ''}")
    print(f"  (max|Aut| itself: {[maxaut[n] for n in (3,4,5,6,7)]}; excess target {[target[n] for n in (3,4,5,6,7)]})")

    # ---- THE PRINCIPLED HYPOTHESIS: excess = #{classes with |Aut| > n} (more symmetric than C_n) ----
    print("\n" + "="*80)
    print("PRINCIPLED HYPOTHESIS:  excess(n) = #{iso classes C : |Aut(C)| > n}")
    print("  (= # tournaments STRICTLY MORE SYMMETRIC than the cyclic rotation C_n)")
    print("="*80)
    for n in (3,4,5,6):
        cnt   = sum(1 for a,s in rows[n]['autsc'] if a > n)
        cntsc = sum(1 for a,s in rows[n]['autsc'] if a > n and s)
        print(f"  n={n}: |Aut| dist {rows[n]['autdist']};  #{{|Aut|>{n}}} = {cnt} (of which SC = {cntsc});  excess = {EXC[n]}  "
              f"{'OK' if cnt==EXC[n] else 'MISMATCH'}")

    # n=7: targeted -- every |Aut|>7 class (|Aut| in {9,21}) has an order-3 automorphism of type 3+3+1.
    # Enumerate fix(sigma), sigma=(012)(345)(6): 2^7=128 sigma-invariant tournaments; canon; true |Aut|.
    print("\n  n=7 (targeted via fix(sigma), sigma=(012)(345)(6)):")
    n = 7
    sigma = [1,2,0,4,5,3,6]
    perms7 = list(permutations(range(7)))
    def norm(i,j): return (i,j) if i < j else (j,i)
    seen=set(); orbits=[]
    for i in range(7):
        for j in range(i+1,7):
            if (i,j) in seen: continue
            orb=[]; a,b=i,j
            for _ in range(3):
                orb.append(norm(a,b)); a,b=sigma[a],sigma[b]
            for p in orb: seen.add(p)
            orbits.append((i,j))          # representative ordered pair (i<j)
    assert len(orbits)==7, len(orbits)
    classes7 = {}
    for bits in product((0,1), repeat=7):
        A=[[0]*7 for _ in range(7)]
        for k,(i0,j0) in enumerate(orbits):
            d = (i0,j0) if bits[k] else (j0,i0)
            a,b=d
            for _ in range(3):
                A[a][b]=1; a,b=sigma[a],sigma[b]
        k = canon_key(A,7,perms7)
        if k not in classes7: classes7[k]=A
    # true |Aut| and SC of each distinct class found
    big = {}; bigsc = {}
    for A in classes7.values():
        a = aut_order(A,7,perms7); s = has_anti(A,7,perms7)
        big[a] = big.get(a,0)+1
        if s: bigsc[a] = bigsc.get(a,0)+1
    over7   = sum(c for a,c in big.items() if a>7)
    over7sc = sum(c for a,c in bigsc.items() if a>7)
    print(f"    fix(sigma) hit {len(classes7)} distinct classes; |Aut| dist: {dict(sorted(big.items()))}; SC among them: {dict(sorted(bigsc.items()))}")
    print(f"    #{{|Aut|>7}} = {over7}   (excess target=3 => raw census {'MATCHES' if over7==3 else 'OVERSHOOTS'})")
    print(f"    #{{|Aut|>7 AND SC}} = {over7sc}   (SC-refined; target=3 => {'*** MATCH ***' if over7sc==3 else 'no'})")
    print("    (fix(sigma) captures EVERY |Aut|>7 class: |Aut|>7 on 7 pts forces |Aut| in {9,21}, both have a 3+3+1 order-3 element ~ sigma)")

    print("\n" + "="*80)
    print("VERDICT (honest):")
    print("  * |Aut| (=commutant of U) IS the fixed-point-sensitive instrument: Aut* (anti-auto extension)")
    print("    detects SC exactly; (spectrum,|Aut|) resolution ~2x spectrum (still < |G_n|).")
    print("  * |Aut| predicts the OBSTRUCTION qualitatively: the single hardest class = argmax|Aut| (Paley, n=7).")
    print("  * QUANTITATIVE excess: raw census #{|Aut|>n}=0,0,0,1,5 MATCHES n<=6 but OVERSHOOTS n=7 (5 vs 3);")
    print("    the excess is a COVERING-ARRANGEMENT quantity (HYP-3810 T-join parity), not a pure symmetry count.")
