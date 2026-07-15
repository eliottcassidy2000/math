# opus-2026-07-15-S319 -- HYP-6950: (1) Milgram coset-constancy + Gauss sums;
# (2) corona widths n = 4,6,8,10 (conjecture floor(n/2)+1); (3) the sedenion
# slice-divergence: E8^2 vs D16+ Sigma=0 slices at norms 2 and 4 (isospectral
# lattices, non-isospectral slices); (4) numeration checks (1001/1024).
from math import comb, gcd
from collections import defaultdict
import itertools, cmath

# ---------- (1) Milgram: coset constancy + Gauss sum
print("(1) MILGRAM: q mod 8 on the all-odd Sigma=0 coset; Gauss sums")
for n in (4, 6, 8, 10):
    # verify constancy on a sample + the Gauss sum of the finite form
    # L = {w in Z^n, Sum w = 0}; coset v0 + 2L with v0 = (1,...,1,1-n)... need
    # all-odd Sigma=0: v0 = (1,1,...,1,-(n-1)) has Sigma = 0 and all odd. q mod 8
    v0 = [1]*(n-1) + [-(n-1)]
    import random
    random.seed(1)
    vals = set()
    for _ in range(2000):
        w = [random.randint(-3, 3) for _ in range(n-1)]
        w.append(-sum(w))
        v = [a + 2*b for a, b in zip(v0, w)]
        vals.add(sum(t*t for t in v) % 8)
    # Gauss sum of the discriminant form on L/2L (rank n-1 over F2):
    # G = Sum_{w in L/2L} i^{q(v0+2w)... use exp(pi i q/4)} over representatives
    G = 0
    for w in itertools.product((0, 1), repeat=n-1):
        wf = list(w) + [-sum(w)]
        v = [a + 2*b for a, b in zip(v0, wf)]
        G += cmath.exp(1j*cmath.pi*sum(t*t for t in v)/4)
    phase = cmath.phase(G)/cmath.pi*4 if abs(G) > 1e-9 else None
    print(f"   n={n}: q mod 8 values on coset: {sorted(vals)} (constant: "
          f"{len(vals)==1}); |G| = {abs(G):.3f}, phase*4/pi = "
          f"{phase if phase is None else round(phase, 4)}")

# ---------- (2) corona widths
print("\n(2) CORONA WIDTHS (even n; conjecture floor(n/2)+1):")
for n in (4, 6, 8, 10):
    shells = defaultdict(set)
    rng = sorted(range(-(n-1), n, 2))
    # enumerate sorted (nondecreasing) d-multisets directly
    def gen(prefix, remaining, lo_idx):
        if remaining == 0:
            if sum(prefix) == 0:
                shells[sum(t*t for t in prefix)].add(tuple(prefix))
            return
        for i in range(lo_idx, len(rng)):
            d = rng[i]
            # prune: bounds on achievable remaining sum
            rem = remaining - 1
            smin = sum(prefix) + d + rem*d
            smax = sum(prefix) + d + rem*rng[-1]
            if smin > 0 or smax < 0: continue
            gen(prefix + [d], rem, i)
    gen([], n, 0)
    def landau_ok(dvec, n=n):
        s = sorted((d + n - 1)//2 for d in dvec)
        return all(sum(s[:k+1]) >= comb(k+1, 2) for k in range(n))
    bites = []
    for x in sorted(shells):
        tot = len(shells[x])
        leg = sum(1 for m in shells[x] if landau_ok(m))
        if leg != tot and leg > 0: bites.append(x)
    ceiling = (n**3 - n)//3
    width = len(bites)
    print(f"   n={n}: corona levels (0<legal<total): {bites}; width = {width} "
          f"(conjecture {n//2 + 1}); ceiling = {ceiling}")

# ---------- (3) the sedenion slice-divergence
print("\n(3) SEDENION RUNG: E8^2 vs D16+ Sigma=0 slices (isospectral lattices):")
def e8_vectors_by_norm(maxnorm):
    # D8+ model; enumerate integer part (small support) + half part
    out = defaultdict(list)
    # integer: even sum, coords in Z: enumerate support patterns up to norm
    # norm 2: (+-1,+-1); norm 4: (+-1)^4, (+-2)
    for pat in [(1,1), (1,1,1,1), (2,)]:
        norm = sum(t*t for t in pat)
        if norm > maxnorm: continue
        from itertools import combinations, product
        for pos in combinations(range(8), len(pat)):
            for signs in product((1,-1), repeat=len(pat)):
                v = [0]*8
                for p, base, sg in zip(pos, pat, signs): v[p] = base*sg
                if sum(v) % 2 == 0: out[norm].append(tuple(v))
    # half: (+-1/2)^8 (norm 2) and (+-3/2, +-1/2^7) (norm 4), even sum
    from itertools import product as pr
    for signs in pr((1,-1), repeat=8):
        v = [s*0.5 for s in signs]
        if sum(v) % 2 == 0: out[2].append(tuple(v))
    if maxnorm >= 4:
        for pos3 in range(8):
            for s3 in (1,-1):
                for signs in pr((1,-1), repeat=7):
                    v = []
                    k = 0
                    for i in range(8):
                        if i == pos3: v.append(s3*1.5)
                        else:
                            v.append(signs[k]*0.5); k += 1
                    if sum(v) % 2 == 0: out[4].append(tuple(v))
    return out

E8 = e8_vectors_by_norm(4)
print(f"   E8 counts: norm2 = {len(E8[2])} (240?), norm4 = {len(E8[4])} (2160?)")
# E8^2 slice: pairs (a,b), |a|^2+|b|^2 = N, Sum(a)+Sum(b) = 0
def slice_count_e8sq(N):
    bysum = {2: defaultdict(int), 4: defaultdict(int), 0: {0: 1}}
    for nm in (2, 4):
        for v in E8[nm]: bysum[nm][sum(v)] += 1
    tot = 0
    for n1 in (0, 2, 4):
        n2 = N - n1
        if n2 not in (0, 2, 4): continue
        for s, c in bysum[n1].items():
            tot += c * bysum[n2].get(-s, 0)
    return tot
def slice_count_d16p(N):
    tot = 0
    if N == 2:
        tot += 16*15   # (+1,-1) pairs ordered positions: 240
    if N == 4:
        tot += comb(16, 8)          # halves (+-1/2)^16 balanced
        tot += comb(16, 2)*comb(14, 2)  # (+1,+1,-1,-1): choose 2 pos +, 2 pos -
    return tot
for N in (2, 4):
    print(f"   norm {N}: E8^2 slice = {slice_count_e8sq(N)}, D16+ slice = "
          f"{slice_count_d16p(N)}  -> DIFFER: {slice_count_e8sq(N) != slice_count_d16p(N)}")
print("   (full-lattice theta series are EQUAL -- the Sigma=0 slice breaks "
      "isospectrality: the degeneration tournaments can see)")

# ---------- (4) numeration checks
print("\n(4) NUMERATION: 1001 and 1024")
print(f"   10^3 mod 1001 = {10**3 % 1001} (= -1 mod 1001: torsion order 2: "
      f"abc x 1001 = abcabc; alternating-block tests for 7, 11, 13)")
print(f"   143 x 7 = {143*7} (the S313 knife-edge equation)")
print(f"   T_5 = 15?? T_5 = {5*6//2}; C(5,2) = {comb(5,2)} = m at n=6; "
      f"m at n=7 = C(6,2) = {comb(6,2)} => tiling space 2^15 at n=7?? "
      f"m = C(n-1,2): n=6: {comb(5,2)}, n=7: {comb(6,2)}")
import math
print(f"   m = 10 at n = 6 (not 7): 2^10 = 1024 tilings at n = 6 -- correction")
cf = []
xv = math.log(10, 2)
a = xv
for _ in range(6):
    cf.append(int(a)); a = 1/(a - int(a))
print(f"   CF of log2(10) = {cf}: convergents 3, 10/3 (2^10 ~ 10^3), 93/28, ...")
print(f"   2^10/10^3 = 1.024: the 3-digit-shift drift per decade = the "
      f"near-dilate offset of the 10/3 convergent")
