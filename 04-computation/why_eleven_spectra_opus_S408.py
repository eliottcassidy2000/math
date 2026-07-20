"""opus-2026-07-20-S408 -- WHY ONLY ELEVEN SPECTRA?

THM-1450 found exactly 11 distinct characteristic polynomials among all 32768 switching
classes of 7-tournaments.  This asks why.

THE CANDIDATE MECHANISM.  Write S for the skew +-1 matrix.  Then
  c5 = e2 = C(7,2) = 21                                    (always)
  c3 = e4 = sum over 4-subsets T of Pf(S_T)^2
  c1 = e6 = sum over 6-subsets T of Pf(S_T)^2
For a 4-subset, Pf = s_ab s_cd - s_ac s_bd + s_ad s_bc, three +-1 terms, so
Pf in {-3,-1,1,3} and Pf^2 in {1,9}.  Hence

  c3 = 35 + 8k,   k = #{4-subsets with |Pf| = 3}          => c3 = 3 (mod 8)

which already forces a congruence.  Check against THM-1450's table:
  115,131,83,99,67 are all = 3 mod 8.  So the congruence holds; but it allows 36 values
  of c3 and we only see a few.  Something stronger is at work.

THE REAL INVARIANT.  Switching at W flips s_ij exactly when |{i,j} cap W| = 1.  Around a
TRIANGLE the number of flipped arcs is 0 or 2 -- never 1 or 3 -- so the PRODUCT
  g(a,b,c) = s_ab * s_bc * s_ca
is a SWITCHING INVARIANT.  That is precisely the "oriented two-graph" / cyclic-triple
indicator of THM-474.  So the triple statistic
  t = #{ triples with g = +1 }
is a switching-class invariant, and it is the natural candidate to control (c3, c1).

TEST: does t determine the spectrum, and does t take exactly 11 values?
"""
import itertools, numpy as np
from collections import defaultdict

n = 7
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
base = {(i, i+1) for i in range(n-1)}
free = [p for p in pairs if p not in base]
tris = list(itertools.combinations(range(n), 3))

def charpoly_coeffs(S):
    """exact integer char poly via Faddeev-LeVerrier"""
    Mk = np.zeros((n, n), dtype=object); coeffs = [1]; Ak = np.eye(n, dtype=object)
    for k in range(1, n+1):
        Ak = S.astype(object) @ (Ak if k == 1 else Mk)
        ck = -int(np.trace(Ak)) // k
        coeffs.append(ck)
        Mk = Ak + ck*np.eye(n, dtype=object)
    return tuple(coeffs)

def triple_stat(S):
    """t = number of triples with s_ab * s_bc * s_ca = +1  (switching invariant)"""
    return sum(1 for (a, b, c) in tris if S[a, b]*S[b, c]*S[c, a] == 1)

spec_by_t = defaultdict(set); t_by_spec = defaultdict(set)
count_t = defaultdict(int); count_spec = defaultdict(int)
sw_inv_fail = 0
for mask in range(1 << len(free)):
    S = np.zeros((n, n), dtype=np.int64)
    for (i, j) in base: S[i, j] = 1; S[j, i] = -1
    for b, (i, j) in enumerate(free):
        v = 1 if (mask >> b) & 1 else -1
        S[i, j] = v; S[j, i] = -v
    co = charpoly_coeffs(S); t = triple_stat(S)
    spec_by_t[t].add(co); t_by_spec[co].add(t)
    count_t[t] += 1; count_spec[co] += 1

print("="*76)
print("(1) DOES THE TRIPLE STATISTIC t DETERMINE THE SPECTRUM?")
print("="*76)
bad = {t: s for t, s in spec_by_t.items() if len(s) > 1}
print(f"   distinct t values : {len(count_t)}   -> {sorted(count_t)}")
print(f"   distinct spectra  : {len(count_spec)}")
print(f"   t values mapping to MORE THAN ONE spectrum: {len(bad)}")
print(f"   spectra arising from MORE THAN ONE t:       "
      f"{len({c for c,s in t_by_spec.items() if len(s)>1})}")
print()
print("   t      c3     c1    #switching classes")
for t in sorted(count_t):
    for co in sorted(spec_by_t[t]):
        print(f"  {t:3d}   {co[4]:5d}  {co[6]:5d}   {count_t[t]:8d}")

print()
print("="*76)
print("(2) THE CONGRUENCE  c3 = 35 + 8k  AND WHAT k COUNTS")
print("="*76)
print("   c3     c3 mod 8    k = (c3-35)/8")
for co in sorted(count_spec, key=lambda c: c[4]):
    c3 = co[4]
    print(f"  {c3:5d}      {c3 % 8}          {(c3-35)/8}")

print()
print("="*76)
print("(3) IS t THE CYCLIC-TRIPLE COUNT, AND WHAT IS ITS RANGE?")
print("="*76)
print("   A 7-tournament has C(7,3) = 35 triples; the number of CYCLIC ones ranges")
print("   0 (transitive) .. 14 (regular/doubly-regular).  Achievable values:")
def cyc_from_scores(sc): return 35 - sum(s*(s-1)//2 for s in sc)
ach = set()
for sc in itertools.combinations_with_replacement(range(7), 7):
    if sum(sc) != 21: continue
    # Landau's condition for a valid score sequence
    ss = sorted(sc)
    if all(sum(ss[:k+1]) >= (k+1)*k//2 for k in range(7)):
        ach.add(cyc_from_scores(ss))
print(f"   achievable cyclic-triple counts: {sorted(ach)}  ({len(ach)} values)")
print(f"   observed t values:               {sorted(count_t)}  ({len(count_t)} values)")
print(f"   observed spectra:                {len(count_spec)}")
