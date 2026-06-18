"""
LRC(14) — RIGOROUS CLOSURE of the COVERING AP family  {1,2,...,12, m}.
angle "scale-invariance-AP-and-k3"  (kind-pasteur-S4-wf)

THEOREM (AP family). Let S = {1,2,...,12, m} with m>=14 (so |S|=13).
  (i)  Covering  <=>  182 | m   (because {1,...,12} has no multiple of 13 and no
       multiple of 14, so m must supply both; LCM(13,14)=182).
  (ii) For EVERY covering S (m=182k, k>=1),  M(S) >= 2/27 > 1/14.
  Hence NO set of the form {1,...,12,m} is an LRC(14) counterexample.

By M-scale-invariance M(cS)=M(S), this immediately closes every set
  {t, 2t, ..., 12t, V}  with  t | V  (it equals c=t times {1,...,12, V/t}).

PROOF = two explicit closed-form GLOBAL WITNESSES tau (each verified by a
finite residue computation, no search):

  WITNESS 1:  tau = 2/27.
     - small part {1,...,12} achieves level exactly 2/27  (12-value check).
     - runner m has ||2m/27|| >= 2/27  unless  m mod 27 in {0,13,14}
       (complete 27-residue table).
     => closes all m with m mod 27 not in {0,13,14}.

  WITNESS 2:  tau = (m/13)/(m+1)   (defined since 13 | m for covering m).
     Write m = 182k, q = m/13 = 14k, tau = 14k/(182k+1).
     - For v=1..12:  v*q = 14vk <= 168k < m+1, so NO wraparound; hence
       ||v*tau|| = min(vq, (m+1)-vq)/(m+1) >= 14k/(182k+1).
     - For runner m:  13q = m ≡ -1 (mod m+1)  =>  mq ≡ -q (mod m+1),
       so ||m*tau|| = q/(m+1) = 14k/(182k+1).
     => level = 14k/(182k+1), which is increasing in k and >= 28/365 (k>=2).
     The covering m with m mod 27 in {0,13,14} are exactly k mod 27 in {0,2,25}
     (since 182 ≡ 20 mod 27); the smallest such k>=1 is k=2, giving 28/365.

Both witness levels exceed 1/14:  28/365 - 1/14 = 27/5110 > 0,  2/27 - 1/14 = 1/378 > 0.
QED.
"""
from fractions import Fraction as F
from math import gcd

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def Mval(S):
    # exact M via complete candidate set (optimum is always at a candidate)
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2))
    return max(min(nrm(x*t) for x in S) for t in C)

def is_cov(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

base = list(range(1, 13))

print("="*68)
print("STEP 0: covering  <=>  182 | m")
print("="*68)
# verify on a range
covm = [m for m in range(14, 4000) if is_cov(base+[m])]
print("  covering m in [14,4000):", covm[:8], "...")
print("  all multiples of 182?:", all(m % 182 == 0 for m in covm),
      " and contains all 182k?:", set(covm) == set(range(182, 4000, 182)))

print()
print("="*68)
print("STEP 1: WITNESS 1  tau = 2/27")
print("="*68)
lv_small = min(nrm(v*F(2,27)) for v in base)
print("  level of {1..12} at tau=2/27:", lv_small, "=", float(lv_small))
unsafe = [r for r in range(27) if nrm(F(2*r,27)) < F(2,27)]
print("  runner-m unsafe residues (m mod 27):", unsafe)
assert lv_small == F(2,27) and unsafe == [0,13,14]
print("  => tau=2/27 is a global witness (level>=2/27) for all m with m mod 27 not in {0,13,14}.")

print()
print("="*68)
print("STEP 2: WITNESS 2  tau = (m/13)/(m+1)  for covering m=182k")
print("="*68)
# verify level formula and the >1/14 bound on the failing residues
def witness2_level(m):
    q = m // 13
    t = F(q, m+1)
    return min(nrm(v*t) for v in base+[m]), t
worst = F(1)
for k in range(1, 600):
    if (182*k) % 27 not in (0,13,14):  # only the failing residues need witness 2
        continue
    m = 182*k
    lv, t = witness2_level(m)
    # closed form predicts 14k/(182k+1)
    assert lv == F(14*k, 182*k+1), (k, lv, F(14*k,182*k+1))
    if lv < worst: worst = lv
print("  closed-form level 14k/(182k+1) VERIFIED for all failing k<600.")
print("  worst (smallest) witness-2 level:", worst, "=", float(worst), "at k=2 (m=364)")
assert worst == F(28,365)

print()
print("="*68)
print("STEP 3: combine — every covering m closed, level > 1/14")
print("="*68)
mn = F(1); arg = None
for k in range(1, 600):
    m = 182*k
    if (182*k) % 27 not in (0,13,14):
        lv = min(nrm(v*F(2,27)) for v in base+[m])   # witness 1
    else:
        lv, _ = witness2_level(m)                    # witness 2
    if lv < mn: mn = lv; arg = m
    assert lv >= F(2,27), (m, lv)   # uniform floor is 2/27 (witness 1's level)
print("  min witness level over all covering m=182k (k<600):", mn, "=", float(mn), "at m=", arg)
print("  (witness 1 floor 2/27 dominates; witness 2 only used on failing residues, gives >=28/365)")
print("  28/365 - 1/14 =", F(28,365)-F(1,14), "> 0 :", F(28,365) > F(1,14))
print("  2/27  - 1/14 =", F(2,27)-F(1,14), "> 0 :", F(2,27) > F(1,14))

print()
print("="*68)
print("STEP 4: cross-check against exact Mval on covering m (independent)")
print("="*68)
bad = []
for k in range(1, 60):
    m = 182*k
    Mv = Mval(base+[m])
    if Mv < F(1,14): bad.append((m, Mv))
print("  exact Mval on covering m=182k, k<60: all M>1/14?", not bad, " violations:", bad)
print("  (witness gives a LOWER bound; true M is >= witness level, here ~1/13)")

print()
print("RESULT: covering AP family {1,...,12,m} fully closed. M >= 2/27 > 1/14. QED.")
print("Scale-invariance => {t,2t,...,12t,V} with t|V also closed.")
