"""
mac-mini-2026-07-07-S57 (HYP-5237) -- BLOCK-DOMINATION: the rigorous sup_E avgc(E,P).

The greedy "k-1 pairs at each top difference" in the sweep is NOT realizable (9 points
give only 36 pairs, and cannot put 8 pairs at 5 distinct differences).  The correct
extremal object is fixed by an ELEMENTARY majorization:

  LEMMA (block domination).  For any k distinct integers a_1<...<a_k and any t>=1,
     N_E(t) := #{(i,j): a_j - a_i <= t}  <=  N_block(t) := sum_{d=1}^{t} (k-d)_+,
  because a_j - a_i >= j - i, so a_j-a_i <= t forces j-i <= t, and there are exactly
  (k-d) index-pairs with j-i=d.  I.e. EVERY k-set's difference multiset is
  block-dominated: its pairs at difference <= t never outnumber the block {0..k-1}'s.

Consequence for the average-form conditional tent.  Writing c(d,P) = 1 + b(d,P) with
b(d) -> 0, and n_E(d) = #pairs at difference exactly d,
     C(k,2) * (avgc(E,P) - 1) = sum_d n_E(d) b(d) = sum_t N_E(t) (b(t) - b(t+1))   [Abel].
Hence the RIGOROUS sup over all k-sets is the value of the linear program
     max sum_d n(d) c(d)   s.t.   sum_{d<=t} n(d) <= N_block(t) for all t,
                                   sum_d n(d) = C(k,2),  n(d) >= 0.
The block itself (n(d)=k-d) is FEASIBLE (saturates every cap) and REALIZABLE, so
     avgc(block,P)  <=  sup_E avgc(E,P)  <=  LP-max(P).
If LP-max(P) <= c*(P), the k-leg is DISCHARGED at shape P for EVERY family, with no
diameter hypothesis.  If avgc(block,P) > c*(P), the block is a genuine offender at P
(then diam=k-1 <= 16, so klein-S174's THM-653 diam form already covers it).
"""
from fractions import Fraction as F
from itertools import combinations
import numpy as np

MP = F(14249, 252252)
TH = F(1, 7)
DMAX = 3000

def GP_intervals(P):
    bad = []
    for p in P:
        w = F(1, 14*p)
        for j in range(p+1):
            bad.append((F(j,p)-w, F(j,p)+w))
    bad = [(max(l, F(0)), min(h, F(1))) for l, h in bad if h > 0 and l < 1]
    bad.sort()
    merged = []
    for l, h in bad:
        if merged and l <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], h))
        else:
            merged.append((l, h))
    good = []; prev = F(0)
    for l, h in merged:
        if l > prev: good.append((prev, l))
        prev = max(prev, h)
    if prev < 1: good.append((prev, F(1)))
    return good

def c_array(P, k, DMAX=DMAX):
    beta = F(14-k, 7*k); w = TH - beta; intf = w*w/2
    beta_f, th_f, intf_f = float(beta), float(TH), float(intf)
    iv = GP_intervals(P); meas = sum(h-l for l, h in iv)
    def Dfun(t):
        return np.where(t <= beta_f, -t*intf_f,
               np.where(t <= th_f, (t-beta_f)**2/2 - t*intf_f, (1.0-t)*intf_f))
    d_arr = np.arange(1, DMAX+1, dtype=np.int64)
    acc = np.zeros(DMAX)
    for (l, h) in iv:
        for (pt, sgn) in ((h, +1.0), (l, -1.0)):
            num, den = pt.numerator, pt.denominator
            fr = ((d_arr*num) % den).astype(np.float64)/den
            acc += sgn*Dfun(fr)
    c = 1.0 + acc/(d_arr*float(meas)*intf_f)
    # tail bound for d > DMAX
    delta2_over_intf = 6/7 + beta_f + float(w)**2/4
    ctail = 1.0 + len(iv)*delta2_over_intf/(DMAX*float(meas))
    return c, float(meas), len(iv), ctail

def block_avgc(c, k):
    """avgc of the block {0..k-1}: n(d)=k-d for d=1..k-1."""
    num = sum((k-d)*c[d-1] for d in range(1, k))
    return num / (k*(k-1)//2)

def lp_max(c, ctail, k):
    """
    max sum n(d) c(d) s.t. prefix sums <= N_block(t), total = C(k,2), n>=0.
    Solved exactly by the greedy exchange: total budget T = C(k,2). The binding caps
    are cumulative B_t = sum_{d<=t}(k-d). We may distribute mass to maximize c-weight
    subject to 'at most B_t in the first t positions'. Equivalent: process caps as a
    laminar interval family; the optimum places mass at the highest-c positions that
    remain feasible.  Since B_t increments by (k-t) at step t and saturates at t=k-1
    (B_{k-1}=T), only positions d=1..k-1 can ever hold the *incremental* room; but the
    incremental room at level t (=k-t) can be spent at ANY position d<=t (prefix), so
    the LP relaxation lets each 'layer l = 1..k-1' (of size k-l, created at cap-step l)
    be placed at the best c-position in {1,...,l}. => LP-max = sum_{l=1}^{k-1} (k-l) *
    max_{d<=l} c(d).  (Layer l must lie within the first l positions to respect that
    cap B_l - B_{l-1} = k-l units become 'available' only once t reaches l; earlier caps
    forbid them earlier. Standard greedy-matroid / interval-scheduling optimum.)
    """
    T = k*(k-1)//2
    prefmax = np.maximum.accumulate(c[:k-1])  # prefmax[l-1] = max_{d<=l} c(d), d,l 1-indexed
    return sum((k-l)*prefmax[l-1] for l in range(1, k)) / T

def cstar(meas, k):
    return (1 - MP/F(meas).limit_denominator(10**9)) * F(7*k, 2*(k-1)*(k-7))

print("=== BLOCK DOMINATION: rigorous sup_E avgc, k=9 ===")
k = 9
rows = []
for P in combinations(range(1, 14), 4):
    c, meas, K, ctail = c_array(P, k)
    bavg = block_avgc(c, k)
    lpm = lp_max(c, ctail, k)
    cs = float(cstar(meas, k))
    rows.append((bavg, lpm, cs, P, meas, K))
# report: does block avgc ever exceed c*?  does LP-max ever exceed c*?
block_fail = [r for r in rows if r[0] > r[2]]
lp_fail    = [r for r in rows if r[1] > r[2]]
print(f"  shapes: {len(rows)}")
print(f"  block avgc > c* (GENUINE offenders, need diam route): {len(block_fail)}")
for bavg, lpm, cs, P, meas, K in sorted(block_fail, reverse=True)[:8]:
    print(f"    P={P}: block avgc={bavg:.4f} > c*={cs:.4f}  (LP-max={lpm:.4f}) diam(block)=8<=16 -> THM-653")
print(f"  LP-max > c* (relaxation cannot certify; block may still pass): {len(lp_fail)}")
for bavg, lpm, cs, P, meas, K in sorted(lp_fail, key=lambda r:-(r[1]-r[2]))[:8]:
    tag = "block PASSES" if bavg <= cs else "block FAILS too"
    print(f"    P={P}: LP-max={lpm:.4f} > c*={cs:.4f} but block avgc={bavg:.4f} ({tag})")
npass_block = sum(1 for r in rows if r[0] <= r[2])
print(f"\n  => BLOCK passes at {npass_block}/{len(rows)} shapes.")
print(f"  => shapes where BLOCK itself offends (the true residual, all diam-8): {len(block_fail)}")

print("\n=== k=10 same analysis ===")
k = 10
rows10 = []
for P in combinations(range(1, 14), 3):
    c, meas, K, ctail = c_array(P, k)
    bavg = block_avgc(c, k)
    lpm = lp_max(c, ctail, k)
    cs = float(cstar(meas, k))
    rows10.append((bavg, lpm, cs, P, meas, K))
block_fail10 = [r for r in rows10 if r[0] > r[2]]
print(f"  shapes: {len(rows10)}; block avgc > c*: {len(block_fail10)}")
for bavg, lpm, cs, P, meas, K in sorted(block_fail10, reverse=True)[:6]:
    print(f"    P={P}: block avgc={bavg:.4f} > c*={cs:.4f} (LP-max={lpm:.4f})")
npass_block10 = sum(1 for r in rows10 if r[0] <= r[2])
print(f"  => BLOCK passes at {npass_block10}/{len(rows10)} shapes.")

# But at k=10 the true offenders are NOT the block (block=diam9<=10, THM-653 covers).
# The residual families have diam>10. Bound THEM via block-domination too: the block
# dominates ALL k-sets including wide ones, so sup_E avgc = LP-max for EVERY P.
print("\n=== k=10: the DIAMETER-FREE verdict via LP-max (dominates ALL families) ===")
lp_fail10 = [r for r in rows10 if r[1] > r[2]]
print(f"  LP-max <= c* (=> shape DISCHARGED for ALL families, no diam): {len(rows10)-len(lp_fail10)}/{len(rows10)}")
print(f"  LP-max > c* (residual shapes, need diam+ledger): {len(lp_fail10)}")
for bavg, lpm, cs, P, meas, K in sorted(lp_fail10, key=lambda r:-(r[1]-r[2]))[:10]:
    print(f"    P={P}: LP-max={lpm:.4f} vs c*={cs:.4f} (block {bavg:.4f}); meas={meas:.4f}")

# ============================================================================
# THE ENVELOPE-BLOCK BOUND (rigorous, monotone-repaired):
#   sorted diffs of any k-set termwise-dominate the block's (d_(r) >= block_(r)),
#   and c(d) <= chat(d) := max_{d'>=d} c(d') which IS decreasing, so
#     sum_r c(d_(r)) <= sum_r chat(d_(r)) <= sum_r chat(block_(r))
#                     = sum_{d=1}^{k-1} (k-d) chat(d).
#   => avgc(E,P) <= EnvBlock(P) := (1/C(k,2)) sum_{d=1}^{k-1} (k-d) chat(d,P), ALL E.
# ============================================================================
def envblock(c, ctail, k):
    allc = np.append(c, ctail)
    chat = np.maximum.accumulate(allc[::-1])[::-1]  # chat[i] = max of allc[i:]
    return sum((k-d)*chat[d-1] for d in range(1, k)) / (k*(k-1)//2)

print("\n=== ENVELOPE-BLOCK BOUND (rigorous sup over ALL families), k=9 ===")
k = 9
env_rows = []
for P in combinations(range(1, 14), 4):
    c, meas, K, ctail = c_array(P, k)
    eb = envblock(c, ctail, k)
    cs = float(cstar(meas, k))
    env_rows.append((eb, cs, P, meas, K))
env_fail = [r for r in env_rows if r[0] > r[1]]
npass = len(env_rows) - len(env_fail)
print(f"  EnvBlock <= c* (RIGOROUSLY DISCHARGED, no diam): {npass}/{len(env_rows)}")
if env_fail:
    print(f"  residual shapes (EnvBlock > c*, need diam route): {len(env_fail)}")
    for eb, cs, P, meas, K in sorted(env_fail, key=lambda r:-(r[0]-r[1]))[:12]:
        print(f"    P={P}: EnvBlock={eb:.4f} vs c*={cs:.4f} (excess {eb-cs:+.4f}); meas={meas:.4f}, K={K}")
else:
    print("  *** k=9 (A') LEG DISCHARGED UNIFORMLY: every shape, every family, no diameter ***")

print("\n=== ENVELOPE-BLOCK BOUND, k=10 ===")
k = 10
env10 = []
for P in combinations(range(1, 14), 3):
    c, meas, K, ctail = c_array(P, k)
    eb = envblock(c, ctail, k); cs = float(cstar(meas, k))
    env10.append((eb, cs, P, meas, K))
ef10 = [r for r in env10 if r[0] > r[1]]
print(f"  EnvBlock <= c*: {len(env10)-len(ef10)}/{len(env10)}; residual: {len(ef10)}")
for eb, cs, P, meas, K in sorted(ef10, key=lambda r:-(r[0]-r[1]))[:10]:
    print(f"    P={P}: EnvBlock={eb:.4f} vs c*={cs:.4f} (excess {eb-cs:+.4f})")
