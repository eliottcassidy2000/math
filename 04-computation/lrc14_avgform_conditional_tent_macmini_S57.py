"""
mac-mini-2026-07-07-S57 (HYP-5237) -- THE AVERAGE-FORM CONDITIONAL TENT:
the k=9 discharge sweep over ALL |P|=4 shapes (and the k=10 verdict).

FRAME (THM-651 + its conditional refinement, kps-S73; my S56 c-table thread;
klein-S174's THM-653 covers k=9 diam<=16, k=10 diam<=10 -- this is the graft
on the CONDITIONAL side klein's NEXT (a) asks for).

THREE INGREDIENTS, each elementary:

(1) AVERAGE FORM (exact, by linearity -- no sup needed):
    E[F 1_{G_P}] = sum_{ordered pairs} int_{G_P} f(frac(dx)) dx
                 = avgc(E,P) * k(k-1) * meas(G_P) * int f,
    where avgc = mean of c(|d|,P) over the C(k,2) unordered pair differences
    (c(-d)=c(d) since G_P is symmetric under x -> 1-x).  Markov with
    min_S F = 1-k*beta then gives
       meas(S cap G_P) <= avgc * (1 - floor_k) * meas(G_P),
    using k(k-1) int f / (1-k beta) = 1 - floor_k exactly (THM-651 optimization).
    DISCHARGE CONDITION (per shape (P,E)):
       avgc(E,P) <= c*(P) := (1 - m_P/meas(G_P)) * 7k / (2(k-1)(k-7)).
    This is strictly weaker than kps's sup-form condition sup_d c <= c*.

(2) CLOSED-FORM c(d,P) VIA THE TENT ANTIDERIVATIVE (no window loop):
    with F(t) = int_0^t f, D(t) = F(t) - t*int f (continuous, D(0)=D(1)=0):
       int_a^b f(frac u) du = (b-a) int f + D(frac b) - D(frac a),
    so, G_P = disjoint union of intervals [l,h],
       c(d,P) = 1 + sum_iv [D(frac(dh)) - D(frac(dl))] / (d * meas * int f).
    frac(d * a/b) = ((d*a) mod b)/b computed in exact integers.

(3) KOKSMA TAIL LEMMA (proved -- "D is bounded"):
       max D = w^2(1-th)/2 = 3w^2/7  (at t = th);
       min D = -(beta w^2/2 + w^4/8) (at t = beta + w^2/2);
       => c(d,P) <= 1 + K_P * Delta2 / (d * meas * int f),
       Delta2/int f = 6/7 + beta + w^2/4,  K_P = #intervals of G_P.
    Covers all d > D0 uniformly.

(4) THE DIFFERENCE-MULTISET CAP: a k-element integer set has at most k-1
    pairs at any fixed difference d (shift-by-d graph is a union of paths).
    => sup_E avgc <= greedy: fill (k-1) pairs at each of the top c-values.
    k=9: 36 pairs = 8+8+8+8+4; k=10: 45 = 9*5.

VERDICT SOUGHT: greedy(P) <= c*(P) for all P => the k-leg is DISCHARGED for
EVERY cluster shape (no diameter condition, no residual).
"""
from fractions import Fraction as F
from itertools import combinations
import numpy as np

MP = F(14249, 252252)   # m_P, THM-530
TH = F(1, 7)
D0 = 2000               # exact-table depth; tail lemma covers d > D0

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

def make_D(beta_f, th_f, intf_f):
    """vectorized D(t) = F(t) - t*intf, piecewise."""
    def D(t):
        out = np.where(t <= beta_f, -t*intf_f,
              np.where(t <= th_f, (t-beta_f)**2/2 - t*intf_f, (1.0-t)*intf_f))
        return out
    return D

def c_values_fast(iv, k, D0=D0):
    """float c(d,P) for d=1..D0 via the antiderivative closed form."""
    beta = F(14-k, 7*k); w = TH - beta; intf = w*w/2
    beta_f, th_f, intf_f = float(beta), float(TH), float(intf)
    Dfun = make_D(beta_f, th_f, intf_f)
    meas = sum(h-l for l, h in iv); meas_f = float(meas)
    d_arr = np.arange(1, D0+1, dtype=np.int64)
    acc = np.zeros(D0)
    for (l, h) in iv:
        for (pt, sgn) in ((h, +1.0), (l, -1.0)):
            num, den = pt.numerator, pt.denominator
            fr = ((d_arr * num) % den).astype(np.float64) / den
            acc += sgn * Dfun(fr)
    c = 1.0 + acc / (d_arr * meas_f * intf_f)
    return c, meas, len(iv), beta, w, intf

def run_level(k, verbose_top=6):
    npairs = k*(k-1)//2
    cap = k-1                          # max pairs at any single difference
    nfull, rem = divmod(npairs, cap)   # greedy fill: nfull full + rem partial
    cstar_mult = F(7*k, 2*(k-1)*(k-7))
    floor_k = 1 - F(2*(k-1)*(k-7), 7*k)
    results = []
    for P in combinations(range(1, 14), 13-k):
        iv = GP_intervals(P)
        c, meas, K, beta, w, intf = c_values_fast(iv, k)
        # tail bound for d > D0 (Koksma lemma)
        delta2_over_intf = 6/7 + float(beta) + float(w)**2/4
        ctail = 1.0 + K * delta2_over_intf / (D0 * float(meas))
        cstar = (1 - MP/meas) * cstar_mult
        allc = np.append(c, ctail)
        top = np.sort(allc)[::-1][:nfull+1]
        greedy = (cap*top[:nfull].sum() + rem*top[nfull]) / npairs
        results.append((float(greedy)/float(cstar), float(greedy), float(cstar),
                        P, float(meas), K, ctail,
                        [(int(np.argmax(c))+1, float(np.max(c)))]))
    results.sort(reverse=True)
    print(f"\n=== k={k}: |P|={13-k} shapes ({len(results)}), floor={floor_k}, "
          f"greedy fill = {cap}x{nfull}+{rem} of {npairs} pairs, D0={D0} ===")
    npass = sum(1 for r in results if r[0] <= 1.0)
    print(f"  PASS (greedy <= c*): {npass}/{len(results)}"
          + ("  => LEVEL DISCHARGED UNIFORMLY (all P, all E, no diam)" if npass == len(results) else ""))
    print(f"  worst {verbose_top} shapes by greedy/c*:")
    for ratio, g, cs, P, meas, K, ctail, topd in results[:verbose_top]:
        print(f"    P={P}: greedy={g:.4f} vs c*={cs:.4f} (ratio {ratio:.4f}); "
              f"meas={meas:.4f}, K={K}, ctail={ctail:.4f}, argmax c: d={topd[0][0]} ({topd[0][1]:.4f})")
    return results

def exact_confirm(P, k, dlist):
    """exact-Fraction c(d,P) at specified d (binding-case confirmation)."""
    beta = F(14-k, 7*k); w = TH - beta; intf = w*w/2
    iv = GP_intervals(P); meas = sum(h-l for l, h in iv)
    def Dex(t):
        if t <= beta: return -t*intf
        if t <= TH:   return (t-beta)**2/2 - t*intf
        return (1-t)*intf
    out = {}
    for d in dlist:
        s = F(0)
        for (l, h) in iv:
            fh = F((d*h.numerator) % h.denominator, h.denominator)
            fl = F((d*l.numerator) % l.denominator, l.denominator)
            s += Dex(fh) - Dex(fl)
        out[d] = 1 + s/(d*meas*intf)
    return out

print("=== AVERAGE-FORM CONDITIONAL TENT: full-shape sweeps ===")
res9 = run_level(9)
res10 = run_level(10)

# exact confirmation at the worst k=9 shape: top-6 d's + c* exact
worst = res9[0]
P9 = worst[3]
iv = GP_intervals(P9)
c, meas, K, beta, w, intf = c_values_fast(iv, 9)
top_d = list(np.argsort(c)[::-1][:6] + 1)
ex = exact_confirm(P9, 9, top_d)
cstar_ex = (1 - MP/meas) * F(63, 32)
print(f"\n=== exact confirmation, worst k=9 shape P={P9} ===")
print(f"  c*(P) exact = {cstar_ex} = {float(cstar_ex):.6f}")
for d in top_d:
    print(f"  c({d}) exact = {ex[d]} = {float(ex[d]):.6f}  (float table: {c[d-1]:.6f})")
gex = (8*sum(sorted(ex.values(), reverse=True)[:4]) + 4*sorted(ex.values(), reverse=True)[4]) / 36
print(f"  greedy from exact top-5: {float(gex):.6f} vs c* {float(cstar_ex):.6f} "
      f"=> {'PASS' if gex <= cstar_ex else 'FAIL'}")

# k=10 residual mapping: true avgc for candidate extremal families vs c*
print("\n=== k=10 residual: true avgc for extremal families (worst shapes) ===")
def true_avgc(E, c_arr, ctail):
    tot = 0.0; n = 0
    E = sorted(E)
    for i in range(len(E)):
        for j in range(i+1, len(E)):
            d = E[j]-E[i]
            tot += c_arr[d-1] if d <= D0 else ctail
            n += 1
    return tot/n
FAMS10 = {
    'block {0..9}':        list(range(10)),
    '2-AP {0,2,..,18}':    list(range(0, 20, 2)),
    'block+1out d11':      list(range(9)) + [11],
    'block+1out d14':      list(range(9)) + [14],
    'block+1out d20':      list(range(9)) + [20],
    'block+1out d30':      list(range(9)) + [30],
    'two-block 5+5 gap12': list(range(5)) + list(range(12, 17)),
    'spread AP*5':         list(range(0, 50, 5)),
}
for ratio, g, cs, P, meas, K, ctail, topd in res10[:3]:
    c, meas_, K_, beta, w, intf = c_values_fast(GP_intervals(P), 10)
    print(f"  P={P} (c* = {cs:.4f}):")
    for name, E in FAMS10.items():
        a = true_avgc(E, c, ctail)
        diam = max(E)-min(E)
        tag = 'PASS' if a <= cs else ('FAIL (diam<=10: THM-653)' if diam <= 10 else 'FAIL *RESIDUAL*')
        print(f"    {name:22s} diam={diam:3d}: avgc = {a:.4f}  {tag}")
