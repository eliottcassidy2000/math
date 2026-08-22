"""ANGLE B2 -- finite certificates for the lemma algebra of the transient
skeleton (amm12592-allR-transient-theorem-boxeph.md). Exact ints only.

T1  coasting: decode_d(2x-1) = ballot vector, admissible, all d; corner -1.
T6  kernel:   decode_{d'}((x^{d-t} q^t)/x) = (1,1) spread at d'=d,
                                             (1,2,1) spread at d'=d+1.
T6  feed:     decode_{d'}(x^j) supported on cells {0,1} for j in {d'-1, d'}.
T7  drain:    one dynamics step maps v x^{d-1} -> (v+2) x^{d'-1}, c0 = -2,
              for even v in [-2(d-1), -2], both d' = d and d' = d+1.
ME  recovery: rule-A closure blocks at R=32,64 satisfy
              sum_{i<=R-2} x^i (Delta_i - (2x-1)) = 2 G_R and Delta_{R-1} = -1.
"""
import sys, os
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from amm12592_allR_family_toolbox_boxeph import (padd, psub, pmul, pshift,
    pscale, qpow, bern_to_poly, poly_to_bern, admissible)
from amm12592_transient_error_dynamics_boxeph import (Em, two_G, run_ruleA,
    exact_profile)

ok_all = True
def check(name, ok):
    global ok_all
    ok_all = ok_all and ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}")

# T1
t1 = True
for d in range(1, 201):
    b = [comb(d-1, t) - (comb(d-1, t-1) if t else 0) for t in range(d+1)]
    t1 &= poly_to_bern([-1, 2], d) == b
    t1 &= admissible(b, d)
    t1 &= bern_to_poly(b, d) == [-1, 2]
    t1 &= bern_to_poly([-comb(d, t) for t in range(d+1)], d) == [-1]
check("T1 ballot decode/admissibility/corner, d = 1..200", t1)

# T6 kernel
t6 = True
for d in range(2, 60):
    for t in range(0, d):          # junk cell t at degree d, t < d (deg >= 1)
        base = pshift(pmul([1], qpow(t)), d - t)      # x^{d-t} q^t
        shifted = base[1:]                            # /x, exact
        # d' = d
        got = poly_to_bern(shifted, d)
        want = [0]*(d+1)
        for s, kcoef in enumerate((1, 1)):
            if t + s <= d: want[t+s] += kcoef
        t6 &= got == want
        # d' = d+1
        got = poly_to_bern(shifted, d+1)
        want = [0]*(d+2)
        for s, kcoef in enumerate((1, 2, 1)):
            if t + s <= d+1: want[t+s] += kcoef
        t6 &= got == want
check("T6 Pascal kernel (1,1)/(1,2,1), d = 2..59, all t < d", t6)

# T6 feed support
tf = True
for dp in range(2, 60):
    for j in (dp - 1, dp):
        cells = poly_to_bern(pshift([1], j), dp)
        tf &= all(v == 0 for t, v in enumerate(cells) if t > dp - j)
check("T6 feed cells subset {0,1} for j in {d'-1, d'}, d' = 2..59", tf)

# T7 drain (hypothesis |v| <= d: at a d-increment row the decode is (v,2v,v)
# and cell 1 with box [-2d, 2] absorbs 2v only when v >= -d; certificate also
# confirms the FAILURE of the naive range: see note Section 8)
t7 = True
for d in range(3, 40):
    for v in range(-2*(d//2), -1, 2):
        e = pshift([v], d - 1)                        # v x^{d-1}
        for dp in (d, d + 1):
            w = poly_to_bern(e[:dp+1], dp)
            c = []
            for t in range(dp + 1):
                lo, hi = -2*comb(dp-1, t), 2*(comb(dp-1, t-1) if t else 0)
                c.append(min(hi, max(lo, w[t])))
            t7 &= c[0] == -2
            rem = psub(e, bern_to_poly(c, dp))
            t7 &= (not rem or rem[0] == 0)
            en = rem[1:] if rem else []
            want = [] if v == -2 else pshift([v + 2], dp - 1)
            t7 &= en == want
check("T7 drain v x^{d-1} -> (v+2) x^{d'-1}, c0 = -2; d = 3..39, even -d<=v<=-2", t7)

# ME recovery
me = True
for R in (32, 64):
    blocks, msg = run_ruleA(R, 0)
    me &= (msg == "CLOSED")
    pr = exact_profile(R, 0)
    S = []
    for i in range(R - 1):
        ci = psub(bern_to_poly(blocks[i], pr[i]), [-1, 2])
        S = padd(S, pshift(ci, i))
    me &= S == two_G(R)
    me &= bern_to_poly(blocks[R-1], pr[R-1]) == [-1]
check("ME master-equation recovery from rule-A closures, R = 32, 64", me)

print("ALL CERTIFICATES PASS" if ok_all else "CERTIFICATE FAILURES PRESENT")
