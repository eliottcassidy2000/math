#!/usr/bin/env python3
"""Audit A, item 2: THM-4601 (ii) / HYP-9240, independent code.

Limit chain for deletion depth D: x* = 1, y* = (2 - 3^D)/3^D.  Both are carried SCALED by P = 3^D
(X = P x, Y = P y are integers; parity of the 2-adic number = parity of the scaled integer since P is odd),
with the scaled Terras map  Z -> Z/2 (Z even),  (3Z + P)/2 (Z odd).   (The session reduces denominators instead.)
Debt k_s = D + #odd(x*) - #odd(y*) over [0, s).  Absorption <=> X == Y and k == 0.
Recorded per D: landing time s0 (first time P | Y), N_D = Y/P, clearing odd steps, entry time se into {1,2}
(or reaching 0), phase at entry, debt at entry and one step later, and the identity k_entry = ceil(se/2) - odd(N_D).

Part R (the finite-j reduction): random residual sources n = 2^K t - 1 with prescribed j = v2(x-1); for D <= min(K-1, 40)
the actual chain state (k_s, e_s), e_s = u_s - 3^k_s v_s (exact Fraction), is compared with the limit chain state
for s <= j+3; the first disagreement time and the actual absorption time are recorded.
"""
import sys, random, math
from fractions import Fraction as Fr

DMAX = 3000   # overridden by argv[1] when run as a script

def limit_chain_scaled(D, cap=10**6):
    P = 3**D
    X, Y = P, 2 - P
    k = D
    s = 0; oddY = 0; s0 = None; N = None; clear_odd = None
    while True:
        if s0 is None and Y % P == 0:
            s0 = s; N = Y // P; clear_odd = oddY
        if Y == 0:
            return dict(D=D, s0=s0, N=N, clear_odd=clear_odd, kind='to0', se=s, k=k)
        if Y == P or Y == 2*P:
            break
        if X == Y and k == 0:                     # cannot happen before entry (X in {P, 2P}); kept as a guard
            return dict(D=D, kind='ABSORBED_EARLY', s=s)
        bx, by = X & 1, Y & 1
        k += bx - by
        oddY += by
        X = (3*X + P) >> 1 if bx else X >> 1
        Y = (3*Y + P) >> 1 if by else Y >> 1
        s += 1
        if s > cap:
            return dict(D=D, kind='CAP')
    se = s
    inphase = (X == Y)
    absorbed = inphase and k == 0
    # one more step and two more steps: debt behaviour in the cycle
    ks = [k]; XX, YY, kk = X, Y, k
    for _ in range(4):
        bx, by = XX & 1, YY & 1
        kk += bx - by
        XX = (3*XX + P) >> 1 if bx else XX >> 1
        YY = (3*YY + P) >> 1 if by else YY >> 1
        ks.append(kk)
    oddN = oddY - clear_odd
    return dict(D=D, s0=s0, N=N, clear_odd=clear_odd, kind='cycle', se=se, inphase=inphase, k=k,
                ks=ks, absorbed=absorbed, oddN=oddN, X_at_entry=X // P, Y_at_entry=Y // P)

def collatz_terras_to_12(N):
    """Terras steps and odd steps of N until it first lies in {1,2} (N >= 1)."""
    s = 0; o = 0
    while N not in (1, 2):
        if N & 1:
            N = (3*N + 1) >> 1; o += 1
        else:
            N >>= 1
        s += 1
    return s, o

def part_B():
    recs = [limit_chain_scaled(D) for D in range(1, DMAX + 1)]
    bad = [r for r in recs if r['kind'] not in ('cycle', 'to0')]
    assert not bad, bad[:3]
    to0 = [r['D'] for r in recs if r['kind'] == 'to0']
    cyc = [r for r in recs if r['kind'] == 'cycle']
    absorbed = [r['D'] for r in cyc if r['absorbed']]
    inph = [r for r in cyc if r['inphase']]
    outph = [r for r in cyc if not r['inphase']]
    # structural checks
    for r in recs:
        assert r['clear_odd'] == r['D'], ("clearing odd steps != D", r['D'])
        assert r['N'] >= 0
        assert r['s0'] >= r['D']
    for r in cyc:
        assert r['N'] >= 1
        assert r['X_at_entry'] in (1, 2) and r['Y_at_entry'] in (1, 2)
        # landing debt = ceil(s0/2)
        # k_entry = ceil(se/2) - oddN (x* odd exactly at even times)
        assert r['k'] == (r['se'] + 1)//2 - r['oddN'], ("entry identity", r['D'])
        sN, oN = collatz_terras_to_12(r['N'])
        assert r['se'] == r['s0'] + sN and r['oddN'] == oN
        if r['inphase']:
            assert all(kk == r['k'] for kk in r['ks']), "in-phase debt constant"
        else:
            assert set(r['ks']) <= {r['k'], r['k'] + 1, r['k'] - 1} and len(set(r['ks'])) == 2, "out-of-phase oscillation"
    kin = [(r['k'], r['D']) for r in inph]
    kout_min = min(min(r['ks']) for r in outph) if outph else None
    print(f"B  limit chains D=1..{DMAX}: absorbed {absorbed}; to0 {to0}; in-phase {len(inph)}; out-of-phase {len(outph)}")
    print(f"   min in-phase debt {min(kin)} (all D with that debt: {[D for k, D in kin if k == min(kin)[0]]}); "
          f"min debt seen by any out-of-phase chain {kout_min}; min debt over all D>=3 (incl. oscillation) "
          f"{min(min(r['ks']) for r in cyc)}")
    print(f"   max landing N_D = {max(r['N'] for r in recs)} (at D = {max(recs, key=lambda r: r['N'])['D']}); "
          f"landing time s0 in [{min(r['s0'] for r in cyc)}, {max(r['s0'] for r in cyc)}]")
    # ratio k/D
    rat = [(r['k']/r['D'], r['D']) for r in inph]
    print(f"   in-phase k/D: min {min(rat)[0]:.3f} (D={min(rat)[1]}), max {max(rat)[0]:.3f} (D={max(rat)[1]})")
    for lo, hi in ((3, 50), (50, 100), (100, 500), (500, 1000), (1000, 2000), (2000, 3001)):
        rr = [x for x, D in rat if lo <= D < hi]
        if rr:
            print(f"     D in [{lo},{hi}): in-phase k/D in [{min(rr):.3f}, {max(rr):.3f}], mean {sum(rr)/len(rr):.3f}")
    viol = [(D, k) for k, D in kin if not (0.75*D <= k <= 1.1*D)]
    print(f"   HYP-9240 'in-phase debts within 0.75D-1.1D': {len(viol)} violations, e.g. {viol[:8]}; "
          f"largest violating D = {max(viol)[0] if viol else None}")
    # entry times and J0(D) = ceil(se/2): end-of-run state independent of J only for J >= J0(D)
    J0 = {r['D']: (r['se'] + 1)//2 for r in cyc}
    print(f"   entry time se(D) for D=3..12: {[(r['D'], r['se']) for r in cyc[:10]]}")
    print(f"   J0(D) = ceil(se/2) for D=3..12: {[(D, J0[D]) for D in range(3, 13)]}; max J0 over D<={DMAX}: {max(J0.values())}")
    # excess of N_D
    ex = []
    for r in cyc:
        sN, oN = collatz_terras_to_12(r['N'])
        ex.append((oN - sN/2, r['N'], r['D']))
    print(f"   odd-step excess of N_D (odd - sigma/2): max {max(ex)[0]} at N={max(ex)[1]}; "
          f"min {min(ex)[0]}; required (landing debt) ceil(s0/2) min over D>=3: {min((r['s0']+1)//2 for r in cyc)}")
    exall = []
    for N in range(1, 881):
        sN, oN = collatz_terras_to_12(N)
        exall.append((oN - sN/2, N))
    print(f"   max excess over all N <= 880: {max(exall)}")
    # first records
    print("   first records (D, s0, N_D, se, inphase, k):",
          [(r['D'], r['s0'], r['N'], r['se'], r.get('inphase'), r['k']) for r in recs[:16]])
    return recs

# ---------------- R: finite-j reduction on actual sources ----------------
def T(x):
    return (3*x + 1) >> 1 if x & 1 else x >> 1

def Tq(q):
    # Terras step on a rational with odd denominator
    if q.numerator & 1:
        return (3*q + 1) / 2
    return q / 2

def make_source(rnd, K, j, extra_bits=300):
    """t odd with x = 2*3^(K-1) t - 1 = 1 + 2^j u, u odd random: 3^(K-1) t = 1 + 2^(j-1) u (mod 2^M)."""
    M = j + extra_bits
    u = rnd.getrandbits(M) | 1
    t = (pow(3, -(K - 1), 1 << M) * (1 + (u << (j - 1)))) % (1 << M)
    assert t & 1
    n = (t << K) - 1
    x = 2*3**(K - 1)*t - 1
    assert (x - 1) % (1 << j) == 0 and ((x - 1) >> j) & 1 == 1
    assert (3**K*t - 1) % 4 == 2                     # first reset letter 2
    return n, t, x

def part_R(nsrc=400):
    rnd = random.Random(424242)
    first_dis = []      # (s_first_disagreement - j)
    absorb_rel = []     # (absorption time - j)
    nchains = 0
    for _ in range(nsrc):
        j = rnd.randint(3, 40)
        K = rnd.randint(2, 45)
        n, t, x = make_source(rnd, K, j)
        for D in range(1, min(K - 1, 40) + 1):
            y = (x + 1)//3**D - 1
            assert (x + 1) % 3**D == 0 and y == 2*3**(K - 1 - D)*t - 1
            # actual chain
            u, v, k = x, y, D
            us, vs = Fr(1), Fr(2 - 3**D, 3**D)
            ks = D
            dis = None; ab = None
            for s in range(j + 250):
                e_act = Fr(u) - Fr(3)**k * v
                if s <= j + 3:
                    e_lim = us - Fr(3)**ks * vs
                    if (k, e_act) != (ks, e_lim) and dis is None:
                        dis = s
                if u == v and k == 0:
                    ab = s; break
                bu, bv = u & 1, v & 1
                k += bu - bv
                u, v = T(u), T(v)
                if s <= j + 3:
                    ks += (1 if us.numerator & 1 else 0) - (1 if vs.numerator & 1 else 0)
                    us, vs = Tq(us), Tq(vs)
            nchains += 1
            if dis is not None:
                first_dis.append(dis - j)
                assert dis > j, ("actual chain disagrees with limit chain at s <= j", K, j, D, dis)
            if ab is not None:
                absorb_rel.append(ab - j)
                assert ab >= j + 1, ("absorbed before time j+1", K, j, D, ab)
    from collections import Counter
    print(f"R  finite-j reduction: {nchains} actual D-chains (400 sources, j in [3,40], K in [2,45], D <= min(K-1,40)):")
    print(f"   state(actual) == state(limit) for all s <= j: yes; first disagreement - j: {sorted(Counter(first_dis).items())[:6]}"
          f" ({len(first_dis)} chains disagree by s = j+3)")
    print(f"   absorbed within j+250: {len(absorb_rel)}; min (absorption time - j) = {min(absorb_rel) if absorb_rel else None}; "
          f"distribution of smallest values: {sorted(Counter(absorb_rel).items())[:8]}")

if __name__ == "__main__":
    DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    part_B()
    part_R()
    print("ALL BARRIER CHECKS PASSED")
