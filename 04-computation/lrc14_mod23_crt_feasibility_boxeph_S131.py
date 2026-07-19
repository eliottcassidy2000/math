#!/usr/bin/env python3
"""
lrc14_mod23_crt_feasibility_boxeph_S131.py  (HYP-7895)

OWNER-DIRECTED (2026-07-19): run the mod-23 CRT feasibility count on the two
depth-minimal rung targets (HYP-7880 lead):

  TARGET A: primitive 12-set C with M(C) = 3/38  (n=12 uniqueness-gap rung)
  TARGET B: primitive 13-set W with M(W) = 4/55  ((1/14,3/41) live-window mediant)

THE STACK (all PROVED necessary conditions; general-p spread lemma = S115/S126):
  * band: at the maximizer t* = m/q0 (gcd(m,q0)=1, forced), every residue
    d = c*m mod q0 lies in [D, q0-D]; two straddlers at exactly +-D (Pinch);
    straddler speeds are units mod q0.
  * covering: duties 7,8,9,10,11,12 (A: M<1/12) / 7..13 (B: M<1/13).
  * spread-or-blocked at p (2/p > M): A: {13,17,19,23,25}; B: {17,19,23,25}.
  * couplings: 38=2*19 (parity<->19: pair {+-1} even-only, {+-2} odd-only, zero
    odd-19-mult-only); 55=5*11 (duty-10/11 pinned to the 1/11- and 1/5-grids);
    25 shares mod-5 with 55.

PART 1: exact counts of feasible residue systems (inclusion-exclusion with
symmetry collapse; 'general' variant + depth-minimal 'k1' variant per v),
cross-checked by Monte Carlo against an independent direct evaluator.
PART 2: exhaustive CRT-guided cap sweeps; every candidate gets an EXACT M
decision (all q <= 2*cap suffice, Pinch), i.e. "no realizer with Vmax <= CAP".
Cross-check: A at cap 26 must be empty (kps-S12).

Pure Python, exact integers. boxeph-2026-07-19-S131.
"""

from math import comb, gcd
from fractions import Fraction
import random

# ---------------------------------------------------------------- utilities

def upairs(p):
    return [j for j in range(1, (p + 1) // 2) if gcd(j, p) == 1 and gcd(p - j, p) == 1]

def pair_of(r, p):
    return min(r % p, (-r) % p)

def spread_count(p, k, hit, haszero):
    """# k-tuples over (Z/p)^k completing 'zero present OR all unit pairs hit',
    given fixed slots already hitting `hit` / providing a zero."""
    if haszero:
        return p ** k
    freeP = [j for j in upairs(p) if j not in hit]
    zero_branch = p ** k - (p - 1) ** k
    cov = sum((-1) ** j * comb(len(freeP), j) * (p - 1 - 2 * j) ** k
              for j in range(len(freeP) + 1))
    return zero_branch + cov

def brute_spread(p, k, hit, haszero):
    tot = 0
    for t in range(p ** k):
        rs = []; x = t
        for _ in range(k):
            rs.append(x % p); x //= p
        z = haszero or any(r == 0 for r in rs)
        got = set(hit) | {pair_of(r, p) for r in rs if r and gcd(r, p) == 1}
        if z or all(j in got for j in upairs(p)):
            tot += 1
    return tot

for (p, k, hit) in [(5, 3, set()), (7, 3, {1}), (25, 2, set()), (25, 2, {1, 2})]:
    assert spread_count(p, k, hit, False) == brute_spread(p, k, hit, False), (p, k, hit)

# ============================================================ TARGET A machinery
# slot coords (c8, c9, c25, r19); parity = c8 % 2; band bans (par, r19) in BAN38.
BAN38 = {(0, 0), (1, 1), (0, 2), (0, 17), (1, 18)}     # parity 0 = even
SLOT_A = 8 * 9 * 25 * 19

def R19(par, zok, a1, a2, g):
    if par == 0:
        return 16 - 2 * a1 - 2 * g
    return 17 - (0 if zok else 1) - 2 * a2 - 2 * g

def C25(par, excl10, zok, j):
    base = (1 if zok else 0) + 4 + (20 - 2 * j)
    if excl10 and par == 0:
        base -= (1 if zok else 0) + 4
    return base

_S89 = {}
def S89(Dfroz, par):
    key = (Dfroz, par)
    if key not in _S89:
        n = 0
        for c8 in range(8):
            if c8 % 2 != par: continue
            for c9 in range(9):
                if 8 in Dfroz and c8 == 0: continue
                if 9 in Dfroz and c9 == 0: continue
                if 12 in Dfroz and c8 % 4 == 0 and c9 % 3 == 0: continue
                n += 1
        _S89[key] = n
    return _S89[key]

def blockA_count(nfree, strad_mode, strad_fix=None):
    """Exact count of free-slot assignments (space SLOT_A^nfree, and for 'general'
    also the two straddler slots' free coords) satisfying: band + duties{8,9,10,12}
    + spread19-sat + spread25-sat.  Straddlers: r19 in pair {+-3} always (never
    zero, always hits P3); 'general': c8 odd free, c9/c25 free; 'k1': all fixed
    (strad_fix = {'serves': duties served, 'hit25': pairs, 'zero25': bool})."""
    duties = (8, 9, 10, 12)
    hit25 = set() if strad_mode == 'general' else set(strad_fix['hit25'])
    z25s = False if strad_mode == 'general' else strad_fix['zero25']
    nf25 = len([j for j in upairs(25) if j not in hit25])
    total = 0
    for Dm in range(16):
        D = frozenset(duties[i] for i in range(4) if Dm >> i & 1)
        sgn_D = (-1) ** len(D)
        excl10 = 10 in D
        if strad_mode == 'k1' and (strad_fix['serves'] & D):
            continue                                    # straddler weight 0
        b19 = [(+1, True, 0, 0, 0), (-1, False, 0, 0, 0)]           # Z: any0 - no0
        b19 += [((-1) ** (a1 + a2 + g) * comb(6, g), False, a1, a2, g)
                for a1 in (0, 1) for a2 in (0, 1) for g in range(7)
                if not (a1 == a2 == g == 0 and False)]
        # NB the C-branch j=(0,0,0) term IS needed (I-E base term)
        if z25s:
            b25 = [(+1, True, 0)]
        else:
            b25 = [(+1, True, 0), (-1, False, 0)]
            b25 += [((-1) ** j * comb(nf25, j), False, j) for j in range(nf25 + 1)]
        for (w19, zok19, a1, a2, g) in b19:
            for (w25, zok25, j25) in b25:
                Wfree = sum(S89(D, par) * C25(par, excl10, zok25, j25)
                            * R19(par, zok19, a1, a2, g) for par in (0, 1))
                if strad_mode == 'general':
                    Ws = 4 * (8 if 9 in D else 9) * C25(1, excl10, zok25, j25)
                    term = Ws * Ws * Wfree ** nfree
                else:
                    term = Wfree ** nfree
                total += sgn_D * w19 * w25 * term
    return total

def directA_sat(slots, strads):
    """independent direct evaluator for MC check. slots = list of (c8,c9,c25,r19)
    free slots; strads likewise (r19 fixed 3/16)."""
    alls = slots + strads
    for (c8, c9, c25, r) in slots:
        if (c8 % 2, r) in BAN38: return False
    for q, f in [(8, lambda s: s[0] == 0), (9, lambda s: s[1] == 0),
                 (10, lambda s: s[0] % 2 == 0 and s[2] % 5 == 0),
                 (12, lambda s: s[0] % 4 == 0 and s[1] % 3 == 0)]:
        if not any(f(s) for s in alls): return False
    z19 = any(s[3] == 0 for s in slots)
    hit19 = {pair_of(s[3], 19) for s in alls if s[3] != 0}
    if not (z19 or all(j in hit19 for j in range(1, 10))): return False
    z25 = any(s[2] == 0 for s in alls)
    hit25 = {pair_of(s[2], 25) for s in alls if s[2] and gcd(s[2], 25) == 1}
    if not (z25 or all(j in hit25 for j in upairs(25))): return False
    return True

def mc_check_A(trials=400000, nfree=3):
    rng = random.Random(38)
    good = 0
    for _ in range(trials):
        slots = [(rng.randrange(8), rng.randrange(9), rng.randrange(25),
                  rng.randrange(19)) for _ in range(nfree)]
        strads = [(2 * rng.randrange(4) + 1, rng.randrange(9), rng.randrange(25), 3),
                  (2 * rng.randrange(4) + 1, rng.randrange(9), rng.randrange(25), 16)]
        if any((s[0] % 2, s[3]) in BAN38 for s in slots):
            continue                    # count only band-passing in numerator? No:
        # we compare against machinery density over the FULL space, so evaluate all
    # redo properly: numerator = satisfies EVERYTHING incl band
    good = 0
    rng = random.Random(38)
    for _ in range(trials):
        slots = [(rng.randrange(8), rng.randrange(9), rng.randrange(25),
                  rng.randrange(19)) for _ in range(nfree)]
        strads = [(2 * rng.randrange(4) + 1, rng.randrange(9), rng.randrange(25), 3),
                  (2 * rng.randrange(4) + 1, rng.randrange(9), rng.randrange(25), 16)]
        if directA_sat(slots, strads):
            good += 1
    exact = blockA_count(nfree, 'general')
    space = SLOT_A ** nfree * 900 ** 2
    mu = exact / space
    sd = (mu * (1 - mu) / trials) ** 0.5
    dev = abs(good / trials - mu) / sd if sd else 0
    return mu, good / trials, dev

# ============================================================ TARGET B machinery
BAN55 = {(0, 0), (1, 1), (2, 2), (3, 3), (2, 8), (3, 9), (4, 10)}   # (d5, d11)
SLOT_B = 8 * 9 * 25 * 11
S_PAIRS = sorted({pair_of(e, 25) for e in range(25) if e % 5 == 4})  # straddler classes
T_PAIRS = [j for j in upairs(25) if j not in S_PAIRS]

def blockB_count(nfree, strad_mode, serves89=(frozenset(), frozenset())):
    """B-block: duties {8,9,10,11,12} + band55 + spread25.  Straddler (d5,d11)
    fixed (4,4)/(1,7); e-digit free (5 values, one per S-pair); c8,c9 free
    ('general') or fixed with served duties serves89[i] ('k1')."""
    duties = (8, 9, 10, 11, 12)
    total = 0
    for Dm in range(32):
        D = frozenset(duties[i] for i in range(5) if Dm >> i & 1)
        sgn_D = (-1) ** len(D)
        if strad_mode == 'k1' and ((serves89[0] & D) or (serves89[1] & D)):
            continue
        b25 = [(+1, True, 0, 0), (-1, False, 0, 0)]
        b25 += [((-1) ** (jS + jT) * comb(5, jS) * comb(5, jT), False, jS, jT)
                for jS in range(6) for jT in range(6)]
        for (w25, zok25, jS, jT) in b25:
            missS, missT = set(S_PAIRS[:jS]), set(T_PAIRS[:jT])
            Wfree = 0
            for par in (0, 1):
                s89 = S89(D & frozenset({8, 9, 12}), par)
                cnt = 0
                for e in range(25):
                    if e == 0:
                        if not zok25: continue
                    elif e % 5 == 0:
                        pass
                    else:
                        pj = pair_of(e, 25)
                        if pj in missS or pj in missT: continue
                    if 10 in D and par == 0 and e % 5 == 0: continue
                    for d11 in range(11):
                        if (e % 5, d11) in BAN55: continue
                        if 11 in D and d11 == 0: continue
                        cnt += 1
                Wfree += s89 * cnt
            ecnt = 5 - jS                      # each missed S-pair kills one digit
            if strad_mode == 'general':
                w89 = S89(D & frozenset({8, 9, 12}), 0) + S89(D & frozenset({8, 9, 12}), 1)
                Ws1 = Ws2 = w89 * ecnt
            else:
                Ws1 = Ws2 = ecnt
            total += sgn_D * w25 * Ws1 * Ws2 * Wfree ** nfree
    return total

def directB_sat(slots, strads):
    """slots = (c8, c9, e25, d11); strads e-digit encoded as full e with e%5 = 4/1,
    d11 = 4/7; c8,c9 free."""
    alls = slots + strads
    for (c8, c9, e, d11) in slots:
        if (e % 5, d11) in BAN55: return False
    for q, f in [(8, lambda s: s[0] == 0), (9, lambda s: s[1] == 0),
                 (10, lambda s: s[0] % 2 == 0 and s[2] % 5 == 0),
                 (11, lambda s: s[3] == 0),
                 (12, lambda s: s[0] % 4 == 0 and s[1] % 3 == 0)]:
        if not any(f(s) for s in alls): return False
    z25 = any(s[2] == 0 for s in alls)
    hit = {pair_of(s[2], 25) for s in alls if s[2] and gcd(s[2], 25) == 1}
    if not (z25 or all(j in hit for j in upairs(25))): return False
    return True

def mc_check_B(trials=400000, nfree=3):
    rng = random.Random(55)
    good = 0
    for _ in range(trials):
        slots = [(rng.randrange(8), rng.randrange(9), rng.randrange(25),
                  rng.randrange(11)) for _ in range(nfree)]
        strads = [(rng.randrange(8), rng.randrange(9), 5 * rng.randrange(5) + 4, 4),
                  (rng.randrange(8), rng.randrange(9), 5 * rng.randrange(5) + 1, 7)]
        if directB_sat(slots, strads):
            good += 1
    exact = blockB_count(nfree, 'general')
    space = SLOT_B ** nfree * (8 * 9 * 5) ** 2
    mu = exact / space
    sd = (mu * (1 - mu) / trials) ** 0.5
    dev = abs(good / trials - mu) / sd if sd else 0
    return mu, good / trials, dev

# ============================================================ direct stack checker

def stack_check(fam, q0, D, duties, spreads):
    band_ms = []
    for m in range(1, q0):
        if gcd(m, q0) != 1: continue
        ds = [(c * m) % q0 for c in fam]
        if all(min(d, q0 - d) >= D for d in ds) and any(min(d, q0 - d) == D for d in ds):
            band_ms.append(m)
    duty_ok = {q: any(c % q == 0 for c in fam) for q in duties}
    sp = {}
    for p in spreads:
        blocked = any(c % p == 0 for c in fam) if p != 25 else any(c % 25 == 0 for c in fam)
        hit = {pair_of(c, p) for c in fam if c % p and gcd(c % p, p) == 1}
        sp[p] = 'blk' if blocked else ('sprd' if all(j in hit for j in upairs(p)) else 'FAIL')
    return band_ms, duty_ok, sp

# ============================================================ PART 2: cap sweep

class LeafLimit(Exception):
    pass

def exact_M(fam, cap2):
    bn, bd = 0, 1
    for q in range(2, cap2 + 1):
        for b in range(1, q // 2 + 1):
            md = q
            for c in fam:
                r = (c * b) % q
                d = min(r, q - r)
                if d < md:
                    md = d
                    if md * bd <= bn * q: break
            if md * bd > bn * q:
                bn, bd = md, q
    g = gcd(bn, bd)
    return bn // g, bd // g

def sweep(q0, D, n, duties, cap, hot, spreads, leaf_limit=9_000_000):
    """Enumerate band+duty+straddle families <= cap; prune on requirement/spread
    reachability (all state packed into two ints); leaf: spread filter, then exact
    beat scan (q <= 2cap decides M by Pinch)."""
    band = {r for r in range(q0) if min(r, q0 - r) >= D}
    qorder = hot + [x for x in range(2, 2 * cap + 1) if x not in hot]
    # requirement bits: one per duty + straddle+/-
    reqbit = {q: 1 << k for k, q in enumerate(duties)}
    BP, BM = 1 << len(duties), 2 << len(duties)
    REQ = (4 << len(duties)) - 1
    # spread packing: per prime an offset; pair bits then a zero bit
    off, width = {}, 0
    for p in sorted(spreads):
        off[p] = width; width += len(upairs(p)) + 1
    FULLSEG = {p: (1 << len(upairs(p))) - 1 for p in spreads}
    SATALL = 0
    survivors, leaves, crtpass, tested = [], 0, 0, 0
    spr = sorted(spreads)
    try:
        for m in [m for m in range(1, q0 // 2 + 1) if gcd(m, q0) == 1]:
            A = [c for c in range(1, cap + 1) if (c * m) % q0 in band]
            if len(A) < n: continue
            if any(all(c % q for c in A) for q in duties): continue
            nA = len(A)
            # per-speed contribution ints
            met = [0] * nA          # requirement bits this speed provides
            sbit = [0] * nA         # spread bits (pair bit or zero bit per prime)
            for i, c in enumerate(A):
                mm = 0
                for q in duties:
                    if c % q == 0: mm |= reqbit[q]
                r = (c * m) % q0
                if r == D: mm |= BP
                if r == q0 - D: mm |= BM
                met[i] = mm
                sb = 0
                for p in spr:
                    if c % p == 0:
                        sb |= 1 << (off[p] + len(upairs(p)))
                    elif gcd(c % p, p) == 1:
                        idx = upairs(p).index(pair_of(c, p))
                        sb |= 1 << (off[p] + idx)
                sbit[i] = sb
            sufmet = [0] * (nA + 1)
            sufsb = [0] * (nA + 1)
            for i in range(nA - 1, -1, -1):
                sufmet[i] = sufmet[i + 1] | met[i]
                sufsb[i] = sufsb[i + 1] | sbit[i]

            def spread_ok(sb):
                for p in spr:
                    seg = sb >> off[p]
                    L = len(upairs(p))
                    if not ((seg >> L) & 1 or (seg & FULLSEG[p]) == FULLSEG[p]):
                        return False
                return True

            stack = []

            def rec(start, got, sb):
                nonlocal leaves, crtpass, tested
                slots = n - len(stack)
                if slots == 0:
                    if REQ & ~got: return
                    leaves += 1
                    if leaves > leaf_limit: raise LeafLimit
                    fam = stack
                    g0 = 0
                    for c in fam: g0 = gcd(g0, c)
                    if g0 != 1: return
                    if not spread_ok(sb): return
                    crtpass += 1
                    tested += 1
                    for q in qorder:
                        beat = False
                        for b in range(1, q // 2 + 1):
                            md = q
                            for c in fam:
                                rr = (c * b) % q
                                d = rr if rr <= q - rr else q - rr
                                if d < md:
                                    md = d
                                    if md * q0 <= D * q: break
                            if md * q0 > D * q:
                                beat = True; break
                        if beat: break
                    else:
                        nu, de = exact_M(list(fam), 2 * cap)
                        if (nu, de) == (D, q0):
                            survivors.append(tuple(fam))
                    return
                if (REQ & ~(got | sufmet[start])): return
                if not spread_ok(sb | sufsb[start]): return
                if nA - start < slots: return
                for i in range(start, nA):
                    stack.append(A[i])
                    rec(i + 1, got | met[i], sb | sbit[i])
                    stack.pop()

            rec(0, 0, 0)
        return sorted(set(survivors)), leaves, crtpass, tested, True
    except LeafLimit:
        return sorted(set(survivors)), leaves, crtpass, tested, False

# ================================================================ main

def main():
    import sys
    part = sys.argv[1] if len(sys.argv) > 1 else 'all'
    print("=" * 100)
    print("MOD-23 CRT FEASIBILITY COUNT -- 3/38 (n=12) and 4/55 (n=13)   [boxeph-S131, HYP-7895]")
    print("=" * 100)

    stradA = [v for v in range(3, 19, 2) if gcd(v, 38) == 1]
    stradB = [v for v in range(4, 28) if gcd(v, 55) == 1]
    print("\nA straddle pairs (v,38-v), v odd unit <19:  %s  (%d)" % (stradA, len(stradA)))
    print("B straddle pairs (v,55-v), v unit <27.5:    %s  (%d)" % (stradB, len(stradB)))

    if part in ('all', 'part1'):
        run_part1()
    if part in ('all', 'part2', 'part2a'):
        run_part2a()
    if part in ('all', 'part2', 'part2b'):
        run_part2b()


def run_part1():
    stradA = [v for v in range(3, 19, 2) if gcd(v, 38) == 1]
    stradB = [v for v in range(4, 28) if gcd(v, 55) == 1]
    print("\n--- structural asserts + MC cross-checks of the exact machinery ---")
    assert blockA_count(0, 'k1', {'serves': set(), 'hit25': set(), 'zero25': False}) == 0
    muA, mcA, devA = mc_check_A()
    print("A-block MC: exact=%.6f  mc=%.6f  dev=%.1f sigma  %s"
          % (muA, mcA, devA, "OK" if devA < 4 else "*** MISMATCH ***"))
    muB, mcB, devB = mc_check_B()
    print("B-block MC: exact=%.6f  mc=%.6f  dev=%.1f sigma  %s"
          % (muB, mcB, devB, "OK" if devB < 4 else "*** MISMATCH ***"))
    assert devA < 4 and devB < 4, "machinery/MC mismatch"

    print("\n--- direct stack check on reference families ---")
    BF38 = [3, 5, 7, 8, 9, 10, 11, 12, 13, 15, 21, 35]
    AP12 = list(range(1, 13))
    for name, fam in [("BF38 (S125 band family)", BF38), ("AP12 {1..12}", AP12)]:
        bms, dok, sp = stack_check(fam, 38, 3, [7, 8, 9, 10, 11, 12], [13, 17, 19, 23, 25])
        print("  %-24s band-m:%-10s duties(7..12):%s  spreads:%s"
              % (name, bms if bms else "NONE", ''.join('1' if dok[q] else '0'
                 for q in [7, 8, 9, 10, 11, 12]), sp))

    # ---------------- PART 1 A
    print("\n--- PART 1, TARGET A = 3/38 ---")
    NA = blockA_count(10, 'general')
    spaceA = SLOT_A ** 10 * 900 ** 2
    outer, ospace = 1, 1
    for p in (13, 17, 23):
        outer *= spread_count(p, 12, set(), False); ospace *= p ** 12
    for q in (7, 11):
        outer *= q ** 12 - (q - 1) ** 12; ospace *= q ** 12
    densA = Fraction(NA, spaceA) * Fraction(outer, ospace)
    print("A general: block=%d  blockdens=%.4e  outerdens=%.4e  TOTAL=%.4e  FEASIBLE=%s"
          % (NA, NA / spaceA, outer / ospace, float(densA), NA * outer > 0))
    print("  certificate ledger (12 slots): band=%.3e sp13=%.3f sp17=%.3f sp23=%.4f sp25=%.4f"
          % ((33 / 38) ** 12,
             spread_count(13, 12, set(), False) / 13 ** 12,
             spread_count(17, 12, set(), False) / 17 ** 12,
             spread_count(23, 12, set(), False) / 23 ** 12,
             spread_count(25, 12, set(), False) / 25 ** 12))
    c23 = sum((-1) ** j * comb(11, j) * (22 - 2 * j) ** 12 for j in range(12)) / 22 ** 12
    c25c = sum((-1) ** j * comb(10, j) * (24 - 2 * j) ** 12 for j in range(11)) / 24 ** 12
    print("  P(all 11 pairs mod23 | no 23-mult, 12 slots) = %.5f (1 in %.0f)" % (c23, 1 / c23))
    print("  P(all 10 pairs mod25 | no 25-mult, 12 slots) = %.5f (1 in %.0f)" % (c25c, 1 / c25c))
    print("  k=1 refinement (straddlers exactly v, 38-v):")
    for v in stradA:
        w = 38 - v
        serves = {9} if (v % 9 == 0 or w % 9 == 0) else set()
        hit25 = {pair_of(x, 25) for x in (v, w) if x % 5}
        z25 = (v % 25 == 0 or w % 25 == 0)
        Nk = blockA_count(10, 'k1', {'serves': serves, 'hit25': hit25, 'zero25': z25})
        ok, osp = 1, 1
        n23eff = None
        for p in (13, 17, 23):
            hz = (v % p == 0 or w % p == 0)
            hp = {pair_of(x, p) for x in (v, w) if x % p}
            cnt = spread_count(p, 10, hp, hz)
            ok *= cnt; osp *= p ** 10
            if p == 23: n23eff = cnt / 23 ** 10
        for q in (7, 11):
            served = (v % q == 0 or w % q == 0)
            ok *= q ** 10 if served else q ** 10 - (q - 1) ** 10
            osp *= q ** 10
        dtot = Fraction(Nk, SLOT_A ** 10) * Fraction(ok, osp)
        print("    v=%2d: block=%.3e  outer=%.3e  TOTAL=%.3e  N23dens=%.4f  strad5mult=%s"
              % (v, Nk / SLOT_A ** 10, ok / osp, float(dtot), n23eff,
                 'Y' if (v % 5 == 0 or w % 5 == 0) else 'n'))

    # ---------------- PART 1 B
    print("\n--- PART 1, TARGET B = 4/55 ---")
    NB = blockB_count(11, 'general')
    spaceB = SLOT_B ** 11 * (8 * 9 * 5) ** 2
    outerB, ospaceB = 1, 1
    for p in (17, 19, 23):
        outerB *= spread_count(p, 13, set(), False); ospaceB *= p ** 13
    for q in (7, 13):
        outerB *= q ** 13 - (q - 1) ** 13; ospaceB *= q ** 13
    densB = Fraction(NB, spaceB) * Fraction(outerB, ospaceB)
    print("B general: block=%d  blockdens=%.4e  outerdens=%.4e  TOTAL=%.4e  FEASIBLE=%s"
          % (NB, NB / spaceB, outerB / ospaceB, float(densB), NB * outerB > 0))
    c23b = sum((-1) ** j * comb(11, j) * (22 - 2 * j) ** 13 for j in range(12)) / 22 ** 13
    print("  P(all 11 pairs mod23 | no 23-mult, 13 slots) = %.5f (1 in %.0f)" % (c23b, 1 / c23b))
    print("  v-collapse: same mod-23 pair at v=%s; same mod-19 pair at v=%s"
          % ([v for v in stradB if (2 * v - 55) % 23 == 0],
             [v for v in stradB if (2 * v - 55) % 19 == 0]))
    print("  k=1 refinement (straddlers exactly v, 55-v):")
    for v in stradB:
        w = 55 - v
        s89 = frozenset(q for q in (8, 9, 12) if v % q == 0), \
              frozenset(q for q in (8, 9, 12) if w % q == 0)
        Nk = blockB_count(11, 'k1', s89)
        ok, osp = 1, 1
        for p in (17, 19, 23):
            hz = (v % p == 0 or w % p == 0)
            hp = {pair_of(x, p) for x in (v, w) if x % p}
            ok *= spread_count(p, 11, hp, hz); osp *= p ** 11
        for q in (7, 13):
            served = (v % q == 0 or w % q == 0)
            ok *= q ** 11 if served else q ** 11 - (q - 1) ** 11
            osp *= q ** 11
        dtot = Fraction(Nk, SLOT_B ** 11) * Fraction(ok, osp)
        print("    v=%2d: block=%.3e  outer=%.3e  TOTAL=%.3e  serves89=%s"
              % (v, Nk / SLOT_B ** 11, ok / osp, float(dtot),
                 sorted(s89[0] | s89[1]) or '-'))


def run_part2a():
    print("\n--- PART 2a: TARGET A sweeps (exact M per candidate) ---")
    spA = [13, 17, 19, 23, 25]; spB = [17, 19, 23, 25]
    hotA = [24, 29, 21, 23, 25, 22, 27, 26, 20]
    sv, lv, cp, ts, done = sweep(38, 3, 12, [7, 8, 9, 10, 11, 12], 26, hotA, spA)
    print("A cap 26 (kps-S12 cross-check): leaves=%d crtpass=%d beatscanned=%d "
          "survivors=%d complete=%s %s"
          % (lv, cp, ts, len(sv), done,
             "== kps-S12 (gap empty in [1,26]) OK" if not sv and done else sv))
    for cap in (34,):
        sv, lv, cp, ts, done = sweep(38, 3, 12, [7, 8, 9, 10, 11, 12], cap, hotA, spA)
        print("A cap %d: leaves=%d crtpass=%d survivors=%d complete=%s %s"
              % (cap, lv, cp, len(sv), done, sv or ""))
        if not done: break

def run_part2b():
    print("\n--- PART 2b: TARGET B sweeps (exact M per candidate) ---")
    spB = [17, 19, 23, 25]
    hotB = [42, 41, 29, 27, 24, 38, 30, 26]
    for cap in (30, 32):
        sv, lv, cp, ts, done = sweep(55, 4, 13, [7, 8, 9, 10, 11, 12, 13], cap, hotB, spB)
        print("B cap %d: leaves=%d crtpass=%d survivors=%d complete=%s %s"
              % (cap, lv, cp, len(sv), done, sv or ""))
        if not done: break

if __name__ == "__main__":
    main()
