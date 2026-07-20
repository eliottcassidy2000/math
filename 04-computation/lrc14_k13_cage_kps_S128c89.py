#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c89 -- HYP-7940: the k=13 CONDITIONAL TIGHTNESS CAGE
(the n=14 analog of opus-S400's k=12 cage HYP-7920), with the multi-terminal-
family complication computed exactly.

CAGE LOGIC (opus-S400's k=12 template, restated for k=13 / n=14):
 (Q) QUANTIZATION: a witness at denominator q certifies M >= ceil(q/14)/q,
     which is STRICTLY above 1/14 whenever 14 does not divide q.  So every
     computational certificate of the S-T sieve lives on a rung above the
     floor; the WEAKEST rung over the sieve's denominator set {l*p} bounds
     the micro-gap the sieve can see.
 (T) TERMINAL FAMILIES: a family is INVISIBLE to the sieve if no sieve
     denominator carries a witness for it.  At k=12 the only terminal class
     was the AP (Prop 4.4).  At k=13 the tight locus already has TWO
     families (AP, GW), and near-floor families can evade by GRID MISS
     (their maximizer denominator is coprime to every sieve denominator and
     the grids are too coarse to capture a witness).  We compute the exact
     capture status of every known ladder family.
 (E) ELIMINATION: conditional on the sieve terminating with terminal classes
     = +-dilates of the T-list mod each caging prime, a family below the
     weakest rung satisfies {+-v_i} == c_p * {+-f_i} mod p for some terminal
     F and unit c_p, at EVERY caging prime p not dividing prod(v) nor the
     needed power sums of F.  Even power sums give the c-free invariants
        R_m(V, F) = S_2(V)^m * S_2m(F) - S_2m(V) * S_2(F)^m == 0 mod p.
     Pigeonhole: >= (usable primes)/T force the SAME F; if their product
     exceeds |R_m| for all m = 2..13, every R_m = 0 exactly, Newton
     determines the multiset {v_i^2}, and primitivity forces V = F.  The
     cage height H0(T) is where the product/T stops exceeding |R_13|.

GATES: (g1) k=12 weakest rung must reproduce opus-S400's 113/1466 = 1/13 +
3/19058 over depths {1,2,4,8} x primes [167,733]; (g2) exact M of every
table family must match canon (evaluator = death-star-S59's, gated); (g3)
the 1/14 families must be terminal BY ARITHMETIC (every rung > 1/14).
"""
import sys, time
from math import log, comb, gcd
sys.path.insert(0, '04-computation')
from lrc14_ladder_realization_crossN_kps_S128c86 import M_exact_arg
from fractions import Fraction as F

def primes_in(lo, hi):
    sieve = bytearray([1])*(hi+1)
    sieve[0:2] = b'\x00\x00'
    for i in range(2, int(hi**0.5)+1):
        if sieve[i]:
            sieve[i*i::i] = b'\x00'*len(sieve[i*i::i])
    return [i for i in range(lo, hi+1) if sieve[i]]

def greedy_primes(start, target_ln):
    ps, s = [], 0.0
    for p in primes_in(start, 5000):
        ps.append(p); s += log(p)
        if s > target_ln: break
    return ps

def rung(q, n):
    # weakest certificate value at denominator q, threshold 1/n
    c = -(-q // n)  # ceil
    return F(c, q)

# ---------------- part 1: quantization ----------------
def part1():
    print("== (1) quantization rungs ==")
    # gate g1: k=12
    P12 = primes_in(167, 733)
    best = None
    for p in P12:
        for l in (1, 2, 4, 8):
            r = rung(l*p, 13)
            if best is None or r < best[0]: best = (r, l, p)
    r, l, p = best
    print(f"  g1 k=12 weakest rung over {{1,2,4,8}}x[167,733]: {r} at l={l},p={p}"
          f"  (opus-S400: 113/1466)  -> {'MATCH' if r == F(113,1466) else 'DIFFERS'}")
    print(f"      margin above 1/13: {r - F(1,13)} (opus: 3/19058)")
    # k=13 target
    lnB13 = 13*(12*log(91) - log(13))
    P13 = greedy_primes(191, lnB13 - log(360360))
    print(f"  k=13 caging primes: {len(P13)} primes {P13[0]}..{P13[-1]}")
    for name, DS in (("level-1 {1}", (1,)),
                     ("S-T style {1,2,4,8}", (1, 2, 4, 8)),
                     ("7-augmented {1,2,4,7,8,14,28,56}", (1, 2, 4, 7, 8, 14, 28, 56))):
        best = None
        for p in P13:
            for l in DS:
                r = rung(l*p, 14)
                if best is None or r < best[0]: best = (r, l, p)
        r, l, p = best
        print(f"  k=13 weakest rung, depths {name}: {r} at l={l},p={p}"
              f"  margin above 1/14 = {r - F(1,14)} ~ {float(r - F(1,14)):.3e}")
    return P13

# ---------------- part 2: terminal-family table ----------------
FAMILIES = [
    ("AP13", list(range(1, 14))),
    ("GW",   list(range(1, 12)) + [13, 24]),
    ("F3-3/41", list(range(1, 12)) + [13, 36]),
    ("K2-2/27", list(range(1, 13)) + [26]),
    ("K3-3/40", list(range(1, 13)) + [39]),
    ("F4-4/53", list(range(1, 12)) + [13, 48]),
    ("K4-4/53", list(range(1, 13)) + [52]),
    ("K5-5/66", list(range(1, 13)) + [65]),
    ("D1a-1/13", list(range(1, 13)) + [14]),
    ("D1b-1/13", list(range(1, 13)) + [15]),
    ("DW-14/183", list(range(1, 13)) + [182]),
]

def has_witness(V, q, need):
    # exists a in [1, q//2] with min_v fold(v*a mod q) >= need
    for a in range(1, q//2 + 1):
        ok = True
        for v in V:
            r = (v*a) % q
            if r > q - r: r = q - r
            if r < need: ok = False; break
        if ok: return a
    return None

def part2(P13, DS=(1, 2, 4, 8)):
    print("\n== (2) terminal-family table (sieve = depths "
          f"{DS} x {len(P13)} primes {P13[0]}..{P13[-1]}) ==")
    qs = sorted(l*p for p in P13 for l in DS)
    terminal = []
    for name, V in FAMILIES:
        M, arg = M_exact_arg(V)
        cap_q = None; cap_a = None; ncap = 0
        t0 = time.time()
        for q in qs:
            c = -(-q // 14)
            if F(c, q) > M: continue          # rung above M: no witness possible
            a = has_witness(V, q, c)
            if a is not None:
                ncap += 1
                if cap_q is None: cap_q, cap_a = q, a
        status = "TERMINAL (sieve-invisible)" if cap_q is None else \
                 f"captured: smallest q={cap_q} (a={cap_a}), {ncap} capturing q's"
        if cap_q is None: terminal.append((name, V, M))
        print(f"  {name:>10} M={str(M):>8} arg={arg}: {status}  [{time.time()-t0:.1f}s]")
    print(f"  => T = {len(terminal)} terminal families: {[t[0] for t in terminal]}")
    return terminal

# ---------------- part 3: elimination arithmetic ----------------
def part3(P13, terminal):
    print("\n== (3) multi-family power-sum elimination: the cage height H0(T) ==")
    LNP = sum(log(p) for p in P13)
    print(f"  total caging ln-product = {LNP:.1f} over {len(P13)} primes")
    # exact usable-prime sets: p must not divide S_2(F) nor S_2m(F), m<=13
    for name, V, M in terminal:
        S = {m: sum(v**(2*m) for v in V) for m in range(1, 14)}
        bad = [p for p in P13 if any(S[m] % p == 0 for m in range(1, 14))]
        print(f"  {name}: caging primes dividing some S_2m(F): {len(bad)} {bad[:6]}")
    Ts = sorted(set([1, 2, 3, len(terminal), len(terminal)+2]))
    print(f"  H0(T) = height where product/T stops forcing R_13 = 0 exactly")
    print(f"  (|R_m| <= S_2(V)^m * S_2m(F) + S_2m(V) * S_2(F)^m; V bounded by H)")
    for T in Ts:
        # self-consistent: loss = caging primes possibly dividing prod(v) (<= 13 log_191 H each)
        H = 10.0
        for _ in range(60):
            loss_ln = min(LNP, 13 * (log(H)/log(191)) * log(max(p for p in P13)))
            avail = (LNP - loss_ln) / T
            # binding invariant m=13: ln|R_13| <= 13*ln(13 H^2) + ln(S_26(F)) + ln 2
            # use the worst terminal F for S_26
            worstS = max(sum(v**26 for v in Vv) for _, Vv, _ in terminal) if terminal else sum(v**26 for v in range(1,14))
            lnR = 13*log(13.0) + 26*log(H) + log(float(worstS)) + log(2.0)
            newH = H * ((avail - lnR) and 1.0)  # placeholder, solve directly below
            # solve avail = 13 ln 13 + 26 ln H + ln worstS + ln 2 for H:
            lnH = (avail - 13*log(13.0) - log(float(worstS)) - log(2.0)) / 26
            newH = max(2.718 ** lnH, 1.0)
            if abs(newH - H) / max(H, 1) < 1e-9: H = newH; break
            H = newH
        print(f"  T = {T}: H0 ~ {H:,.0f}  (micro-gap emptiness + rigidity below the "
              f"weakest rung, for max speed <= H0, CONDITIONAL on the sieve run)")

# ---------------- part 4: the 12-even branch ----------------
def part4():
    print("\n== (4) the 12-of-13-even branch (k=12's 11-even analog), spot check ==")
    import random
    random.seed(89)
    worst = None; n = 0
    for _ in range(300):
        W = random.sample(range(1, 60), 12)
        odd = random.choice(range(1, 120, 2))
        V = sorted(set([2*w for w in W] + [odd]))
        if len(V) != 13: continue
        n += 1
        M, arg = M_exact_arg(V)
        if worst is None or M < worst[0]: worst = (M, V, arg)
    M, V, arg = worst
    print(f"  {n} random 12-even+1-odd families: min M = {M} ~ {float(M):.4f} at {V[:5]}...")
    print(f"  vs 1/13 = {float(F(1,13)):.4f}: branch floor {'HOLDS >= 1/13' if M >= F(1,13) else 'BELOW 1/13 -- needs care'}")

if __name__ == "__main__":
    P13 = part1()
    terminal = part2(P13)
    part3(P13, terminal)
    part4()

# ---------------- part 5: the J-separator two-stage cage ----------------
def part5(P13):
    """Repair the T=2 collapse: a c-free invariant separating AP from GW
    congruence classes per prime turns the pigeonhole into a two-stage
    elimination.  Stage 1: one of the two J-identities is forced EXACTLY
    (heights up to e^{(LNP/2 - C_J)/8}).  Stage 2: an exact J-identity
    de-classes the other family's primes except divisors of the fixed
    separator integer, restoring a near-full-product T=1 tower."""
    print("\n== (5) the J-separator two-stage cage ==")
    AP = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
    def S(Vv, m): return sum(v**(2*m) for v in Vv)
    # candidate c-free invariants (degree-balanced ratios), pick separators
    from fractions import Fraction as Fr
    cands = [
        ("S2*S6/S4^2",  lambda Vv: Fr(S(Vv,1)*S(Vv,3), S(Vv,2)**2)),
        ("S2^2*S8/S4^3? deg mismatch skip", None),
        ("S2*S10/S4*S8", lambda Vv: Fr(S(Vv,1)*S(Vv,5), S(Vv,2)*S(Vv,4))),
        ("S4*S12/S8^2",  lambda Vv: Fr(S(Vv,2)*S(Vv,6), S(Vv,4)**2)),
    ]
    LNP = sum(log(p) for p in P13)
    for name, f in cands:
        if f is None: continue
        jA, jG = f(AP), f(GW)
        sep = jA != jG
        print(f"  J = {name}: J(AP) = {jA}, J(GW) = {jG}  -> separator: {sep}")
        if not sep: continue
        # separator integer: cross-multiplied difference
        Dj = abs(jA.numerator*jG.denominator - jG.numerator*jA.denominator)
        bad = [p for p in P13 if Dj % p == 0]
        print(f"     separator integer D_J: {Dj}  (ln = {log(Dj):.1f}); caging primes dividing D_J: {bad if bad else 'NONE'}")
        # stage-1 height: pigeonhole class (>= (LNP-loss)/2 of the ln-product)
        # forces its J-identity exactly while (LNP-loss)/2 > ln|N_F| with
        # N_F = S2(V)S6(V)*d_F - S4(V)^2*n_F, |N_F| <= 169 H^8 (n_F + d_F):
        # C_J = ln(169*(n_F+d_F)), worst family binds.
        CJ = max(log(169.0*(jA.numerator + jA.denominator)),
                 log(169.0*(jG.numerator + jG.denominator)))
        # closed-form: (LNP - c*x)/2 - CJ = 8x, c = 13*ln(859)/ln(191)
        c = 13*log(859)/log(191)
        x1 = (LNP/2 - CJ)/(8 + c/2)
        H1 = 2.718 ** max(x1, 0.0)
        print(f"     stage-1 forcing height: H1 ~ {H1:,.3e}  (below this, one J-identity is EXACT; C_J={CJ:.1f}, loss-slope {c:.1f}/lnH)")
        # stage-2: winning class keeps all primes except divisors of D_J (and its own S-divisors)
        # -> T=1 tower with LNP - (those few) product
        lost = sum(log(p) for p in bad)
        H = 10.0
        AP26 = sum(v**26 for v in AP); GW26 = sum(v**26 for v in GW)
        worstS = max(AP26, GW26)
        for _ in range(60):
            loss_ln = min(LNP, 13*(log(H)/log(191))*log(859)) + lost
            avail = LNP - loss_ln
            lnH = (avail - 13*log(13.0) - log(float(worstS)) - log(2.0))/26
            newH = 2.718 ** lnH
            if abs(newH - H)/max(H, 1) < 1e-9: H = newH; break
            H = newH
        print(f"     stage-2 tower height: H0 ~ {H:,.0f}  (near-T=1 strength restored)")
        print(f"     => TWO-STAGE CAGE HEIGHT = min(stage1, stage2) ~ {min(H1, H):,.0f}")
        break  # first separator suffices

if __name__ == "__main__" and "--part5" in sys.argv:
    P13 = greedy_primes(191, 13*(12*log(91) - log(13)) - log(360360))
    part5(P13)
