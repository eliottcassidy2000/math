#!/usr/bin/env python3
"""Lane G2 (death-star-2026-07-30-coinC2): the 2-adic coupled functional Phi.

Frame: THM-2966 spine normal form / HYP-9061 sec 2e. Depth law
d_m = floor(gamma*m) + D0; row-m cells on anti-diagonal A_m = m + d_m + 1;
0-side cell (m,k): (z,o) = (m+d_m-k, k+1), box binom(d_m,k);
1-side cell (m,k): (z,o) = (k+1, m+d_m-k), same box.
Doubled deficit delta = 2w - binom(d_m,k), w in [0, binom]: |delta| <= box,
delta == box (mod 2).  D_M(p) = (1/2) sum delta p^z q^o and exactly

    D_M = S_M - 1/2 + (p^{M+1}+q^{M+1})/2        (telescoping identity)

Certificate biases (context CONFIRMED):
  p_A = 1285/2181,        b_A - a_A =  896 = 2^7 * 7        (s_A=7, u_A=7)
  p_B = 8847357/11821757, b_B - a_B = 2974400 = 2^6 * 46475 (s_B=6, u_B odd)

N_X := 2 D_M(p_X) b_X^{A_M} = sum delta a^z (b-a)^o b^{A_M-z-o}  (integer).
Coupled functional (Psi := 2*Phi so everything is an integer):

    Psi = x N_A / 2^{s_A} + y N_B / 2^{s_B},   x = 1, y odd,
    y solves  u_A a_A^{A-1} + y u_B a_B^{A-1} == 0  (mod 2^7)   [canonical]

Claims tested (all arithmetic exact int/Fraction):
  T1  (D1) Psi mod 128 (= Phi mod 64) is choice-independent across schemes.
  T2  per-cell shift spectrum: min v_2(shift change in Psi) = 7 exactly;
      binding cells identified; no odd y (mod 256) does better  => K = 6 cap.
  T3  forced value == boundary closed form
      Psi_bnd = sum_X wgt_X (a^{M+1} + (b-a)^{M+1} - b^{M+1}) b^{d_M} / 2^{s_X}
      (= Psi of the all-w=0 scheme); tabulate mod 128; periodicity probe.
  T4  (D3) timing scan: levels with A_M = 2^r - 1, and deep-corner sublevels
      d_M = 2^rho - 1, for gamma in {1/2, 3/5, 2457/6592}, D0 = 0..8.
  T5  wall inequalities: env_X(M) = (a^{M+1} + (b-a)^{M+1}) b^{d_M}
      satisfies env >= 2^{s*A_M} >= 2^{s*(M+1)} >> 2^{s+1}; exact for M<=40.
  T6  weighted sums over M (single bias B): forcing invariance mod 2^{7+j}
      costs |c_{M0}| >= 2^j exactly (triangular system) -- dollar-for-dollar;
      log2(envelope) - forced-modulus is j-independent.
"""

import math
import random
from fractions import Fraction
from math import comb, gcd

MOD7 = 1 << 7  # 128


# ---------------------------------------------------------------- bias data
def mk_bias(a, b):
    assert 0 < a < b and gcd(a, b) == 1 and a % 2 == 1 and b % 2 == 1
    c = b - a
    s = (c & -c).bit_length() - 1
    u = c >> s
    assert u % 2 == 1 and c == u << s
    return dict(a=a, b=b, c=c, s=s, u=u)


BIAS = {"A": mk_bias(1285, 2181), "B": mk_bias(8847357, 11821757)}
assert BIAS["A"]["s"] == 7 and BIAS["A"]["u"] == 7
assert BIAS["B"]["s"] == 6 and BIAS["B"]["u"] == 46475


def v2(n):
    """2-adic valuation; None = +infinity (n = 0)."""
    if n == 0:
        return None
    return (n & -n).bit_length() - 1


def exact_div(n, twopow):
    assert n % twopow == 0
    return n // twopow


# ---------------------------------------------------------------- ledger
def depth(gamma, D0, m):
    return math.floor(gamma * m) + D0


def A_level(gamma, D0, M):
    return M + depth(gamma, D0, M) + 1


def cells(gamma, D0, M):
    """List of (m, side, k, box, z, o) for rows 1..M, both sides."""
    out = []
    for m in range(1, M + 1):
        dm = depth(gamma, D0, m)
        for k in range(dm + 1):
            B = comb(dm, k)
            out.append((m, 0, k, B, m + dm - k, k + 1))
            out.append((m, 1, k, B, k + 1, m + dm - k))
    return out


def scheme_zero(cs):
    return [-B for (_, _, _, B, _, _) in cs]  # all w = 0 (valid box point)


def scheme_full(cs):
    return [B for (_, _, _, B, _, _) in cs]  # all w = box


def scheme_canonical(cs):
    return [B & 1 for (_, _, _, B, _, _) in cs]  # delta = box mod 2


def scheme_random(cs, rng):
    return [2 * rng.randint(0, B) - B for (_, _, _, B, _, _) in cs]


def N_eval(cs, deltas, X, A):
    a, b, c = BIAS[X]["a"], BIAS[X]["b"], BIAS[X]["c"]
    tot = 0
    for (_, _, _, _, z, o), d in zip(cs, deltas):
        if d:
            tot += d * a**z * c**o * b ** (A - z - o)
    return tot


def weights(A, mod=MOD7):
    """x = 1; y = unique odd solution mod `mod` of the unit congruence."""
    tA = (BIAS["A"]["u"] * pow(BIAS["A"]["a"], A - 1, mod)) % mod
    tB = (BIAS["B"]["u"] * pow(BIAS["B"]["a"], A - 1, mod)) % mod
    y = (-tA * pow(tB, -1, mod)) % mod
    assert y % 2 == 1
    return 1, y


def Psi_of(cs, deltas, gamma, D0, M):
    A = A_level(gamma, D0, M)
    x, y = weights(A)
    NA = N_eval(cs, deltas, "A", A)
    NB = N_eval(cs, deltas, "B", A)
    # per-term v2 >= s*o >= s, so these divisions are exact:
    return x * exact_div(NA, 1 << 7) + y * exact_div(NB, 1 << 6), NA, NB


# ================================================================ T1
def t1_invariance():
    print("=" * 72)
    print("T1: (D1) choice-independence of Psi mod 128  (= Phi mod 64)")
    rng = random.Random(20260730)
    cases = [
        (Fraction(1, 2), 0, 8),
        (Fraction(1, 2), 0, 11),
        (Fraction(1, 2), 0, 14),
        (Fraction(3, 5), 0, 10),
        (Fraction(3, 5), 0, 14),
        (Fraction(1, 2), 4, 7),  # D3 deep-corner timing level (r=4, rho=3)
    ]
    for gamma, D0, M in cases:
        cs = cells(gamma, D0, M)
        schemes = [scheme_zero(cs), scheme_full(cs), scheme_canonical(cs)]
        schemes += [scheme_random(cs, rng) for _ in range(30)]
        psis, minv, sA_ok, sB_ok = [], [], True, True
        for dl in schemes:
            ps, NA, NB = Psi_of(cs, dl, gamma, D0, M)
            psis.append(ps)
            assert ps % 2 == 0, "Psi should always be even (Phi integer)"
            # single-bias forced classes: N_X mod 2^{s_X+1} scheme-independent
            sA_ok &= (NA - N_eval(cs, schemes[0], "A", A_level(gamma, D0, M))) % (1 << 8) == 0
            sB_ok &= (NB - N_eval(cs, schemes[0], "B", A_level(gamma, D0, M))) % (1 << 7) == 0
            # v2(N_X) = s_X exactly iff M even (c_1 = M+1 mod 2)
            if M % 2 == 0:
                assert v2(NA) == 7 and v2(NB) == 6
            else:
                assert (v2(NA) is None or v2(NA) >= 8) and (v2(NB) is None or v2(NB) >= 7)
        base = psis[0]
        diffs = [p - base for p in psis[1:]]
        assert all(d % 128 == 0 for d in diffs), "INVARIANCE mod 128 FAILS"
        nz = [v2(d) for d in diffs if d != 0]
        minv = min(nz) if nz else "inf"
        print(
            f"  gamma={gamma} D0={D0} M={M}: 33 schemes, Psi mod 128 = "
            f"{base % 128:3d} for ALL; min v2(Psi-diff) = {minv} "
            f"(7 = sharp); single-bias mod 2^(s+1) invariance: A={sA_ok} B={sB_ok}"
        )
    print("  T1 PASS: Phi mod 64 is choice-independent; single-bias forced")
    print("  info is exactly N_X mod 2^(s_X+1), and v2(N_X)=s_X iff M even.")


# ================================================================ T2
def shift_changes(gamma, D0, M, x, y):
    """Exact change in Psi for a +2 shift of each cell's delta."""
    A = A_level(gamma, D0, M)
    out = []
    for m, side, k, B, z, o in cells(gamma, D0, M):
        chg = 0
        for X, wgt in (("A", x), ("B", y)):
            a, b, c, s = (BIAS[X][t] for t in "abcs")
            term = a**z * c**o * b ** (A - z - o)
            chg += 2 * wgt * exact_div(term, 1 << s)
        out.append(((m, side, k, z, o), chg))
    return out


def t2_shift_spectrum():
    print("=" * 72)
    print("T2: per-cell shift spectrum and the K = 6 cap")
    gamma, D0, M = Fraction(1, 2), 0, 12
    A = A_level(gamma, D0, M)
    x, y = weights(A)
    sc = shift_changes(gamma, D0, M, x, y)
    by_o = {}
    for (m, side, k, z, o), chg in sc:
        key = "o=1" if o == 1 else ("o=2" if o == 2 else "o>=3")
        by_o.setdefault(key, []).append((v2(chg), (m, side, k, z, o)))
    for key in ("o=1", "o=2", "o>=3"):
        vs = sorted(by_o[key])
        print(f"  {key:5s}: min v2(shift) = {vs[0][0]} at cell (m,side,k,z,o)={vs[0][1]}; "
              f"count={len(vs)}")
    allmin = min(v for v, _ in sum(by_o.values(), []))
    assert allmin == 7, "expected sharp min shift valuation 7 in Psi units"
    # binding o=1 cells: which z-parity class?
    o1 = [(v, cell) for v, cell in by_o["o=1"]]
    for v, (m, side, k, z, o) in sorted(o1)[:4]:
        print(f"    o=1 cell z={z} (A-1-z={A - 1 - z}, parity {(A - 1 - z) % 2}): v2={v}")
    # exhaustive scan over odd y mod 256 (x=1 WLOG: odd unit scaling of (x,y)
    # multiplies every shift change by an odd unit, preserving all v2):
    best = -1
    for yy in range(1, 256, 2):
        mn = min(v2(chg) for _, chg in shift_changes(gamma, D0, M, 1, yy))
        best = max(best, mn)
    print(f"  max over odd y mod 256 of min v2(shift) = {best}  (cap: cannot exceed 7)")
    assert best == 7
    print("  T2 PASS: K = 6 exactly (Psi mod 128).  Binding: side-B o=2 cells")
    print("  (v2 = 1 + 2 s_B - s_B = 7) AND the odd z-parity class at o=1")
    print("  (v2((b_B-a_B)) = 6 exactly => z-dependence mod 2^7).  The side-A")
    print("  z-independence mod 2^7 does NOT lift K to 7: side B binds twice.")


# ================================================================ T3
def Psi_boundary(gamma, D0, M):
    """Closed form: Psi of the all-w=0 scheme (S_M = 0), boundary terms only."""
    A = A_level(gamma, D0, M)
    dM = depth(gamma, D0, M)
    x, y = weights(A)
    tot = 0
    for X, wgt in (("A", x), ("B", y)):
        a, b, c, s = (BIAS[X][t] for t in "abcs")
        num = (a ** (M + 1) + c ** (M + 1) - b ** (M + 1)) * b**dM
        tot += wgt * exact_div(num, 1 << s)
    return tot


def Psi_boundary_mod128(gamma, D0, M):
    """Same value mod 128, via modular arithmetic (fast for large M)."""
    A = A_level(gamma, D0, M)
    dM = depth(gamma, D0, M)
    x, y = weights(A)
    tot = 0
    for X, wgt in (("A", x), ("B", y)):
        a, b, c, s = (BIAS[X][t] for t in "abcs")
        mod = 1 << (s + 7)
        num = (pow(a, M + 1, mod) + pow(c, M + 1, mod) - pow(b, M + 1, mod)) * pow(b, dM, mod) % mod
        assert num % (1 << s) == 0
        tot += wgt * (num >> s)
    return tot % 128


def t3_forced_value():
    print("=" * 72)
    print("T3: forced value of Psi mod 128 = boundary closed form")
    rng = random.Random(4711)
    for gamma, D0 in [(Fraction(1, 2), 0), (Fraction(3, 5), 0), (Fraction(1, 2), 4)]:
        for M in (1, 2, 5, 9, 14):
            cs = cells(gamma, D0, M)
            ps0, _, _ = Psi_of(cs, scheme_zero(cs), gamma, D0, M)
            assert ps0 == Psi_boundary(gamma, D0, M), "closed form must be EXACT"
            psr, _, _ = Psi_of(cs, scheme_random(cs, rng), gamma, D0, M)
            assert (psr - ps0) % 128 == 0
            assert ps0 % 128 == Psi_boundary_mod128(gamma, D0, M)
    print("  closed form == Psi(all-w=0) exactly; random schemes agree mod 128. PASS")
    gammas = [
        ("1/2   ", Fraction(1, 2), 0),
        ("3/5   ", Fraction(3, 5), 0),
        ("cert  ", Fraction(2457, 6592), 0),
        ("1/2+D4", Fraction(1, 2), 4),
    ]
    print("  forced Psi mod 128 (Phi mod 64 = value/2), M = 1..40:")
    for label, g, D0 in gammas:
        row = [Psi_boundary_mod128(g, D0, M) for M in range(1, 41)]
        nz = sum(1 for r in row if r != 0)
        print(f"    {label}: {row}")
        print(f"        nonzero at {nz}/40 levels")
    # ONE-BIT closed form (hand derivation, M >= 2):
    #   a_A == b_A mod 2^7  =>  Sigma_A == (M+1) b_A^M  (mod 2^7);
    #   a_B^{A-1} == b_B^{A-1} (1 - 64 u_B (A-1))  (mod 2^7);
    #   Sigma_B == (M+1) b_B^M + 64*odd*(M(M+1)/2)  (mod 2^7);
    #   unit congruence => u_A b_A^{A-1} + y u_B b_B^{A-1} == 64 (A-1) (mod 128);
    #   Psi == -(M+1)*64*(A-1) + 64*M(M+1)/2
    #       == 64 * [ (M+1) d_M + floor((M+1)/2) ]   (mod 128)
    # (using (M+1)(A-1) == (M+1)(M+d_M) == (M+1) d_M mod 2 since M(M+1) even).
    ok = True
    for g, D0 in [(Fraction(1, 2), 0), (Fraction(3, 5), 0), (Fraction(2457, 6592), 0),
                  (Fraction(1, 2), 4), (Fraction(3, 5), 3), (Fraction(9, 10), 2)]:
        for M in range(2, 3000):
            bit = ((M + 1) * depth(g, D0, M) + (M + 1) // 2) % 2
            if Psi_boundary_mod128(g, D0, M) != 64 * bit:
                ok = False
    print(f"  ONE-BIT law  Psi_forced == 64*[(M+1) d_M + floor((M+1)/2)] mod 128,")
    print(f"  M in [2,3000), six depth laws: {ok}")
    assert ok
    print("  (M=1 is special: the extra term (b_B-a_B)^2 b_B^d/2^6 == 64 flips the bit.)")
    print("  => the coupled invariant carries exactly ONE BIT beyond triviality, and")
    print("     it depends on the depth law ONLY through d_M mod 2: pure boundary-")
    print("     term bookkeeping, completely blind to the band/interior structure.")


# ================================================================ T4
def is_pow2(n):
    return n >= 1 and (n & (n - 1)) == 0


def t4_timing():
    print("=" * 72)
    print("T4: (D3) timing scan: A_M = 2^r - 1 levels; deep corner d_M = 2^rho - 1")
    for gname, g in [("1/2", Fraction(1, 2)), ("3/5", Fraction(3, 5)),
                     ("2457/6592", Fraction(2457, 6592))]:
        for D0 in range(0, 9):
            hits, deep = [], []
            for M in range(1, 8193):
                A = A_level(g, D0, M)
                if is_pow2(A + 1):
                    hits.append(M)
                    dM = depth(g, D0, M)
                    if is_pow2(dM + 1):
                        deep.append((M, (A + 1).bit_length() - 1, (dM + 1).bit_length() - 1))
            deepstr = ""
            if deep:
                vals = [(M, r, rho, Psi_boundary_mod128(g, D0, M)) for (M, r, rho) in deep]
                deepstr = f"; DEEP-CORNER hits (M,r,rho,Psi mod 128): {vals}"
            print(f"  gamma={gname:9s} D0={D0}: {len(hits):3d} all-odd-band levels "
                  f"(first {hits[:4]}){deepstr}")


# ================================================================ T5
def t5_wall():
    print("=" * 72)
    print("T5: the product-formula wall, exact inequalities (M = 1..40)")
    gamma, D0 = Fraction(1, 2), 0
    print("  need for contradiction: forced modulus 2^K > envelope;"
          " available K <= s+1 (+O(1) coupling)")
    for X in ("A", "B"):
        a, b, c, s = (BIAS[X][t] for t in "abcs")
        min_log_gap = None
        for M in range(1, 41):
            dM = depth(gamma, D0, M)
            A = A_level(gamma, D0, M)
            env = (a ** (M + 1) + c ** (M + 1)) * b**dM  # = (p^{M+1}+q^{M+1}) b^{A_M}
            # exact wall chain: env >= c^{M+1} b^{d_M} >= c^{A_M} >= 2^{s A_M}
            assert env > c ** (M + 1) * b**dM - 1 and c ** (M + 1) * b**dM >= c**A >= 1 << (s * A)
            assert env >= 1 << (s * A)          # kills EVERY per-cell forced 2^0..2^{sA}
            assert env >= 1 << (s * (M + 1))    # a fortiori
            assert env > 1 << (s + 1)           # actual forced modulus: hopeless
            gap = env.bit_length() - 1 - (s + 1)
            min_log_gap = gap if min_log_gap is None else min(min_log_gap, gap)
            if M in (1, 10, 20, 40):
                print(f"    bias {X}, M={M:2d}: log2(env) ~ {env.bit_length()-1:5d}, "
                      f"available K <= {s+1};  log2 gap >= {gap}")
        print(f"    bias {X}: min over M<=40 of [log2(env) - (s+1)] = {min_log_gap}")
    print("  wall chain PROVED for every bias p=a/b, every M >= 1:")
    print("    envelope = (a^(M+1)+(b-a)^(M+1)) b^(d_M) >= (b-a)^(A_M) >= 2^(s A_M)")
    print("    >= 2^(s(M+1)) > 2^(s+1) * 2^(sM - 1);   coupled envelopes ADD.")
    print("  T5 PASS: no single-evaluation (or difference: envelopes add) variant")
    print("  can reach 'nonzero forced integer below its envelope'.")


# ================================================================ T6
def t6_weighted_sums():
    print("=" * 72)
    print("T6: weighted sums over M (bias B): dollar-for-dollar")
    g, D0, M0 = Fraction(1, 2), 0, 6
    a, b, c, s = (BIAS["B"][t] for t in "abcs")
    As = {m: A_level(g, D0, m) for m in range(1, M0 + 1)}
    for j in (0, 4, 8, 12):
        # triangular system: T_m = sum_{M>=m} c_M b^{A_M - A_m} == 0 mod 2^j,
        # minimal-|.| back-substitution; T_{M0} = c_{M0}  =>  2^j | c_{M0}.
        cvec = {M0: 1 << j}
        for m in range(M0 - 1, 0, -1):
            resid = sum(cvec[M] * b ** (As[M] - As[m]) for M in cvec) % (1 << j) if j else 0
            cm = (-resid) % (1 << j) if j else 0
            if j and cm > 1 << (j - 1):
                cm -= 1 << j
            cvec[m] = cm
        Ts = {m: sum(cvec[M] * b ** (As[M] - As[m]) for M in cvec if M >= m)
              for m in range(1, M0 + 1)}
        assert all(T % (1 << j) == 0 for T in Ts.values())
        # forced modulus of Psi_c = sum c_M N_M(bias B): min over cells of
        # v2(shift) = 1 + s*o + v2(T_m); verify by direct shift computation
        minshift = None
        for m, side, k, B, z, o in cells(g, D0, M0):
            if Ts[m] == 0:
                continue
            chg = 2 * a**z * c**o * b ** (As[m] - z - o) * Ts[m]
            minshift = v2(chg) if minshift is None else min(minshift, v2(chg))
        K = minshift  # invariance modulus exponent
        env = sum(abs(cvec[M]) * (a ** (M + 1) + c ** (M + 1)) * b ** (As[M] - M - 1)
                  for M in cvec)
        print(f"    j={j:2d}: c = {[cvec[m] for m in range(1, M0+1)]}, forced K = {K} "
              f"(= 1+s+j = {1+s+j}), log2(env/2^K) = {env.bit_length()-1-K}")
        assert K == 1 + s + j
    print("  T6 PASS: modulus gain j costs |c_M0| = 2^j exactly; the envelope/2^K")
    print("  gap is j-independent (b^(A) factor untouched). Weighted sums lose.")


if __name__ == "__main__":
    t1_invariance()
    t2_shift_spectrum()
    t3_forced_value()
    t4_timing()
    t5_wall()
    t6_weighted_sums()
    print("=" * 72)
    print("ALL LANE-G2 CHECKS PASS")
