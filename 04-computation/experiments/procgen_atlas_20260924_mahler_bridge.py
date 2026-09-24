#!/usr/bin/env python3
"""
procgen_atlas_20260924_mahler_bridge.py

Implication-atlas lane (session collatz-procgen-20260922, 2026-09-24).
The two-place bridge between Collatz no-divergence on constant-weight block classes and
Mahler-type fractional-part problems for the rational ratio 3^a/2^L ("Theorem M" of the note),
and the swap-principle trichotomy SP(rho).

Notation.  T(x) = x/2 (x even), (3x+1)/2 (x odd).  A block B is a binary word of length L
(parity letters, B[0] first) with a ones.  On the cylinder of B (one class mod 2^L),
T^L(x) = (3^a x + R_B)/2^L, R_() = 0, R_{z e} = m_e R_z + e 2^{|z|} (m_0 = 1, m_1 = 3).
For A a set of >= 2 blocks of the same length L and weight a with p = 3^a > q = 2^L:
    alpha = p/q,  rho = q/p,  I_A = [min R_B, max R_B]/(p - q).

Sections:
  B1  identity checks: T^L on cylinders; cylinder of B = {x : q | p x + R_B}.
  B2  census of instances (L <= 12): minimal-spread pairs, |I_A|, alpha*|I_A|, FLP comparison.
  B3  the bridge on finite prefixes, exact rationals: xi_n alpha^j = x_j + t_j^{(n)}, t_j in I_A.
  B4  decoupled Mahler-type problem Z(alpha, I_A): exact survivor search over integer parts <= X.
  B5  swap principle SP(rho): the boundary rho = 2/3 is FALSE (Mahler's map, M-orbit of 1).
  B6  two-place orthogonality for Mahler's map: every positive integer has a 2-adic M-word; the
      real tails f_n = (1/3) sum_j r_{n+j} (2/3)^j are the fractional parts {xi (3/2)^n}.
  B7  independent check of the swap identity (Bernstein's formula vs the block formula, mod 2^N).
  B8  Theorem M': for adjacent pairs (R_B' = R_B + 1) and 3^a > 2^(L+1) the Mahler-type statement is
      EQUIVALENT to no-divergence on {B,B'}^N; census for L <= 20 and a finite-horizon equality check.
"""
from fractions import Fraction as Fr
from itertools import combinations
import math, random

def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)

def T(x):
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2

def RB(B):
    R = 0
    for i, e in enumerate(B):
        R = (3 if e else 1) * R + (e << i)
    return R

def blocks(L, a):
    out = []
    for ones in combinations(range(L), a):
        B = [0] * L
        for i in ones:
            B[i] = 1
        out.append(tuple(B))
    return out

def word(B):
    return "".join(map(str, B))

# ---------------------------------------------------------------- B1
def B1():
    hdr("B1  block identities (exact)")
    random.seed(20260924)
    checked = 0
    for L in range(1, 11):
        for a in range(0, L + 1):
            p, q = 3 ** a, 2 ** L
            for B in blocks(L, a):
                R = RB(B)
                # cylinder class: x = -R p^{-1} mod q
                c = (-R * pow(p, -1, q)) % q
                for _ in range(3):
                    x = c + q * random.randrange(0, 10 ** 12)
                    y = x
                    wd = []
                    for _ in range(L):
                        wd.append(y & 1)
                        y = T(y)
                    assert tuple(wd) == B, (L, a, B, x)
                    assert y * q == p * x + R
                    checked += 1
    print(f"  T^L(x) = (3^a x + R_B)/2^L and 'x in cylinder(B) <=> 2^L | 3^a x + R_B': {checked} checks PASS")
    print("  (for a constant-weight alphabet the divisibility alone forces the parity letters, so an")
    print("   integer sequence x_{j+1} = (p x_j + R_{c_j})/q with c_j in A IS a T-orbit with word c.)")

# ---------------------------------------------------------------- B2
def B2():
    hdr("B2  census of constant-weight instances: 1 <= a <= L-1, 3^a > 2^L, L <= 12")
    print("  columns: L a  alpha=3^a/2^L  best pair (min spread)  R-values  |I_A|  alpha*|I_A|  "
          "FLP needs |I|<1/p?  status")
    best = []
    for L in range(2, 13):
        for a in range(1, L):
            p, q = 3 ** a, 2 ** L
            if p <= q:
                continue
            Bs = blocks(L, a)
            Rs = sorted((RB(B), B) for B in Bs)
            # minimal spread among pairs = minimal gap between consecutive R-values
            gaps = [(Rs[i + 1][0] - Rs[i][0], Rs[i], Rs[i + 1]) for i in range(len(Rs) - 1)]
            dR, (r1, b1), (r2, b2) = min(gaps)
            I_len = Fr(dR, p - q)
            aI = Fr(p, q) * I_len
            flp = I_len < Fr(1, p)
            if I_len >= 1:
                st = "vacuous (|I|>=1)"
            elif aI >= 1:
                st = "non-vacuous, heuristically FALSE (alpha|I|>=1)"
            else:
                st = "non-vacuous, heuristically TRUE (alpha|I|<1)"
            best.append((L, a, float(aI)))
            print(f"  {L:2d} {a:2d}  {p}/{q}={p/q:9.4f}  {word(b1)},{word(b2)}  R={r1},{r2}  "
                  f"|I|={float(I_len):.5f}  a|I|={float(aI):.4f}  FLP:{flp}  {st}")
    print("  In every row |I_A| >= 1/(p-q) > 1/p, so the Flatto-Lagarias-Pollington spread bound 1/p")
    print("  never excludes an instance (PROVED: dR >= 1).")
    # the Y3 pair
    L, a = 10, 9
    p, q = 3 ** 9, 2 ** 10
    B = tuple([1] * 9 + [0])
    Bp = tuple([1] * 8 + [0, 1])
    print(f"  Y3 blocks: {word(B)} (R={RB(B)}), {word(Bp)} (R={RB(Bp)}), dR={RB(Bp)-RB(B)}, "
          f"|I|={RB(Bp)-RB(B)}/{p-q}={ (RB(Bp)-RB(B))/(p-q):.6f}, alpha|I|={(p/q)*(RB(Bp)-RB(B))/(p-q):.4f}")
    # single-zero blocks: R(z) - R(z+1) = 3^{L-2-z} 2^z (zero at position z vs z+1)
    for L in range(3, 12):
        a = L - 1
        Rz = [RB(tuple(1 if i != z else 0 for i in range(L))) for z in range(L)]
        for z in range(L - 1):
            assert Rz[z] - Rz[z + 1] == 3 ** (L - 2 - z) * 2 ** z
    print("  single-zero blocks: R(zero at z) - R(zero at z+1) = 3^(L-2-z) 2^z (checked L<=11); the minimum 2^(L-2)")
    print("  is attained by the Y3-type pair 1^(L-2)01, 1^(L-1)0.")

# ---------------------------------------------------------------- B3
def B3():
    hdr("B3  the bridge on finite prefixes (exact rationals)")
    inst = [(4, 3, (1, 1, 1, 0), (1, 1, 0, 1)), (5, 4, (1, 1, 1, 1, 0), (1, 1, 1, 0, 1)),
            (10, 9, tuple([1] * 9 + [0]), tuple([1] * 8 + [0, 1]))]
    for (L, a, B, Bp) in inst:
        p, q = 3 ** a, 2 ** L
        RBv, RBpv = RB(B), RB(Bp)
        lo, hi = Fr(min(RBv, RBpv), p - q), Fr(max(RBv, RBpv), p - q)
        alpha, rho = Fr(p, q), Fr(q, p)
        # integers x <= X whose T-word starts with n blocks from {B,Bp}: build classes mod q^n
        # recursively: classes C_n = {x mod q^n}.
        cls = [0]
        n_max = 6 if L <= 5 else 3
        counts = []
        for n in range(1, n_max + 1):
            new = []
            mod_prev = q ** (n - 1)
            for c in cls:
                for t in range(q):
                    x = c + mod_prev * t
                    # check block n-1 of x's word is in {B,Bp}
                    y = x
                    for _ in range((n - 1) * L):
                        y = T(y)
                    wd = []
                    for _ in range(L):
                        wd.append(y & 1); y = T(y)
                    if tuple(wd) in (B, Bp):
                        new.append(x)
            cls = new
            counts.append(len(cls))
        # verify the identity for the smallest positive representative of each class (and a shifted one)
        ok = 0
        for c in cls[:200]:
            for x0 in (c if c > 0 else q ** n_max, c + 7 * q ** n_max):
                xs = [x0]; cs = []
                y = x0
                for j in range(n_max):
                    wd = []
                    for _ in range(L):
                        wd.append(y & 1); y = T(y)
                    cs.append(RBv if tuple(wd) == B else RBpv)
                    xs.append(y)
                xi_n = x0 + Fr(1, p) * sum(cs[i] * rho ** i for i in range(n_max))
                for j in range(n_max):
                    t_j = Fr(1, p) * sum(cs[i] * rho ** (i - j) for i in range(j, n_max))
                    assert xi_n * alpha ** j == xs[j] + t_j
                    # the finite tail lies in [R_min, R_max]/p * (1 - rho^{n-j})/(1-rho) subset I_A
                    assert lo * (1 - rho ** (n_max - j)) <= t_j <= hi
                ok += 1
        print(f"  L={L} a={a} A={{{word(B)},{word(Bp)}}}: classes mod 2^(L n) with n A-blocks: "
              f"{counts} (ratio 2/2^L per block);  identity xi alpha^j = x_j + t_j checked on {ok} integers")

# ---------------------------------------------------------------- B4
def survivors(alpha, lo, hi, X, depth):
    """For each integer x0 in [1,X], the set of xi in [x0+lo, x0+hi] with xi alpha^j in Z+[lo,hi]
    for j=0..depth.  Returns the list of death depths (depth+1 = survived)."""
    deaths = []
    for x0 in range(1, X + 1):
        ivs = [(x0 + lo, x0 + hi)]
        d = 0
        aj = Fr(1)
        for j in range(1, depth + 1):
            aj *= alpha
            new = []
            for (u, v) in ivs:
                U, V = u * aj, v * aj
                n_lo = math.floor(U - hi)
                n_hi = math.floor(V - lo) + 1
                for n in range(n_lo, n_hi + 1):
                    a2 = max(U, n + lo); b2 = min(V, n + hi)
                    if a2 <= b2:
                        new.append((a2 / aj, b2 / aj))
            ivs = new
            if not ivs:
                break
            d = j
        deaths.append(d if ivs == [] else depth + 1)
    return deaths

def B4():
    hdr("B4  decoupled Mahler-type problem Z(alpha, I_A): exact survivor search")
    print("  Z(alpha,I): no xi > 0 with xi alpha^j in Z + I for all j >= 0.  Theorem M: Z(alpha,I_A) => no positive")
    print("  integer has a T-parity vector eventually in A^N.  Survivors below are xi with integer part x0 <= X.")
    inst = [(4, 3, (1, 1, 1, 0), (1, 1, 0, 1), 20000, 60),
            (5, 4, (1, 1, 1, 1, 0), (1, 1, 1, 0, 1), 20000, 60),
            (7, 5, (1, 0, 1, 1, 0, 1, 1), (0, 1, 1, 1, 1, 1, 0), 20000, 40),
            (10, 9, tuple([1] * 9 + [0]), tuple([1] * 8 + [0, 1]), 20000, 30)]
    for (L, a, B, Bp, X, depth) in inst:
        p, q = 3 ** a, 2 ** L
        lo, hi = Fr(min(RB(B), RB(Bp)), p - q), Fr(max(RB(B), RB(Bp)), p - q)
        alpha = Fr(p, q)
        dd = survivors(alpha, lo, hi, X, depth)
        hist = {}
        for d in dd:
            hist[d] = hist.get(d, 0) + 1
        mx = max(dd)
        print(f"  L={L} a={a} alpha={p}/{q} I=[{lo},{hi}] (|I|={float(hi-lo):.4f}, alpha|I|={float(alpha*(hi-lo)):.4f}):")
        print(f"     x0 <= {X}: death-depth histogram {dict(sorted(hist.items()))}; max depth reached {mx}"
              + ("  (SURVIVORS to the horizon!)" if mx > depth else ""))

# ---------------------------------------------------------------- B5
def B5():
    hdr("B5  swap principle SP(rho): Theta_S(rho) = sum_{s in S} rho^s in Q_2 is irrational for every")
    print("      non-eventually-periodic S.  rho = 2^L/3^a: PROVED for 3^a < 2^L (Lemma L); implied by PC")
    print("      for 2^L < 3^a, a < L (open, HARD); FALSE at a = L (rho = 2/3) by Mahler's map:")
    N = 4000
    # Mahler map M(g) = ceil(3g/2); odd positions D of the M-orbit of 1
    g = 1
    D = []
    for n in range(N):
        if g & 1:
            D.append(n)
        g = (3 * g + 1) // 2 if g & 1 else 3 * g // 2
    mod = 2 ** N
    inv3 = pow(3, -1, mod)
    s = 0
    for d in D:
        s = (s + pow(2, d, mod) * pow(inv3, d, mod)) % mod
    print(f"  D = odd positions of the M-orbit of 1 (first 20: {D[:20]}), non-periodic (the orbit is increasing).")
    print(f"  Theta_D(2/3) = sum_(d in D) (2/3)^d == -3 (mod 2^{N}): {(s + 3) % mod == 0}")
    assert (s + 3) % mod == 0
    print("  So SP(2/3) is FALSE: a non-periodic 0/1 series in base 2/3 has the rational 2-adic value -3.")
    print("  (Every positive integer g gives such a witness: Theta_{D_g}(2/3) = -3g, THM-2228's carry identity.)")

# ---------------------------------------------------------------- B6
def B6():
    hdr("B6  two-place orthogonality for Mahler's map (real tails of the same series)")
    H = 400
    rows = []
    for g0 in range(1, 13):
        g = g0; r = []
        gs = []
        for n in range(H):
            gs.append(g)
            r.append(g & 1)
            g = (3 * g + 1) // 2 if g & 1 else 3 * g // 2
        # xi = lim g_n (2/3)^n ; f_n = (1/3) sum_j r_{n+j} (2/3)^j  (truncated at H)
        xi = Fr(gs[H - 1]) * Fr(2, 3) ** (H - 1)
        bad = 0
        maxdev = 0.0
        for n in range(40):
            f = sum(Fr(r[n + j], 3) * Fr(2, 3) ** j for j in range(H - n - 1))
            frac = xi * Fr(3, 2) ** n - gs[n]
            maxdev = max(maxdev, abs(float(frac - f)))
            if f >= Fr(1, 2):
                bad += 1
        rows.append((g0, float(xi), bad, maxdev))
    for g0, xi, bad, md in rows:
        print(f"  g0={g0:3d}: xi={xi:.10f}  #n<40 with f_n >= 1/2: {bad:2d}  max|{{xi(3/2)^n}} - f_n| = {md:.1e}")
    print("  Every positive integer has a 2-adic M-word (no 2-adic constraint); the Z-number condition is the")
    print("  real tail condition f_n < 1/2 for all n.  For T-swap orbits it is the opposite: the real tails lie in")
    print("  I_A automatically and the whole difficulty is 2-adic integrality (B3).")

def B7():
    hdr("B7  independent check of the swap identity behind Proposition SP (Bernstein formula, mod 2^N)")
    import random
    random.seed(7)
    N = 3000
    mod = 2 ** N
    ok = 0
    for (L, a) in ((3, 1), (4, 3), (5, 2), (10, 9), (7, 5)):
        Bs = blocks(L, a)
        for trial in range(3):
            B, Bp = random.sample(Bs, 2)
            nblk = N // L + 2
            S = [random.random() < 0.3 for _ in range(nblk)]
            word = []
            for j in range(nblk):
                word.extend(Bp if S[j] else B)
            # Bernstein: Phi(w) = - sum_l 2^{d_l} 3^{-l}, d_l = position of the l-th one (l >= 1)
            inv3 = pow(3, -1, mod)
            phi = 0
            l = 0
            for pos, e in enumerate(word):
                if pos >= N:
                    break
                if e:
                    l += 1
                    phi = (phi - pow(2, pos, mod) * pow(inv3, l, mod)) % mod
            # swap identity: -(1/3^a)[R_B/(1-rho) + (R_B' - R_B) Theta_S(rho)], rho = 2^L/3^a
            p3 = pow(3, a, mod); ip3 = pow(p3, -1, mod)
            rho = (pow(2, L, mod) * ip3) % mod
            one_minus = (1 - rho) % mod
            theta = 0
            for j in range(nblk):
                if S[j] and j * L < N:
                    theta = (theta + pow(rho, j, mod)) % mod
            rhs = (-ip3 * (RB(B) * pow(one_minus, -1, mod) + (RB(Bp) - RB(B)) * theta)) % mod
            assert phi == rhs, (L, a)
            ok += 1
    print(f"  Phi_T(Y_S) == -(1/3^a)[R_B/(1-rho) + (R_B'-R_B) Theta_S(rho)]  (mod 2^{N}) on {ok} random (L,a,B,B',S): PASS")

def B8():
    hdr("B8  Theorem M': exact equivalence for adjacent pairs (R_B' = R_B + 1) with alpha = 3^a/2^L > 2")
    print("  If R_B' = R_B + 1 and p > 2q, then for I = [R_B, R_B+1]/(p-q) the integer e_j = q x_{j+1} - p x_j = p t_j - q t_{j+1}")
    print("  lies in [R_B - q/(p-q), R_B + 1 + q/(p-q)], whose only integers are R_B, R_B + 1: the real confinement")
    print("  xi alpha^j in Z + I FORCES the digits, hence the Collatz blocks.  So Z(alpha, I) <=> T1 on {B,B'}^N.")
    from itertools import combinations
    rows = []
    for L in range(2, 21):
        for a in range(1, L):
            p, q = 3 ** a, 2 ** L
            if p <= 2 * q:
                continue
            byR = {}
            for ones in combinations(range(L), a):
                B = [0] * L
                for i in ones:
                    B[i] = 1
                byR.setdefault(RB(tuple(B)), []).append(tuple(B))
            Rs = sorted(byR)
            pairs = [(r, byR[r][0], byR[r + 1][0]) for r in Rs if r + 1 in byR]
            if pairs:
                r, B, Bp = pairs[0]
                rows.append((L, a, p, q, len(pairs), r, B, Bp))
    for (L, a, p, q, npairs, r, B, Bp) in rows:
        lo = r - Fr(q, p - q); hi = r + 1 + Fr(q, p - q)
        ints = [e for e in range(math.floor(lo), math.ceil(hi) + 1) if lo <= e <= hi]
        assert ints == [r, r + 1]
        print(f"  L={L:2d} a={a:2d} alpha={p}/{q}={p/q:.4f}: {npairs} adjacent pair(s), e.g. {word(B)} (R={r}), {word(Bp)} (R={r+1});"
              f" I = [{r}, {r+1}]/{p-q}; |I|*p = {p/(p-q):.4f} (FLP needs < 1)")
    # finite-horizon equality for the (10,7) instance
    L, a = 10, 7
    p, q = 3 ** a, 2 ** L
    D = p - q
    B = (0, 1, 1, 1, 1, 0, 1, 1, 1, 0); Bp = (1, 1, 0, 1, 1, 0, 0, 1, 1, 1)
    R = RB(B)
    assert R == 4726 and RB(Bp) == 4727
    lo, hi = Fr(R, D), Fr(R + 1, D)
    alpha = Fr(p, q)
    X = 2 ** 22
    # S_1 by an exact integer test: xi in [(x0 D + R)/D, (x0 D + R + 1)/D] and xi*alpha in [n + R/D, n + (R+1)/D]
    S1 = []
    for x0 in range(1, X + 1):
        num_lo = (x0 * D + R) * p - (R + 1) * q      # n >= num_lo/(D q)
        num_hi = (x0 * D + R + 1) * p - R * q        # n <= num_hi/(D q)
        if (num_hi // (D * q)) >= -((-num_lo) // (D * q)):
            S1.append(x0)
    def surv(x0, d):
        ivs = [(x0 + lo, x0 + hi)]
        aj = Fr(1)
        for j in range(1, d + 1):
            aj *= alpha
            new = []
            for (u, v) in ivs:
                U, V = u * aj, v * aj
                for n in range(math.floor(U - hi), math.floor(V - lo) + 2):
                    a2 = max(U, n + lo); b2 = min(V, n + hi)
                    if a2 <= b2:
                        new.append((a2 / aj, b2 / aj))
            ivs = new
            if not ivs:
                return False
        return True
    def blocks_ok(x0, d):
        y = x0
        for _ in range(d):
            wd = []
            for _ in range(L):
                wd.append(y & 1); y = T(y)
            if tuple(wd) not in (B, Bp):
                return False
        return True
    C1 = [x0 for x0 in range(1, X + 1) if blocks_ok(x0, 1)]
    print(f"  (10,7), x0 <= 2^22: d=1: #S_1 (integer test) = {len(S1)}, #C_1 (T-iteration) = {len(C1)}, equal: {S1 == C1}")
    assert S1 == C1
    for d in (2, 3):
        Sd = [x0 for x0 in S1 if surv(x0, d)]
        Cd = [x0 for x0 in C1 if blocks_ok(x0, d)]
        print(f"  (10,7), x0 <= 2^22: d={d}: S_d = {Sd}  C_d = {Cd}  equal: {Sd == Cd}")
        assert Sd == Cd
    print("  So, for these explicit pairs, 'no divergent Collatz orbit with parity vector eventually in {B,B'}^N' is")
    print("  EQUIVALENT to the Mahler-type statement 'no xi > 0 has xi (3^a/2^L)^j in Z + [R_B, R_B+1]/(3^a-2^L) for all j'.")

if __name__ == "__main__":
    B1(); B2(); B3(); B4(); B5(); B6(); B7(); B8()
    print("\nDONE mahler_bridge")
