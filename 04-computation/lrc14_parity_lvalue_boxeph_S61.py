#!/usr/bin/env python3
"""
THE PARITY LAW AND THE L-VALUE WEIGHTS (boxeph-2026-07-17-S61)
Owner directive: prove the mod-7 vanishing; evaluate the W-hat_g character sums.

LEM-032 referee, five parts, on the frame-cross of LEM-031:
  cross(w) = sum_{a != 0} W(a) X(aw),  W(a) = (pi^2/P^2) csc^2(pi a/P),
  X(m) = |sum_e S_e(m)|^2 - sum_e |S_e(m)|^2,   w in (Z/P)*,
  c-hat(chi) = (1/phi(P)) sum_w cross(w) conj(chi(w)).

(A) THE PARITY LAW: W(-a) = W(a) and X(-m) = X(m), so cross(-w) = cross(w),
    hence c-hat(chi) = chi(-1) c-hat(chi): EVERY ODD CHARACTER CARRIES ZERO
    MASS -- and per-factor, W-hat_g(chi) = X-hat_g(chi) = 0 too (double kill).
    DECODE OF S60: the quadratic (Legendre) character mod p is odd iff
    p == 3 (mod 4).  7 == 3 == 3 (mod 4)  => the S60 vanishings at mod 7 and
    mod 3 are PARITY zeros, NOT the seven-section structure; 5 == 1 (mod 4)
    => mod 5 survives.  Sharpened: the EVEN mod-7 characters (the two cubics)
    are predicted to carry mass -- measured here.  (S60's "section symmetry"
    conjecture is CORRECTED by this part.)
(B) THE SUPPORT LAW: the class-g term of c-hat(chi) vanishes unless
    cond(chi) | P/g  (fiber orthogonality of U_P -> U_{P/g}).
(C) THE TWISTED JORDAN LEMMA: for chi a character mod q,
      T_q(chi) := sum_{u in (Z/q)*} chi(u) csc^2(pi u/q)
                =  0                                   (chi odd)
                =  J_2(q)/3                            (chi trivial: THM-892 (C*))
                =  (2 q^2 / pi^2) L_q(2, chi)          (chi even; L_q = the
                       L-series with Euler factors at p | q removed)
                =  2 tau(chi) B_{2, chi-bar}           (chi even PRIMITIVE mod q)
    Proof route: the tent-kernel inversion (THM-892 (K)) turns csc^2 into the
    quadratic Bernoulli moment; the Gauss sum carries the twist; the partial-
    fraction face csc^2(pi x) = (1/pi^2) sum_n (x+n)^{-2} gives the L-value.
    Three-route battery over q in {5,7,9,12,13,15,20,36,63,180}, ALL characters.
(D) THE CLOSED-FORM WEIGHT SIDE: W-hat_g(chi) = (2/g^2) L_{P/g}(2, chi).
    The factorization law becomes
      c-hat(chi) = sum_{g | P, g < P, cond(chi) | P/g}
                     (2/g^2) L_{P/g}(2, chi) * X-hat_g(chi),
    weight side PURE L-VALUES (no csc^2 sums).  Referee: reproduce the S60
    spectrum (Legendre mod 5 = -8.3403 on the balanced cluster) and the new
    cubic mod-7 masses from L-values alone.
(E) THE EVEN-CHARACTER CENSUS (three clusters): all characters of (Z/P)*;
    every odd mass == 0 (mass instances of (A)); Parseval: sum of nontrivial
    masses == frame variance; conjugation symmetry; ranked co-resonant list
    by conductor.

Conventions (fixed by the exactness of the factorization referee):
  W-hat_g(chi) = sum_{u in (Z/q)*} W(g u) chi(u)          (UNnormalized, chi plain)
  X-hat_g(chi) = (1/phi(q)) sum_{u in (Z/q)*} X(g u) conj(chi(u)),   q = P/g.

Pure Python.  Reuses S26's owner_data (S25 endpoint conventions).
"""

import sys
from math import gcd, lcm, pi, sin, sqrt
from fractions import Fraction as Fr
from itertools import product
import cmath

sys.path.insert(0, '04-computation')
from lrc14_general_resonance_law_boxeph_S26 import owner_data


# ---------------------------------------------------------------- utilities

def csc2(x):
    """csc^2(pi x) for x in (0,1)."""
    s = sin(pi * x)
    return 1.0 / (s * s)


def factorize(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def divisors(n):
    ds = [1]
    for p, a in factorize(n).items():
        ds = [x * p ** b for x in ds for b in range(a + 1)]
    return sorted(ds)


def J2(q):
    r = q * q
    for p in factorize(q):
        r = r // (p * p) * (p * p - 1)
    return r


def prim_root(pp, p):
    ph = pp - pp // p
    for g in range(2, pp):
        if gcd(g, pp) != 1:
            continue
        x, k = g % pp, 1
        while x != 1:
            x = x * g % pp
            k += 1
        if k == ph:
            return g
    raise RuntimeError(f"no primitive root mod {pp}")


def hurwitz2(x, N0=200):
    """zeta(2, x) by direct sum + Euler-Maclaurin tail; abs err << 1e-13."""
    s = 0.0
    for n in range(N0):
        t = n + x
        s += 1.0 / (t * t)
    y = N0 + x
    s += 1.0 / y + 0.5 / y ** 2 + 1.0 / (6 * y ** 3) - 1.0 / (30 * y ** 5) \
        + 1.0 / (42 * y ** 7)
    return s


# ---------------------------------------------------------------- characters

class UGroup:
    """(Z/N)* as a product of cyclic components with dlog tables.
    Supports N with 2-part in {1, 2, 4} (all moduli in this project)."""

    def __init__(self, N):
        self.N = N
        self.comps = []                      # (pp, p, gen, order, dlog)
        for p in sorted(factorize(N)):
            a = factorize(N)[p]
            pp = p ** a
            if pp == 2:
                continue                     # U(2) trivial
            assert not (p == 2 and a >= 3), "2-part >= 8 unsupported"
            g = 3 if pp == 4 else prim_root(pp, p)
            order = pp - pp // p
            dlog, x = {}, 1
            for k in range(order):
                dlog[x] = k
                x = x * g % pp
            self.comps.append((pp, p, g, order, dlog))
        self.orders = [c[3] for c in self.comps]
        self.units = [w for w in range(1, max(N, 2)) if gcd(w, N) == 1]
        self.phi = len(self.units)
        self.vec = {w: tuple(c[4][w % c[0]] for c in self.comps)
                    for w in self.units}
        self.roots = [[cmath.exp(2j * pi * k / o) for k in range(o)]
                      for o in self.orders]

    def chars(self):
        if not self.comps:
            return [()]
        return list(product(*[range(o) for o in self.orders]))

    def chi(self, js, w):
        w %= self.N
        if self.N > 1 and gcd(w, self.N) != 1:
            return 0j
        if not self.comps:
            return 1 + 0j
        v = self.vec[w]
        z = 1 + 0j
        for i in range(len(js)):
            if js[i]:
                z *= self.roots[i][(js[i] * v[i]) % self.orders[i]]
        return z

    def parity(self, js):
        if self.N <= 2:
            return 1
        return round(self.chi(js, self.N - 1).real)

    def order_of(self, js):
        o = 1
        for i, j in enumerate(js):
            o = lcm(o, self.orders[i] // gcd(self.orders[i], j) if j else 1)
        return o

    def conductor(self, js):
        cond = 1
        for i, (pp, p, g, order, dlog) in enumerate(self.comps):
            j = js[i]
            if j == 0:
                continue
            b = 1
            while True:
                f = p ** b
                if all((j * dlog[u]) % order == 0
                       for u in dlog if u % f == 1):
                    break
                b += 1
            cond *= p ** b
        return cond

    def chistar(self, js, a, f):
        """the primitive character (conductor f) attached to js, at a."""
        a %= f
        if f > 1 and gcd(a, f) != 1:
            return 0j
        while gcd(a, self.N) != 1 or a == 0:
            a += f
        return self.chi(js, a)


def L2_primitive(G, js, f):
    """L(2, chi*) for the primitive character of conductor f attached to js."""
    tot = 0j
    for a in range(1, f + 1):
        c = G.chistar(js, a, f)
        if c != 0:
            tot += c * hurwitz2(a / f)
    return tot / (f * f)


def euler_restrict(G, js, f, q, Lstar):
    """L_q(2, chi) = L(2, chi*) prod_{p | q, p not| f} (1 - chi*(p)/p^2)."""
    L = Lstar
    for p in factorize(q):
        if f % p != 0:
            L *= (1 - G.chistar(js, p, f) / (p * p))
    return L


def gauss_sum(G, js, f):
    return sum(G.chistar(js, a, f) * cmath.exp(2j * pi * a / f)
               for a in range(1, f + 1))


def B2chi(G, js, f):
    """generalized Bernoulli B_{2,chi*} = f sum_a chi*(a) B_2(a/f)."""
    return f * sum(G.chistar(js, a, f) * ((a / f) ** 2 - a / f + 1 / 6.0)
                   for a in range(1, f + 1))


# ---------------------------------------------------------------- cluster side

def cluster_spectra(E, s):
    P, M, data = owner_data(E, s)
    owners = sorted(data)
    root = [cmath.exp(2j * pi * t / P) for t in range(P)]
    Se = {}
    for e in owners:
        d = data[e]
        arr = [0j] * P
        for sg, q in zip(d["sgn"], d["pos"]):
            for m in range(P):
                arr[m] += sg * root[(m * q) % P]
        Se[e] = arr
    X = [0.0] * P
    for m in range(P):
        tot = 0j
        diag = 0.0
        for e in owners:
            z = Se[e][m]
            tot += z
            diag += z.real * z.real + z.imag * z.imag
        X[m] = (tot.real * tot.real + tot.imag * tot.imag) - diag
    W = [0.0] + [(pi * pi / (P * P)) * csc2(a / P) for a in range(1, P)]
    return P, M, owners, X, W


def cross_all(P, X, W, units):
    out = {}
    for w in units:
        c = 0.0
        aw = 0
        for a in range(1, P):
            aw += w
            if aw >= P:
                aw -= P * (aw // P)
            c += W[a] * X[aw]
        out[w] = c
    return out


def chat_direct(G, js, cross):
    return sum(cross[w] * G.chi(js, w).conjugate() for w in G.units) / G.phi


def chat_factorized(G, js, P, X):
    """sum over classes g of (2/g^2) L_{P/g}(2,chi) X-hat_g(chi); L-values only."""
    f = G.conductor(js)
    Lstar = L2_primitive(G, js, f)
    tot = 0j
    for g in divisors(P):
        if g == P:
            continue
        q = P // g
        if q % f != 0:
            continue
        Uq = [u for u in range(1, q + 1) if gcd(u, q) == 1]
        xh = 0j
        for u in Uq:
            xh += X[(g * u) % P] * G.chistar(js, u, f).conjugate()
        xh /= len(Uq)
        Wg = (2.0 / (g * g)) * euler_restrict(G, js, f, q, Lstar)
        tot += Wg * xh
    return tot


def class_contribution(G, js, g, P, X, W, cross_units):
    """direct (1/phi) sum_w conj(chi(w)) c_g(w), c_g(w) = sum_u W(gu) X(guw)."""
    q = P // g
    Uq = [u for u in range(1, q + 1) if gcd(u, q) == 1]
    tot = 0j
    for w in G.units:
        cg = 0.0
        for u in Uq:
            cg += W[(g * u) % P] * X[(g * u * w) % P]
        tot += cg * G.chi(js, w).conjugate()
    return tot / G.phi


# ---------------------------------------------------------------- PART C

def part_C():
    print("=" * 78)
    print("PART C -- THE TWISTED JORDAN LEMMA: T_q(chi) = sum chi(u) csc^2(pi u/q)")
    print("  faces: odd -> 0; trivial -> J_2(q)/3; even -> (2q^2/pi^2) L_q(2,chi);")
    print("         even primitive -> 2 tau(chi) B_{2,chi-bar}")
    battery = [5, 7, 9, 12, 13, 15, 20, 36, 63, 180]
    worst = {"odd": 0.0, "triv": 0.0, "L": 0.0, "GB": 0.0}
    counts = {"odd": 0, "triv": 0, "L": 0, "GB": 0}
    for q in battery:
        G = UGroup(q)
        rows = []
        for js in G.chars():
            par = G.parity(js)
            f = G.conductor(js)
            T = sum(G.chi(js, u) * csc2(u / q) for u in G.units)
            if par == -1:
                worst["odd"] = max(worst["odd"], abs(T))
                counts["odd"] += 1
                continue
            if f == 1:
                err = abs(T - J2(q) / 3)
                worst["triv"] = max(worst["triv"], err)
                counts["triv"] += 1
                rows.append((f, G.order_of(js), T.real, "J2/3"))
                continue
            Lstar = L2_primitive(G, js, f)
            Lq = euler_restrict(G, js, f, q, Lstar)
            face_L = (2 * q * q / (pi * pi)) * Lq
            errL = abs(T - face_L) / (1 + abs(T))
            worst["L"] = max(worst["L"], errL)
            counts["L"] += 1
            tag = "L"
            if f == q:                       # primitive: Gauss x Bernoulli face
                face_GB = 2 * gauss_sum(G, js, f) * B2chi(G, js, f).conjugate()
                errGB = abs(T - face_GB) / (1 + abs(T))
                worst["GB"] = max(worst["GB"], errGB)
                counts["GB"] += 1
                tag = "L+GB"
            rows.append((f, G.order_of(js), T, tag))
        ev = sum(1 for js in G.chars() if G.parity(js) == 1)
        od = G.phi - ev
        print(f"  q={q:>4}: phi={G.phi:>3} chars ({ev} even, {od} odd); "
              f"J_2(q)/3 = {J2(q) // 3 if J2(q) % 3 == 0 else J2(q) / 3}")
        if q <= 15:
            for f, o, T, tag in rows:
                Ts = (f"{T:.6f}" if isinstance(T, float)
                      else f"{T.real:+.6f}{T.imag:+.6f}i")
                print(f"        even chi: cond={f:>2} ord={o:>2}  T = {Ts}   [{tag}]")
    # exact-rational face at q = 5 and q = 13 (even primitive quadratic)
    print("  exact faces (quadratic, q == 1 mod 4): B_{2,chi} exact, tau = sqrt(q):")
    for q in (5, 13):
        G = UGroup(q)
        qjs = None
        for js in G.chars():
            if G.order_of(js) == 2 and G.parity(js) == 1:
                qjs = js
        B = Fr(0)
        for a in range(1, q):
            c = round(G.chistar(qjs, a, q).real)
            B += c * (Fr(a, q) ** 2 - Fr(a, q) + Fr(1, 6))
        B *= q
        T = sum(G.chi(qjs, u) * csc2(u / q) for u in G.units).real
        pred = 2 * sqrt(q) * float(B)
        print(f"    q={q}: B_2,chi = {B} exact; T = {T:.10f} vs 2 sqrt(q) B = "
              f"{pred:.10f}  (match {abs(T - pred) < 1e-9})")
    print(f"  BATTERY: odd chars {counts['odd']}: max|T| = {worst['odd']:.2e} "
          f"(all == 0); trivial {counts['triv']}: max err = {worst['triv']:.2e}; "
          f"L-face {counts['L']}: max rel err = {worst['L']:.2e}; "
          f"Gauss-Bernoulli {counts['GB']}: max rel err = {worst['GB']:.2e}")
    assert worst["odd"] < 1e-9 and worst["triv"] < 1e-8
    assert worst["L"] < 1e-9 and worst["GB"] < 1e-9
    print("  PART C: ALL FOUR FACES EXACT.")


# ---------------------------------------------------------------- PARTS A/B/D/E

def run_cluster(E, s, name, ncf=6):
    print("=" * 78)
    print(f"CLUSTER {name}: E={E}, s={s}")
    P, M, owners, X, W = cluster_spectra(E, s)
    G = UGroup(P)
    print(f"  P = {P}, M = {M}, owners = {owners}, phi(P) = {G.phi}")
    cross = cross_all(P, X, W, G.units)
    mean = sum(cross.values()) / G.phi
    var = sum(c * c for c in cross.values()) / G.phi - mean * mean

    # ---- full character spectrum
    spec = {}
    for js in G.chars():
        spec[js] = chat_direct(G, js, cross)

    # (A) parity law
    odd_max = max(abs(spec[js]) for js in spec if G.parity(js) == -1)
    scale = max(abs(c) for c in cross.values())
    print(f"  (A) PARITY LAW: {sum(1 for js in spec if G.parity(js) == -1)} odd "
          f"characters, max |c-hat| = {odd_max:.2e} (scale {scale:.0f}) -> ALL ZERO")
    assert odd_max < 1e-6 * max(scale, 1.0)
    # conjugation symmetry
    conj_err = max(abs(spec[js] -
                       spec[tuple((-j) % o for j, o in zip(js, G.orders))]
                       .conjugate()) for js in spec)
    print(f"      conjugation symmetry c-hat(chi-bar) = conj(c-hat(chi)): "
          f"max err {conj_err:.2e}")
    assert conj_err < 1e-7 * max(scale, 1.0)

    # decode table: the named small-conductor characters
    print("  (A) decode of the small conductors (parity is the whole story):")
    for js in sorted(spec, key=lambda t: (G.conductor(t), G.order_of(t))):
        f, o, par = G.conductor(js), G.order_of(js), G.parity(js)
        if f in (3, 4, 5, 7) and o in (2, 3, 4):
            z = spec[js]
            zs = f"{z.real:+9.4f}{z.imag:+.4f}i" if abs(z.imag) > 1e-8 \
                else f"{z.real:+9.4f}"
            veh = "ZERO (odd: parity law)" if par == -1 else "carries mass (even)"
            print(f"      cond={f} ord={o} {'even' if par == 1 else 'odd '}: "
                  f"c-hat = {zs}   {veh}")

    # (B) support law: two spot checks with cond not dividing P/g
    checks = []
    for js in spec:
        f = G.conductor(js)
        if f <= 1 or G.parity(js) == -1:
            continue
        for g in divisors(P):
            if g in (1, P):
                continue
            if (P // g) % f != 0:
                checks.append((js, g))
                break
        if len(checks) >= 2:
            break
    for js, g in checks:
        cc = class_contribution(G, js, g, P, X, W, cross)
        print(f"  (B) SUPPORT LAW: cond={G.conductor(js)} char, class g={g} "
              f"(q={P // g}): direct class term = {abs(cc):.2e} -> ZERO")
        assert abs(cc) < 1e-6 * max(scale, 1.0)

    # (E) census: even masses, Parseval, co-resonant list
    triv = tuple(0 for _ in G.orders)
    masses = [(abs(spec[js]) ** 2, js) for js in spec if js != triv]
    mass_sum = sum(m for m, _ in masses)
    print(f"  (E) mean(cross) = c-hat(chi_0) = {mean:.4f} (the LEM-030 baseline row)")
    print(f"      Parseval: sum of nontrivial masses = {mass_sum:.1f} vs "
          f"frame variance = {var:.1f} (rel err "
          f"{abs(mass_sum - var) / max(var, 1e-9):.1e})")
    assert abs(mass_sum - var) < 1e-6 * max(var, 1.0)
    masses.sort(reverse=True)
    print("      TOP MASSES (the co-resonant list):")
    for m, js in masses[:10]:
        z = spec[js]
        zs = f"{z.real:+9.3f}{z.imag:+7.3f}i" if abs(z.imag) > 1e-7 \
            else f"{z.real:+9.3f}        "
        print(f"        cond={G.conductor(js):>4} ord={G.order_of(js):>2} "
              f"mass={m:10.2f}  c-hat = {zs}")
    bycond = {}
    for m, js in masses:
        bycond[G.conductor(js)] = bycond.get(G.conductor(js), 0.0) + m
    tot = sum(bycond.values())
    top = sorted(bycond.items(), key=lambda kv: -kv[1])[:6]
    print("      mass by conductor: " + ", ".join(
        f"{f}: {v / tot * 100:.1f}%" for f, v in top if v > 1e-9))

    # (D) closed-form weight side on selected characters (top masses + named)
    sel = [js for _, js in masses[:ncf] if abs(spec[js]) > 1e-8]
    named = [js for js in spec
             if G.conductor(js) in (5, 7) and G.order_of(js) in (2, 3)
             and G.parity(js) == 1]
    for js in named:
        if js not in sel:
            sel.append(js)
    worstD = 0.0
    print("  (D) CLOSED-FORM WEIGHT SIDE: c-hat(chi) from PURE L-VALUES:")
    for js in sel:
        direct = spec[js]
        fact = chat_factorized(G, js, P, X)
        err = abs(direct - fact) / (1 + abs(direct))
        worstD = max(worstD, err)
        z = f"{direct.real:+.4f}" + (f"{direct.imag:+.4f}i"
                                     if abs(direct.imag) > 1e-7 else "")
        print(f"      cond={G.conductor(js):>4} ord={G.order_of(js):>2}: "
              f"direct {z:>20}; sum_g (2/g^2) L_(P/g)(2,chi) Xhat_g = "
              f"{fact.real:+.4f}{fact.imag:+.4f}i  (rel err {err:.1e})")
    print(f"      worst rel err over {len(sel)} characters: {worstD:.2e}")
    assert worstD < 1e-7
    return P, G, spec, mean, var


if __name__ == "__main__":
    print("THE PARITY LAW AND THE L-VALUE WEIGHTS -- LEM-032 referee (boxeph S61)")
    part_C()
    run_cluster([12, 15, 20, 21, 28, 30, 35], 0, "balanced (S60 flagship)")
    run_cluster([1, 2, 3, 4, 5, 36, 60], 0, "two-owner")
    run_cluster([1, 2, 3, 4, 5, 6, 60], 0, "family {1..6,60}")
    print("=" * 78)
    print("done")
