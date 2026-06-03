#!/usr/bin/env python3
"""
claudebox-2026-06-02-S579 : the apex-lift certificate sheaf (HYP-2101)

Goal. Recast the S559/HYP-2063 polynomial-method corrector (Sungkawichai-
Trakulthongchai Prop 4.1 at k+1 = 2q, the n=14 wall) as a *sheaf of loneliness
certificates* on a projective line of corrector-directions over F_q, and make
the "apex obstruction" a precise geometric statement: the apex runner is the
UNIQUE local certificate whose support is the whole line (not a point), so it
kills every global section. Then test the "lift": adjoining the r/p (mod-p) time
freedom of t = s/(2q) + r/p turns that whole-line section back into a proper
hyperplane (codim 1), restoring transversality at the apex.

Everything here re-encodes / cross-checks S559 (must agree) and then adds the
geometric organization + the explicit lift computation.

Definitions (n = 2q, runners i = 1..2q-1, tuple v in Z_{2q}^{2q-1}):
  CRT reduction (S559, PROVED): a strict-interior corrector
      exists  <=>  exists s',r' in F_q^x : s' w_i + r' c_i != f_i  for all i,
  where w_i = v_i mod q, c_i = i mod q, f_i = 0 if (v_i+i) even else q-1.

  CERTIFICATE LINE L_i  := { (s,r) in A^2(F_q) : s w_i + r c_i = f_i }   (forbidden)
  CERTIFICATE LOCUS     := T \ U_i (L_i cap T),  T = (F_q^x)^2  (correctors)
  CERTIFICATE SHEAF C_v := extension-by-zero of the constant sheaf on the locus;
                           H^0(C_v) = correctors. Corrector exists <=> H^0 != 0.

  PARITY-MATCHED case (all f_i = 0): every L_i passes through the origin, so it
  is the single projective point [slope] rho_i = -c_i / w_i in P^1(F_q). The
  certificate sheaf lives on P^1; non-apex runner i deletes the one point rho_i;
  the apex i=q has c_q = 0 so rho_q is the slope-0 point, and when w_q = 0 (apex
  speed divisible by q, e.g. tight tuple) its "line" degenerates to 0 = 0, i.e.
  the WHOLE plane is forbidden -> the unique whole-line section.
"""

import itertools
from math import gcd

def units_mod(N):
    return [x for x in range(1, N) if gcd(x, N) == 1]

# ---------------------------------------------------------------------------
# A. Brute full-Z_{2q} corrector vs CRT-reduced mod-q corrector (reproduce S559)
# ---------------------------------------------------------------------------
def corrector_full(v, q):
    """Exists units s,r mod 2q with (s v_i + r i) mod 2q in {1,..,2q-2} for all i?"""
    N = 2 * q
    U = units_mod(N)
    idx = list(range(1, N))  # runners 1..2q-1
    for s in U:
        for r in U:
            ok = True
            for i in idx:
                land = (s * v[i - 1] + r * i) % N
                if land == 0 or land == N - 1:
                    ok = False
                    break
            if ok:
                return True
    return False

def corrector_crt(v, q, exclude_apex=False):
    """S559 mod-q reduction: exists s',r' in F_q^x with s' w_i + r' c_i != f_i all i."""
    N = 2 * q
    Uq = units_mod(q)
    idx = list(range(1, N))
    w = {i: v[i - 1] % q for i in idx}
    c = {i: i % q for i in idx}
    f = {i: (0 if (v[i - 1] + i) % 2 == 0 else q - 1) for i in idx}
    for s in Uq:
        for r in Uq:
            ok = True
            for i in idx:
                if exclude_apex and i == q:
                    continue
                if (s * w[i] + r * c[i]) % q == f[i]:
                    ok = False
                    break
            if ok:
                return True
    return False

# ---------------------------------------------------------------------------
# B. Certificate sheaf on P^1 : forbidden slopes, apex degeneration, counts
# ---------------------------------------------------------------------------
def certificate_locus(v, q, exclude_apex=False):
    """Return (#correctors over T=(F_q^x)^2, apex_is_whole_plane, forbidden_slopes)."""
    N = 2 * q
    Uq = units_mod(q)
    idx = list(range(1, N))
    w = {i: v[i - 1] % q for i in idx}
    c = {i: i % q for i in idx}
    f = {i: (0 if (v[i - 1] + i) % 2 == 0 else q - 1) for i in idx}

    # apex degeneration test (line 0 = 0 forbids everything)
    apex_whole = (c[q] % q == 0 and w[q] % q == 0 and f[q] % q == 0)

    # forbidden slopes for parity-matched lines through origin (f_i = 0, w_i != 0)
    forbidden_slopes = set()
    parity_matched = all(f[i] == 0 for i in idx)
    if parity_matched:
        for i in idx:
            if exclude_apex and i == q:
                continue
            if w[i] % q != 0:
                forbidden_slopes.add((-c[i] * pow(w[i], q - 2, q)) % q)  # rho_i = -c/w
            # w_i == 0 with f_i=0 -> whole-plane forbidder (the apex when tight)

    # honest count over the torus
    count = 0
    for s in Uq:
        for r in Uq:
            ok = True
            for i in idx:
                if exclude_apex and i == q:
                    continue
                if (s * w[i] + r * c[i]) % q == f[i]:
                    ok = False
                    break
            if ok:
                count += 1
    return count, apex_whole, sorted(forbidden_slopes), parity_matched

# ---------------------------------------------------------------------------
# C. The lift : adjoin mod-p freedom; does the apex line regain codim 1?
# ---------------------------------------------------------------------------
def is_prime(p):
    if p < 2:
        return False
    d = 2
    while d * d <= p:
        if p % d == 0:
            return False
        d += 1
    return True

def smallest_prime_not_dividing(N):
    p = 2
    while True:
        if is_prime(p) and N % p != 0:
            return p
        p += 1

def apex_forbidden_fraction(v, q, lifted):
    """
    Fraction of corrector-parameter space forbidden by the APEX runner alone.
      unlifted: parameters (s,r) in (F_q^x)^2, apex constraint s w_q + r c_q = f_q.
      lifted  : t = s/(2q) + r/p ; over modulus M = 2q*p the apex speed mod p is
                q mod p (a unit since p does not divide q), so the apex constraint
                acquires a mod-p slot and is no longer parameter-free.
    Returns (fraction_forbidden, detail).
    """
    N = 2 * q
    if not lifted:
        Uq = units_mod(q)
        wq, cq = v[q - 1] % q, q % q
        fq = 0 if (v[q - 1] + q) % 2 == 0 else q - 1
        bad = sum(1 for s in Uq for r in Uq if (s * wq + r * cq) % q == fq)
        return bad / (len(Uq) ** 2), f"apex line over F_q: s*{wq}+r*{cq}={fq}"
    else:
        p = smallest_prime_not_dividing(N)
        M = N * p
        UM = units_mod(M)
        # apex runner q lands at (s*v_q + r_tight*q) mod M ; tight corrector r=s.
        # Forbidden = lands in {0, M-1} (strict-interior failure) for the apex only.
        vq = v[q - 1]
        bad = sum(1 for s in UM if (s * vq + s * q) % M in (0, M - 1))
        return bad / len(UM), f"lifted mod {M}=2*{q}*{p}: apex speed mod p = {q % p}"

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def tight_tuple(q):
    return list(range(1, 2 * q))  # (1,2,...,2q-1)

def main():
    print("=" * 72)
    print("S579  apex-lift certificate sheaf  (HYP-2101)  --  building on S559")
    print("=" * 72)

    for q in (3, 5, 7, 11):
        N = 2 * q
        print(f"\n########## q={q}  (n={N}={2}*{q}) ##########")
        v = tight_tuple(q)

        # --- A. cross-check S559: full vs CRT reduction must agree, both should
        #        fail for the tight tuple (apex obstruction), succeed apex-excluded
        full = corrector_full(v, q)
        crt = corrector_crt(v, q)
        crt_excl = corrector_crt(v, q, exclude_apex=True)
        print(f"[A] tight tuple corrector  full-Z_{N}: {full}   CRT mod-{q}: {crt}"
              f"   (must match: {full == crt})")
        print(f"    apex-excluded corrector (CRT): {crt_excl}  "
              f"(S559 tight-tuple repair => expect True)")

        # --- B. sheaf picture
        cnt, apex_whole, slopes, pm = certificate_locus(v, q)
        cnt_x, _, slopes_x, _ = certificate_locus(v, q, exclude_apex=True)
        print(f"[B] certificate locus  |H^0|(all runners) = {cnt}   "
              f"parity-matched={pm}")
        print(f"    apex section is the WHOLE plane (degenerate): {apex_whole}")
        print(f"    |H^0|(apex excluded) = {cnt_x}   forbidden slopes (non-apex) "
              f"rho_i = {slopes_x}")
        Fqx = units_mod(q)
        print(f"    F_q^x = {Fqx} ; slopes cover F_q^x: "
              f"{set(slopes_x) >= set(Fqx)}  "
              f"(S559 residual: tight collapses slopes to a single value)")

        # --- C. the lift
        fr0, d0 = apex_forbidden_fraction(v, q, lifted=False)
        fr1, d1 = apex_forbidden_fraction(v, q, lifted=True)
        print(f"[C] apex forbidden-fraction  UNLIFTED = {fr0:.3f}   ({d0})")
        print(f"    apex forbidden-fraction  LIFTED   = {fr1:.3f}   ({d1})")
        print(f"    lift restores codim-1 (fraction drops from 1): {fr0==1.0 and fr1<1.0}")

    # --- random tuple sanity: CRT reduction == full check (S559 '0 mismatches')
    print("\n" + "=" * 72)
    print("Cross-check on random tuples: CRT reduction vs full-Z_{2q} (0 mismatches expected)")
    print("=" * 72)
    import random
    random.seed(579)
    mism = 0
    trials = 0
    for q in (3, 5, 7):
        N = 2 * q
        for _ in range(400):
            v = [random.randrange(1, N) for _ in range(N - 1)]
            if corrector_full(v, q) != corrector_crt(v, q):
                mism += 1
            trials += 1
    print(f"trials={trials}  mismatches={mism}  "
          f"(S559 reduction theorem confirmed: {mism == 0})")

if __name__ == "__main__":
    main()
