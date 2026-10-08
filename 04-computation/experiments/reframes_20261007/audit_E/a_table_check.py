#!/usr/bin/env python3
"""Audit E, task A: re-derive and check the px+1 pair-chain table and every line of THM-4606 Proof 1.

Independent of the author's scripts. Three layers:
 (A1) The table itself, checked against DIRECT 2-adic arithmetic: y is a random 2-adic integer (mod 2^M), u = p^k y + e
     with k in [-6, 6] and e in Z_(2) cap Q (denominators: powers of p times other odd numbers, so e need not lie in
     Z[1/p]); one T_p step is applied to u and to v = y directly, and u' = p^k' v' + e' is checked mod 2^(M-1) for the
     (k', e') the table predicts. All four (sigma, beta) branches are hit.
 (A2) Every displayed f' formula of Proof 1 (k >= 1, k <= -1, k = 0), with f = e p^-max(k,0), checked as an exact rational
     identity on random states (also for non-integer e), plus: lambda in {1/2, p/2}, sharp max |eta|, the claimed bound
     1/2 + 1/(2p), bijectivity of the multiplier assignment off departures, lambda = 1/2 for both coins at departures,
     and the k-move rule k' = k + 1 - 2 beta at every flip.
 (A3) The chain run alongside actual big-integer orbits (y random 1500-bit integer, u_0 = y + e_0 or p^k0 y + e_0,
     k0 >= 0, and k0 < 0 handled 2-adically): state identity u_n = p^k_n v_n + e_n, k_n - k_0 = (odd steps of u) - (odd
     steps of v), and merge u_n = v_n <=> state (0,0) (on these samples).
 (A4) Skeleton: empirical fairness/independence of beta at flip times, split by the sign of k and by p (sanity only;
     the proof is optional skipping).
"""
from fractions import Fraction as Fr
import random, math, sys

rnd = random.Random(20261008)
M = 400                       # 2-adic precision for (A1)
MOD = 1 << M

def inv2(b, mod=MOD):         # inverse of an odd integer mod 2^M
    return pow(b, -1, mod)

def to2adic(x, mod=MOD):      # rational with odd denominator -> residue mod 2^M
    x = Fr(x)
    assert x.denominator % 2 == 1
    return (x.numerator * inv2(x.denominator, mod)) % mod

def T_mod(p, x, mod):         # one T_p step on a residue x mod 'mod' (=2^m); result valid mod mod/2
    if x % 2 == 0:
        return (x // 2) % (mod // 2)
    return ((p * x + 1) // 2) % (mod // 2)

def par(e):                   # 2-adic parity of a rational with odd denominator
    e = Fr(e)
    return (e.numerator * pow(e.denominator, -1, 2)) % 2   # = numerator mod 2 since denominator odd

def table(p, k, e, b):        # THM-4606 Setting table (the object under audit)
    s = par(e)
    P = Fr(p)
    if s == 0 and b == 0: return k, e / 2
    if s == 0 and b == 1: return k, (p * e + 1 - P ** k) / 2
    if s == 1 and b == 0: return k + 1, (p * e + 1) / 2
    return k - 1, (e - P ** (k - 1)) / 2

def rand_e(p, k):
    # random e in Z_(2) cap Q: numerator / (p^j * odd), j random (covers Z[1/p] and beyond)
    num = rnd.randint(-10**6, 10**6)
    j = rnd.randint(0, 3) + max(0, -k)
    odd = rnd.choice([1, 1, 1, 3, 5, 7, 9, 11, 13, 15, 21, 25, 27])
    return Fr(num, p ** j * odd)

# ---------------- (A1) table vs direct 2-adic arithmetic ----------------
def check_A1(ntrials=60000):
    hits = {}
    bad = 0
    for _ in range(ntrials):
        p = rnd.choice([1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 63, 101])
        k = rnd.randint(-6, 6)
        e = rand_e(p, k)
        y = rnd.getrandbits(M)
        pk = to2adic(Fr(p) ** k)
        u = (pk * y + to2adic(e)) % MOD
        v = y
        s, b = (u - v) % 2, v % 2      # NB: sigma computed from u - v... check consistency below
        # u mod 2 should be beta xor sigma(e)
        if (u % 2) != (b ^ par(e)):
            bad += 1; print("parity mismatch", p, k, e); continue
        u1, v1 = T_mod(p, u, MOD), T_mod(p, v, MOD)
        k1, e1 = table(p, k, e, b)
        mod1 = MOD // 2
        rhs = (to2adic(Fr(p) ** k1, mod1) * v1 + to2adic(e1, mod1)) % mod1
        key = (par(e), b)
        hits[key] = hits.get(key, 0) + 1
        if rhs != u1:
            bad += 1
            if bad < 5: print("TABLE MISMATCH", p, k, e, b, k1, e1)
    print(f"(A1) table vs direct 2-adic step, {ntrials} random (p,k,e,y) incl. p=1, k in [-6,6], e in Z_(2) cap Q:"
          f" mismatches {bad}; branch counts {dict(sorted(hits.items()))}")
    return bad == 0

# ---------------- (A2) Proof 1 formulas ----------------
def proof1_claim(p, k, e, s, b):
    """Return (k', f') as displayed in THM-4606 Proof 1 (f = e p^-max(k,0))."""
    P = Fr(p)
    f = e / P ** max(k, 0)
    if k >= 1:
        if (s, b) == (0, 0): return k, f / 2
        if (s, b) == (0, 1): return k, P / 2 * f - (1 - P ** (-k)) / 2
        if (s, b) == (1, 0): return k + 1, f / 2 + 1 / (2 * P ** (k + 1))
        return k - 1, P / 2 * f - Fr(1, 2)
    if k <= -1:
        if (s, b) == (0, 0): return k, e / 2
        if (s, b) == (0, 1): return k, P / 2 * e + (1 - P ** k) / 2
        if (s, b) == (1, 0): return k + 1, P / 2 * e + Fr(1, 2)
        return k - 1, e / 2 - P ** (k - 1) / 2
    # k == 0
    if (s, b) == (0, 0): return 0, e / 2
    if (s, b) == (0, 1): return 0, P * e / 2
    if (s, b) == (1, 0): return 1, e / 2 + 1 / (2 * P)
    return -1, e / 2 - 1 / (2 * P)

def lam_of(p, k, s, b):
    """multiplier lambda assigned by Proof 1 to the branch"""
    P = Fr(p)
    if s == 0: return Fr(1, 2) if b == 0 else P / 2
    if k == 0: return Fr(1, 2)                       # departure: both coins 1/2
    if k >= 1: return Fr(1, 2) if b == 0 else P / 2  # away (up) pays 1/2, toward 0 (down) costs p/2
    return P / 2 if b == 0 else Fr(1, 2)             # k <= -1: toward 0 (up) costs p/2

def check_A2(ntrials=200000):
    bad_formula = bad_bij = bad_dep = bad_move = 0
    max_eta = {}
    sup_eta = Fr(0)
    for t in range(ntrials):
        p = rnd.choice([3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 25, 27, 31, 33, 45, 63, 101])
        k = rnd.randint(-8, 8)
        if t % 7 == 0: k = 0
        e = rand_e(p, k)
        if t % 11 == 0: e = Fr(rnd.randint(-50, 50), p ** max(0, -k))   # small states too
        s = par(e)
        P = Fr(p)
        f = e / P ** max(k, 0)
        lams = []
        for b in (0, 1):
            k1, e1 = table(p, k, e, b)
            f1 = e1 / P ** max(k1, 0)
            kc, fc = proof1_claim(p, k, e, s, b)
            if (kc, fc) != (k1, f1):
                bad_formula += 1
                if bad_formula < 5: print("FORMULA MISMATCH", p, k, e, s, b, (k1, f1), (kc, fc))
            if s == 1 and k1 != k + 1 - 2 * b: bad_move += 1
            lam = lam_of(p, k, s, b)
            eta = f1 - lam * f
            lams.append(lam)
            key = (('k>=1' if k >= 1 else 'k<=-1' if k <= -1 else 'k=0'), s, b)
            if abs(eta) > max_eta.get(key, -1): max_eta[key] = abs(eta)
            # the claimed bound
            if abs(eta) > Fr(1, 2) + 1 / (2 * P): sup_eta = max(sup_eta, abs(eta))
        if s == 1 and k == 0:
            if lams != [Fr(1, 2), Fr(1, 2)]: bad_dep += 1
        else:
            if sorted(lams) != sorted([Fr(1, 2), Fr(p, 2)]): bad_bij += 1
    print(f"(A2) Proof 1 displayed formulas vs table on {ntrials} random states (|k|<=8, e incl. non-Z[1/p]):"
          f" formula mismatches {bad_formula}; flip move k'=k+1-2beta violations {bad_move};"
          f" departures with lambda != 1/2 {bad_dep}; non-departure states without one 1/2 and one p/2 {bad_bij};"
          f" |eta| > 1/2 + 1/(2p) found: {'none' if sup_eta == 0 else sup_eta}")
    print("(A2) sup |eta| observed per branch (key = (k-range, sigma, beta)); compare analytic sup:")
    for key in sorted(max_eta):
        print(f"      {key}: {float(max_eta[key]):.6f}")
    return bad_formula == bad_bij == bad_dep == bad_move == 0 and sup_eta == 0

def analytic_eta():
    # exact sup of |eta| per branch over all k in the region (closed forms):
    print("(A2) analytic |eta| per branch: k>=1: (0,0) 0; (0,1) (1-p^-k)/2 < 1/2; (1,0) 1/(2p^(k+1)) <= 1/(2p^2); (1,1) 1/2 exactly."
          " k<=-1: (0,0) 0; (0,1) (1-p^k)/2 < 1/2; (1,0) 1/2 exactly; (1,1) p^(k-1)/2 <= 1/(2p^2). k=0: (0,·) 0; flips 1/(2p)."
          " So max|eta| = 1/2 (attained), which is <= 1/2 + 1/(2p): the stated bound is TRUE but not sharp.")

# ---------------- (A3) chain vs actual big-integer orbits ----------------
def Tint(p, x):
    return x >> 1 if x % 2 == 0 else (p * x + 1) >> 1

def chain_int_step(p, k, E, b):
    """Integer-coordinate chain, E = p^max(0,-k) * e (derived independently in the report)."""
    s = E & 1
    if s == 0:
        if b == 0: return k, E >> 1
        if k >= 0: return k, (p * E + 1 - p ** k) >> 1
        return k, (p * E + p ** (-k) - 1) >> 1
    if b == 0:
        if k >= 0: return k + 1, (p * E + 1) >> 1
        return k + 1, (E + p ** (-k - 1)) >> 1
    if k >= 1: return k - 1, (E - p ** (k - 1)) >> 1
    return k - 1, (p * E - 1) >> 1          # k <= 0 (k = 0 departure down and k <= -1 flip down)

def check_A3(nsamp=3000, nsteps=300):
    bad_rel = bad_odd = bad_merge = 0
    merges = 0
    for i in range(nsamp):
        p = rnd.choice([3, 5, 7, 9, 11, 13])
        k0 = rnd.choice([0, 0, 0, 1, 2, 3])
        E0 = rnd.choice([1, -1, 2, 3, -3, 5, 7, rnd.randint(-1000, 1000)])
        if k0 == 0 and E0 == 0: E0 = 1
        y = rnd.getrandbits(1500) | (1 << 1499)
        v = y
        u = p ** k0 * y + E0
        k, E = k0, E0
        ou = ov = 0
        for n in range(nsteps):
            b = v & 1
            ou += u & 1; ov += v & 1
            k, E = chain_int_step(p, k, E, b)
            u, v = Tint(p, u), Tint(p, v)
            # relation: k >= 0: u = p^k v + E ; k < 0: p^|k| u = v + E
            ok = (u == p ** k * v + E) if k >= 0 else (p ** (-k) * u == v + E)
            if not ok:
                bad_rel += 1; print("RELATION FAIL", p, k0, E0, n); break
            if k - k0 != ou - ov: bad_odd += 1
            if (u == v) != (k == 0 and E == 0):
                bad_merge += 1; print("MERGE/ABSORB MISMATCH", p, k0, E0, n, k, E)
            if k == 0 and E == 0:
                merges += 1; break
    print(f"(A3) chain (integer coords) vs actual big-integer orbits: {nsamp} pairs x <= {nsteps} steps, p in 3..13,"
          f" k0 in 0..3: relation failures {bad_rel}, odd-count identity failures {bad_odd},"
          f" (u_n = v_n) != absorbed mismatches {bad_merge}; absorbed pairs {merges}")
    return bad_rel == bad_odd == bad_merge == 0

def check_A3_negk(ntrials=4000, nsteps=60):
    """k0 < 0: u_0 = p^k0 y + e0 is a 2-adic integer; work mod 2^Mbig and compare with the integer-coordinate chain."""
    Mbig = 800; mod = 1 << Mbig
    bad = 0
    for _ in range(ntrials):
        p = rnd.choice([3, 5, 7, 9, 11])
        k0 = rnd.randint(-5, -1)
        E0 = rnd.randint(-500, 500)
        y = rnd.getrandbits(Mbig)
        u = (pow(p, -(-k0), mod) * y + E0 * pow(p, k0, mod)) % mod   # p^k0 y + E0/p^|k0|
        v = y
        k, E = k0, E0
        m = mod
        for n in range(nsteps):
            b = v & 1
            k, E = chain_int_step(p, k, E, b)
            u, v = T_mod(p, u, m), T_mod(p, v, m)
            m //= 2
            if k >= 0: rhs = (pow(p, k, m) * v + E) % m
            else: rhs = (pow(p, k, m) * v + E * pow(p, k, m)) % m
            if rhs != u % m:
                bad += 1; print("NEG-K RELATION FAIL", p, k0, E0, n); break
    print(f"(A3') chain with k0 < 0 vs 2-adic orbits mod 2^{Mbig}: {ntrials} pairs x {nsteps} steps: failures {bad}")
    return bad == 0

# ---------------- (A4) skeleton sanity ----------------
def check_A4(npaths=1500, nsteps=1000):
    """beta is read off the ACTUAL orbit v_n = T^n(y) of a random 1200-bit integer y (parities of v_n, n < 1200, are
    exactly uniform i.i.d. under a uniform y mod 2^1200), not drawn from the RNG."""
    import collections
    for p in (5, 7, 9):
        up = collections.Counter(); tot = collections.Counter()
        pairs = collections.Counter()
        for _ in range(npaths):
            k, E = 0, 1
            prev = None
            v = rnd.getrandbits(1200)
            for n in range(nsteps):
                b = v & 1
                v = Tint(p, v)
                s = E & 1
                if s == 1:
                    cls = 'k>0' if k > 0 else 'k<0' if k < 0 else 'k=0'
                    tot[cls] += 1; up[cls] += (b == 0)
                    if prev is not None: pairs[(prev, b)] += 1
                    prev = b
                k, E = chain_int_step(p, k, E, b)
                if k == 0 and E == 0: break
                if abs(E) > 10**60 * p ** abs(k): break
        tot_pairs = sum(pairs.values())
        corr = (pairs[(0, 0)] + pairs[(1, 1)] - pairs[(0, 1)] - pairs[(1, 0)]) / max(1, tot_pairs)
        print(f"(A4) p={p}: P(up | flip) by sign of k: " +
              ", ".join(f"{c}: {up[c]/tot[c]:.4f} (n={tot[c]})" for c in ('k<0', 'k=0', 'k>0') if tot[c]) +
              f"; lag-1 correlation of successive flip coins {corr:+.4f} (n={tot_pairs})")

if __name__ == '__main__':
    ok = True
    ok &= check_A1()
    ok &= check_A2()
    analytic_eta()
    ok &= check_A3()
    ok &= check_A3_negk()
    check_A4()
    print("ALL A CHECKS PASS" if ok else "SOME A CHECK FAILED")
