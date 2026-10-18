"""S40 audit: THM-4666 statements 1-4 and the FINITE-EXACT counts of statement 7, with independent code.

Checks
 (1) excess of the natural order of Gamma_D: direct pair count (|D| <= 2000) and the distance formula, against
     h(D) from (a) Cohen-style reduced-form count, (b) the analytic theta/erfc series for L(1, chi_D);
     multiplier orders k -> a k + s (all a, several s) and all cyclic shifts by direct pair count (|D| <= 230);
     all fundamental D >= -20000, a fresh random sample of fundamental D in [-10^6, -20000), all primes p = 3 mod 4
     below 2*10^5 and a random sample of primes in (2*10^5, 10^6].
 (2) half window sum_(n <= (p-1)/2) (n/p) = (2 - (2/p)) h(-p), p > 3.
 (3) arcs at a vertex = split primes: roots of x^2 - x + (p+1)/4 mod q by brute force against (q/p).
 (4) h(-p) = 1 list; Euler's lucky numbers by brute force; Rabinowitsch correspondence.
Usage: python audit_thm4666_identities.py
"""
import sys, os, math, json, random, time
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_common import sieve, legendre, h_forms, fundamental_neg, char_table, L1_theta, paley_chi

OUT = {}


@njit(cache=True)
def excess_pairs(chi):
    n = chi.shape[0]
    s = 0
    for i in range(n):
        for j in range(i + 1, n):
            s += chi[(j - i) % n]
    return s


@njit(cache=True)
def excess_formula(chi):
    n = chi.shape[0]
    s = 0
    for d in range(1, n):
        s += (n - d) * chi[d]
    return s


@njit(cache=True)
def excess_affine(chi, a, s):
    n = chi.shape[0]
    tot = 0
    for k in range(n):
        u = (a * k + s) % n
        for l in range(k + 1, n):
            v = (a * l + s) % n
            tot += chi[(v - u) % n]
    return tot


@njit(cache=True)
def walk_area_and_half(chi):
    p = chi.shape[0]
    s = 0
    area = 0
    half = 0
    for t in range(1, p):
        s += chi[t]
        area += s
        if t == (p - 1) // 2:
            half = s
    return area, half, s


def wD(D):
    return 6 if D == -3 else (4 if D == -4 else 2)


def check_disc(D, direct=False):
    chi = char_table(D).astype(np.int64)
    n = -D
    hF = h_forms(D)
    L = L1_theta(chi, n)
    hL = wD(D) * math.sqrt(n) * L / (2 * math.pi)
    ex = excess_formula(chi)
    ok = (ex * wD(D) == 2 * n * hF) and abs(hL - hF) < 1e-6
    exd = excess_pairs(chi) if direct else None
    if direct:
        ok = ok and (exd == ex)
    return ok, hF, hL, ex


def main():
    t0 = time.time()
    # ---------- fundamental discriminants D >= -20000 ----------
    fund = [-n for n in range(3, 20001) if fundamental_neg(-n)]
    bad = []
    maxdev = 0.0
    for D in fund:
        ok, hF, hL, ex = check_disc(D, direct=(-D <= 2000))
        maxdev = max(maxdev, abs(hL - hF))
        if not ok:
            bad.append(D)
    print('fundamental D in [-20000, -3]: %d discriminants; identity/analytic/direct failures: %d; max |h_L - h_forms| = %.2e'
          % (len(fund), len(bad), maxdev), flush=True)
    OUT['fundamental_20000'] = dict(count=len(fund), failures=bad, max_analytic_dev=maxdev)
    # ---------- fresh random sample of fundamental D in [-10^6, -20000) ----------
    random.seed(20261017)
    sample = []
    while len(sample) < 1500:
        n = random.randrange(20001, 1000001)
        if fundamental_neg(-n):
            sample.append(-n)
    bad = []
    maxdev = 0.0
    for D in sample:
        ok, hF, hL, ex = check_disc(D)
        maxdev = max(maxdev, abs(hL - hF))
        if not ok:
            bad.append(D)
    print('random fundamental D in [-10^6, -20000): %d; failures %d; max |h_L - h_forms| = %.2e' % (len(sample), len(bad), maxdev),
          flush=True)
    OUT['fundamental_random_1e6'] = dict(count=len(sample), failures=bad, max_analytic_dev=maxdev)
    # ---------- multiplier orders and cyclic shifts (direct pair counts) ----------
    badm = 0
    nchk = 0
    for D in [d for d in fund if -d <= 230]:
        chi = char_table(D).astype(np.int64)
        n = -D
        nat = excess_formula(chi)
        for a in range(1, n):
            if math.gcd(a, n) != 1:
                continue
            for s in (0, 1, n // 2, n - 1):
                e = excess_affine(chi, a, s)
                nchk += 1
                if e != chi[a] * nat:
                    badm += 1
        for s in range(n):
            e = excess_affine(chi, 1, s)
            nchk += 1
            if e != nat:
                badm += 1
    print('multiplier/shift orders, fundamental |D| <= 230: %d orders checked, failures %d' % (nchk, badm), flush=True)
    OUT['affine_orders'] = dict(checked=nchk, failures=badm)
    # ---------- Paley primes p = 3 mod 4 below 2*10^5 ----------
    S = sieve(1000000)
    P = [p for p in range(3, 200001, 4) if S[p]]
    bad1 = []; bad2 = []; badchi = []; one = []; area_bad = []
    hs = {}
    maxdev = 0.0
    for p in P:
        chi = paley_chi(p)
        if p < 3000 or p % 97 == 3:   # cross-check the squares table against Euler's criterion on a sample
            for a in range(1, p, max(1, p // 50)):
                if chi[a] != legendre(a, p):
                    badchi.append(p)
                    break
        hF = h_forms(-p)
        L = L1_theta(chi, p)
        hL = wD(-p) * math.sqrt(p) * L / (2 * math.pi)
        maxdev = max(maxdev, abs(hL - hF))
        hs[p] = hF
        ex = excess_formula(chi)
        area, half, last = walk_area_and_half(chi)
        if ex * wD(-p) != 2 * p * hF or abs(hL - hF) > 1e-6 or area != ex or last != 0:
            bad1.append(p)
        if p > 3:
            two = 1 if p % 8 == 7 else -1
            if half != (2 - two) * hF:
                bad2.append(p)
            if area != p * hF:
                area_bad.append(p)
        if hF == 1:
            one.append(p)
    print('Paley primes p = 3 mod 4, p <= 2e5: %d primes (incl. 3); statement-1 failures %d; half-window failures %d; '
          'area != p h failures %d; chi table vs Euler failures %d; max |h_L - h_forms| = %.2e'
          % (len(P), len(bad1), len(bad2), len(area_bad), len(badchi), maxdev), flush=True)
    print('h(-p) = 1 for p = 3 mod 4 prime <= 2e5:', one, flush=True)
    OUT['paley'] = dict(primes=len(P), stmt1_failures=bad1, half_failures=bad2, area_failures=area_bad, h1=one)
    # random primes in (2e5, 1e6]
    big = [p for p in range(200003, 1000001, 4) if S[p]]
    random.seed(4666)
    smp = random.sample(big, 400)
    badb = []
    for p in smp:
        chi = paley_chi(p)
        hF = h_forms(-p)
        hL = math.sqrt(p) * L1_theta(chi, p) / math.pi
        ex = excess_formula(chi)
        area, half, last = walk_area_and_half(chi)
        two = 1 if p % 8 == 7 else -1
        if ex != p * hF or abs(hL - hF) > 1e-6 or half != (2 - two) * hF:
            badb.append(p)
    print('random primes p = 3 mod 4 in (2e5, 1e6]: %d checked, failures %d' % (len(smp), len(badb)), flush=True)
    OUT['paley_random_big'] = dict(checked=len(smp), failures=badb)
    # ---------- statement 3: split primes vs arcs ----------
    qs = [q for q in range(2, 3000) if S[q]]
    bad3 = 0; n3 = 0
    for p in [p for p in P if p < 1200]:
        c = (p + 1) // 4
        for q in qs:
            if q == p:
                continue
            # roots of x^2 - x + c mod q (minimal polynomial of (1 + sqrt(-p))/2)
            xs = np.arange(q, dtype=np.int64)
            roots = int(np.count_nonzero((xs * xs - xs + c) % q == 0))
            split = (roots == 2)
            inert = (roots == 0)
            assert split or inert, (p, q, roots)
            arc = (pow(q % p, (p - 1) // 2, p) == 1)      # 0 -> q mod p in P_p
            n3 += 1
            if split != arc:
                bad3 += 1
    print('statement 3: %d (p, q) pairs (p < 1200, q < 3000 incl. q = 2), split <=> arc failures %d' % (n3, bad3), flush=True)
    OUT['split_vs_arc'] = dict(pairs=n3, failures=bad3)
    # ---------- statement 4: Euler's lucky numbers and Rabinowitsch ----------
    lucky = []
    for q in range(2, 991):
        if all(S[m * m + m + q] for m in range(0, q - 1)):
            lucky.append(q)
    rab = [q for q in range(1, 50001) if (4 * q - 1) in hs and hs[4 * q - 1] == 1]
    print('Euler lucky numbers q in [2, 990] (n^2+n+q prime for 0 <= n <= q-2):', lucky, '; q in [1, 50000] with 4q-1 prime and h(1-4q) = 1:', rab,
          flush=True)
    OUT['lucky'] = lucky
    OUT['rabinowitsch'] = rab
    print('elapsed %.1f s' % (time.time() - t0))
    json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'audit_thm4666_identities.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
