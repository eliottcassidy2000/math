#!/usr/bin/env python3
"""Primes on the Collatz Fermat and Mersenne lines: a calibration of the random model behind "no more Fermat primes"
(opus-2026-10-26-S58).

Lines (exact Syracuse runs, THM-4793 (F) and THM-4556):
  Fermat-type line of 2^m + 1:  x_j = 3^j 2^(m-2j) + 1,  j = 1 .. floor((m-1)/2)   (valuation-two run)
  Mersenne-type line of 2^p - 1: y_k = 3^k 2^(p-k) - 1,  k = 1 .. p-1              (valuation-one run)
  genuine Fermat lines: m = 2^e.
Prediction (mode 'predict', run BEFORE any primality test): sieve every value by all primes q <= Q = 2^16 using the
exact congruences q | x_j  <=>  (3/4)^j = -2^(-m) (mod q),  q | y_k  <=>  (3/2)^k = 2^(-p) (mod q); a value below 2^32
is classified exactly and excluded from the random part.  For a survivor x >= 2^32 the model probability is
P(x) = prod_{q <= Q} (1 - 1/q)^(-1) / ln x  (the same model that gives P(F_n prime) ~ 2/ln F_n before sieving).
Expected count E = sum P, variance V = sum P(1-P); pre-registered window E +- 2.5 sqrt(V).
Mode 'test': count probable primes (sympy.isprime = BPSW) among the survivors >= 2^32.
Post-audit of the pre-registration (L-P1 failed low): modes 'predict2'/'test2' add the algebraic factorisations that
the small-prime sieve cannot see -- 3^j 2^b + 1 = w^r + 1 when an odd r divides gcd(j, b), 4u^4 + 1 (Aurifeuille)
when j = 0 mod 4 and b = 2 mod 4, and 3^k 2^c - 1 = w^r - 1 when r = gcd(k, c) > 1 -- and are run on NEW ranges
(out of sample), predicted and committed before testing.
Usage: python collatz_line_primes_20261026.py predict|test|posthoc|predict2|test2
"""
import math as _m


def algebraic_composite(kind, n, j):
    """True if the line value has an algebraic (or Aurifeuillean) factorisation with both factors > 1."""
    if kind == 'F':
        g = _m.gcd(j, n)                   # gcd(j, n - 2j) = gcd(j, n)
        while g % 2 == 0:
            g //= 2
        if g > 1:
            return True                    # w^r + 1, r odd
        return j % 4 == 0 and n % 4 == 2   # 4u^4 + 1
    return _m.gcd(j, n) > 1                # w^r - 1
import sys, os, json, math, time
import numpy as np
from sympy import primerange, isprime

HERE = os.path.dirname(os.path.abspath(__file__))
Q = 1 << 16
SMALL = 1 << 32
FAMILIES = {
    'fermat_type_m8_1200': ('F', list(range(8, 1201))),
    'mersenne_type_p8_400': ('M', list(range(8, 401))),
    'genuine_fermat_e3_12': ('F', [2 ** e for e in range(3, 13)]),
}


def line_values(kind, n):
    if kind == 'F':
        return [(j, 3 ** j * 2 ** (n - 2 * j) + 1) for j in range(1, (n - 1) // 2 + 1)]
    return [(k, 3 ** k * 2 ** (n - k) - 1) for k in range(1, n)]


def sieve_family(kind, ns, primes):
    """returns dict n -> set of indices j killed by some prime q <= Q (q < value)."""
    jmax = max(((n - 1) // 2 if kind == 'F' else n - 1) for n in ns)
    killed = {n: set() for n in ns}
    for q in primes:
        if q < 5:
            continue
        t = 3 * pow(4, -1, q) % q if kind == 'F' else 3 * pow(2, -1, q) % q
        # positions of t^j, j = 1..jmax (store first occurrence j0 and period o)
        seen = {}
        v = 1
        o = None
        for j in range(1, jmax + 1):
            v = v * t % q
            if v in seen:
                o = j - seen[v]
                break
            seen[v] = j
        for n in ns:
            c = (-pow(2, -n, q)) % q if kind == 'F' else pow(2, -n, q)
            if c in seen:
                j0 = seen[c]
                J = (n - 1) // 2 if kind == 'F' else n - 1
                step = o if o is not None else jmax + 1
                j = j0
                while j <= J:
                    killed[n].add(j)
                    j += step
    return killed


FAMILIES2 = {
    'fermat_type_m1201_1500': ('F', list(range(1201, 1501))),
    'mersenne_type_p401_520': ('M', list(range(401, 521))),
}


def predict_family(kind, ns, primes, mertens, algebraic):
    killed = sieve_family(kind, ns, primes)
    E = V = 0.0
    surv = alg = 0
    for n in ns:
        for j, x in line_values(kind, n):
            if x < SMALL or j in killed[n]:
                continue
            if algebraic and algebraic_composite(kind, n, j):
                alg += 1
                continue
            surv += 1
            p = mertens / (x.bit_length() * math.log(2))
            E += p
            V += p * (1 - p)
    return E, math.sqrt(V), surv, alg, killed


if __name__ == '__main__':
    mode = sys.argv[1]
    primes = list(primerange(2, Q + 1))
    mertens = 1.0
    for q in primes:
        mertens *= 1 / (1 - 1 / q)
    pred_path = os.path.join(HERE, 'collatz_line_primes_prediction.json')
    if mode in ('posthoc', 'predict2', 'test2'):
        if mode == 'posthoc':
            res = json.load(open(os.path.join(HERE, 'collatz_line_primes_result.json')))
            out = {}
            for name, (kind, ns) in FAMILIES.items():
                E, sd, surv, alg, _ = predict_family(kind, ns, primes, mertens, True)
                obs = res[name]['observed']
                out[name] = {'E_corrected': E, 'sd': sd, 'observed': obs, 'z': (obs - E) / sd,
                             'survivors_removed_as_algebraic': alg}
                print(f'POST-HOC {name}: corrected E = {E:.2f} (sd {sd:.2f}); observed {obs}; z = {(obs - E) / sd:+.2f}; '
                      f'{alg} sieve survivors removed as algebraic composites', flush=True)
            json.dump(out, open(os.path.join(HERE, 'collatz_line_primes_posthoc.json'), 'w'), indent=1)
        elif mode == 'predict2':
            out = {}
            for name, (kind, ns) in FAMILIES2.items():
                E, sd, surv, alg, _ = predict_family(kind, ns, primes, mertens, True)
                E0, sd0, _, _, _ = predict_family(kind, ns, primes, mertens, False)
                out[name] = {'E_corrected': E, 'sd': sd, 'window': [E - 2.5 * sd, E + 2.5 * sd],
                             'E_uncorrected': E0, 'survivors': surv, 'algebraic_removed': alg}
                print(f'{name}: corrected E = {E:.2f} (sd {sd:.2f}), window [{E - 2.5 * sd:.1f}, {E + 2.5 * sd:.1f}]; '
                      f'uncorrected {E0:.2f}', flush=True)
            json.dump(out, open(os.path.join(HERE, 'collatz_line_primes_prediction2.json'), 'w'), indent=1)
        else:
            pred = json.load(open(os.path.join(HERE, 'collatz_line_primes_prediction2.json')))
            out = {}
            for name, (kind, ns) in FAMILIES2.items():
                killed = sieve_family(kind, ns, primes)
                cnt = 0
                for n in ns:
                    for j, x in line_values(kind, n):
                        if x < SMALL or j in killed[n]:
                            continue
                        if isprime(x):
                            assert not algebraic_composite(kind, n, j)
                            cnt += 1
                f = pred[name]
                ok = f['window'][0] <= cnt <= f['window'][1]
                out[name] = {'observed': cnt, **f, 'pass': ok}
                print(f'{name}: observed {cnt}; corrected prediction {f["E_corrected"]:.2f}, window '
                      f'[{f["window"][0]:.1f}, {f["window"][1]:.1f}] -> {"PASS" if ok else "FAIL"}; uncorrected '
                      f'{f["E_uncorrected"]:.2f}', flush=True)
            json.dump(out, open(os.path.join(HERE, 'collatz_line_primes_result2.json'), 'w'), indent=1)
        sys.exit(0)
    if mode == 'predict':
        out = {'Q': Q, 'mertens_factor': mertens, 'families': {}}
        t0 = time.time()
        for name, (kind, ns) in FAMILIES.items():
            killed = sieve_family(kind, ns, primes)
            E = V = 0.0
            exact_small = 0
            survivors = 0
            per_line = {}
            for n in ns:
                e_line = 0.0
                for j, x in line_values(kind, n):
                    if x < SMALL:
                        exact_small += int(isprime(x))       # deterministic, excluded from the random part
                        continue
                    if j in killed[n]:
                        continue
                    survivors += 1
                    p = mertens / (x.bit_length() * math.log(2))
                    E += p
                    V += p * (1 - p)
                    e_line += p
                per_line[n] = e_line
            sd = math.sqrt(V)
            out['families'][name] = {'E': E, 'sd': sd, 'window': [E - 2.5 * sd, E + 2.5 * sd],
                                     'survivors_ge_2^32': survivors, 'exact_small_primes_lt_2^32': exact_small,
                                     'mean_E_per_line': E / len(ns)}
            if name == 'genuine_fermat_e3_12':
                out['families'][name]['E_per_line'] = per_line
            print(f'{name}: E = {E:.2f}, sd = {sd:.2f}, window [{E - 2.5 * sd:.1f}, {E + 2.5 * sd:.1f}], '
                  f'{survivors} survivors >= 2^32, {exact_small} exact small primes ({time.time() - t0:.0f}s)', flush=True)
        json.dump(out, open(pred_path, 'w'), indent=1)
    else:
        pred = json.load(open(pred_path))
        res = {}
        t0 = time.time()
        for name, (kind, ns) in FAMILIES.items():
            killed = sieve_family(kind, ns, primes)
            cnt = 0
            per_line = {}
            for n in ns:
                c = 0
                for j, x in line_values(kind, n):
                    if x < SMALL or j in killed[n]:
                        continue
                    if isprime(x):
                        c += 1
                per_line[n] = c
                cnt += c
            f = pred['families'][name]
            ok = f['window'][0] <= cnt <= f['window'][1]
            res[name] = {'observed': cnt, 'E': f['E'], 'window': f['window'], 'pass': ok}
            if name == 'genuine_fermat_e3_12':
                res[name]['per_line'] = per_line
            print(f'{name}: observed {cnt} probable primes; predicted {f["E"]:.2f}, window '
                  f'[{f["window"][0]:.1f}, {f["window"][1]:.1f}] -> {"PASS" if ok else "FAIL"} ({time.time() - t0:.0f}s)', flush=True)
        json.dump(res, open(os.path.join(HERE, 'collatz_line_primes_result.json'), 'w'), indent=1)
