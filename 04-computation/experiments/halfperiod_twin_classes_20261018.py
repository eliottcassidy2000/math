#!/usr/bin/env python3
"""Half-period twin classes among generated binomial sums (S43, generalizing THM-4669/THM-4670).

Generate u(n) = sum_k prod of 3 or 4 binomial factors from a library L (as in clausen_twin_search_20261018.py, plus a
few more factors). Two sums u, v are 'half-period twins' if u((p-1)/2) == eps(p) v((p-1)/2) (mod p^2) for all primes
5 <= p <= PMAX with eps in {+1, (-1)^n}, but u and v are different sequences (differ at some n <= 12). We screen on a
signature (residues at the first primes), confirm on all primes <= PMAX, and print the nontrivial classes, flagging
sequences that are identically equal (rewritings) and classes whose values vanish (trivially congruent to 0).
Usage: python halfperiod_twin_classes_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import os, sys, json, itertools
from math import comb


def C(a, b):
    if b < 0 or a < 0 or b > a:
        return 0
    return comb(a, b)


LIB = {
    'C(n,k)': lambda n, k: C(n, k),
    'C(n+k,k)': lambda n, k: C(n + k, k),
    'C(2k,k)': lambda n, k: C(2 * k, k),
    'C(2k,n)': lambda n, k: C(2 * k, n),
    'C(2n-2k,n-k)': lambda n, k: C(2 * n - 2 * k, n - k),
    'C(k,n-k)': lambda n, k: C(k, n - k),
    'C(n,2k)': lambda n, k: C(n, 2 * k),
    'C(2k,n-k)': lambda n, k: C(2 * k, n - k),
    'C(n+k,2k)': lambda n, k: C(n + k, 2 * k),
    'C(2n,n+k)': lambda n, k: C(2 * n, n + k),
}


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    primes = [p for p in range(5, PMAX + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    screen = primes[:10]
    names = list(LIB)
    seqs = {}
    for r in (3, 4):
        for combo in itertools.combinations_with_replacement(names, r):
            fs = [LIB[c] for c in combo]
            def u(n, fs=fs):
                tot = 0
                for k in range(n + 1):
                    v = 1
                    for f in fs:
                        v *= f(n, k)
                        if v == 0:
                            break
                    tot += v
                return tot
            seqs[' '.join(combo)] = u
    # small-n values to detect identical sequences
    small = {name: tuple(u(n) for n in range(13)) for name, u in seqs.items()}
    # merge identical sequences
    reps = {}
    for name, sv in small.items():
        reps.setdefault(sv, []).append(name)
    distinct = {names_[0]: (sv, names_) for sv, names_ in reps.items() if any(sv)}
    print('generated %d sums, %d distinct nonzero sequences' % (len(seqs), len(distinct)), flush=True)
    # signatures at the half-period, normalized for the two sign conventions
    sig = {}
    for name in distinct:
        u = seqs[name]
        vals = []
        for p in screen:
            n = (p - 1) // 2
            vals.append(u(n) % (p * p))
        sig[name] = vals
    # group: key = residues with eps = +1 and with eps = (-1)^n applied
    def key(vals, eps):
        out = []
        for v, p in zip(vals, screen):
            n = (p - 1) // 2
            e = 1 if eps == '+' else (-1) ** n
            out.append((e * v) % (p * p))
        return tuple(out)
    buckets = {}
    for name, vals in sig.items():
        buckets.setdefault(key(vals, '+'), []).append((name, '+'))
        buckets.setdefault(key(vals, '-'), []).append((name, '(-1)^n'))
    classes = []
    seen = set()
    for k, members in buckets.items():
        nm = sorted(set(m for m, e in members))
        if len(nm) < 2 or all(x == 0 for x in k):
            continue
        # confirm pairwise against the first member on all primes
        base, beps = members[0]
        conf = []
        for m, e in members[1:]:
            if m == base:
                continue
            ok = True
            for p in primes:
                n = (p - 1) // 2
                eb = 1 if beps == '+' else (-1) ** n
                em = 1 if e == '+' else (-1) ** n
                if (eb * seqs[base](n) - em * seqs[m](n)) % (p * p):
                    ok = False
                    break
            if ok:
                conf.append((m, e))
        if conf:
            cl = tuple(sorted([base] + [m for m, e in conf]))
            if cl in seen:
                continue
            seen.add(cl)
            classes.append(dict(base=base, base_eps=beps, members=[dict(seq=m, eps=e) for m, e in conf],
                                aliases={x: distinct[x][1] for x in cl}))
    print('nontrivial half-period twin classes (confirmed on all primes 5..%d):' % PMAX)
    for c in classes:
        print('  [%s] %s' % (c['base_eps'], c['base']))
        for m in c['members']:
            print('      ~ [%s] %s' % (m['eps'], m['seq']))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(classes, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'halfperiod_twin_classes_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
