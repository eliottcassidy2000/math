#!/usr/bin/env python3
"""Audit check (a) of THM-4710's Lemma with the independent engine vir_pbw.py.

Claim: for n >= 1 and a in V^c with top PBW length r, L_n a has length <= r and its length-r component is
D_n sigma_r(a), where D_n is the derivation x_j -> (j+n+2) x_(j-n) (j >= n), 0 (j < n), x_j <-> L_(-j-2); no c.

Tests:
  1. engine sanity: the commutator [L_m, L_n] = (m-n) L_(m+n) + (c/12)(m^3-m) delta on random vectors (m, n in -4..4);
  2. for random vectors (random PBW monomials, random integer coefficients) at weights 8..16, n = 1..5, and three
     values of c (including a symbolic-looking pair to test c-independence): top-length part of L_n a equals
     D_n(top part of a), and the length of L_n a never exceeds r;
  3. the central term: the length-(r-1) part of L_n a does depend on c for suitable a (so c really is present one level
     down), while the length-r part never does.
Output: lemma_check.out (printed) and lemma_check.json.
"""
import random, json, os
from fractions import Fraction
from vir_pbw import VirVacuum, partitions_ge2

HERE = os.path.dirname(os.path.abspath(__file__))


def length_part(vec, r):
    return {m: v for m, v in vec.items() if len(m) == r}


def top_length(vec):
    return max((len(m) for m in vec), default=-1)


def Dn_pbw(vec, n):
    """D_n on symbols written as PBW-monomial multisets: part k (x_(k-2)) with k - n >= 2 -> (k + n) * part (k - n)."""
    out = {}
    for m, a in vec.items():
        lst = list(m)
        for i, k in enumerate(lst):
            if k - n >= 2:
                new = tuple(sorted(lst[:i] + [k - n] + lst[i + 1:], reverse=True))
                out[new] = out.get(new, 0) + a * (k + n)
    return {m: v for m, v in out.items() if v}


def sub(u, v):
    out = dict(u)
    for m, a in v.items():
        out[m] = out.get(m, 0) - a
    return {m: x for m, x in out.items() if x}


def main():
    rng = random.Random(20261018)
    cs = [Fraction(7, 10), Fraction(-22, 5), Fraction(1234, 577)]
    engines = [VirVacuum(c) for c in cs]
    report = {'commutator_tests': 0, 'lemma_tests': 0, 'lemma_failures': 0, 'length_increase': 0,
              'top_c_dependence': 0, 'lower_c_dependence_seen': 0}
    # 1. commutator sanity
    for V in engines:
        for trial in range(40):
            W = rng.randint(4, 12)
            monos = list(partitions_ge2(W))
            a = {m: rng.randint(-5, 5) for m in rng.sample(monos, min(len(monos), 4))}
            a = {m: x for m, x in a.items() if x}
            m_, n_ = rng.randint(-4, 4), rng.randint(-4, 4)
            lhs = sub(V.apply(m_, V.apply(n_, a)), V.apply(n_, V.apply(m_, a)))
            rhs = {k: (m_ - n_) * x for k, x in V.apply(m_ + n_, a).items()}
            if m_ + n_ == 0:
                for k, x in a.items():
                    rhs[k] = rhs.get(k, 0) + V.c12 * (m_ ** 3 - m_) * x
            rhs = {k: x for k, x in rhs.items() if x}
            assert lhs == rhs, (m_, n_, a)
            report['commutator_tests'] += 1
    # 2. lemma on random vectors
    for trial in range(300):
        W = rng.randint(8, 16)
        monos = list(partitions_ge2(W))
        k = min(len(monos), rng.randint(2, 8))
        chosen = rng.sample(monos, k)
        # make sure at least some monomials contain small parts (so central terms / L_-1 creations occur)
        a = {m: rng.randint(-9, 9) for m in chosen}
        a = {m: x for m, x in a.items() if x}
        if not a:
            continue
        r = top_length(a)
        sig = length_part(a, r)
        for n in range(1, 6):
            tops = []
            lowers = []
            for V in engines:
                b = V.apply(n, a)
                rb = top_length(b)
                if rb > r:
                    report['length_increase'] += 1
                tops.append(length_part(b, r))
                lowers.append(length_part(b, r - 1))
            pred = Dn_pbw(sig, n)
            report['lemma_tests'] += 1
            if any(t != pred for t in tops):
                report['lemma_failures'] += 1
                print('FAIL', a, n, tops[0], pred)
            if not all(t == tops[0] for t in tops):
                report['top_c_dependence'] += 1
            if not all(l == lowers[0] for l in lowers):
                report['lower_c_dependence_seen'] += 1
    # 3. explicit central-term examples: L_n L_-n L_-2^s |0>
    examples = []
    for n in (2, 3, 4):
        for s in (0, 1, 2):
            mono = tuple(sorted([n] + [2] * s, reverse=True)) if n >= 2 else None
            a = {mono: 1}
            res = []
            for V in engines:
                b = V.apply(n, a)
                res.append({str(m): str(v) for m, v in b.items()})
            examples.append({'n': n, 'vector': str(mono), 'L_n a at c=7/10, -22/5, 1234/577': res})
    report['central_examples'] = examples
    # 4. stress test: bigger vectors (weights 10..20, 12 monomials), n = 1..6, two more c values
    rng2 = random.Random(5)
    Vs = [VirVacuum(Fraction(7, 10)), VirVacuum(Fraction(-3, 5))]
    nz = tot = bad = 0
    for t in range(200):
        W = rng2.randint(10, 20)
        monos = list(partitions_ge2(W))
        a = {m: rng2.randint(1, 9) for m in rng2.sample(monos, min(len(monos), 12))}
        r = top_length(a)
        sig = length_part(a, r)
        for n in range(1, 7):
            pred = Dn_pbw(sig, n)
            for V in Vs:
                b = V.apply(n, a)
                tot += 1
                nz += bool(pred)
                if length_part(b, r) != pred or top_length(b) > r:
                    bad += 1
    report['stress_tests'] = tot
    report['stress_nonzero_predictions'] = nz
    report['stress_failures'] = bad
    print(json.dumps({k: v for k, v in report.items() if k != 'central_examples'}, indent=1))
    for e in examples[:4]:
        print(e)
    ok = (report['lemma_failures'] == 0 and report['length_increase'] == 0 and report['top_c_dependence'] == 0
          and report['stress_failures'] == 0)
    print('LEMMA CHECK', 'PASSED' if ok else 'FAILED')
    report['passed'] = ok
    json.dump(report, open(os.path.join(HERE, 'lemma_check.json'), 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
