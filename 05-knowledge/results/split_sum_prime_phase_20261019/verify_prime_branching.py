#!/usr/bin/env python3
"""Exact finite controls for the split-sum, phase, and binary-balancing theorems.

Standard library only. No assert statements: every check survives python -O.
The computations support the accompanying proofs; they do not establish
infinitude or exhaust any open prime-generating process.
"""
from collections import Counter
from itertools import combinations_with_replacement, product
from math import isqrt, prod
from pathlib import Path
import hashlib
import json

CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def prime(n):
    return type(n) is int and n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def primes_to(n):
    sieve = bytearray(b"\1") * (n + 1)
    sieve[:2] = b"\0\0"
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            sieve[p*p::p] = b"\0" * ((n-p*p)//p + 1)
    return [p for p in range(2, n + 1) if sieve[p]]


def legendre(a, q):
    z = pow(a % q, (q-1)//2, q)
    return -1 if z == q-1 else z


def divisors(P):
    D = [1]
    for p in P:
        D += [d*p for d in D]
    return D


def sign_class(a, q):
    a %= q
    return min(a, q-a)


def signed_counts(P, q):
    if type(q) is not int or not prime(q) or q == 2:
        raise ValueError("q must be an odd prime")
    if any(type(p) is not int or not prime(p) or p == q for p in P) or len(set(P)) != len(P):
        raise ValueError("P must contain distinct primes other than q")
    D = {1: 1}
    for p in P:
        E = Counter(D)
        for x, count in D.items():
            E[sign_class(x*p, q)] += count
        D = dict(E)
    return [D.get(x, 0) for x in range(1, (q+1)//2)]


def residue_witnesses(P, q):
    D = {1: 1}
    for p in P:
        for x, d in tuple(D.items()):
            D.setdefault(x*p % q, d*p)
    return D


def split_witness(P, q):
    N = prod(P)
    return next((d for x, d in residue_witnesses(P, q).items() if (x*x+N) % q == 0), None)


def subset_counts_cyclic(logs, h):
    C = [0]*h
    C[0] = 1
    for t in logs:
        C = [C[x] + C[(x-t) % h] for x in range(h)]
    return C


def binary_balance_predicate(logs, h):
    if h & (h-1):
        return False
    s = h.bit_length()-1
    valuations = {(t & -t).bit_length()-1 for t in logs if t % h}
    return set(range(s)) <= valuations


def primitive_root(q):
    n = q-1
    factors = [p for p in primes_to(n) if n % p == 0]
    return next(g for g in range(2, q) if all(pow(g, n//p, q) != 1 for p in factors))


def logarithms(q, g):
    logs, x = {}, 1
    for j in range(q-1):
        logs[x] = j
        x = x*g % q
    check(len(logs) == q-1 and x == 1, "primitive-root permutation")
    return logs


def main():
    out = {"scope": "Finite controls, not resolution of a prime-generation conjecture"}
    cyclic_cases = 0
    for h in range(1, 17):
        for k in range(6):
            for ts in combinations_with_replacement(range(h), k):
                C = subset_counts_cyclic(ts, h)
                check((len(set(C)) == 1) == binary_balance_predicate(ts, h), "cyclic balancing iff")
                cyclic_cases += 1
    out["cyclic_multisets_h_le_16_length_le_5"] = cyclic_cases

    # Independent literal divisor enumeration versus quotient dynamic programming.
    count_cases = 0
    base = [2, 3, 5, 7, 11, 13]
    for q in primes_to(101)[1:]:
        pool = [p for p in base if p != q]
        for bits in product((0, 1), repeat=len(pool)):
            P = [p for p, b in zip(pool, bits) if b]
            N, D = prod(P), divisors(P)
            C = signed_counts(P, q)
            literal = Counter(sign_class(d, q) for d in D)
            check(C == [literal[x] for x in range(1, (q+1)//2)], "signed DP versus literal")
            Z = sum((d + N//d) % q == 0 for d in D)
            roots = [x for x in range(1, q) if (x*x+N) % q == 0]
            predicted = 0 if not roots else C[sign_class(roots[0], q)-1]
            check(Z == predicted, "phase plus root-class multiplicity")
            if P:
                check(Z % 2 == 0, "complement-pair parity")
            if min(C) > 0:
                check((Z > 0) == (legendre(-N, q) == 1), "saturated criterion")
            count_cases += 1
    out["literal_divisor_cases"] = count_cases

    # Full primorial coverage is CITED in the proof; these are bounded controls.
    coverage = []
    for q in primes_to(499)[1:]:
        P = primes_to(q-1)
        W = residue_witnesses(P, q)
        check(len({sign_class(x, q) for x in W}) == (q-1)//2, "primorial signed support")
        missing = sorted(set(range(1, q))-W.keys())
        check(missing == ([4] if q in (5, 7) else []), "Booker-Pomerance small control")
        coverage.append({"q": q, "missing_unsigned": missing})
    out["primorial_coverage"] = coverage

    # Exact Chua greedy steps: test every unused prime below the selected one,
    # with subset-residue witnesses, not just a quadratic-residue assumption.
    P, N, rows, ps = [], 1, [], primes_to(5000)
    for index in range(80):
        least_missing = next(p for p in ps if p not in P)
        candidate_rows = []
        selected = None
        for q in ps:
            if q in P:
                continue
            if q != 2 and legendre(-N, q) != 1:
                candidate_rows.append([q, "quadratic obstruction"])
                continue
            d = split_witness(P, q)
            if d is not None:
                selected = q
                break
            candidate_rows.append([q, "root class absent"])
        check(selected is not None, "finite greedy search found next prime")
        check(N % d == 0 and (d+N//d) % selected == 0, "actual split witness")
        if least_missing != 2:
            eligible = legendre(-N, least_missing) == 1
            check((selected == least_missing) == eligible, "least-missing equivalence")
        else:
            check(selected == 2, "initial parity step")
        rows.append({"index": index+1, "selected": selected, "least_missing": least_missing,
                     "source": str(N), "divisor": str(d), "rejections": candidate_rows})
        P.append(selected)
        N *= selected
    for before, after in zip(rows, rows[1:]):
        q, r = before["least_missing"], before["selected"]
        if q != 2 and r != q:
            check((after["selected"] == q) == (legendre(r, q) == -1), "one-step nonresidue trigger")
    out["greedy_split_sum_steps"] = rows

    # Same complete kernel, same factor count, same smallest absent q,
    # opposite candidate eligibility. No randomness or Fourier rounding.
    hostile = []
    for q, A, B in [(5, [2,3,7], [2,3,11]), (7, [2,3,5,11], [2,3,5,17])]:
        CA, CB = signed_counts(A, q), signed_counts(B, q)
        check(CA == CB and min(CA) > 0, "identical complete signed kernels")
        ZA = sum((d+prod(A)//d) % q == 0 for d in divisors(A))
        ZB = sum((d+prod(B)//d) % q == 0 for d in divisors(B))
        check(ZA == 0 and ZB > 0, "kernel forgets target phase")
        hostile.append({"q":q, "A":A, "B":B, "kernel":CA, "Z_A":ZA, "Z_B":ZB})
    check(legendre(-2,11) == 1 and split_witness([2],11) is None, "Legendre alone without saturation")
    out["same_kernel_phase_hostiles"] = hostile

    fermat = []
    for q in (3,5,17,257,65537):
        h, g = (q-1)//2, primitive_root(q)
        s, logs = h.bit_length()-1, logarithms(q,g)
        layers = [[] for _ in range(s+1)]
        P = primes_to(q-1)
        for p in P:
            t = logs[p] % h
            v = (t & -t).bit_length()-1 if t else s
            layers[v].append(p)
        missing = [v for v in range(s) if not layers[v]]
        row = {"q":q,"generator":g,"layer_counts":[len(x) for x in layers],
               "missing_layers":missing,"layer_first_primes":[x[:4] for x in layers]}
        if q <= 257:
            check(len(set(signed_counts(P,q))) == 1, "small Fermat primorial balance")
        if q == 65537:
            check(missing == [14], "65537 unique missing top layer")
            top = [x for x in range(1,q) if logs[x] % h == h//2]
            check(top == [256,65281], "missing layer residues")
            check(256 == 2**8 and 65281 == 97*673, "composite top-layer representatives")
            repair = next(p for p in primes_to(327941) if p % q in top)
            check(repair == 327941, "least repair prime")
            basis = [layers[v][0] for v in range(14)] + [repair]
            check(len(basis) == s and len(set(basis)) == s, "minimal binary basis")
            C0 = signed_counts(basis,q)
            check(set(C0) == {1}, "one divisor per signed class")
            check(split_witness(basis,q) is None, "minimal balanced basis cannot generate q")
            extra = next(p for p in layers[0] if p not in basis)
            generating = basis + [extra]
            C1 = signed_counts(generating,q)
            check(set(C1) == {2}, "phase-corrected balanced basis")
            N1 = prod(generating)
            hits = sorted(d for d in divisors(generating) if (d+N1//d) % q == 0)
            check(len(hits) == 2 and hits[0]*hits[1] == N1, "unique unordered split")
            row.update({"least_repair_prime":repair,"minimal_basis":basis,
                        "phase_correcting_prime":extra,"generating_source":str(N1),
                        "split_divisors":[str(x) for x in hits],
                        "split_quotient":str((hits[0]+hits[1])//q)})
        fermat.append(row)
    out["Fermat_layer_certificates"] = fermat

    for P,q in [([2],True),([2],9),([True],5),([2,2],5),([5],5),([4],5)]:
        try:
            signed_counts(P,q)
        except ValueError:
            check(True,"malformed input rejected")
        else:
            check(False,"malformed input accepted")
    out["always_active_checks"] = CHECKS
    payload = json.dumps(out, indent=2, sort_keys=True)+"\n"
    Path(__file__).with_name("prime_branching_receipt.json").write_text(payload, encoding="utf-8")
    print(f"PASS: {CHECKS} always-active exact checks")
    print("receipt sha256:", hashlib.sha256(payload.encode()).hexdigest())
    print("first 25 greedy terms:", [r["selected"] for r in rows[:25]])
    print("65537 certificate:",json.dumps(fermat[-1],sort_keys=True))


if __name__ == "__main__":
    main()
