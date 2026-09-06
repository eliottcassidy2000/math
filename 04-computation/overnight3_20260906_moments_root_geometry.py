"""Exact controls for complete-channel real-rootedness and the OPEN two-row sign.

This is a finite certificate bank, not the proof of the all-exponent theorem.
Every gate stays active under python -O. No numerical root approximations are
used: rational isolating intervals and rational interval Horner certify signs.
"""
from collections import Counter, defaultdict
from functools import reduce
from hashlib import sha256
from math import comb, factorial, gcd
from pathlib import Path
import json
import sympy as s

z = s.Symbol("z")
GATES = 0


def need(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def params(n, A, B, h, r, t):
    N, d = B*h+r, B-A
    a = A*N+B*t
    g = n+N+t
    b, c = g*A-a, g*B-a
    need(0 < A < B and gcd(A, B) == 1, "reduced directions")
    need(0 <= r < B and 0 <= t < A and h >= 1, "canonical residues")
    need(n > 0 and 0 < b < c and gcd(a, g) == 1, "primitive genuine support")
    return N, d, a, b, c, g


def row(p, k):
    n, A, B, h, r, t = p
    N, d, a, b, c, g = params(*p)
    low, high = -(k*t//A), k*h+k*r//B
    weights = {}
    for j in range(low, high+1):
        v = (k*n+d*j, k*N-B*j, k*t+A*j)
        need(min(v) >= 0 and sum(v) == k*g, "row multiplicity")
        need(-a*v[0]+b*v[1]+c*v[2] == 0, "row balance")
        weights[j] = factorial(k*g)//reduce(lambda x, y: x*y, map(factorial, v), 1)
    content = reduce(gcd, weights.values())
    poly = s.Poly(sum((v//content)*z**(j-low) for j, v in weights.items()), z)
    return poly, low, weights


def line(p, q, U, V):
    low, high = -(V//p), U//q
    return {j: comb(U+V+(p-q)*j, V+p*j) for j in range(low, high+1)}


def check_factorization(p, k):
    n, A, B, h, r, t = p
    N, d, a, b, c, g = params(*p)
    ell = k*t//A
    nk, Nk, tk = k*n-d*ell, k*N+B*ell, k*t-A*ell
    need(nk > 0 and 0 <= tk < A, "canonical higher-row shift")
    alpha = line(A, A, k*g-tk, tk)
    beta = line(d, B, Nk, nk)
    need(min(alpha) == 0, "alpha complete lower support")
    L = -min(beta)
    product = {j: alpha.get(j, 0)*beta.get(j, 0)
               for j in range(max(alpha)+1) if alpha.get(j, 0)*beta.get(j, 0)}
    P, low, weights = row(p, k)
    need({j-ell: v for j, v in product.items()} == weights, "full shifted Hadamard identity")
    # Individual factors are checked from their FULL finite supports.
    for seq in (alpha, beta):
        first = min(seq)
        R = s.Poly(sum(v*z**(j-first) for j, v in seq.items()), z)
        need(R.count_roots(-s.oo, 0) == R.sqf_part().degree(), "complete binomial line negative roots")
    return L


def interval_value(poly, left, right):
    lower = upper = s.Rational(0)
    for coeff in poly.all_coeffs():
        products = (lower*left, lower*right, upper*left, upper*right)
        lower, upper = min(products)+coeff, max(products)+coeff
    return lower, upper


def exact_signs(P, Q, low):
    need(s.gcd(P, Q).degree() == 0, "no first/second common root in finite bank")
    R = Q.rem(P)
    certificates = []
    for (left, right), multiplicity in P.intervals(eps=s.Rational(1, 10**8)):
        while right >= 0 and left != right:
            left, right = P.refine_root(left, right, eps=(right-left)/1000)
        need(multiplicity == 1 and right < 0, "finite simple negative first roots")
        for attempt in range(150):
            lower, upper = interval_value(R, left, right)
            if lower > 0 or upper < 0:
                sign = 1 if lower > 0 else -1
                raw_sign = sign*(-1 if low % 2 else 1)
                certificates.append({"interval": [str(left), str(right)],
                                     "remainder_enclosure": [str(lower), str(upper)],
                                     "compressed_sign": sign, "raw_laurent_sign": raw_sign})
                break
            need(left != right, "rational root must have nonzero exact remainder")
            left, right = P.refine_root(left, right, eps=(right-left)/1000)
        else:
            raise RuntimeError("unresolved exact enclosure after bounded refinement")
    need(len(certificates) == P.degree(), "every first root covered")
    return certificates


def literal_rows(p):
    """Independent ordered-word dynamic program in charge and c2-count."""
    n, A, B, h, r, t = p
    N, d, a, b, c, g = params(*p)
    counts = {(0, 0): 1}
    for mass in range(1, 2*g+1):
        nxt = defaultdict(int)
        for (charge, c2_count), multiplicity in counts.items():
            for value, delta in ((-a, 0), (b, 0), (c, 1)):
                nxt[charge+value, c2_count+delta] += multiplicity
        counts = nxt
        if mass in (g, 2*g):
            k = mass//g
            raw = {count: value for (charge, count), value in counts.items() if charge == 0}
            weights = row(p, k)[2]
            need(raw == {k*t+A*j: value for j, value in weights.items()}, "literal word row agreement")


def bank():
    # Explicit finite universe: each listed (A,B), degree h, corner residue,
    # and the first two valid n>=1 under b>0 and primitive-charge filters.
    output = set()
    for A, B in ((1, 2), (2, 3), (3, 5), (4, 7), (5, 9), (5, 12)):
        for h in (2, 3, 4, 6, 8):
            for r, t in {(0, 0), (B-1, 0), (0, A-1), (B-1, A-1)}:
                a = A*(B*h+r)+B*t
                valid = []
                for n in range(1, 500):
                    g = n+B*h+r+t
                    if n*A-(B-A)*t > 0 and gcd(a, g) == 1:
                        valid.append(n)
                    if len(valid) == 2:
                        break
                need(len(valid) == 2, "finite bank completion")
                output.update((n, A, B, h, r, t) for n in valid)
    # Explicit truncation witness and the large-coefficient midpoint trap.
    output.update(((2, 1, 2, 2, 1, 0), (16, 5, 9, 3, 6, 4)))
    return sorted(output)


def main():
    source = Path(__file__)
    print("source_sha256", sha256(source.read_bytes()).hexdigest(), flush=True)
    print("status FINITE-EXACT; analytical theorem in companion report", flush=True)
    hostile = s.Poly(21+20*z+5*z*z, z)
    need(s.discriminant(hostile) == -20, "truncation hostile")
    witness = (2, 1, 2, 2, 1, 0)
    need(row(witness, 1)[2] == {0: 21, 1: 140, 2: 105}, "hostile genuine multinomial row")
    print("truncated_beta 21+20*z+5*z^2; discriminant -20", flush=True)
    print("full_beta 1+8*z+21*z^2+20*z^3+5*z^4; P 21+140*z+105*z^2", flush=True)
    # g=1 is NOT the first support-return mass for (-1,1,2).
    scalar_counts = {0: 1}
    scalar_returns = []
    for mass in range(1, 7):
        nxt = defaultdict(int)
        for charge, count in scalar_counts.items():
            for step in (-1, 1, 2):
                nxt[charge+step] += count
        scalar_counts = nxt
        scalar_returns.append(scalar_counts.get(0, 0))
    need(scalar_returns == [0, 2, 3, 6, 20, 35], "nonsemigroup first-return mass firewall")
    print("nonsemigroup_hostile (-1,1,2), g=1, actual_first_support_mass=2; CT rows", scalar_returns, flush=True)
    cases, signs, carry_counts, degrees = [], Counter(), Counter(), Counter()
    universe = bank()
    print("universe_supports", len(universe), flush=True)
    for index, p in enumerate(universe, 1):
        P, _, first = row(p, 1)
        Q, low, second = row(p, 2)
        need(P.count_roots(-s.oo, 0) == P.degree(), "first row real negative roots")
        need(Q.count_roots(-s.oo, 0) == Q.degree(), "second row real negative roots")
        certs = exact_signs(P, Q, low)
        signs.update(c["raw_laurent_sign"] for c in certs)
        n, A, B, h, r, t = p
        carry = (2*t//A, 2*r//B)
        carry_counts[carry] += 1
        degrees[h] += 1
        # Check full factor sequences for one representative per direction
        # and carry signature; the identity itself is checked for every case.
        if h == 2:
            check_factorization(p, 1)
            check_factorization(p, 2)
        else:
            for k, weights in ((1, first), (2, second)):
                N, d, a, b, c, g = params(*p)
                for j, weight in weights.items():
                    need(comb(k*g, k*t+A*j)*comb(k*(n+N)-A*j, k*n+d*j) == weight,
                         "unshifted full weight identity")
        cases.append({"parameters_n_A_B_h_r_t": list(p),
                      "support": list(params(*p)[2:5]), "g": params(*p)[5],
                      "first_coefficients_ascending": [int(P.nth(j)) for j in range(P.degree()+1)],
                      "second_coefficients_ascending": [int(Q.nth(j)) for j in range(Q.degree()+1)],
                      "second_laurent_low": low, "certificates": certs})
        if index % 40 == 0:
            print("completed", index, flush=True)
    for p in ((1, 1, 2, 2, 0, 0), witness, (1, 2, 3, 1, 2, 1)):
        literal_rows(p)
    print("literal_ordered_word_controls 3 supports, 6 raw rows", flush=True)
    print("first_degrees", sorted(degrees.items()), flush=True)
    print("carry_signature_counts", sorted(carry_counts.items()), flush=True)
    print("raw_second_sign_counts_at_first_roots", sorted(signs.items()), flush=True)
    print("all_roots_simple_in_finite_bank", True, flush=True)
    print("OPEN all-exponent second-row signs and pairwise root separation", flush=True)
    payload = {"scope": "FINITE-EXACT", "source_sha256": sha256(source.read_bytes()).hexdigest(),
               "cases": cases, "active_gates": GATES}
    certificate = source.with_name(source.stem+"_certificates.json")
    certificate.write_text(json.dumps(payload, sort_keys=True, indent=2)+"\n", encoding="utf-8", newline="\n")
    print("certificate_sha256", sha256(certificate.read_bytes()).hexdigest(), flush=True)
    print("active_gates", GATES, flush=True)
    print("PASS", flush=True)


if __name__ == "__main__":
    main()
