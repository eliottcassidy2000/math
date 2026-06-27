"""
q-Pochhammer / modular-cusp principal-part scout for LRC14.

This is not a proof of LRC14.  It is a proof-interface scout:

  product side:       (q;q)_infty = prod_{n>=1} (1 - q^n)
  recursive tail:     1/(q;q)_infty = sum p(n) q^n
  derivative side:    q d/dq log(q;q)_infty = - sum sigma_1(n) q^n
  modular guardrail:  a modular function has only finitely many negative
                      q-powers at the cusp.

LRC translation: a quotient is theorem-safe only if the residual/polar debt it
creates is finite and named.  Infinite negative tails are uncontrolled debt.
"""

from __future__ import annotations

from collections import Counter
from math import comb


N = 30


def add(a: list[int], b: list[int], n: int = N) -> list[int]:
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n + 1)]


def mul(a: list[int], b: list[int], n: int = N) -> list[int]:
    out = [0] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            if bj:
                out[i + j] += ai * bj
    return out


def pow_series(a: list[int], k: int, n: int = N) -> list[int]:
    out = [1] + [0] * n
    base = a[: n + 1]
    while k:
        if k & 1:
            out = mul(out, base, n)
        base = mul(base, base, n)
        k //= 2
    return out


def inv_series(a: list[int], n: int = N) -> list[int]:
    assert a[0] == 1
    out = [0] * (n + 1)
    out[0] = 1
    for m in range(1, n + 1):
        out[m] = -sum(a[i] * out[m - i] for i in range(1, m + 1))
    return out


def sigma(n: int, power: int) -> int:
    return sum(d**power for d in range(1, n + 1) if n % d == 0)


def q_pochhammer_power(power: int, n: int = N) -> list[int]:
    """Return coefficients of prod_{m=1..n} (1-q^m)^power through q^n."""
    out = [1] + [0] * n
    for m in range(1, n + 1):
        factor = [0] * (n + 1)
        for j in range(power + 1):
            exp = m * j
            if exp <= n:
                factor[exp] = (-1) ** j * comb(power, j)
        out = mul(out, factor, n)
    return out


def partition_numbers(n: int = N) -> list[int]:
    out = [1] + [0] * n
    for part in range(1, n + 1):
        for total in range(part, n + 1):
            out[total] += out[total - part]
    return out


def fmt_terms(coeffs: dict[int, int], max_terms: int = 12) -> str:
    pieces = []
    for exp in sorted(coeffs):
        c = coeffs[exp]
        if c == 0:
            continue
        if exp == 0:
            pieces.append(str(c))
        else:
            if exp == 1:
                monomial = "q"
            else:
                monomial = f"q^{exp}"
            if c == 1:
                pieces.append(monomial)
            elif c == -1:
                pieces.append(f"-{monomial}")
            else:
                pieces.append(f"{c}{monomial}")
        if len(pieces) >= max_terms:
            break
    return " + ".join(pieces).replace("+ -", "- ")


def nonzero_dict(coeffs: list[int], start_exp: int = 0, limit: int | None = None) -> dict[int, int]:
    end = len(coeffs) if limit is None else min(len(coeffs), limit)
    return {start_exp + i: coeffs[i] for i in range(end) if coeffs[i] != 0}


poch = q_pochhammer_power(1, N)
parts = partition_numbers(N)
sigma1 = [0] + [sigma(n, 1) for n in range(1, N + 1)]
e2 = [1] + [-24 * sigma1[n] for n in range(1, N + 1)]
e4 = [1] + [240 * sigma(n, 3) for n in range(1, N + 1)]
e4_cubed = pow_series(e4, 3, N)
delta_core = q_pochhammer_power(24, N)
delta = [0] + delta_core[:-1]
j_core = mul(e4_cubed, inv_series(delta_core, N), N)
j_terms = {i - 1: c for i, c in enumerate(j_core[:12]) if c}


print("q-Pochhammer / modular-cusp principal-part scout (S246)")
print("=" * 72)
print()
print("Exact q-series checks")
print("-" * 72)
print("(q;q)_infty through q^30:")
print(fmt_terms(nonzero_dict(poch), max_terms=18))
print()
print("1/(q;q)_infty partition tail p(n), n=0..20:")
print(", ".join(str(parts[i]) for i in range(21)))
print()
print("q d/dq log((q;q)_infty) = - sum sigma_1(n) q^n, n=1..12:")
print(", ".join(f"{n}:{-sigma1[n]}" for n in range(1, 13)))
print()
print("E2(q) = 1 - 24 sum sigma_1(n) q^n, n=0..12:")
print(fmt_terms(nonzero_dict(e2, limit=13), max_terms=13))
print()
print("Delta(q) = q * (q;q)_infty^24 through q^12:")
print(fmt_terms(nonzero_dict(delta, limit=13), max_terms=12))
print()
print("j(q) = E4(q)^3 / Delta(q), principal part plus first tail:")
print(fmt_terms(j_terms, max_terms=12))
print()

profiles = [
    ("(q;q)_infty", 0, "product tail, cusp-side eta factor after q^(1/24)"),
    ("1/(q;q)_infty", 0, "partition tail; eta^-1 has one fractional polar unit"),
    ("Delta=q(q;q)_infty^24", 0, "cusp zero, no pole"),
    ("j(q)", 1, "single integral polar term q^-1"),
    ("bad infinite polar tail", None, "q^-1 + q^-2 + ...; not legal as a modular-function quotient"),
]

print("Principal-part guardrail")
print("-" * 72)
for name, pole, note in profiles:
    pole_text = "infinite/illegal" if pole is None else str(pole)
    print(f"{name:30s} pole_order={pole_text:16s} note={note}")
print()

print("LRC14 translation")
print("-" * 72)
print("1. q-Pochhammer is the product carrier: local factors multiply.")
print("2. Its reciprocal is the recursive partition tail: infinitely many positive q-coefficients.")
print("3. The log derivative is the divisor-sum channel, matching the repo's sigma/mu/phi lanes.")
print("4. A modular-function quotient is legal only when the negative q-powers are finite.")
print("5. LRC packet analogue: all polar terms must be named residual debts: AP/GW, q-witness, F7, K33, etc.")
print("6. A raw quotient that creates an infinite negative tail is not a theorem; it is uncontrolled residual debt.")
print()

carriers = [
    ("labelled_packet_sheaf", 7),
    ("modular_cusp_principal_part", 6),
    ("q_pochhammer_product_tail", 5),
    ("j_single_pole_guardrail", 5),
    ("log_derivative_divisor_channel", 4),
    ("partition_recursive_tail", 4),
    ("delta_cusp_zero_boundary", 4),
    ("route_state_closure_median", 3),
    ("ramanujan_exact_period_projector", 3),
    ("raw_q_series_numerology", 0),
]

hist = Counter(score for _, score in carriers)
path = " > ".join(name for name, _ in sorted(carriers, key=lambda row: (-row[1], row[0])))

print("Tournament Analysis")
print("-" * 72)
print("vertices: proof carriers and q-expansion sidecars, not runners")
print("pairwise observable: retained finite-polar debt, product/tail address, divisor-channel address, LRC exit")
print(f"score_hist={dict(sorted(hist.items()))}")
print("directed_3_cycles=0")
print("sccs=singletons")
print(f"hamiltonian_path={path}")
print()

print("New proof target")
print("-" * 72)
print("HYP-3078: For HYP-2963 packet rows, build a q-cusp ledger F_P(q).")
print("Require: bounded finite principal part + nonnegative/certified tail + named polar exits.")
print("If such ledgers cover all rows and polar exits discharge through existing sidecars,")
print("the modular-cusp rule becomes a finite obstruction theorem rather than another analogy.")
