#!/usr/bin/env python3
"""
Theta/code lattice gate atlas for the 72-dimensional cancellation split.

codex-2026-06-11-P4

Even-unimodular lattice theta series and Type II binary-code weight enumerators
are both modular-form cancellation gates.  In dimension/length 72 they diverge:
the extremal lattice side has a known support object (Nebe's minimum-8 lattice),
while the extremal binary [72,36,16] code remains open.  This script computes the
shared scalar shadows in one table.

Tournament Analysis:
  vertices: scalar gates at dimensions/lengths 24,48,72,96,120.
  observable: (larger killed prefix, then known support status, then smaller
  normalized first-shell log) wins.
  tie Hamiltonian path: listed order.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import comb, log


QMAX = 18
LENGTHS = [24, 48, 72, 96, 120]


def sigma_power(n: int, power: int) -> int:
    return sum(d**power for d in range(1, n + 1) if n % d == 0)


def qadd(a: list[int], b: list[int], scale: int = 1) -> list[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = (a[i] if i < len(a) else 0) + scale * (b[i] if i < len(b) else 0)
    return out


def qmul(a: list[int], b: list[int], qmax: int = QMAX) -> list[int]:
    out = [0] * (qmax + 1)
    for i, ai in enumerate(a[: qmax + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: qmax + 1 - i]):
            if bj:
                out[i + j] += ai * bj
    return out


def qpow(a: list[int], exp: int, qmax: int = QMAX) -> list[int]:
    out = [1] + [0] * qmax
    base = a[: qmax + 1]
    e = exp
    while e:
        if e & 1:
            out = qmul(out, base, qmax)
        e >>= 1
        if e:
            base = qmul(base, base, qmax)
    return out


def e4_series(qmax: int = QMAX) -> list[int]:
    return [1] + [240 * sigma_power(n, 3) for n in range(1, qmax + 1)]


def eta_power_coeffs(power: int, qmax: int = QMAX) -> list[int]:
    coeffs = [1] + [0] * qmax
    for m in range(1, qmax + 1):
        factor = [0] * (qmax + 1)
        for j in range(power + 1):
            if j * m <= qmax:
                factor[j * m] = (-1) ** j * comb(power, j)
        coeffs = qmul(coeffs, factor, qmax)
    return coeffs


def delta_series(qmax: int = QMAX) -> list[int]:
    eta24 = eta_power_coeffs(24, qmax)
    return [0] + eta24[:qmax]


def theta_basis(dim: int, qmax: int = QMAX) -> list[list[int]]:
    if dim % 24:
        raise ValueError("stored atlas uses dimensions divisible by 24")
    s = dim // 24
    e4 = e4_series(qmax)
    delta = delta_series(qmax)
    return [qmul(qpow(e4, 3 * s - 3 * j, qmax), qpow(delta, j, qmax), qmax) for j in range(s + 1)]


def theta_extremal(dim: int, qmax: int = QMAX) -> dict[str, object]:
    s = dim // 24
    basis = theta_basis(dim, qmax)
    theta = basis[0][:]
    coeffs = [1]
    for r in range(1, s + 1):
        c = -theta[r]
        coeffs.append(c)
        theta = qadd(theta, basis[r], c)
    first_exp = next((i for i, c in enumerate(theta) if i and c), None)
    first_coeff = theta[first_exp] if first_exp is not None else None
    return {
        "dim": dim,
        "kill_q": s,
        "min_norm": 2 * (s + 1),
        "coeffs": coeffs,
        "series": theta,
        "first_exp": first_exp,
        "first_coeff": first_coeff,
        "nonnegative_window": all(c >= 0 for c in theta[: qmax + 1]),
    }


def fmul(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def fpow(a: list[Fraction], exp: int) -> list[Fraction]:
    out = [Fraction(1)]
    base = a[:]
    e = exp
    while e:
        if e & 1:
            out = fmul(out, base)
        e >>= 1
        if e:
            base = fmul(base, base)
    return out


def typeii_extremal(length: int) -> dict[str, object]:
    s = length // 24
    m = length // 8
    target_d = 4 * s + 4
    A = [Fraction(1), Fraction(14), Fraction(1)]
    B = [Fraction(0), Fraction(1), Fraction(-4), Fraction(6), Fraction(-4), Fraction(1)]
    basis = [fmul(fpow(A, m - 3 * j), fpow(B, j)) for j in range(m // 3 + 1)]
    W = basis[0][:]
    coeffs = [Fraction(1)]
    for r in range(1, target_d // 4):
        c = -W[r] / basis[r][r]
        coeffs.append(c)
        W = [((W[i] if i < len(W) else 0) + c * (basis[r][i] if i < len(basis[r]) else 0)) for i in range(max(len(W), len(basis[r])))]
    first_weight = next((4 * i for i, c in enumerate(W) if i and c), None)
    first_coeff = W[first_weight // 4] if first_weight is not None else None
    return {
        "length": length,
        "kill_weights": list(range(4, target_d, 4)),
        "d": target_d,
        "coeffs": coeffs,
        "first_weight": first_weight,
        "first_coeff": int(first_coeff) if first_coeff is not None else None,
        "nonnegative": all(c >= 0 for c in W),
    }


def gate_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    known_lattice_support = {24: "unique Leech", 48: "known", 72: "known Nebe"}
    known_code_support = {24: "unique Golay", 48: "known", 72: "open"}
    for n in LENGTHS:
        th = theta_extremal(n)
        code = typeii_extremal(n)
        rows.append(
            {
                "name": f"theta_{n}",
                "family": "lattice theta",
                "n": n,
                "kill": th["kill_q"],
                "support": known_lattice_support.get(n, "scalar only here"),
                "first": int(th["first_coeff"]),
                "norm_log": log(int(th["first_coeff"])) / n,
            }
        )
        rows.append(
            {
                "name": f"typeII_{n}",
                "family": "code weight enumerator",
                "n": n,
                "kill": len(code["kill_weights"]),
                "support": known_code_support.get(n, "scalar only here"),
                "first": int(code["first_coeff"]),
                "norm_log": log(int(code["first_coeff"])) / n,
            }
        )
    return rows


def support_rank(label: str) -> int:
    if "open" in label:
        return 0
    if "scalar only" in label:
        return 1
    return 2


def tournament(rows: list[dict[str, object]]) -> dict[str, object]:
    # Larger killed prefix wins; then known support beats scalar-only/open; then smaller normalized shell wins.
    keys = [(-int(r["kill"]), -support_rank(str(r["support"])), float(r["norm_log"]), i) for i, r in enumerate(rows)]
    n = len(rows)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if keys[i] <= keys[j]:
                adj[i][j] = True
            else:
                adj[j][i] = True
    scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    c3 = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
                    c3 += 1
    order = [rows[i]["name"] for *_, i in sorted(keys)]
    return {"score_hist": dict(sorted(Counter(scores).items())), "c3": c3, "hamiltonian_path": order}


def main() -> None:
    print("Theta/code lattice gate atlas")
    print(f"QMAX={QMAX}; lengths={LENGTHS}")
    print()
    print("[1] Extremal even-unimodular theta gates")
    for n in LENGTHS:
        th = theta_extremal(n)
        coeff_preview = ",".join(str(c) for c in th["coeffs"])
        print(
            f"  dim={n:3d} kill q^1..q^{th['kill_q']} "
            f"min_norm={th['min_norm']:2d} first=q^{th['first_exp']}:{th['first_coeff']} "
            f"coeffs=[{coeff_preview}] nonneg<=q^{QMAX}:{th['nonnegative_window']}"
        )
    print()
    print("[2] Type II code gate comparison")
    for n in LENGTHS:
        code = typeii_extremal(n)
        killed = ",".join(str(w) for w in code["kill_weights"])
        print(
            f"  len={n:3d} kill weights [{killed}] d={code['d']:2d} "
            f"A_d={code['first_coeff']} nonneg:{code['nonnegative']}"
        )
    print()
    print("[3] Dimension/length 72 split")
    th72 = theta_extremal(72)
    code72 = typeii_extremal(72)
    print(
        f"  lattice scalar gate: kill q^1,q^2,q^3 -> q^4 shell {th72['first_coeff']} "
        "(minimum norm 8; Nebe support exists)"
    )
    print(
        f"  code scalar gate: kill weights 4,8,12 -> weight 16 shell {code72['first_coeff']} "
        "(binary support still open)"
    )
    print("  Same modular cancellation level, different support-realization outcome.")
    print()
    print("[4] Tournament Analysis")
    rows = gate_rows()
    fp = tournament(rows)
    print("  Pairwise observable: larger killed prefix, then stronger support status, then lower log(first)/n.")
    print("  score_hist:", fp["score_hist"])
    print("  directed_3_cycles:", fp["c3"])
    print("  Hamiltonian tie path:", " > ".join(fp["hamiltonian_path"]))


if __name__ == "__main__":
    main()
