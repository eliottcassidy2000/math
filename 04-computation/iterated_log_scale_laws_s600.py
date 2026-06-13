#!/usr/bin/env python3
"""S600: an atlas for iterated-log inequality shapes.

This is not a theorem prover.  It records the scale-calculus behind many
Tao-style logarithmic improvements:

  product over scales (1 - saving_j) ~= exp(-sum saving_j).

Harmonic savings on ordinary scales give log; harmonic savings after one
scale-compression give log log; after two compressions give log log log.
The printed table compares the classic lonely-runner Tao-shape surplus with a
few candidate "own" surplus templates that would require additional structure.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, log


def L(x: float, depth: int) -> float:
    for _ in range(depth):
        x = log(x)
    return x


def safe_L(x: float, depth: int) -> float:
    y = x
    for _ in range(depth):
        if y <= 1:
            return float("nan")
        y = log(y)
    return y


def tao_shape(k: float) -> float:
    return log(k) / (k * k * log(log(k)) ** 2)


def third_log_tax_shape(k: float) -> float:
    return tao_shape(k) / max(log(log(log(k))), 1.0)


def third_log_dividend_shape(k: float) -> float:
    return tao_shape(k) * max(log(log(log(k))), 1.0)


def entropy_rank_shape(k: float, rank: float) -> float:
    """A rank-discounted template: rank extra orthogonal scale channels."""
    return log(k) / (k * k * (log(log(k)) + rank * max(log(log(log(k))), 0.0)) ** 2)


def fmt_float(x: float) -> str:
    if x != x:
        return "       -"
    return f"{x:8.3f}"


def scale_product(start: int, stop: int, alpha: float, mode: str) -> float:
    """Survival after multiplying (1-saving_j) over a scale ladder."""
    out = 1.0
    for j in range(start, stop + 1):
        if mode == "log":
            saving = alpha / j
        elif mode == "loglog":
            saving = alpha / (j * log(j))
        elif mode == "logloglog":
            saving = alpha / (j * log(j) * log(log(j)))
        else:
            raise ValueError(mode)
        saving = min(saving, 0.49)
        out *= 1.0 - saving
    return out


@dataclass(frozen=True)
class Template:
    name: str
    formula: str
    meaning: str


def main() -> None:
    print("S600 iterated-log scale-law atlas")
    print()
    print("ABSTRACTION")
    print("  If a proof has savings eps_j over scale levels, then residual")
    print("  survival is product_j(1-eps_j) <= exp(-sum_j eps_j).")
    print("  log/loglog/logloglog are bookkeeping for how many independent")
    print("  scale currencies survive after one, two, or three compressions.")
    print()

    templates = [
        Template(
            "scale_harmonic",
            "prod_{j<=J}(1-a/j) <= J^{-a+o(1)}",
            "ordinary dyadic scale dividend",
        ),
        Template(
            "compressed_harmonic",
            "prod_{j<=J}(1-a/(j log j)) <= (log J)^{-a+o(1)}",
            "one compressed scale ladder; loglog when J is log N",
        ),
        Template(
            "meta_compressed_harmonic",
            "prod_{j<=J}(1-a/(j log j loglog j)) <= (loglog J)^{-a+o(1)}",
            "two compressed ladders; logloglog when J is log N",
        ),
        Template(
            "Tao_shape_unit",
            "log k / (k^2 (loglog k)^2)",
            "second-moment overlap dividend spread over compressed scales",
        ),
        Template(
            "third_log_tax",
            "Tao_shape / logloglog k",
            "one extra unresolved meta-scale tax",
        ),
        Template(
            "third_log_dividend",
            "Tao_shape * logloglog k",
            "possible if meta-scale blocks are independently profitable",
        ),
    ]
    print("TEMPLATES")
    for t in templates:
        print(f"  {t.name:<26} {t.formula:<62} {t.meaning}")
    print()

    print("NUMERIC ORIENTATION")
    print(" k        L1        L2       L3      Tao      /L3      *L3      rank1")
    for k in (13, 50, 10**3, 10**6, 10**12):
        kf = float(k)
        l1 = safe_L(kf, 1)
        l2 = safe_L(kf, 2)
        l3 = safe_L(kf, 3)
        rank1 = entropy_rank_shape(kf, 1.0)
        print(
            f"{k:<8g} {fmt_float(l1)} {fmt_float(l2)} {fmt_float(l3)} "
            f"{tao_shape(kf):8.3e} {third_log_tax_shape(kf):8.3e} "
            f"{third_log_dividend_shape(kf):8.3e} {rank1:8.3e}"
        )
    print()

    print("PRODUCT SURVIVAL EXAMPLES")
    print(" J      ordinary       compressed     meta-compressed")
    for J in (20, 100, 1000, 10000):
        ordinary = scale_product(3, J, 0.5, "log")
        compressed = scale_product(3, J, 0.5, "loglog")
        meta = scale_product(4, J, 0.5, "logloglog")
        print(f"{J:<6d} {ordinary:12.5e} {compressed:12.5e} {meta:16.5e}")
    print()

    print("OWN INEQUALITY TEMPLATES")
    print("1. Meta-scale dividend template:")
    print("   If a second-moment proof has one profitable block in each")
    print("   logloglog k independent meta-scale class, try surplus")
    print("   >= c log k logloglog k / (k^2 (loglog k)^2).")
    print("2. Rank-tax template:")
    print("   If r unresolved orthogonal scale channels must all be paid, replace")
    print("   loglog k by loglog k + r logloglog k in the denominator.")
    print("3. Residual-core product template:")
    print("   If a sieve-covered residual exports alpha/(j log j) frontier mass")
    print("   on compressed scale j, its unresolved mass falls like (log J)^(-alpha).")
    print("4. Helly-scale template:")
    print("   If every live automaton state of size M has a witness among the first")
    print("   H(M) determinant rows, then global CRT search can be replaced by")
    print("   a sum over H(M); loglog factors come from harmonic H(M), not from time.")
    print()

    print("ASSUMPTION CHALLENGE")
    print("  Vertices considered: integers, primes, denominator scales, prime blocks,")
    print("  proof obligations, determinant rows, residual packets, and entropy channels.")
    print("  Chosen vertices: scale currencies.  This preserves where logarithmic")
    print("  losses are paid and destroys the original object labels.")


if __name__ == "__main__":
    main()
