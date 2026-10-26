"""Count-only cross-checks and rigorous coverage budgets; no high-order atlas.

Uses the inherited audited reverse-valuation engine, not a new Burnside solver.
The saved family census is independent of that engine.
"""
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT / "05-knowledge/results/lca_ports_reverse_20261025/counting/reverse.py"
REFERENCE = ROOT / "04-computation/b000568_v2.txt"


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def guaranteed_binary_decay(numerator, denominator):
    """Largest integer b with ratio <= 2^-b, evaluated without floating point."""
    b = denominator.bit_length() - numerator.bit_length()
    if b >= 0 and numerator << b > denominator:
        b -= 1
    return b


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "counts-normal.json")
    args = parser.parse_args()
    spec = spec_from_file_location("inherited_reverse_family_audit", SOURCE)
    engine = module_from_spec(spec)
    sys.modules[spec.name] = engine
    spec.loader.exec_module(engine)
    reference = {}
    for line in REFERENCE.read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            n, a = line.split()
            reference[int(n)] = int(a)
    computed, checks = {}, []
    # Small checks connect literal families to count arithmetic; large controls
    # exercise consecutive orders needed by the exact coverage bound.
    for n in list(range(10)) + [50, 76, 99, 100, 127, 128]:
        dag = engine.OrdinaryDAG(n)
        value, stats = dag.exact(audit_original=True)
        need(value == reference[n], "repository reference disagreement")
        normalized, residue_stats = dag.residue(32)
        need(normalized == (value >> dag.B) % (1 << 32), "normalized residue mismatch")
        computed[n] = value
        checks.append({"order": n, "digits": len(str(value)),
                       "sha256_decimal": hashlib.sha256(str(value).encode()).hexdigest(),
                       "baseline": dag.B, "valuation": (value & -value).bit_length()-1,
                       "normalized_mod_2_32": normalized, "gauge_coordinates": stats["coordinate_checks"],
                       "structure": dag.structure(), "residue_updates": residue_stats["modular_updates"]})
        print(f"exact count and normalized residue PASS n={n}", flush=True)
    census = json.loads((HERE / "order9-normal.json").read_text())
    for row in census["census"]:
        n = row["order"]
        need(row["classes"] == computed[n], "literal versus count-only comparison")
        for m in range(1, n+1):
            coverage = sum(v for q, v in row["by_core_order"].items() if int(q) <= m)
            need(coverage <= comb(n-1, m-1)*computed[m], "finite core cutoff coverage bound")
    budgets = []
    for n in [100, 128]:
        numerator = (n-1)*computed[n-1]
        denominator = computed[n]
        ratio = Fraction(numerator, denominator)
        budgets.append({"order": n, "core_cutoff": n-1,
                        "upper_fraction_numerator": str(ratio.numerator),
                        "upper_fraction_denominator": str(ratio.denominator),
                        "upper_bound_2_to_minus": guaranteed_binary_decay(numerator, denominator),
                        "interpretation": "all families with core order below n occupy at most this fraction; remainder are new pair-free cores"})
    # Exact asymptotic hostility: a fixed finite core bank has polynomial
    # response, while A(n)>=2^(n(n-1)/2)/n!; no fixed bank covers all orders.
    result = {"status": "FINITE-EXACT; counts only beyond order9", "checks": checks,
              "coverage_budgets": budgets,
              "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "inherited_engine_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "reference_sha256": hashlib.sha256(REFERENCE.read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("PASS: count checks, reachable gauges, residues, finite coverage bounds, high-order budgets")


if __name__ == "__main__":
    main()
