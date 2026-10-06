"""Compose the proved run-family floor with the arithmetic-shell selector.

Production takes a source, not moment values or a supplied ROOT word. Sources
outside this proper family return a typed open obligation, never a zero floor.
"""
from fractions import Fraction as F
import json
import collatz_run_block_family_20261005 as runs
import collatz_refinement_energy_dual_20261005 as shells


def compile_measurement(source):
    receipt = runs.receipt(source)
    if receipt is None:
        return {"source":source,"status":"unresolved_by_this_family"}
    epsilon = receipt["leaf_floor"]
    degree = shells.conditional_degree(epsilon)
    error = shells.error_bound(degree)
    coefficients = shells.minorant(degree)
    return {
        "source":source,"status":"proved_on_run_family",
        "leaf":receipt["leaf"],"source_index":receipt["index"],
        "expense_ceiling":receipt["expense_ceiling"],
        "odd_step_deadline":receipt["odd_step_deadline"],
        "atom_floor":epsilon,"polynomial_degree":degree,
        "linear_form_coefficients":coefficients,
        "uniform_selector_error":error,
        "exact_readout_lower":epsilon-error,
        "simple_exact_readout_lower":3*epsilon/4,
        "conditional_moment_width":epsilon/2048,
        "conditional_noisy_readout_lower":epsilon/2,
        "direct_coefficient_norm":sum(map(abs,coefficients),F(0)),
        "original_selector_polynomial_degree":receipt["degree"]+1,
        "moment_values_consumed":0,"supplied_ROOT_words_consumed":0}


def main():
    checks = 0
    def need(condition, why):
        nonlocal checks
        checks += 1
        if not condition:
            raise ValueError(why)

    # These are neither inputs nor fallback oracles for the composed compiler.
    old_export, old_literal = runs.export_word, runs.prior.literal_root
    def forbidden(*args, **kwargs):
        raise RuntimeError("ROOT-word oracle was consumed by production")
    runs.export_word = runs.prior.literal_root = forbidden
    try:
        rows = [compile_measurement(n) for n in (5,23,35,739,7,27)]
    finally:
        runs.export_word, runs.prior.literal_root = old_export, old_literal

    for row in rows[:4]:
        e,d = row["atom_floor"],row["polynomial_degree"]
        need(row["status"] == "proved_on_run_family", "family status retained")
        need(d >= 3 and d % 2 == 1, "odd arithmetic selector")
        need(row["uniform_selector_error"] <= e/4, "paid selector error")
        need(d == 3 or shells.error_bound(d-2) > e/4, "least permitted degree")
        need(row["direct_coefficient_norm"] < 512, "uniform conditioning")
        need(row["exact_readout_lower"] >= 3*e/4 > 0, "positive exact readout")
        need(row["exact_readout_lower"]-512*row["conditional_moment_width"] >= e/2,
             "truthful interval budget leaves a positive noisy readout")
        need(sum(row["linear_form_coefficients"]) == 1, "target normalization")
        need(d < row["original_selector_polynomial_degree"], "strict degree reduction")
        need(row["moment_values_consumed"] == row["supplied_ROOT_words_consumed"] == 0,
             "input provenance")
    for row in rows[4:]:
        need(row == {"source":row["source"],"status":"unresolved_by_this_family"},
             "uncovered source has no manufactured floor")
    for invalid in (True,1.0,0,-3,4,"23"):
        try:
            compile_measurement(invalid)
        except (TypeError,ValueError):
            need(True,"invalid source rejected")
        else:
            need(False,"invalid source accepted")
    need(rows[1]["leaf"] == 15 and rows[1]["atom_floor"] == F(1,5148)
         and rows[1]["polynomial_degree"] == 9
         and rows[1]["simple_exact_readout_lower"] == F(1,6864),
         "source23 explicit all-height instance")
    compact = []
    for row in rows:
        compact.append({k:v for k,v in row.items() if k not in
                        ("linear_form_coefficients","exact_readout_lower","direct_coefficient_norm")})
    print(json.dumps({"checks":checks,"plans":compact,
        "scope":"Family theorem supplies the atom floor; arithmetic support supplies the selector. Universal coverage remains OPEN."},
        indent=2,sort_keys=True,default=str))


if __name__ == "__main__":
    main()
