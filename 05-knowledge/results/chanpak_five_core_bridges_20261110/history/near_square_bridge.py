"""Exact near-square and fixed-chart consecutive-count controls."""
from functools import lru_cache
from math import isqrt
from pathlib import Path
import hashlib
import json

from audit_compiler_selector import coefficients, construct, formula, g, path_data
from check_five_cores import CORES, check, paths_dp, ready_word, rooted_profiles


@lru_cache(None)
def old_source_f(n):
    if n <= 2:
        return 1
    if n % 2 == 0:
        return old_source_f(n//2)
    m = (n+1)//2
    return old_source_f(m)+m-(m % 2)


def main():
    check(all(old_source_f(n) == g(n-1) for n in range(1, 4097)), "g shift mismatch")
    cap = 1000001
    witnesses = {}
    for q in range(cap):
        value = g(2*q)
        if value <= cap:
            witnesses.setdefault(value, 2*q)
    misses = [(n, n*n+1) for n in range(2, isqrt(cap-1)+1, 2)
              if n*n+1 not in witnesses]
    check(misses[0] == (10,101) and len(misses) == 231, "square probe mismatch")
    ready101 = ready_word(["11", "01", "00", "00", "10"])
    check(paths_dp(ready101) == 101, "inherited C6 source does not give101")
    profile101 = rooted_profiles(ready101)[0]
    check(profile101["ready_gap"] == 6 and profile101["joint_delta"] == 6,
          "wrong neutral defect at101")
    specs = [(3, 5, "left", 62), (4, 5, "left", 64), (3, 14, "right", 66)]
    records = []
    checks = 0
    for core_index, out, side, constant in specs:
        rows = CORES[core_index]
        z = coefficients(rows, out, path_data(rows))
        h, d, ell, i, o, s, e = z
        linear, correction = (i, s) if side == "left" else (o, e)
        check(d == 9 and 2*linear == 82 and h-d+linear+correction == constant,
              "symbolic consecutive-family coefficient mismatch")
        examples = []
        for k in range(16):
            a, b = (2*k+1, 0) if side == "left" else (0, 2*k+1)
            child = construct(rows, out, a, b)
            expected = 9*g(k)+82*k+constant
            check(paths_dp(child) == expected, "literal adjacent-family count mismatch")
            examples.append({"k": k, "H": expected, "order": len(child)})
            checks += 1
        records.append({"core_index": core_index, "out_mask": out, "varying_side": side,
                        "address": "2k+1", "other_address": 0,
                        "coefficients_h_d_L_i_o_s_e": list(z), "constant": constant,
                        "examples": examples})
    five101 = None
    for core_index, rows in enumerate(CORES):
        data = path_data(rows)
        for out in range(1 << len(rows)):
            z = coefficients(rows, out, data)
            for a in range(21):
                for b in range(21):
                    if formula(z, a, b) == 101:
                        child = construct(rows, out, a, b)
                        check(paths_dp(child) == 101, "five-bank101 witness failed")
                        five101 = {"core_index": core_index, "out_mask": out, "A": a, "B": b,
                                   "child_rows": list(child), "order": len(child)}
                        break
                if five101:
                    break
            if five101:
                break
        if five101:
            break
    check(five101 is not None, "no positive five-bank101 control found")
    report = {
        "status": "FINITE-EXACT controls plus symbolic fixed-chart coefficient identities",
        "single_apex_universe": "Every g-value<=1,000,001 via normalized even input2q, 0<=q<1,000,001",
        "normalization_proof": "Delete trailing ones; g(2q)>=q+1 makes the q cap exhaustive",
        "even_square_parameters": "2<=n<=1000, n even", "targets": 500,
        "single_apex_missing_count": len(misses), "single_apex_misses": misses,
        "first_20_controls": [{"n": n, "target": n*n+1, "even_address": witnesses.get(n*n+1)}
                              for n in range(2,42,2)],
        "old_source_function_shift_checks": 4096,
        "ready101": {"word": ["11","01","00","00","10"], "rows": list(ready101),
                     "root_profile": profile101},
        "five_core101": five101,
        "consecutive_family": records, "literal_consecutive_checks": checks,
        "scope": "No all-square or prime infinitude assertion; no arbitrary lower-anchor selection",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name("near_square_receipt.json").write_text(
        json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"missing": len(misses), "first_missing": misses[0],
                      "ready101": profile101, "five_core101": five101,
                      "consecutive_checks": checks}, indent=2))


if __name__ == "__main__":
    main()
