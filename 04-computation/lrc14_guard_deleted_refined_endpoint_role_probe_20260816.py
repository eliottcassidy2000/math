#!/usr/bin/env python3
"""Exact guard-packet-row-refined 13^3 endpoint-role theorem component.

Only the H-labelled row is removed from the algebraic owner-packet span.
The endpoint patterns retain their literal Boolean ``guard_safe`` factors.
"""

import ast
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIDGE_PATH = ROOT / "04-computation/lrc_half_twist_relation_current_bridge_thm3479.py"
BRIDGE_SHA = "ad2a620cdc238f28e3384698b2c612f38cdf2566bd56b76d1cbabcc03107ec0b"
EXPECTED_COARSE_GAMMA_SHA = "afcdb043eb1bf8095c313473a3d3bdcf4ce027f86b01b5f3cecbc7c87e6484b3"
EXPECTED_COARSE_TARGET_SHA = "726423df1b9e1c93b356966e5c3c386669e2f6b19da0bf8818606204eb2e9ee5"
EXPECTED_GAMMA_SHA = "1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682"
EXPECTED_CHART_SHA = "b7d8c2c9860e4f1aa542b1c85fdb7b65cf4985aba5a81a84ff3a324834d51c51"
EXPECTED_SEMANTIC_SHA = "f06a258b85212ba981de4c18027147a8875c67fc35d1f55d9128c17d746994a9"
P = 13
ROLE_CLASSES = {
    "c1": (0, 0, 0), "c2": (1, 0, 0), "c3": (0, 1, 0),
    "H": (1, 0, 1), "q2": (1, 12, 0), "q3": (1, 0, 0),
    "q4": (1, 0, 0), "q5": (1, 0, 0),
}
EDGES = (
    (0, 3), (0, 4), (0, 5), (1, 2), (1, 4), (1, 7), (2, 4),
    (2, 7), (3, 4), (3, 5), (4, 5), (4, 6), (4, 7),
)
HUB, LEAF = 4, 6
WINGS = ((0, 3, 5), (1, 2, 7))
BLOCKERS = ("c1", "c2", "c3")
UNITS = ("q2", "q3", "q4")


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def load_module():
    actual = hashlib.sha256(BRIDGE_PATH.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual == BRIDGE_SHA, (actual, BRIDGE_SHA))
    spec = importlib.util.spec_from_file_location("thm3479_refined", BRIDGE_PATH)
    require(spec is not None and spec.loader is not None, "bad module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_module()
CTX = None


def context():
    global CTX
    if CTX is None:
        word = M.to_current(M.U_FULL_REL)
        t_den = 182 * M.lcm_tuple(word)
        nn = M.R_DILATION * t_den
        p, root, prime_factors, lucas_witnesses = M.FULL_EMBEDDINGS[0]
        M.verify_lucas_certificate(p, prime_factors, lucas_witnesses,
                                   "U_full refined prime")
        M.verify_embedding(p, root, nn, M.FULL_NN_FACTORS, "U_full refined")
        zero = (0,) * 9
        q_intervals = M.fast_build_set(word, t_den, M.PATTERN_QA, zero)
        q_starts = [left for left, _ in q_intervals]
        embeddings = ((p, root),)
        tabs = M.fast_make_tabs(q_intervals, M.X_FREQUENCY, nn, embeddings)
        CTX = word, t_den, nn, p, root, q_intervals, q_starts, embeddings, tabs
    return CTX


def worker(alpha):
    word, t_den, nn, p, root, q_intervals, q_starts, embeddings, tabs = context()
    rows, interval_count = [], 0
    for beta in range(P):
        for tau in range(P):
            ell = (tau, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
            e_intervals = M.fast_build_set(word, t_den, M.PATTERN_E, ell)
            interval_count += len(e_intervals)
            ax, _ = M.fast_x_sweep(
                e_intervals, q_intervals, q_starts, M.X_FREQUENCY,
                t_den, nn, embeddings, tabs,
            )
            by = M.fast_endpoint_sum(
                e_intervals, -(M.X_FREQUENCY + word[M.TARGET_B]), nn, embeddings
            )
            phase = pow(root, beta * (nn // P) % nn, p)
            rows.append(phase * ax[0] % p * by[0] % p)
    return alpha, interval_count, tuple(rows)


def digest_integers(values):
    return hashlib.sha256(",".join(str(v) for v in values).encode("ascii")).hexdigest()


def digest_json(value):
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(body).hexdigest()


def inverse_value(gamma, q, p, zeta, normalize=True):
    total = 0
    qa, qb, qt = q
    powers = tuple(pow(zeta, k, p) for k in range(P))
    index = 0
    for alpha in range(P):
        for beta in range(P):
            for tau in range(P):
                exponent = -(alpha * qa + beta * qb + tau * qt) % P
                total = (total + gamma[index] * powers[exponent]) % p
                index += 1
    return total * (pow(P**3, -1, p) if normalize else 1) % p


def coarse_inverse(gamma0, q, p, zeta, normalize=True):
    total, index = 0, 0
    powers = tuple(pow(zeta, k, p) for k in range(P))
    for alpha in range(P):
        for beta in range(P):
            total = (total + gamma0[index] * powers[-(alpha*q[0] + beta*q[1]) % P]) % p
            index += 1
    return total * (pow(P**2, -1, p) if normalize else 1) % p


def permutations3(values):
    a, b, c = values
    return ((a,b,c),(a,c,b),(b,a,c),(b,c,a),(c,a,b),(c,b,a))


def role_charts():
    charts = []
    for swap in (0, 1):
        for blockers in permutations3(BLOCKERS):
            for units in permutations3(UNITS):
                chart = {HUB: "H", LEAF: "q5"}
                chart.update(zip(WINGS[swap], blockers))
                chart.update(zip(WINGS[1-swap], units))
                charts.append(tuple(chart[i] for i in range(8)))
    answer = tuple(sorted(set(charts)))
    require(len(answer) == 72, len(answer))
    return answer


def det_mod(matrix, p):
    rows = [list(row) for row in matrix]
    answer = 1
    for column in range(len(rows)):
        pivot = next((r for r in range(column, len(rows)) if rows[r][column] % p), None)
        if pivot is None:
            return 0
        if pivot != column:
            rows[pivot], rows[column] = rows[column], rows[pivot]
            answer = -answer
        value = rows[column][column] % p
        answer = answer * value % p
        inverse = pow(value, -1, p)
        for r in range(column + 1, len(rows)):
            factor = rows[r][column] * inverse % p
            for c in range(column, len(rows)):
                rows[r][c] = (rows[r][c] - factor * rows[column][c]) % p
    return answer % p


def k4_factor(values, chart, vertices, p):
    positions = {v: i for i, v in enumerate(vertices)}
    lap = [[0] * 4 for _ in range(4)]
    for left, right in EDGES:
        if left not in positions or right not in positions:
            continue
        weight = (values[chart[left]] - values[chart[right]]) % p
        i, j = positions[left], positions[right]
        lap[i][i] += weight; lap[j][j] += weight
        lap[i][j] -= weight; lap[j][i] -= weight
    return det_mod([row[:-1] for row in lap[:-1]], p)


def chart_census(values, p):
    rows = []
    for chart in role_charts():
        bridge = (values["H"] - values["q5"]) % p
        left = k4_factor(values, chart, WINGS[0] + (HUB,), p)
        right = k4_factor(values, chart, WINGS[1] + (HUB,), p)
        rows.append((chart, bridge, left, right, bridge * left % p * right % p))
    return tuple(rows)


def main():
    word, t_den, nn, p, root, q_intervals, _, _, _ = context()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(a for a, _, _ in chunks) == tuple(range(P)), "alpha order")
    gamma = tuple(value for _, _, rows in chunks for value in rows)
    require(len(gamma) == P**3, len(gamma))
    require(digest_integers(gamma) == EXPECTED_GAMMA_SHA, digest_integers(gamma))
    gamma0 = tuple(gamma[(alpha*P + beta)*P] for alpha in range(P) for beta in range(P))
    require(digest_integers(gamma0) == EXPECTED_COARSE_GAMMA_SHA, digest_integers(gamma0))
    require((gamma0[0], gamma0[1]) == (248447851579556601771, 510954897124935772821), "coarse controls")
    zeta = pow(root, nn // P, p)
    require(pow(zeta, P, p) == 1 and zeta != 1, "zeta order")
    coarse_targets_unnormalized = tuple(
        coarse_inverse(gamma0, (a,b), p, zeta, False)
        for a in range(P) for b in range(P)
    )
    require(digest_integers(coarse_targets_unnormalized) == EXPECTED_COARSE_TARGET_SHA,
            digest_integers(coarse_targets_unnormalized))
    distinct_classes = tuple(sorted(set(ROLE_CLASSES.values())))
    refined_values = {q: inverse_value(gamma, q, p, zeta) for q in distinct_classes}
    values = {label: refined_values[q] for label, q in ROLE_CLASSES.items()}
    coarse_bases = tuple(sorted(set((q[0],q[1]) for q in distinct_classes)))
    fibre_rows = []
    for base in coarse_bases:
        refined_fibre = tuple(inverse_value(gamma, base + (tau,), p, zeta) for tau in range(P))
        coarse = coarse_inverse(gamma0, base, p, zeta)
        require(sum(refined_fibre) % p == coarse, (base, sum(refined_fibre) % p, coarse))
        fibre_rows.append((base, coarse, refined_fibre, sum(refined_fibre) % p))
    charts = chart_census(values, p)
    bridge_zeros = sum(row[1] == 0 for row in charts)
    left_zeros = sum(row[2] == 0 for row in charts)
    right_zeros = sum(row[3] == 0 for row in charts)
    product_zeros = sum(row[4] == 0 for row in charts)
    require((bridge_zeros, left_zeros, right_zeros, product_zeros) == (0, 0, 0, 0),
            (bridge_zeros, left_zeros, right_zeros, product_zeros))
    require(digest_json(charts) == EXPECTED_CHART_SHA, digest_json(charts))
    flat = chart_census({label: 1 for label in ROLE_CLASSES}, p)
    require(all(row[1:] == (0,0,0,0) for row in flat), "flat hostile")
    record = (
        BRIDGE_SHA, p, root, t_den, nn, len(q_intervals),
        tuple(count for _, count, _ in chunks), digest_integers(gamma),
        digest_integers(gamma0), digest_integers(coarse_targets_unnormalized),
        tuple((q, refined_values[q]) for q in distinct_classes), tuple(fibre_rows),
        charts, (bridge_zeros,left_zeros,right_zeros,product_zeros),
    )
    require(digest_json(record) == EXPECTED_SEMANTIC_SHA, digest_json(record))
    source = Path(__file__).read_text(encoding="utf-8")
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
            "assert node")
    print("LRC REFINED ENDPOINT ROLE PROBE U_FULL")
    print("status=FINITE-EXACT theorem component; THM-3479 PROVED STRUCTURAL + INDEPENDENTLY AUDITED; LRC(14) OPEN")
    print(f"dependency={BRIDGE_PATH.name}:{BRIDGE_SHA}")
    print(f"embedding=(prime={p},root={root},order={nn}); primitive_order_certified=True")
    print(f"universe=gamma(alpha*v1+beta*v2+tau*e_H), F13^3, size={len(gamma)}; H owner-packet row removed, Boolean guard_safe retained")
    print(f"interval_counts_by_alpha={tuple(count for _,count,_ in chunks)} total={sum(count for _,count,_ in chunks)}")
    print(f"gamma_sha256={digest_integers(gamma)}")
    print(f"tau0_coarse_recovery=(gamma_sha256={digest_integers(gamma0)},target_sha256={digest_integers(coarse_targets_unnormalized)})")
    print(f"refined_role_classes={tuple((label,ROLE_CLASSES[label]) for label in ROLE_CLASSES)}")
    print(f"refined_role_values={tuple((q,refined_values[q]) for q in distinct_classes)}")
    print(f"fibre_sum_recovery={tuple(fibre_rows)}")
    print(f"charts=72 factor_zero_counts=(bridge,left_K4,right_K4,product)={(bridge_zeros,left_zeros,right_zeros,product_zeros)}")
    print(f"factor_histograms=(bridge={tuple(sorted(Counter(r[1] for r in charts).items()))},left={tuple(sorted(Counter(r[2] for r in charts).items()))},right={tuple(sorted(Counter(r[3] for r in charts).items()))})")
    print(f"chart_record_sha256={digest_json(charts)}")
    print("flat_response_hostile=all role potentials 1; all 72 factor quadruples (0,0,0,0)")
    print(f"semantic_schema=(dependency,embedding,scales,q_count,interval_counts,gamma_digest,coarse_digests,role_values,fibre_rows,chart_rows,zero_counts)")
    print(f"semantic_sha256={digest_json(record)}")
    print("scope=unrestricted mod13 endpoint aggregates only; not grouped C(a;X,m), all-unit B(q), physical current, ancestry, bispectrum, scalar exclusion, or LRC(14)")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
