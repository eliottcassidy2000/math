#!/usr/bin/env python3
"""Complete nauty/THM-4131 exact audit for prospective THM-4135."""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
import importlib.util
import multiprocessing as mp
import os
import subprocess


BASE_PATH = os.path.join(
    os.path.dirname(__file__),
    "tournament_strong_centrality_through_order_eight_thm4131.py",
)
EXPECTED_BASE_SHA256 = "6b195b0379d1ae3e5d215aa1c495f7180daeecae189df86269d07ef855867881"
EXPECTED_SEMANTIC = "0b3d2c65723e0ecc78cbf02d1735320794a6bd6e5f7c3a371ee94432e19c5b49"
EXPECTED_STRONG_COUNT = 178133
EXPECTED_ALL_COUNT = 191536


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def file_digest(path):
    hasher = sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


require(file_digest(BASE_PATH) == EXPECTED_BASE_SHA256, "THM-4131 evaluator hash")
SPEC = importlib.util.spec_from_file_location("thm4131", BASE_PATH)
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)


def audit_label(label):
    adjacency = BASE.decode_bits(label, 9)
    require(BASE.is_strong(adjacency), "gentourng strong record")
    record = BASE.analyze(adjacency)
    central = {-1, 1}
    best_central = max(
        layer["coset_floor"] for layer in record["layers"]
        if layer["t"] in central
    )
    best_outer = max(
        layer["coset_floor"] for layer in record["layers"]
        if layer["t"] not in central
    )
    return (
        label,
        BASE.digest(BASE.profile_signature(record)),
        BASE.compact(record),
        record["normalized_tilt"],
        record["theta"] == 0,
        record["rational_t"],
        record["coset_t"],
        record["actual_t"],
        set(record["rational_t"]) <= central,
        set(record["coset_t"]) <= central,
        not central.isdisjoint(record["actual_t"]),
        best_central - best_outer,
    )


def raw_universe(executable, strong):
    command = [executable, "-q"]
    if strong:
        command.append("-c")
    command.append("9")
    process = subprocess.run(command, check=False, capture_output=True, text=True)
    require(process.returncode == 0 and not process.stderr.strip(), "gentourng universe")
    lines = tuple(line.strip() for line in process.stdout.splitlines() if line.strip())
    return lines, sha256(("\n".join(lines) + "\n").encode()).hexdigest()


def main():
    executable = BASE.find_gentourng()
    all_labels, all_digest = raw_universe(executable, False)
    require(len(all_labels) == EXPECTED_ALL_COUNT, "all-class count")
    require(all_digest == "4f7d6c43cfed87e1e5293dc751736efe2d7bc1554946cdc83f4026a575fbbbf8",
            "all-class raw digest")

    strong_labels, strong_digest = raw_universe(executable, True)
    require(len(strong_labels) == EXPECTED_STRONG_COUNT, "strong-class count")
    require(strong_digest == "3073d5ecf5f34345aa5f35c349c51b35f4c244e687db25c16627fcb602b019a1",
            "strong-class raw digest")

    profile_hasher = sha256()
    rational_hist = Counter()
    coset_hist = Counter()
    actual_hist = Counter()
    rational_fail = coset_fail = actual_fail = 0
    theta_zero = coset_reorders = 0
    min_margin = None
    worst_rho = None
    worst = []
    workers = int(os.environ.get("THM4135_WORKERS", min(8, max(1, (os.cpu_count() or 2) - 1))))
    with mp.Pool(workers) as pool:
        results = pool.imap(audit_label, strong_labels, chunksize=64)
        for result in results:
            (label, profile_digest, packet, rho, is_zero, rational_t, coset_t,
             actual_t, rational_only, coset_only, actual_has_central, margin) = result
            profile_hasher.update(profile_digest.encode())
            profile_hasher.update(b"\n")
            rational_hist[rational_t] += 1
            coset_hist[coset_t] += 1
            actual_hist[actual_t] += 1
            rational_fail += not rational_only
            coset_fail += not coset_only
            actual_fail += not actual_has_central
            theta_zero += is_zero
            coset_reorders += coset_t != rational_t
            min_margin = margin if min_margin is None else min(min_margin, margin)
            if worst_rho is None or rho > worst_rho:
                worst_rho = rho
                worst = [(label, packet)]
            elif rho == worst_rho:
                worst.append((label, packet))

    ordered_profile_digest = profile_hasher.hexdigest()
    require(ordered_profile_digest ==
            "f000c999cab374ffcdc03cffaa08675f5e97152340a511c3d8de745230141746",
            "ordered primary profile digest")
    require((rational_fail, coset_fail, actual_fail) == (0, 0, 3248), "failure counts")
    require(theta_zero == 2661 and coset_reorders == 10772, "zero/reorder counts")
    require(min_margin == 90, "strict central coset margin")
    require(worst_rho == Fraction(16430, 44663) and len(worst) == 2, "worst tilt")
    require(tuple(label for label, _ in worst) == (
        "111111011100111110111111111111111101",
        "101111101111111111111110011101111111",
    ), "worst labels")

    summary = {
        "theorem": "THM-4135",
        "all_count": len(all_labels),
        "strong_count": len(strong_labels),
        "all_raw_sha256": all_digest,
        "strong_raw_sha256": strong_digest,
        "ordered_profile_sha256": ordered_profile_digest,
        "rational_fail": rational_fail,
        "coset_fail": coset_fail,
        "actual_fail": actual_fail,
        "theta_zero": theta_zero,
        "coset_reorders": coset_reorders,
        "minimum_strict_coset_margin": min_margin,
        "worst_rho": BASE.pair(worst_rho),
        "worst_labels": tuple(label for label, _ in worst),
        "worst_packets": tuple(packet for _, packet in worst),
        "rational_histogram": tuple(sorted(rational_hist.items())),
        "coset_histogram": tuple(sorted(coset_hist.items())),
        "actual_histogram": tuple(sorted(actual_hist.items())),
        "scope": (
            "complete strong order-nine rational and exact-coset centrality; actual maxima "
            "may be noncentral; THM-4133 refutes all-order centrality at order twelve"
        ),
    }
    semantic = BASE.digest(summary)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")
    print("status=PASS")
    print("theorem=THM-4135 complete strong order-nine centrality")
    print(f"class_counts=all:{summary['all_count']};strong:{summary['strong_count']}")
    print(f"raw_sha256=all:{all_digest};strong:{strong_digest}")
    print(f"ordered_profile_sha256={ordered_profile_digest}")
    print(f"failures=rational:{rational_fail};coset:{coset_fail};actual:{actual_fail}")
    print(f"theta_zero={theta_zero};coset_reorders={coset_reorders};minimum_strict_coset_margin={min_margin}")
    print(f"worst_rho={summary['worst_rho']};multiplicity={len(worst)}")
    print(f"worst_labels={summary['worst_labels']}")
    print(f"worst_packets={summary['worst_packets']}")
    print(f"rational_histogram={summary['rational_histogram']}")
    print(f"coset_histogram={summary['coset_histogram']}")
    print(f"actual_histogram={summary['actual_histogram']}")
    print("scope=complete n=9 positive row; actual maxima separate; all-order centrality refuted by THM-4133")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
