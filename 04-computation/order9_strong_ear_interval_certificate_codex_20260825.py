#!/usr/bin/env python3
"""Construct explicit strong order-nine witnesses for every odd H in [85,2881].

This is a certificate generator, not an exhaustive order-nine census.  It
generates the 6,880 order-eight tournament isomorphism classes, retains the
strong parents, and evaluates all 254 nonconstant one-vertex cuts by an exact
ear formula.

For a parent T on V and a new vertex e, put

    I = {a : a -> e},    O = {b : e -> b}.

Let E[S,a] count Hamiltonian paths of T[S] ending at a and B[S,b]
Hamiltonian paths of T[S] beginning at b.  Define

    Q[a,b] = sum_{S: a in S, b notin S} E[S,a] B[V-S,b].

Splitting a child Hamiltonian path at e gives the disjoint exact formula

    H(T+e) = sum_{a in I} E[V,a]
           + sum_{b in O} B[V,b]
           + sum_{a in I,b in O} Q[a,b].

The first two sums count paths ending/starting at e; Q counts the internal-ear
paths.  If T is strong and I,O are nonempty, T+e is strong.  Thus every row
of the emitted table is a compact constructive certificate.

Mask convention (both parent and child): lexicographic pairs i<j, with bit 1
meaning i->j.  A cut bit b=1 means the new vertex e beats b.
"""

from hashlib import sha256
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from h_spectrum_n9_exhaustive_monad_s6 import generate_level
from h21_finite_check_v2_monad_s4 import beats_from_canon


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "05-knowledge/results/order9_strong_ear_interval_certificates_codex_20260825.tsv"
)
TARGETS = tuple(range(85, 2882, 2))


def require(condition, label):
    if not condition:
        raise RuntimeError("CHECK FAILED: " + label)


def edge(adj, i, j):
    return bool((adj[i] >> j) & 1)


def is_strong(adj):
    n = len(adj)
    full = (1 << n) - 1
    seen = 1
    frontier = 1
    while frontier:
        nxt = 0
        bits = frontier
        while bits:
            bit = bits & -bits
            bits ^= bit
            nxt |= adj[bit.bit_length() - 1]
        nxt &= full ^ seen
        seen |= nxt
        frontier = nxt
    if seen != full:
        return False
    reverse = [0] * n
    for i in range(n):
        bits = adj[i]
        while bits:
            bit = bits & -bits
            bits ^= bit
            reverse[bit.bit_length() - 1] |= 1 << i
    seen = 1
    frontier = 1
    while frontier:
        nxt = 0
        bits = frontier
        while bits:
            bit = bits & -bits
            bits ^= bit
            nxt |= reverse[bit.bit_length() - 1]
        nxt &= full ^ seen
        seen |= nxt
        frontier = nxt
    return seen == full


def endpoint_dp(adj):
    """Return induced-subtournament path counts by final and initial vertex."""
    n = len(adj)
    size = 1 << n
    end = [[0] * n for _ in range(size)]
    begin = [[0] * n for _ in range(size)]
    for v in range(n):
        end[1 << v][v] = 1
        begin[1 << v][v] = 1
    for mask in range(1, size):
        if mask & (mask - 1) == 0:
            continue
        bits = mask
        while bits:
            bit = bits & -bits
            bits ^= bit
            v = bit.bit_length() - 1
            rest = mask ^ bit
            end[mask][v] = sum(
                end[rest][u]
                for u in range(n)
                if (rest >> u) & 1 and edge(adj, u, v)
            )
            begin[mask][v] = sum(
                begin[rest][u]
                for u in range(n)
                if (rest >> u) & 1 and edge(adj, v, u)
            )
    return end, begin


def ear_kernel(adj):
    n = len(adj)
    full = (1 << n) - 1
    end, begin = endpoint_dp(adj)
    kernel = [[0] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            total = 0
            for subset in range(1, full):
                if (subset >> a) & 1 and not ((subset >> b) & 1):
                    total += end[subset][a] * begin[full ^ subset][b]
            kernel[a][b] = total
    return end[full], begin[full], kernel


def ear_value(end_full, begin_full, kernel, cut):
    n = len(end_full)
    full = (1 << n) - 1
    incoming = full ^ cut
    terminal = sum(end_full[a] for a in range(n) if (incoming >> a) & 1)
    initial = sum(begin_full[b] for b in range(n) if (cut >> b) & 1)
    internal = sum(
        kernel[a][b]
        for a in range(n)
        if (incoming >> a) & 1
        for b in range(n)
        if (cut >> b) & 1
    )
    return terminal, initial, internal


def mask_from_adjacency(adj):
    n = len(adj)
    code = 0
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if edge(adj, i, j):
                code |= 1 << bit
            bit += 1
    return code


def extend(adj, cut):
    n = len(adj)
    out = list(adj) + [0]
    for i in range(n):
        if (cut >> i) & 1:
            out[n] |= 1 << i
        else:
            out[i] |= 1 << n
    return tuple(out)


def h_count_direct(adj):
    """Held--Karp control used only on calibration cuts in this generator."""
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v, count in enumerate(dp[mask]):
            if not count:
                continue
            available = adj[v] & ~mask
            while available:
                bit = available & -available
                available ^= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return sum(dp[-1])


def calibration():
    # C3: bits for 0->1, 1->2, 2->0.
    c3 = (1 << 1, 1 << 2, 1 << 0)
    end, begin, kernel = ear_kernel(c3)
    rows = []
    for cut in range(1, 7):
        parts = ear_value(end, begin, kernel, cut)
        child = extend(c3, cut)
        direct = h_count_direct(child)
        require(sum(parts) == direct, "C3 all-cut ear calibration")
        require(is_strong(child), "C3 nonconstant-cut strong calibration")
        rows.append((cut, direct))
    require(rows == [(1, 5), (2, 5), (3, 5), (4, 5), (5, 5), (6, 5)],
            "C3 calibration transcript")
    return tuple(rows)


def main():
    print("ORDER-NINE STRONG-EAR INTERVAL CERTIFICATE")
    print("formula_calibration_C3=" + repr(calibration()))

    generation_stream = StringIO()
    with redirect_stdout(generation_stream):
        class_codes = generate_level(8)
    require(class_codes is not None and len(class_codes) == 6880,
            "A000568(8) isomorphism-class count")
    require("MISMATCH" not in generation_stream.getvalue(),
            "canonical-augmentation count transcript")
    found = {}
    strong_parents = 0
    cuts_scanned = 0
    formula_spot_checks = 0
    h613_parent_classes = 0
    h613_parent_cuts = 0
    h613_to_h623_ears = 0

    for parent_code in sorted(class_codes):
        parent = tuple(beats_from_canon(parent_code, 8))
        if not is_strong(parent):
            continue
        strong_parents += 1
        end_full, begin_full, kernel = ear_kernel(parent)
        parent_h = sum(end_full)
        if parent_h == 613:
            h613_parent_classes += 1
        for cut in range(1, 255):
            cuts_scanned += 1
            terminal, initial, internal = ear_value(
                end_full, begin_full, kernel, cut
            )
            value = terminal + initial + internal
            if parent_h == 613:
                h613_parent_cuts += 1
                h613_to_h623_ears += value == 623
            if formula_spot_checks < 512 and (cut % 31 == parent_code % 31):
                require(
                    value == h_count_direct(extend(parent, cut)),
                    "deterministic formula spot check",
                )
                formula_spot_checks += 1
            if value in TARGETS and value not in found:
                child = extend(parent, cut)
                child_code = mask_from_adjacency(child)
                found[value] = (
                    parent_code,
                    cut,
                    child_code,
                    parent_h,
                    terminal,
                    initial,
                    internal,
                )

    require(formula_spot_checks == 512, "formula spot-check count")
    require(h613_parent_classes == 2, "two strong order-eight H=613 classes")
    require(h613_parent_cuts == 2 * 254, "complete H=613-parent cut universe")
    require(h613_to_h623_ears == 0, "no nonconstant H=613 to H=623 ear")
    missing = tuple(value for value in TARGETS if value not in found)
    require(not missing, "complete odd interval 85..2881")

    header = (
        "H\tparent_code\tcut\tchild_code\tH_parent\tterminal"
        "\tinitial\tinternal\n"
    )
    rows = [header]
    for value in TARGETS:
        fields = (value,) + found[value]
        rows.append("\t".join(map(str, fields)) + "\n")
    payload = "".join(rows)
    CERTIFICATE.write_text(payload)
    digest = sha256(payload.encode()).hexdigest()

    print("order8_classes=6880")
    print("strong_order8_parents=" + str(strong_parents))
    print("nonconstant_cuts_scanned=" + str(cuts_scanned))
    print("formula_spot_checks=" + str(formula_spot_checks))
    print("H613_parent_classes=" + str(h613_parent_classes))
    print("H613_parent_nonconstant_cuts=" + str(h613_parent_cuts))
    print("H613_to_H623_ears=" + str(h613_to_h623_ears))
    print("certificate_rows=" + str(len(TARGETS)))
    print("interval=all odd H in [85,2881]")
    print("certificate_sha256=" + digest)
    print("certificate_path=" + str(CERTIFICATE.relative_to(ROOT)))
    for value in (85, 613, 623, 2881):
        print("witness_" + str(value) + "=" + repr(found[value]))
    print("conclusion=FINITE-EXACT explicit strong order-nine ear witnesses")


if __name__ == "__main__":
    main()
