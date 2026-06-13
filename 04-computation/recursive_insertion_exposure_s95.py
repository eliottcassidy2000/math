#!/usr/bin/env python3
"""
recursive_insertion_exposure_s95.py

Vertex insertion becomes simple if the recursive state carries the right
boundary statistic.

For a tournament T and a new vertex v, write sig[x]=1 if v -> x and sig[x]=0
if x -> v.  For a Hamiltonian path P in T, v can be inserted:
  - at the start if sig[first(P)] = 1
  - at the end if sig[last(P)] = 0
  - between consecutive a,b in P if sig[a]=0 and sig[b]=1

The tempting simplification "insert v into Hamiltonian paths of T" is false:
v can bridge an old ordered pair a,b even when a -> b is not an edge of T.

The correct boundary statistic is the bridge-exposure matrix
  Q_T[a,b] = # old-vertex permutations where a is immediately followed by b
             and every other consecutive step is a valid edge of T.

Together with start/end counts of Hamiltonian paths in T:
  H(T+v) =
    sum_{b: v -> b} HP_start[b]
  + sum_{a: a -> v} HP_end[a]
  + sum_{a: a -> v, b: v -> b} Q_T[a,b].

This script verifies the formula exhaustively for n<=5 bases and all insertion
signatures, then measures how much state is needed to predict all insertions.
"""

from collections import Counter, defaultdict
from itertools import permutations


def bits_to_adj(bits, n):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


def add_vertex(A, sig):
    n = len(A)
    B = [row[:] + [0] for row in A]
    B.append([0] * (n + 1))
    for x, bit in enumerate(sig):
        if bit:
            B[n][x] = 1
        else:
            B[x][n] = 1
    return B


def hamiltonian_paths_list(A):
    n = len(A)
    paths = []
    for p in permutations(range(n)):
        if all(A[p[i]][p[i + 1]] for i in range(n - 1)):
            paths.append(p)
    return paths


def hp_exposure_matrix(A):
    n = len(A)
    E = [[0] * n for _ in range(n)]
    paths = hamiltonian_paths_list(A)
    for p in paths:
        for i in range(n - 1):
            E[p[i]][p[i + 1]] += 1
    return len(paths), E


def boundary_state(A):
    """Return (start_count, end_count, bridge_exposure Q)."""
    n = len(A)
    starts = [0] * n
    ends = [0] * n
    for p in hamiltonian_paths_list(A):
        starts[p[0]] += 1
        ends[p[-1]] += 1

    Q = [[0] * n for _ in range(n)]
    for p in permutations(range(n)):
        for k in range(n - 1):
            ok = True
            for i in range(n - 1):
                if i == k:
                    continue
                if not A[p[i]][p[i + 1]]:
                    ok = False
                    break
            if ok:
                Q[p[k]][p[k + 1]] += 1
    return starts, ends, Q


def endpoint_matrix_from_paths(A):
    n = len(A)
    M = [[0] * n for _ in range(n)]
    for p in hamiltonian_paths_list(A):
        M[p[0]][p[-1]] += 1
    return M


def insertion_formula(starts, ends, Q, sig):
    n = len(sig)
    total = 0
    for b in range(n):
        if sig[b]:
            total += starts[b]
    for a in range(n):
        if not sig[a]:
            total += ends[a]
    for a in range(n):
        if sig[a]:
            continue
        for b in range(n):
            if sig[b]:
                total += Q[a][b]
    return total


def response_vector(A):
    n = len(A)
    H, E = hp_exposure_matrix(A)
    starts, ends, Q = boundary_state(A)
    values = []
    values_formula = []
    for mask in range(1 << n):
        sig = tuple((mask >> i) & 1 for i in range(n))
        values.append(len(hamiltonian_paths_list(add_vertex(A, sig))))
        values_formula.append(insertion_formula(starts, ends, Q, sig))
    return tuple(values), tuple(values_formula), H, E, starts, ends, Q


def flatten_matrix(M):
    return tuple(x for row in M for x in row)


def score_sequence(A):
    return tuple(sorted((sum(row) for row in A), reverse=True))


def first_collision(records, key_name):
    buckets = defaultdict(list)
    for rec in records:
        buckets[rec[key_name]].append(rec)
    for key, rows in buckets.items():
        responses = {r["R"] for r in rows}
        if len(responses) > 1:
            a, b = rows[0], next(r for r in rows if r["R"] != rows[0]["R"])
            return key, a, b
    return None


def main():
    print("=" * 78)
    print("RECURSIVE INSERTION EXPOSURE S95")
    print("=" * 78)
    print()
    print("Formula: H(T+v) = starts(sig=1) + ends(sig=0) + bridge cut Q(sig=0, sig=1)")
    print("where Q[a,b] counts old-vertex spines with only the a _ b slot left to repair.")

    all_summaries = []
    for n in range(2, 6):
        records = []
        failures = []
        for bits in range(1 << (n * (n - 1) // 2)):
            A = bits_to_adj(bits, n)
            R, R_formula, H, E, starts, ends, Q = response_vector(A)
            if R != R_formula:
                failures.append((bits, R, R_formula))
            records.append({
                "bits": bits,
                "H": H,
                "score": score_sequence(A),
                "endpoint_M": flatten_matrix(endpoint_matrix_from_paths(A)),
                "hp_exposure_E": flatten_matrix(E),
                "boundary_Q": (tuple(starts), tuple(ends), flatten_matrix(Q)),
                "R": R,
            })

        distinct = {
            "H": len({r["H"] for r in records}),
            "score": len({r["score"] for r in records}),
            "endpoint_M": len({r["endpoint_M"] for r in records}),
            "hp_exposure_E": len({r["hp_exposure_E"] for r in records}),
            "boundary_Q": len({r["boundary_Q"] for r in records}),
            "insertion_response": len({r["R"] for r in records}),
        }
        all_summaries.append((n, distinct))

        print("\n" + "-" * 78)
        print(f"n={n}: tournaments={len(records)}, formula failures={len(failures)}")
        for name, val in distinct.items():
            print(f"  distinct {name}: {val}")

        for key_name in ["H", "score", "endpoint_M", "hp_exposure_E"]:
            col = first_collision(records, key_name)
            if col is None:
                print(f"  {key_name} predicts insertion response for n={n}")
            else:
                key, a, b = col
                print(f"  {key_name} is insufficient: key={key}")
                print(f"    bits {a['bits']} response starts {a['R'][:min(8, len(a['R']))]}")
                print(f"    bits {b['bits']} response starts {b['R'][:min(8, len(b['R']))]}")

        # By construction boundary_Q determines R; check whether R ever splits a Q-bucket.
        q_collision = first_collision(records, "boundary_Q")
        print(f"  boundary_Q determines response: {q_collision is None}")

        # Aggregate identities.
        sum_E_identity_ok = True
        for rec in records:
            if sum(rec["hp_exposure_E"]) != (n - 1) * rec["H"]:
                sum_E_identity_ok = False
                break
        print(f"  Σ_ab E[a,b] = (n-1)H for every T: {sum_E_identity_ok}")

    print("\n" + "=" * 78)
    print("SEQUENCE TABLE")
    print("=" * 78)
    for label in ["H", "score", "endpoint_M", "hp_exposure_E", "boundary_Q", "insertion_response"]:
        seq = [d[label] for _, d in all_summaries]
        print(f"  distinct {label}: {seq}")

    print("\n" + "=" * 78)
    print("CLARIFICATION")
    print("=" * 78)
    print("""
Endpoint data is a terminal statistic: it remembers where Hamiltonian paths
start and end.  HP exposure remembers valid consecutive arcs inside old
Hamiltonian paths; in this small range it is already an extremely strong
fingerprint.  But the direct local insertion formula uses the larger boundary
object Q, because insertion can also repair invalid old adjacencies.

The recursive simplification is:

  carry Q_T = (start counts, end counts, bridge-exposure matrix)
  instead of only H(T) or endpoint M(T).

Then every insertion signature is evaluated by one cut sum across the bipartition
sig^{-1}(0) | sig^{-1}(1), plus the two boundary terms.  The correction is
literally the number of broken Hamiltonian spines that the new vertex repairs.
""")


if __name__ == "__main__":
    main()
