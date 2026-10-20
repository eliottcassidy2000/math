"""Independent finite audit of the incoming THM-4730 odd-circulant argument.

Run: python 05-knowledge/results/weighted_core_transport_20261020/incoming_audit.py
No external dependencies. All checks remain active under python -O.

The actual predicate tests every tournament edge under positional reversal.
The constant-partner-sum proxy tests only affine reflections. Neither this
bounded verification nor a failed elementary repair proves/refutes B1.
"""

import ast
from itertools import combinations, product
from pathlib import Path


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def valid_connection_set(n, connection):
    return (
        0 not in connection
        and len(connection) == (n - 1) // 2
        and all(((d in connection) != ((n - d) in connection)) for d in range(1, n))
    )


def valid_path(n, connection, path):
    return sorted(path) == list(range(n)) and all(
        (v - u) % n in connection for u, v in zip(path, path[1:])
    )


def reversal(path):
    n = len(path)
    result = [None] * n
    for i, vertex in enumerate(path):
        result[vertex] = path[n - 1 - i]
    return result


def actual_palindromic(n, connection, path):
    reverse = reversal(path)
    return all(
        (((v - u) % n in connection) == ((reverse[u] - reverse[v]) % n in connection))
        for u, v in combinations(range(n), 2)
    )


def partner_sums(n, path):
    return sorted({(path[i] + path[n - 1 - i]) % n for i in range(n)})


def three_interval_path(n, a, b):
    x, y = a - 1, b - 1
    return list(range(x)) + list(range(x + y, n)) + list(range(x, x + y))


def all_connection_sets(n):
    for choices in product((0, 1), repeat=(n - 3) // 2):
        yield {1} | {
            d if choice == 0 else n - d
            for d, choice in zip(range(2, (n + 1) // 2), choices)
        }


def known_arc_witness(n, a, b, path):
    """Find a contradiction using only the three already specified offsets."""
    known = {1, a, b}
    negative = {n - d for d in known}
    reverse = reversal(path)
    for u, v in combinations(range(n), 2):
        original = (v - u) % n
        image = (reverse[v] - reverse[u]) % n
        if (original in known and image in known) or (
            original in negative and image in negative
        ):
            return u, v, original, image
    return None


def main():
    n, connection = 9, {1, 3, 4, 7}
    path = [0, 1, 2, 3, 7, 5, 6, 4, 8]
    need(valid_connection_set(n, connection), "Invalid nonaffine witness tournament")
    need(valid_path(n, connection, path), "Invalid nonaffine witness path")
    need(actual_palindromic(n, connection, path), "Witness reversal is not anti")
    need(partner_sums(n, path) == [5, 8], "Unexpected witness partner sums")
    # Load only the repaired predicate, avoiding the legacy script's main census.
    source = (Path(__file__).resolve().parents[3] / "04-computation" / "experiments"
              / "metagraph_20261020" / "atom_alphabet_proofcheck.py")
    tree = ast.parse(source.read_text(encoding="utf-8"))
    function = next(node for node in tree.body
                    if isinstance(node, ast.FunctionDef) and node.name == "palindromic")
    namespace = {}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), "exec"), namespace)
    need(namespace["palindromic"](n, path, lambda u, v: (v - u) % n in connection),
         "Repaired legacy predicate still rejects the nonaffine palindrome")
    print("PROVED FINITE WITNESS: the constant-partner-sum proxy is not necessary.")
    print(f"n={n}; S={sorted(connection)}; P={path}")
    print(f"reversal={reversal(path)}; partner_sums={partner_sums(n, path)}")
    print("Repaired legacy predicate accepts the nonaffine palindromic witness.")

    total_sets = total_paths = descending_cases = 0
    for n in range(7, 26, 2):
        sets_here = paths_here = 0
        for connection in all_connection_sets(n):
            need(valid_connection_set(n, connection), f"Invalid S at n={n}")
            sets_here += 1
            for a, b in combinations(sorted(connection - {1}), 2):
                if a + b > n + 1:
                    continue
                path = three_interval_path(n, a, b)
                need(valid_path(n, connection, path), f"B1 path invalid: {n,a,b}")
                need(
                    not actual_palindromic(n, connection, path),
                    f"B1 counterexample: n={n}, S={sorted(connection)}, a={a}, b={b}",
                )
                paths_here += 1
            small = connection.intersection(range(2, (n + 1) // 2))
            if not small:
                path = [0, 1] + [
                    vertex
                    for j in range(1, (n - 3) // 2 + 1)
                    for vertex in (n - 2 * j, n - 2 * j + 1)
                ] + [2]
                need(valid_path(n, connection, path), f"Descending path invalid: n={n}")
                need(
                    not actual_palindromic(n, connection, path),
                    f"Descending path counterexample: n={n}",
                )
                descending_cases += 1
        total_sets += sets_here
        total_paths += paths_here
        print(f"FINITE-EXACT n={n}: sets={sets_here}; actual-anti B1 checks={paths_here}")
    need(total_sets == 4092, "Unexpected connection-set universe")
    need(total_paths == 103423, "Unexpected eligible-path universe")
    need(descending_cases == 10, "Unexpected descending-case universe")
    print(
        f"FINITE-EXACT total: {total_sets} sets; {total_paths} B1 paths; "
        f"{descending_cases} descending paths; no actual-anti counterexample."
    )

    tested = 0
    first_miss = None
    for n in range(7, 52, 2):
        for a, b in combinations(range(2, n - 1), 2):
            # These offsets can coexist with 1 in an oriented circulant.
            if a + b > n + 1 or a + b == n:
                continue
            path = three_interval_path(n, a, b)
            tested += 1
            if known_arc_witness(n, a, b, path) is None:
                first_miss = n, a, b, path
                break
        if first_miss is not None:
            break
    need(first_miss == (11, 2, 7, [0, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6]),
         "Unexpected first known-offset witness miss")
    need(tested == 16, "Unexpected minimal-witness search prefix")
    n, a, b, path = first_miss
    completions = [s for s in all_connection_sets(n) if {1, a, b} <= s]
    need(len(completions) == 4, "Unexpected completion count")
    need(all(valid_path(n, s, path) for s in completions), "Invalid completed path")
    need(all(not actual_palindromic(n, s, path) for s in completions),
         "Actual specialized-path counterexample among completions")
    print(f"ONE-PAIR REPAIR MISS after {tested} parameter cases: n={n}, a={a}, b={b}.")
    print(f"P={path}; no contradiction from only {{1,a,b}} or its negative.")
    print("All 4 full connection-set completions still fail actual palindromicity.")
    print("SCOPE: proof gap established; specialized B1 and alphabet classification not refuted.")
    print("PASS")


if __name__ == "__main__":
    main()
