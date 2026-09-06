"""Exact native creation minors and odd-order complement characters.

Universe: complete Eulerian orbit sets n=3,...,7; all nonempty cycle-length
subsets. No new larger census or unbounded kernel-saturation assertion.
"""
from hashlib import sha256
from itertools import combinations
import json
import sympy as sp
from sympy.polys.matrices import DomainMatrix
import eulerian_boolean_kernel_overnight_hexagon_sep05 as native

gates = 0


def check(condition, label):
    global gates
    if not condition:
        raise RuntimeError(label)
    gates += 1


def rank(matrix):
    if not matrix.rows or not matrix.cols:
        return 0
    return len(DomainMatrix.from_Matrix(matrix).convert_to(sp.QQ).rref()[1])


def kernel(matrix):
    return DomainMatrix.from_Matrix(matrix).convert_to(sp.QQ).nullspace().to_Matrix().T


def components(mask, n, index):
    adjacency = [set() for _ in range(n)]
    for edge, bit in index.items():
        if mask >> bit & 1:
            a, b = edge
            adjacency[a].add(b)
            adjacency[b].add(a)
    isolates = [v for v in range(n) if not adjacency[v]]
    unseen = set(range(n)) - set(isolates)
    cycles = {}
    while unseen:
        start = min(unseen)
        reached = {start}
        todo = [start]
        while todo:
            v = todo.pop()
            for w in adjacency[v] - reached:
                reached.add(w)
                todo.append(w)
        unseen -= reached
        if all(len(adjacency[v]) == 2 for v in reached):
            length = len(reached)
            cycles[length] = cycles.get(length, 0) + 1
    return isolates, cycles


def character_basis(parity, complement, sign, q):
    columns = []
    for i in parity:
        j = complement[i]
        if j < i:
            continue
        if j == i:
            if sign == 1:
                columns.append(sp.eye(q)[:, i])
        else:
            columns.append(sp.eye(q)[:, i] + sign * sp.eye(q)[:, j])
    return sp.Matrix.hstack(*columns) if columns else sp.zeros(q, 0)


def main():
    q_counts = {0: 1, 1: 1, 2: 1, 3: 2, 4: 3, 5: 7, 6: 16, 7: 54}
    records = []
    family_count = 0
    weighted_two_hostile = False
    for n in range(3, 8):
        reps, lookup, sizes, index = native.generator_orbits(n)
        q = len(reps)
        check(q == q_counts[n], "complete inherited orbit count")
        info = [components(F, n, index) for F in reps]
        cycles = {k: native.old.simple_cycles(n, k, index) for k in range(3, n + 1)}
        layers = {k: sp.Matrix(native.old.weighted_operator(reps, lookup, cycles[k]))
                  for k in cycles}
        local_records = []
        for length_count in range(1, n - 1):
            for lengths in combinations(range(3, n + 1), length_count):
                family_count += 1
                k = max(lengths)
                raw = sum((layers[j] for j in lengths), sp.zeros(q))
                boolean = sp.Matrix(native.old.boolean_support(raw.tolist()))
                check(boolean == boolean.T, "native Boolean undirected support")
                domain = sorted((i for i in range(q) if len(info[i][0]) >= k),
                                key=lambda i: (-reps[i].bit_count(), -info[i][1].get(k, 0), i))
                check(len(domain) == q_counts[n - k], "complete isolate stratum")
                targets = []
                for i in domain:
                    vertices = info[i][0][:k]
                    added = sum(1 << index[tuple(sorted((vertices[j], vertices[(j + 1) % k])))]
                                for j in range(k))
                    H = reps[i] ^ added
                    targets.append(lookup[H])
                    check(reps[lookup[H]].bit_count() == reps[i].bit_count() + k,
                          "cycle creation edge count")
                    for ell in lengths:
                        for toggle in cycles[ell]:
                            result = H ^ toggle
                            check(result.bit_count() >= reps[i].bit_count() + k - ell,
                                  "literal toggle lower edge bound")
                            if result.bit_count() == reps[i].bit_count():
                                check(ell == k and (H & toggle) == toggle,
                                      "equal-edge terms are full maximum-cycle deletions")
                                j = lookup[result]
                                if j != i:
                                    check(info[j][1].get(k, 0) > info[i][1].get(k, 0),
                                          "noncomponent deletion increases isolated-cycle count")
                check(len(set(targets)) == len(targets), "creation injective on isomorphism classes")
                U = boolean.extract(targets, domain)
                raw_U = raw.extract(targets, domain)
                for a, i in enumerate(domain):
                    check(U[a, a] == 1 and raw_U[a, a] == info[i][1].get(k, 0) + 1,
                          "Boolean unit versus exact raw component multiplicity")
                    for b in range(a + 1, len(domain)):
                        check(U[a, b] == raw_U[a, b] == 0, "ordered creation minor lower triangular")
                    if n == 6 and k == 3 and reps[i].bit_count() == 3:
                        check(raw_U[a, a] == 2 and U[a, a] == 1, "dyadic raw-pivot hostile")
                        weighted_two_hostile = True
                check(U.det() == 1 and all(x.q == 1 for x in U.inv()),
                      "native creation block unimodular")
                local_records.append((lengths, len(domain), str(raw_U.det())))

        B = sp.Matrix(native.old.boolean_support(layers[3].tolist()))
        even = [i for i, F in enumerate(reps) if F.bit_count() % 2 == 0]
        odd = [i for i, F in enumerate(reps) if F.bit_count() % 2]
        C = B.extract(odd, even)
        domain = sorted((i for i in even if len(info[i][0]) >= 3),
                        key=lambda i: (-reps[i].bit_count(), -info[i][1].get(3, 0), i))
        targets = []
        for i in domain:
            vertices = info[i][0][:3]
            triangle = sum(1 << index[tuple(sorted(e))] for e in combinations(vertices, 2))
            targets.append(lookup[reps[i] ^ triangle])
        free_even = [i for i in even if i not in domain]
        free_odd = [i for i in odd if i not in targets]
        U = B.extract(targets, domain)
        V = B.extract(targets, free_even)
        W = B.extract(free_odd, domain)
        Z = B.extract(free_odd, free_even)
        S = Z - W * U.inv() * V
        check(all(x.q == 1 for x in S), "integer native Schur complement")
        check(rank(C) == len(domain) + rank(S), "unconditional Schur rank identity")
        small_kernel = kernel(S)
        lifted = sp.zeros(len(even), small_kernel.cols)
        first = -U.inv() * V * small_kernel
        for j, i in enumerate(domain):
            lifted[even.index(i), :] = first[j, :]
        for j, i in enumerate(free_even):
            lifted[even.index(i), :] = small_kernel[j, :]
        check(C * lifted == sp.zeros(len(odd), lifted.cols), "actual native Schur kernel lift")
        check(rank(lifted) == len(even) - rank(C), "complete even-side kernel parametrization")

        row = {"n": n, "q": q, "cycle_families": len(local_records),
               "creation_cases": local_records, "triangle_integer_pivots": len(domain),
               "triangle_crossrank": rank(C), "native_even_nullity": lifted.cols}
        if n % 2:
            full = (1 << len(index)) - 1
            complement = [lookup[F ^ full] for F in reps]
            J = sp.zeros(q)
            for i, j in enumerate(complement):
                J[j, i] = 1
            check(J * B == B * J, "actual primal complement commutes with Boolean adjacency")
            if n % 4 == 1:
                delta = len(even) - len(odd)
                fixed = [i for i, j in enumerate(complement) if i == j]
                k = (n - 1) // 4
                check(len(fixed) == delta and all(reps[i].bit_count() == k * n for i in fixed),
                      "fixed primal count and exact half-edge parity")
                sector_rows = []
                for sign in (-1, 1):
                    E = character_basis(even, complement, sign, q)
                    O = character_basis(odd, complement, sign, q)
                    check(J * E == sign * E and J * O == sign * O,
                          "literal complement character bases")
                    deficit = E.cols - O.cols
                    check(deficit == (delta if sign == (-1) ** k else 0),
                          "all forced index in one complement character")
                    check(J * B * E == sign * B * E, "native character preservation")
                    nullity = E.cols - rank(B * E)
                    check(nullity >= deficit, "unconditional native character kernel lower bound")
                    sector_rows.append((sign, E.cols, O.cols, nullity))
                row["complement_sectors"] = sector_rows
                four = next(i for i, F in enumerate(reps) if F.bit_count() == 4)
                six = next(i for i, F in enumerate(reps) if F.bit_count() == 6)
                witness = sp.eye(q)[:, four] - sp.eye(q)[:, six]
                check(B * witness == sp.zeros(q, 1) and J * witness == -witness,
                      "n5 complement-twin native kernel")
                check(max(len(info[four][0]), len(info[six][0])) <= 2,
                      "n5 kernel reaches nonisolated core")
            else:
                check(all(complement[i] in odd for i in even), "odd-complement parity exchange")
                folded = B.extract([complement[i] for i in even], even)
                check(folded == folded.T, "balanced crossblock symmetric after complement identification")
        if n == 4:
            check((1 << len(index)) - 1 not in lookup, "even-order complement is not an Eulerian map")
        if n == 7:
            check(abs(C.det(method="domain-ge")) == 1076,
                  "local unimodular pivots do not imply global unimodularity or mod2 invertibility")
        records.append(row)
        print("NATIVE STRUCTURE", json.dumps({key: value for key, value in row.items()
                                             if key != "creation_cases"}, sort_keys=True))
    check(weighted_two_hostile, "raw multiplicity hostile retained")
    print("COMPLETE CYCLE-FAMILY CASES", family_count)
    print("SEMANTIC SHA256", sha256(json.dumps(records, sort_keys=True).encode()).hexdigest())
    print("PASS EXPLICIT GATES", gates)
    print("PROVED creation minors and complement-character index; FINITE-EXACT controls n3..7")
    print("OPEN full native nullity, opposite-character kernels, and Laplacian gap")


if __name__ == "__main__":
    main()
