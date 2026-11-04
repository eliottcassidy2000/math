"""Exact odd-module and complete three-path response controls, stdlib only."""
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def tournament(n, mask):
    rows = [0] * n
    for k, (u, v) in enumerate(combinations(range(n), 2)):
        if mask >> k & 1: rows[u] |= 1 << v
        else: rows[v] |= 1 << u
    return rows


def mask_of(rows):
    return sum(1 << k for k, (u, v) in enumerate(combinations(range(len(rows)), 2))
               if rows[u] >> v & 1)


def converse(rows):
    full = (1 << len(rows)) - 1
    return [full ^ row ^ (1 << u) for u, row in enumerate(rows)]


def cyclic(n):
    return [sum(1 << ((u + d) % n) for d in range(1, (n + 1) // 2)) for u in range(n)]


def endpoints(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for u in range(n): dp[1 << u][u] = 1
    for s in range(1, 1 << n):
        for u in range(n):
            value = dp[s][u]
            if not value: continue
            targets = rows[u] & ~s
            while targets:
                bit = targets & -targets
                targets -= bit
                dp[s | bit][bit.bit_length() - 1] += value
    return dp


def h(rows):
    return sum(endpoints(rows)[-1]) if rows else 1


def perm_h(rows):
    return sum(all(rows[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in permutations(range(len(rows))))


def cover_matrix(rows):
    n = len(rows)
    e, b = endpoints(rows), endpoints(converse(rows))
    c = [[0] * n for _ in rows]
    full = (1 << n) - 1
    for s in range(1, full):
        for u, eu in enumerate(e[s]):
            if not eu: continue
            for v, bv in enumerate(b[full ^ s]):
                if bv: c[u][v] += eu * bv
    return c


def extend(rows, masks):
    """Bit u=1 means old u -> new port; new ports are transitively ordered."""
    n = len(rows)
    out = rows[:] + [0] * len(masks)
    for j, mask in enumerate(masks):
        for u in range(n):
            if mask >> u & 1: out[u] |= 1 << (n + j)
            else: out[n + j] |= 1 << u
        for k in range(j + 1, len(masks)): out[n + j] |= 1 << (n + k)
    return out


def inflate(quotient, modules):
    blocks, start = [], 0
    for module in modules:
        blocks.append(list(range(start, start + len(module))))
        start += len(module)
    out = [0] * start
    for i, module in enumerate(modules):
        for a, u in enumerate(blocks[i]):
            for b, v in enumerate(blocks[i]):
                if module[a] >> b & 1: out[u] |= 1 << v
        for j in range(len(modules)):
            if quotient[i] >> j & 1:
                for u in blocks[i]:
                    for v in blocks[j]: out[u] |= 1 << v
    return out, blocks


def one_module(quotient, module):
    return inflate(quotient, [module] + [[0]] * (len(quotient) - 1))[0]


def skeletons(quotient):
    m = len(quotient) - 1
    c = [0] * (m + 2)
    for order in permutations(range(1, m + 1)):
        for gaps in range(1, 1 << (m + 1)):
            valid = True
            for j in range(m + 1):
                if gaps >> j & 1:
                    if j and not (quotient[order[j - 1]] & 1): valid = False
                    if j < m and not (quotient[0] >> order[j] & 1): valid = False
                elif 0 < j < m and not (quotient[order[j - 1]] >> order[j] & 1):
                    valid = False
            if valid: c[gaps.bit_count()] += 1
    require(c[1] == h(quotient) and c[1] % 2 == 1, "c1 quotient parity")
    return c


def path_cover_profile(rows):
    n = len(rows)
    all_h = [sum(row) for row in endpoints(rows)]
    all_h[0] = 1
    q = [[0] * (n + 2) for _ in range(1 << n)]
    q[0][0] = 1
    for s in range(1, 1 << n):
        sub = s
        while sub:
            for k in range(1, s.bit_count() + 1):
                q[s][k] += all_h[sub] * q[s ^ sub][k - 1]
            sub = (sub - 1) & s
    return q[-1]


def literal_profile_derivatives(rows):
    n = len(rows)
    q, total_gradient = [0] * (n + 2), [0] * (n + 1)
    for p in permutations(range(n)):
        bad = sum(1 << j for j, (u, v) in enumerate(zip(p, p[1:])) if not (rows[u] >> v & 1))
        for cuts in range(1 << (n - 1)):
            k = cuts.bit_count() + 1
            missing = (bad & ~cuts).bit_count()
            if missing == 0:
                q[k] += 1
                total_gradient[k] += n - k
            elif missing == 1:
                total_gradient[k] += 1
    require(q == path_cover_profile(rows), "profile permutation/DP mismatch")
    for k in range(1, n + 1):
        require(total_gradient[k] == k * q[k + 1], "ordered derivative multiplicity")
    return q


def module_controls():
    totals = {"profile_derivative_cores": 0, "ambient_pairs": 0,
              "gradient_transpose_entries": 0, "higher_run_contexts": 0}
    qdata = {}
    for r in (3, 5):
        for mask in range(1 << (r * (r - 1) // 2)):
            qdata[r, mask] = literal_profile_derivatives(tournament(r, mask))
            totals["profile_derivative_cores"] += 1
    selected5 = sorted({0, 1, 12, 123, 511, 729, 1023, mask_of(cyclic(5))})
    cases = []
    for qn in range(1, 5):
        for qm in range(1 << (qn * (qn - 1) // 2)):
            quotient = tournament(qn, qm)
            coefficients = skeletons(quotient)
            module_masks3 = range(8)
            module_masks5 = range(1024) if qn <= 2 else selected5
            cases.extend((3, mm, quotient, coefficients, False) for mm in module_masks3)
            cases.extend((5, mm, quotient, coefficients, False) for mm in module_masks5)
    # Alternating uniform contexts with 4 and6 exterior vertices expose
    # k=3 andk=4 coefficients, rather than only the leading two runs.
    for runs in (2, 3):
        qn = 2 * runs + 1
        quotient = [0] * qn
        for u in range(1, qn):
            for v in range(u + 1, qn): quotient[u] |= 1 << v
            if u % 2: quotient[0] |= 1 << u
            else: quotient[u] |= 1
        coefficients = skeletons(quotient)
        for r, mask_list in ((3, range(8)), (5, selected5)):
            cases.extend((r, mm, quotient, coefficients, True) for mm in mask_list)
    for r, mm, quotient, coefficients, higher in cases:
        module = tournament(r, mm)
        first = one_module(quotient, module)
        second = one_module(quotient, converse(module))
        c, cp = cover_matrix(first), cover_matrix(second)
        q = qdata[r, mm]
        prediction = sum(coefficients[k] * (2 * (r - k) * q[k] - k * q[k + 1])
                         for k in range(1, min(len(coefficients), r + 1)))
        actual = sum(c[u][v] - c[v][u] for u in range(r) for v in range(r) if module[u] >> v & 1)
        require(actual == prediction and actual % 4 == 2, ("module current", r, mm, quotient))
        require(any(c[u][v] != cp[u][v] for u in range(r) for v in range(r)), "odd module invisible")
        for u in range(r):
            for v in range(r):
                require(cp[u][v] == c[v][u], "internal transpose")
                totals["gradient_transpose_entries"] += 1
        totals["ambient_pairs"] += 1
        totals["higher_run_contexts"] += higher
    totals["all_three_vertex_module_quotient_orders"] = [1, 4]
    totals["all_five_vertex_module_quotient_orders"] = [1, 2]
    totals["selected_five_vertex_masks_for_quotients_3_4_and_deep_probes"] = selected5
    return totals


def multi_module_controls():
    families = variants = 0
    for qm in range(8):
        quotient = tournament(3, qm)
        for sizes in ((3, 3, 1), (3, 5, 1), (3, 3, 3)):
            modules = [cyclic(r) if r > 1 else [0] for r in sizes]
            active = [i for i, r in enumerate(sizes) if r > 1]
            baseline, blocks = inflate(quotient, modules)
            base_c = cover_matrix(baseline)
            signatures = set()
            for choice in range(1 << len(active)):
                changed = modules[:]
                for bit, i in enumerate(active):
                    if choice >> bit & 1: changed[i] = converse(changed[i])
                rows, _ = inflate(quotient, changed)
                c = cover_matrix(rows)
                detector = []
                for bit, i in enumerate(active):
                    block = blocks[i]
                    reverse_it = choice >> bit & 1
                    for u in block:
                        for v in block:
                            require(c[u][v] == (base_c[v][u] if reverse_it else base_c[u][v]), "multimodule transfer")
                    detector.append(sum(c[u][v] - c[v][u] for u in block for v in block if baseline[u] >> v & 1) % 8)
                signatures.add(tuple(detector))
                variants += 1
            require(len(signatures) == 1 << len(active), "mod8 reversal bits lost")
            families += 1
    return {"families": families, "variants": variants, "quotients": "all8 labelled order3", "module_orders": [[3, 3, 1], [3, 5, 1], [3, 3, 3]]}


def addp(*polys):
    result = {}
    for poly in polys:
        for monomial, value in poly.items(): result[monomial] = result.get(monomial, 0) + value
    return {key: value for key, value in result.items() if value}


def mulp(a, b):
    result = {}
    for s, av in a.items():
        for t, bv in b.items():
            if not s & t: result[s | t] = result.get(s | t, 0) + av * bv
    return result


def path_matrix(rows):
    n = len(rows)
    g = [[{} for _ in rows] for _ in rows]
    for start in range(n):
        def visit(last, support):
            cell = g[start][last]
            cell[support] = cell.get(support, 0) + 1
            for nxt in range(n):
                if not support >> nxt & 1 and rows[last] >> nxt & 1:
                    visit(nxt, support | 1 << nxt)
        visit(start, 1 << start)
    e = [addp(*(g[i][j] for i in range(n))) for j in range(n)]
    s = [addp(*row) for row in g]
    return g, e, s


def squarefree_response(g, e, s, a, c):
    n = len(g); full = (1 << n) - 1; b, d = full ^ a, full ^ c
    lin = lambda vec, mask: addp(*(vec[i] for i in range(n) if mask >> i & 1))
    middle = lambda left, right: addp(*(g[j][k] for j in range(n) for k in range(n)
                                      if left >> j & 1 and right >> k & 1))
    first = mulp(mulp(addp({0: 1}, lin(e, a)), addp({0: 1}, middle(b, c))), addp({0: 1}, lin(s, d)))
    second = mulp(mulp(addp({0: 1}, lin(e, c)), middle(d, a)), addp({0: 1}, lin(s, b)))
    return first.get(full, 0) + second.get(full, 0)


def literal_three_tensor(rows):
    n = len(rows); result = {}
    for p in permutations(range(n)):
        for left, right in combinations(range(1, n), 2):
            if all(rows[p[j]] >> p[j + 1] & 1 for j in range(n - 1) if j + 1 not in (left, right)):
                key = (p[left - 1], p[left], p[right - 1], p[right])
                result[key] = result.get(key, 0) + 1
    return result


def submasks(mask):
    sub = mask
    while True:
        yield sub
        if sub == 0: return
        sub = (sub - 1) & mask


def response_controls():
    counts = {"cores": 0, "full_response_entries": 0, "independent_permutation_counts": 0,
              "three_tensor_entries": 0, "middle_singleton_nonzero_entries": 0,
              "quartic_coefficients": 0}
    for n in range(5):
        for mask in range(1 << (n * (n - 1) // 2)):
            rows = tournament(n, mask)
            g, e, s = path_matrix(rows)
            full = (1 << n) - 1
            table = [[h(extend(rows, [a, c])) for c in range(1 << n)] for a in range(1 << n)]
            for a in range(1 << n):
                for c in range(1 << n):
                    require(squarefree_response(g, e, s, a, c) == table[a][c], ("full joint transfer", n, mask, a, c))
                    counts["full_response_entries"] += 1
                    if n <= 3:
                        require(perm_h(extend(rows, [a, c])) == table[a][c], "independent full count")
                        counts["independent_permutation_counts"] += 1
            tensor = literal_three_tensor(rows)
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        for l in range(n):
                            value = mulp(mulp(e[i], g[j][k]), s[l]).get(full, 0)
                            require(value == tensor.get((i, j, k, l), 0), "joint three-path tensor")
                            counts["three_tensor_entries"] += 1
                            if j == k and value: counts["middle_singleton_nonzero_entries"] += 1
            pairs = list(combinations(range(n), 2))
            for pair_a in pairs:
                for pair_c in pairs:
                    a = sum(1 << u for u in pair_a); c = sum(1 << u for u in pair_c)
                    coefficient = sum((-1) ** (i.bit_count() + j.bit_count()) * table[i][j]
                                      for i in submasks(a) for j in submasks(c))
                    predicted = 0
                    for i, j in (pair_a, pair_a[::-1]):
                        for k, l in (pair_c, pair_c[::-1]): predicted += tensor.get((i, j, k, l), 0)
                    for i, j in (pair_c, pair_c[::-1]):
                        for k, l in (pair_a, pair_a[::-1]): predicted += tensor.get((i, j, k, l), 0)
                    require(coefficient == predicted and coefficient >= 0, "quartic three-path coefficient")
                    common = set(pair_a) & set(pair_c)
                    if len(common) == 2: require(coefficient == 0, "quartic diagonal")
                    elif len(common) == 1:
                        v = next(iter(common)); u = next(iter(set(pair_a) - common)); w = next(iter(set(pair_c) - common))
                        labels = [z for z in range(n) if z != v]
                        card = [sum(1 << j for j, z in enumerate(labels) if rows[t] >> z & 1) for t in labels]
                        cc = cover_matrix(card); iu, iw = labels.index(u), labels.index(w)
                        require(coefficient == cc[iu][iw] + cc[iw][iu] and coefficient > 0, "overlapping quartic pair")
                    elif n == 4: require(coefficient == 4, "order4 disjoint quartic")
                    counts["quartic_coefficients"] += 1
            counts["cores"] += 1
    return counts


def actual_query_current(rows, module_size):
    potential = [h(extend(rows, [1 << u, 0])) - h(extend(rows, [1 << u])) for u in range(module_size)]
    total = 0
    for u in range(module_size):
        for v in range(module_size):
            if rows[u] >> v & 1:
                contrast = h(extend(rows, [1 << u, 1 << v])) - h(extend(rows, [1 << v, 1 << u]))
                total += potential[u] - potential[v] - contrast
    return total


def named_controls():
    query_examples = []
    quotient = tournament(3, 5)
    for r in (3, 5):
        first = one_module(quotient, cyclic(r)); second = one_module(quotient, converse(cyclic(r)))
        c = cover_matrix(first)
        predicted = sum(c[u][v] - c[v][u] for u in range(r) for v in range(r) if first[u] >> v & 1)
        got = actual_query_current(first, r)
        require(got == predicted and got % 4 == 2, "actual query current")
        # Detector uses the first core's fixed reference directions in both states.
        potentials = [h(extend(second, [1 << u, 0])) - h(extend(second, [1 << u])) for u in range(r)]
        other = sum(potentials[u] - potentials[v] - h(extend(second, [1 << u, 1 << v]))
                    + h(extend(second, [1 << v, 1 << u]))
                    for u in range(r) for v in range(r) if first[u] >> v & 1)
        require(other == -got and (got - other) % 8 == 4, "mod8 actual-query separation")
        query_examples.append({"module_order": r, "core_order": len(first), "reference_currents": [got, other]})
    first = one_module(quotient, cyclic(9)); second = one_module(quotient, converse(cyclic(9)))
    full = (1 << len(first)) - 1
    values = {"base": [h(first), h(second)],
              "one_port": [h(extend(t, [full ^ 1])) for t in (first, second)],
              "two_port": [h(extend(t, [full ^ 1, full ^ 2])) for t in (first, second)]}
    require(values == {"base": [52875, 52875], "one_port": [133903, 133903], "two_port": [414055, 417121]}, "eleven-vertex witness")
    hostile = tournament(5, 0)
    c = cover_matrix(hostile)
    subset = (0, 1, 4)
    hostile_current = sum(c[u][v] - c[v][u] for u in subset for v in subset if hostile[u] >> v & 1)
    require(hostile_current == -12, "odd nonmodule hostile")
    require(any(len({bool(hostile[v] >> u & 1) for u in subset}) > 1
                for v in range(5) if v not in subset), "hostile is a module")
    return {"query_current_controls": query_examples,
            "odd_nonmodule_hostile": {"rows": hostile, "subset": list(subset), "current": hostile_current},
            "eleven_vertex": {"first": first, "second": second, **values}}


def main():
    source = Path(__file__)
    digest = sha256(source.read_bytes()).hexdigest()
    result = {"status": "PASS", "source_sha256": digest,
              "odd_module": module_controls(), "independent_modules": multi_module_controls(),
              "full_joint_response": response_controls(), "named_controls": named_controls()}
    require(sha256(source.read_bytes()).hexdigest() == digest, "source changed during run")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__": main()
