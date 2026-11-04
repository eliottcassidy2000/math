#!/usr/bin/env python3
"""Exact standard-library checks for constant-region port-cycle accounting."""
import itertools
import json
import hashlib
from pathlib import Path


def check(ok, *context):
    if not ok:
        raise RuntimeError(context)


def gate_edges(pins):
    return {(p, g) for g, ps in pins.items() for p in ps if p in pins}


def components(vertices, edges):
    unseen = set(vertices)
    adjacency = {v: set() for v in unseen}
    for u, v in edges:
        adjacency[u].add(v)
        adjacency[v].add(u)
    out = []
    while unseen:
        root = min(unseen)
        stack = [root]
        block = set()
        unseen.remove(root)
        while stack:
            u = stack.pop()
            block.add(u)
            for v in adjacency[u] & unseen:
                unseen.remove(v)
                stack.append(v)
        out.append(block)
    return out


def beta(vertices, edges):
    return len(edges) - len(vertices) + len(components(vertices, edges))


def phi(pins):
    degrees = {g: 0 for g in pins}
    for p, g in gate_edges(pins):
        degrees[p] += 1
    return sum(max(d - 1, 0) for d in degrees.values())


def ancestors(pins, output):
    live = set()
    stack = [output]
    while stack:
        g = stack.pop()
        if g in pins and g not in live:
            live.add(g)
            stack.extend(pins[g])
    return live


def gf2_rank(vectors):
    basis = {}
    for value in vectors:
        while value:
            p = value.bit_length() - 1
            if p in basis:
                value ^= basis[p]
            else:
                basis[p] = value
                break
    return len(basis)


def cycle_basis(vertices, edges):
    """Fundamental binary cycle vectors, with one bit per original edge."""
    adjacency = {v: [] for v in vertices}
    result = []
    for i, (u, v) in enumerate(edges):
        stack = [(u, 0)]
        seen = {u}
        path = None
        while stack:
            w, mask = stack.pop()
            if w == v:
                path = mask
                break
            for x, j in adjacency[w]:
                if x not in seen:
                    seen.add(x)
                    stack.append((x, mask ^ (1 << j)))
        if path is None:
            adjacency[u].append((v, i))
            adjacency[v].append((u, i))
        else:
            result.append(path ^ (1 << i))
    check(len(result) == beta(vertices, edges), "cycle basis size")
    check(gf2_rank(result) == len(result), "cycle basis independence")
    return result


def minor_valid(original_vertices, original_edges, pins, branch):
    occupied = set()
    for g in pins:
        block = branch[g]
        check(g in block, "surviving gate absent from branch")
        check(not occupied.intersection(block), "branch intersection")
        check(block <= set(original_vertices), "branch outside universe")
        occupied.update(block)
        es = {(u, v) for u, v in original_edges if u in block and v in block}
        check(len(components(block, es)) == 1, "disconnected branch", block)
    for p, g in gate_edges(pins):
        check(any((u in branch[p] and v in branch[g]) or
                  (v in branch[p] and u in branch[g]) for u, v in original_edges),
              "missing branch crossing", p, g)


def remove_region(pins, output, C, B):
    ps = {g: tuple(-999 if p in C else p for p in pair)
          for g, pair in pins.items() if g not in C}
    branch = {g: {g} for g in ps}
    for g in sorted(B):
        rest = set(ps[g]) - {-999}
        check(len(rest) == 1, "nonunary boundary", pins, C, g)
        p = next(iter(rest))
        if p in ps:
            branch[p].update(branch[g])
        del branch[g]
        del ps[g]
        ps = {h: tuple(p if q == g else q for q in pair) for h, pair in ps.items()}
        if output == g:
            output = p
    live = ancestors(ps, output)
    ps = {g: pair for g, pair in ps.items() if g in live}
    branch = {g: block for g, block in branch.items() if g in live}
    return ps, output, branch


def region_data(pins, output, C, edges=None, cycles=None, full_check=True):
    if edges is None:
        edges = sorted(gate_edges(pins))
    if cycles is None:
        cycles = cycle_basis(set(pins), edges)
    outside = set(pins) - C
    es_out = {(u, v) for u, v in edges if u in outside and v in outside}
    comps = components(outside, es_out)
    kappa = len(comps)
    incoming = [(u, v) for u, v in edges if u not in C and v in C]
    a = len(incoming)
    deg = {g: 0 for g in pins}
    for u, v in edges:
        deg[u] += 1
    E = sum(deg[g] - 1 for g in C)
    Q = E + a - kappa + 1
    r = sum(not(set(pins[g]) & C) for g in C)
    f = sum(len(set(pins[g]) & C) == 1 for g in C)
    B = {g for g in outside if set(pins[g]) & C}
    lam = sum(2 - sum(p in pins for p in set(pins[g])) for g in C)
    check(all(len(set(pins[g]) & C) == 1 for g in B), "closure boundary")
    b0 = beta(set(pins), edges)
    b_out = beta(outside, es_out)
    check(b0 == phi(pins), "live beta/Phi")
    check(b0 - b_out == Q, "port loss", pins, C)
    check(Q >= E >= 0, "nonnegative port correction", pins, C)
    check(len(C) + len(B) == 2*r + f + E, "inherited certificate")
    check(len(C) + len(B) == Q + lam + kappa - 1, "capacity identity")
    check(Q <= len(C) + len(B), "charge above certificate")
    equality = True
    for block in comps:
        count = sum(u in block for u, v in incoming)
        equality &= (count == (0 if output in block else 1))
    check((Q == E) == equality, "equality characterization")
    # The quotient map keeps physical incident-edge IDs, including parallels.
    incident = sum(1 << i for i, (u, v) in enumerate(edges) if u in C or v in C)
    image = [z & incident for z in cycles]
    check(gf2_rank(image) == Q, "quotient cycle rank")
    check(len(cycles) - gf2_rank(image) == b_out, "quotient kernel")
    quotient_label = {g: (0, g) for g in C}
    for i, block in enumerate(comps):
        quotient_label.update({g: (1, i) for g in block})
    for z in image:
        odd = set()
        for i, (u, v) in enumerate(edges):
            if z >> i & 1:
                for t in (quotient_label[u], quotient_label[v]):
                    if t in odd:
                        odd.remove(t)
                    else:
                        odd.add(t)
        check(not odd, "image not cycle")
    if full_check:
        final, out, branch = remove_region(pins, output, C, B)
        minor_valid(outside, es_out, final, branch)
        check(phi(pins) - phi(final) >= Q, "sequential budget")
        check(len(pins) - len(final) >= len(C) + len(B), "deleted gate certificate")
    return {"r": r, "f": f, "E": E, "a": a, "kappa": kappa,
            "Q": Q, "lambda": lam, "C": sorted(C), "B": sorted(B)}


def enumerate_topologies(max_gates, repeated):
    rows = []
    for n in range(1, max_gates + 1):
        pools = []
        for g in range(n):
            sources = range(-2, g)
            pairs = itertools.combinations_with_replacement(sources, 2) if repeated else itertools.combinations(sources, 2)
            pools.append(tuple(pairs))
        totals = {"gates": n, "all_topologies": 0, "live_topologies": 0,
                  "closed_regions": 0, "positive_port_corrections": 0,
                  "E_zero_Q_positive": 0, "whole_certificate_consumed": 0}
        for pairs in itertools.product(*pools):
            totals["all_topologies"] += 1
            pins = dict(enumerate(pairs))
            if ancestors(pins, n - 1) != set(pins):
                continue
            totals["live_topologies"] += 1
            edges = sorted(gate_edges(pins))
            cycles = cycle_basis(set(pins), edges)
            for mask in range(1, 1 << (n - 1)):
                C = {g for g in range(n - 1) if mask >> g & 1}
                if any(set(ps) <= C and g not in C for g, ps in pins.items()):
                    continue
                totals["closed_regions"] += 1
                row = region_data(pins, n - 1, C, edges, cycles)
                totals["positive_port_corrections"] += (row["Q"] > row["E"])
                totals["E_zero_Q_positive"] += (row["E"] == 0 and row["Q"] > 0)
                totals["whole_certificate_consumed"] += (row["Q"] == len(C) + len(row["B"]))
        rows.append(totals)
    return rows


def table_value(table, x, y):
    return (table >> (2*x+y)) & 1


def truth_vectors(pins, tables, restriction, n_inputs):
    values = {p: [] for p in range(-n_inputs, 0)}
    values.update({g: [] for g in pins})
    for bits in itertools.product((0, 1), repeat=n_inputs):
        if any(v is not None and bits[i] != v for i, v in enumerate(restriction)):
            continue
        point = dict(zip(range(-n_inputs, 0), bits))
        for g, (p, q) in pins.items():
            point[g] = table_value(tables[g], point[p], point[q])
        for name in values:
            values[name].append(point[name])
    return {k: tuple(v) for k, v in values.items()}


def truth_universe():
    out = {"inputs": 3, "gates_through": 2, "live_circuits": 0,
           "proper_coordinate_restrictions": 0, "nonconstant_output_cases": 0,
           "nonempty_constant_region_cases": 0}
    for n in (1, 2):
        pools = [tuple(itertools.combinations_with_replacement(range(-3, g), 2)) for g in range(n)]
        for pairs in itertools.product(*pools):
            pins = dict(enumerate(pairs))
            if ancestors(pins, n-1) != set(pins):
                continue
            for tables in itertools.product(range(16), repeat=n):
                out["live_circuits"] += 1
                for rho in itertools.product((None, 0, 1), repeat=3):
                    if rho == (None, None, None):
                        continue
                    out["proper_coordinate_restrictions"] += 1
                    values = truth_vectors(pins, tables, rho, 3)
                    if len(set(values[n-1])) != 2:
                        continue
                    out["nonconstant_output_cases"] += 1
                    C = {g for g in pins if len(set(values[g])) == 1}
                    if C:
                        out["nonempty_constant_region_cases"] += 1
                        region_data(pins, n-1, C)
    return out


def xor_grid(m):
    pins = {}
    forms = {}
    source = -1
    L, T, V = {}, {}, {}
    # Two independent primary inputs per boundary gate.
    for target in (L, T):
        for i in range(1, m+1):
            g = len(pins)
            pins[g] = (source, source-1)
            forms[source] = 1 << (-source-1)
            forms[source-1] = 1 << (-source)
            forms[g] = forms[source] ^ forms[source-1]
            source -= 2
            target[i] = g
    for i in range(1, m+1):
        for j in range(1, m+1):
            p = T[j] if i == 1 else V[i-1, j]
            q = L[i] if j == 1 else V[i, j-1]
            g = len(pins)
            pins[g] = (p, q)
            forms[g] = forms[p] ^ forms[q]
            check(forms[g] and forms[g] != forms[p] and forms[g] != forms[q], "grid irreducibility")
            V[i, j] = g
    check(len(set(forms[g] for g in pins)) == len(pins), "grid distinct gate forms")
    check(ancestors(pins, V[m,m]) == set(pins), "grid liveness")
    check(len(pins) == m*m+2*m and len(gate_edges(pins)) == 2*m*m, "grid sizes")
    check(phi(pins) == (m-1)**2, "grid Phi")
    return pins, V[m,m], V[1,1]


def shield(core, output, p, k):
    result = {}
    old_to_new = {}
    A, chain = [], []
    next_source = min(q for ps in core.values() for q in ps) - 1
    for g, pair in core.items():
        mapped = tuple((chain[-1] if q == p and chain else old_to_new[q]) if q in old_to_new else q for q in pair)
        new = len(result)
        result[new] = mapped
        old_to_new[g] = new
        if g == p:
            prev = new
            for i in range(k):
                ai = len(result)
                result[ai] = (next_source, next_source-1)
                next_source -= 2
                gi = len(result)
                result[gi] = (prev, ai)
                A.append(ai)
                chain.append(gi)
                prev = gi
    return result, old_to_new[output], old_to_new, A, chain


def shielding_controls():
    rows = []
    for m in range(2, 10):
        core, output, p = xor_grid(m)
        for k in (1, 2, 3, 8, 16, 32):
            current, out, mapping, A, chain = shield(core, output, p, k)
            check(len(current) == len(core)+2*k and phi(current) == phi(core), "shield initial size/Phi")
            check(ancestors(current, out) == set(current), "shield liveness")
            total_Q = 0
            charged_boundary = 0
            trace = []
            for ai, gi in reversed(list(zip(A, chain))):
                row = region_data(current, out, {ai})
                check(row["B"] == [gi] and row["Q"] == row["E"] == row["a"] == 0 and row["kappa"] == 1, "shield charge")
                charged_boundary += sum(u == gi for u, v in gate_edges(current)) - 1
                n0 = len(current)
                current, out, branch = remove_region(current, out, {ai}, {gi})
                check(len(current) == n0-2 and phi(current) == phi(core), "shield exact step")
                total_Q += row["Q"]
                trace.append({"constant": ai, "boundary": gi, "Q": row["Q"], "remaining_gates": len(current)})
            restored = {mapping[g]: tuple(mapping.get(q,q) for q in pair) for g,pair in core.items()}
            check(current == restored and out == mapping[output], "core not exactly restored")
            check(total_Q == 0 and charged_boundary == k, "shield total charge")
            rows.append({"m": m, "k": k, "core_gates": len(core), "Phi": phi(core),
                         "deleted_gates": 2*k, "Q_sum": total_Q,
                         "boundary_fanout_sum": charged_boundary,
                         "trace": trace if m == 2 and k == 3 else None})
    return rows


def hostile():
    # Inputs a,b,c,y are -4,-3,-2,-1; AND=8, XOR=6 in the declared convention.
    pins = {0:(-3,-2), 1:(0,-4), 2:(1,-1), 3:(0,2)}
    row = region_data(pins, 3, {1})
    check(row == {"r":1,"f":0,"E":0,"a":1,"kappa":1,"Q":1,"lambda":1,"C":[1],"B":[2]}, "hostile inventory")
    final, out, _ = remove_region(pins, 3, {1}, {2})
    check(final == {0:(-3,-2),3:(0,-1)}, "hostile final topology")
    old = truth_vectors(pins, {0:8,1:8,2:6,3:8}, (0,None,None,None), 4)
    new = truth_vectors(final, {0:8,3:8}, (0,None,None,None), 4)
    check(old[3] == new[3], "hostile Boolean equivalence")
    for g, pair in final.items():
        check(len(set(new[g])) == 2, "hostile surviving constant")
        for p in pair:
            check(new[g] != new[p] and new[g] != tuple(1-x for x in new[p]), "hostile surviving literal")
    return {"inventory":row, "initial_Phi":phi(pins),"final_Phi":phi(final),
            "deleted_gates":[1,2],"output_vector_a_zero":list(old[3]),
            "final_gates":{str(g):list(ps) for g,ps in final.items()}}


def main():
    report = {"status":"PASS", "scope":"Finite checks corroborating the displayed all-size proofs; no lower-bound claim.",
              "distinct_pin_topologies": enumerate_topologies(6, False),
              "repeated_pin_topologies": enumerate_topologies(5, True),
              "truth_table_universe":truth_universe(),
              "four_gate_hostile":hostile(),
              "grid_shielding":shielding_controls()}
    rendered = json.dumps(report, sort_keys=True, indent=2) + "\n"
    path = Path(__file__).with_name("RECEIPT.json")
    path.write_text(rendered)
    print(rendered)
    print("RECEIPT_SHA256=" + hashlib.sha256(rendered.encode()).hexdigest())

if __name__ == "__main__":
    main()
