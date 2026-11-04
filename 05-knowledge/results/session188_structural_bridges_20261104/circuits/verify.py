"""Finite regression checks for the constant-region incidence identity.

No third-party dependencies. This corroborates, rather than replaces, PROOF.md.
"""
from itertools import combinations, combinations_with_replacement, product
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def live(parents, inputs):
    seen = {inputs + len(parents) - 1}
    for v in range(inputs + len(parents) - 1, inputs - 1, -1):
        if v in seen:
            seen.update(parents[v - inputs])
    return set(range(inputs, inputs + len(parents))) <= seen


def incidence(parents, inputs, C, check_components=True):
    gates = set(range(inputs, inputs + len(parents)))
    require(C and C <= gates, "invalid constant region")
    require(inputs + len(parents) - 1 not in C, "constant output")
    children = {v: set() for v in gates}
    for i, ps in enumerate(parents):
        for p in set(ps):
            if p in gates:
                children[p].add(inputs + i)
    indeg = {v: len(set(parents[v-inputs]) & C) for v in gates}
    require(all(not set(ps) <= C or inputs+i in C
                for i, ps in enumerate(parents)), "region not closed")
    B = {v for v in gates - C if indeg[v]}
    require(all(indeg[v] == 1 for v in B), "boundary has two constant parents")
    require(all(children[v] for v in C), "dead constant vertex")
    r = sum(indeg[v] == 0 for v in C)
    f = sum(indeg[v] == 1 for v in C)
    E = sum(len(children[v])-1 for v in C)
    actual = len(C | B)
    require(actual == 2*r+f+E, "surplus identity failed")
    equality = actual == 2*r
    require(equality == (f == E == 0), "equality algebra failed")
    if equality and check_components:
        # A separate undirected traversal validates the claimed merge forest.
        unseen = set(C)
        boundaries = set()
        while unseen:
            root = next(iter(unseen))
            component, todo = set(), [root]
            while todo:
                v = todo.pop()
                if v in component:
                    continue
                component.add(v)
                todo.extend((set(parents[v-inputs]) | children[v]) & C)
            unseen -= component
            leaves = sum(indeg[v] == 0 for v in component)
            internal_edges = sum(indeg[v] for v in component)
            exits = {(v, w) for v in component for w in children[v] - C}
            require(len(component) == 2*leaves-1, "not a full binary tree")
            require(internal_edges == len(component)-1, "component not a tree")
            require(len(exits) == 1, "component has nonunique boundary")
            boundary = next(iter(exits))[1]
            require(boundary not in boundaries, "two trees share boundary")
            boundaries.add(boundary)
        require(boundaries == B, "boundary/component correspondence")
    return r, f, E, len(B), equality


def enumerate_topologies(repeated, max_gates):
    rows = []
    choose = combinations_with_replacement if repeated else combinations
    for n in range(1, max_gates+1):
        choices = [list(choose(range(i+2), 2)) for i in range(n)]
        total = kept = regions = equalities = forcing = branching = 0
        for parents in product(*choices):
            total += 1
            if not live(parents, 2):
                continue
            kept += 1
            for mask in range(1, 1 << (n-1)):
                C = {i+2 for i in range(n-1) if mask >> i & 1}
                if any(set(ps) <= C and i+2 not in C for i, ps in enumerate(parents)):
                    continue
                r, f, E, b, eq = incidence(parents, 2, C)
                regions += 1
                equalities += eq
                forcing += f > 0
                branching += E > 0
        rows.append({"gates": n, "all_topologies": total,
                     "live_topologies": kept, "closed_constant_regions": regions,
                     "merge_forest_equalities": equalities,
                     "regions_with_one_parent_constants": forcing,
                     "regions_with_fanout_surplus": branching})
    return rows


def truth_gate(a, b, table, full):
    value = 0
    for x, y in product((0, 1), repeat=2):
        if table >> (2*x+y) & 1:
            value |= (a if x else full ^ a) & (b if y else full ^ b)
    return value


def eval_table(table, x, y):
    return table >> (2*x+y) & 1


def normalized(parents, tables, fixed, actual, mask, input_truths, full):
    n = len(fixed)
    aliases = [("c", v) if v is not None else ("w", i, 0)
               for i, v in enumerate(fixed)]
    new = []
    for i, (ps, table) in enumerate(zip(parents, tables)):
        restricted = actual[i] & mask
        if restricted in (0, mask):
            aliases.append(("c", int(restricted == mask)))
            continue
        args = [aliases[p] for p in ps]
        roots = sorted({a[1] for a in args if a[0] == "w"})
        outputs = []
        for bits in product((0, 1), repeat=len(roots)):
            env = dict(zip(roots, bits))
            values = [a[1] if a[0] == "c" else env[a[1]] ^ a[2] for a in args]
            outputs.append(eval_table(table, *values))
        if len(set(outputs)) == 1:
            aliases.append(("c", outputs[0]))
        elif len(roots) == 1:
            aliases.append(("w", roots[0], outputs[0]))
        else:
            # The complete 4-value table is retained, including asymmetric gates.
            require(len(roots) == 2, "too many normalized parents")
            if outputs[0] == outputs[1] and outputs[2] == outputs[3]:
                aliases.append(("w", roots[0], outputs[0]))
            elif outputs[0] == outputs[2] and outputs[1] == outputs[3]:
                aliases.append(("w", roots[1], outputs[0]))
            else:
                newtable = sum(v << j for j, v in enumerate(outputs))
                aliases.append(("w", n+len(new), 0))
                new.append((tuple(roots), newtable))
    out = aliases[-1]
    require(out[0] == "w", "expected nonconstant output")
    ancestor, todo = set(), [out[1]]
    while todo:
        v = todo.pop()
        if v < n or v in ancestor:
            continue
        ancestor.add(v)
        todo.extend(new[v-n][0])
    signed = len(ancestor)
    strict = signed + int(out[1] < n and out[2] == 1)
    values = list(input_truths)
    for ps, table in new:
        values.append(truth_gate(values[ps[0]], values[ps[1]], table, full))
    result = values[out[1]] ^ (full if out[2] else 0)
    require((result ^ actual[-1]) & mask == 0, "normalization changed function")
    return signed, strict


def semantic_two_gate_census():
    n = 3
    full = (1 << (1 << n))-1
    inputs = [sum(((a >> i) & 1) << a for a in range(1 << n)) for i in range(n)]
    restrictions = []
    for fixed in product((None, 0, 1), repeat=n):
        if all(v is None for v in fixed) or all(v is not None for v in fixed):
            continue
        mask = sum(1 << a for a in range(1 << n)
                   if all(v is None or (a >> i) & 1 == v for i, v in enumerate(fixed)))
        restrictions.append((fixed, mask))
    circuits = nonconstant_cases = constant_region_cases = endpoint_losses = 0
    rows = []
    for s in (1, 2):
        before = circuits
        choices = [list(combinations_with_replacement(range(n+i), 2)) for i in range(s)]
        for parents in product(*choices):
            if not live(parents, n):
                continue
            for tables in product(range(16), repeat=s):
                circuits += 1
                values = list(inputs)
                for ps, table in zip(parents, tables):
                    values.append(truth_gate(values[ps[0]], values[ps[1]], table, full))
                actual = values[n:]
                for fixed, mask in restrictions:
                    if actual[-1] & mask in (0, mask):
                        continue
                    nonconstant_cases += 1
                    C = {n+i for i, val in enumerate(actual) if val & mask in (0, mask)}
                    if not C:
                        continue
                    constant_region_cases += 1
                    r, f, E, b, eq = incidence(parents, n, C)
                    D = 2*r+f+E
                    signed, strict = normalized(parents, tables, fixed, actual, mask, inputs, full)
                    require(s-signed >= D, "signed deletion count")
                    require(s-strict >= D-1, "strict endpoint bound")
                    endpoint_losses += s-strict < D
        rows.append({"gates": s, "live_circuit_presentations": circuits-before})
    require(endpoint_losses > 0, "terminal NOT hostile not detected")
    return {"universe": "all 16 gate tables, 3 inputs, unordered repeated pins, at most 2 gates",
            "coordinate_restrictions_per_circuit": len(restrictions),
            "circuits": circuits, "nonconstant_restriction_cases": nonconstant_cases,
            "constant_region_cases": constant_region_cases,
            "strict_terminal_correction_cases": endpoint_losses, "rows": rows}


def controls():
    tests = {}
    # Three external inputs a,b,x. Constant A, then a unary-forced AND, then XOR x.
    parents = ((0, 1), (3, 0), (4, 2))
    tests["comparable_chain"] = incidence(parents, 3, {3, 4})[:4]
    require(tests["comparable_chain"] == (1, 1, 0, 1), "chain control")
    # Two constant leaves merge into one constant, then meet x.
    parents = ((0, 1), (2, 3), (5, 6), (7, 4))
    tests["merge_equality"] = incidence(parents, 5, {5, 6, 7})[:4]
    require(tests["merge_equality"] == (2, 0, 0, 1), "merge control")
    # One constant fans out to two boundaries, which meet at the output.
    parents = ((0, 1), (3, 0), (3, 2), (4, 5))
    tests["fork_surplus"] = incidence(parents, 3, {3})[:4]
    require(tests["fork_surplus"] == (1, 0, 1, 2), "fork control")
    # Repeated pins make a genuine one-parent constant and cannot escape closure.
    parents = ((0, 1), (2, 2), (3, 0))
    tests["repeated_pin_constant"] = incidence(parents, 2, {2, 3})[:4]
    require(tests["repeated_pin_constant"] == (1, 1, 0, 1), "repeated control")
    # Full merge trees attain 2r for every displayed size.
    for r in range(1, 25):
        n = 2*r+1
        ps = [(2*i, 2*i+1) for i in range(r)]
        roots = list(range(n, n+r))
        while len(roots) > 1:
            u, v = roots.pop(0), roots.pop(0)
            roots.append(n+len(ps))
            ps.append((u, v))
        C = set(range(n, n+len(ps)))
        ps.append((roots[0], n-1))
        rr, f, E, b, eq = incidence(ps, n, C)
        require((rr, f, E, b, eq) == (r, 0, 0, 1, True), "sharp family")
    tests["sharp_merge_family_r"] = [1, 24]
    # Strict B2 hostile: (a AND b) NOR x, restricted a=0.
    ps, ts = ((0, 1), (3, 2)), (8, 1)
    n, full = 3, 255
    ins = [sum(((a >> i) & 1) << a for a in range(8)) for i in range(3)]
    vals = list(ins)
    for p, t in zip(ps, ts):
        vals.append(truth_gate(vals[p[0]], vals[p[1]], t, full))
    mask = sum(1 << a for a in range(8) if a & 1 == 0)
    signed, strict = normalized(ps, ts, (0, None, None), vals[3:], mask, ins, full)
    require((signed, strict) == (0, 1), "strict terminal hostile")
    tests["terminal_not_sizes"] = {"before": 2, "signed_after": signed, "strict_after": strict}
    # A nonliteral surviving output absorbs signs without a terminal charge.
    # ((a AND b) XOR x) XOR y becomes x XOR y when a=0.
    ps, ts = ((0, 1), (4, 2), (5, 3)), (8, 6, 6)
    n, full = 4, 65535
    ins = [sum(((a >> i) & 1) << a for a in range(16)) for i in range(n)]
    vals = list(ins)
    for p, t in zip(ps, ts):
        vals.append(truth_gate(vals[p[0]], vals[p[1]], t, full))
    mask = sum(1 << a for a in range(16) if a & 1 == 0)
    signed, strict = normalized(ps, ts, (0, None, None, None), vals[n:], mask, ins, full)
    require((signed, strict) == (1, 1), "two-essential-input control")
    tests["two_essential_input_sizes"] = {"before": 3, "signed_after": signed, "strict_after": strict}
    return tests


def main():
    report = {"status": "PASS", "claim": "constant-region source/forcing/fanout identity",
              "distinct_pin_topologies": enumerate_topologies(False, 6),
              "repeated_pin_topologies": enumerate_topologies(True, 5),
              "all_table_semantic_census": semantic_two_gate_census(),
              "hostile_and_sharp_controls": controls()}
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
