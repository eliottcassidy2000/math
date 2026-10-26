#!/usr/bin/env python3
"""Portable exact controls. Standard library; no assertions; deterministic JSON.

This verifies finite consequences, not the universal proofs in REPORT.md.
"""
import argparse
from collections import Counter
from itertools import combinations, permutations, product
import json
from math import factorial, isqrt


def need(condition, label):
    if not condition:
        raise ValueError(label)


def edge(a, b):
    return tuple(sorted((a, b)))


def components(vertices, edges):
    parent = {v: v for v in vertices}
    def root(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    for a, b in edges:
        parent[root(a)] = root(b)
    groups = {}
    for v in vertices:
        groups.setdefault(root(v), []).append(v)
    return sorted(sorted(g) for g in groups.values())


def forest(vertices, edges):
    return len(edges) == len(vertices) - len(components(vertices, edges))


def path_graph(vertices, edges):
    degree = Counter(v for e in edges for v in e)
    return (len(components(vertices, edges)) == 1 and
            sorted(degree.get(v, 0) for v in vertices) == [1, 1] + [2] * (len(vertices)-2))


def matchings(vertices):
    if not vertices:
        yield ()
        return
    a = vertices[0]
    for i, b in enumerate(vertices[1:], 1):
        for rest in matchings(vertices[1:i] + vertices[i+1:]):
            yield tuple(sorted((edge(a, b),) + rest))


def distinguishing_context(p, q):
    """Returns disjoint connector pairs: p+connectors path, q+connectors cyclic."""
    a, b = min(set(q)-set(p))
    mate = {a: b for e in p for a, b in (e, e[::-1])}
    used = {a, mate[a], b, mate[b]}
    connectors = [edge(a, b)]
    tail = mate[b]
    for x, y in p:
        if x not in used:
            connectors.append(edge(tail, x))
            used.update((x, y))
            tail = y
    return tuple(connectors)


def expanded_union(pairing, connectors):
    n = 2*len(pairing)
    result = list(pairing)
    for i, (a, b) in enumerate(connectors):
        result.extend((edge(a, n+i), edge(b, n+i)))
    return list(range(n+len(connectors))), result


def check_matching_carrier():
    rows = []
    total_contexts = 0
    for r in range(1, 6):
        vertices = tuple(range(2*r))
        bank = list(matchings(vertices))
        expected = factorial(2*r)//(2**r*factorial(r))
        need(len(bank) == expected, "matching count")
        pairs = tuple(combinations(vertices, 2))
        profiles = {tuple(e not in p for e in pairs) for p in bank}
        need(len(profiles) == expected, "guard profiles must distinguish")
        checked = 0
        if r <= 4:
            for p in bank:
                for q in bank:
                    if p == q:
                        continue
                    connectors = distinguishing_context(p, q)
                    v, pe = expanded_union(p, connectors)
                    _, qe = expanded_union(q, connectors)
                    need(len(connectors) == r-1, "context connector count")
                    need(path_graph(v, pe), "positive final-path context")
                    need(not forest(v, qe), "negative context has cycle")
                    checked += 1
        total_contexts += checked
        rows.append({"r": r, "states": len(bank), "guard_profiles": len(profiles),
                     "ordered_final_context_checks": checked})
    p = ((0, 1), (2, 3))
    q = ((0, 3), (1, 2))
    c = ((1, 2),)
    v, pe = expanded_union(p, c)
    _, qe = expanded_union(q, c)
    need(path_graph(v, pe) and not path_graph(v, qe), "same-marginal hostile")
    return {"rows": rows, "total_contexts": total_contexts,
            "hostile": {"p": p, "q": q, "connector": [1, 4, 2],
                        "positive_components": components(v, pe),
                        "negative_components": components(v, qe)}}


def square_edges(n):
    return {edge(a, b) for a in range(1, n+1) for b in range(a+1, n+1)
            if isqrt(a+b)**2 == a+b}


def hamilton_paths(vertices, edges):
    vertices = tuple(vertices)
    adj = {v: set() for v in vertices}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    if len(vertices) == 1:
        return [vertices]
    if len(components(vertices, edges)) != 1:
        return []
    leaves = [v for v in vertices if len(adj[v]) == 1]
    if len(leaves) > 2:
        return []
    starts = leaves[:1] if leaves else vertices
    found = set()
    def visit(path, remaining):
        if not remaining:
            p = tuple(path)
            found.add(min(p, p[::-1]))
            return
        for w in sorted(adj[path[-1]] & remaining):
            visit(path+[w], remaining-{w})
    for s in starts:
        visit([s], set(vertices)-{s})
    return sorted(found)


def check_square_sum():
    expected = {1: 1, 15: 1, 16: 1, 17: 1, 23: 3, 25: 10}
    paths = {}
    for n in range(1, 26):
        paths[n] = hamilton_paths(range(1, n+1), square_edges(n))
        need(len(paths[n]) == expected.get(n, 0), "Q_n path count at "+str(n))
    hostile = next(p for p in paths[23] if 22 in (p[0], p[-1]))
    need(sum(22 in e for e in square_edges(23)) == 2, "degree-two endpoint")
    need(edge(14, 22) not in {edge(a,b) for a,b in zip(hostile, hostile[1:])},
         "degree-two endpoint omits one edge")
    blocks = [{2,7,9,16,18}, {4,5,11,12}, {1,3,6,8,10,17}]
    contracted = {v: i for i, block in enumerate(blocks) for v in block}
    next_id = 3
    for v in range(1,26):
        if v not in contracted:
            contracted[v] = next_id
            next_id += 1
    qedges = [edge(contracted[a],contracted[b]) for a,b in sorted(square_edges(25))
              if contracted[a] != contracted[b]]
    core = set(range(3, next_id))
    core_edges = [e for e in qedges if all(x in core for x in e)]
    need(next_id == 13 and len(qedges) == 20, "Q25 exact multigraph contraction")
    need(len(components(core, core_edges)) == 6, "Q25 quotient cut certificate")
    used = Counter(e for p in paths[25] for e in (edge(a,b) for a,b in zip(p,p[1:])))
    need(used[edge(1,15)] == 0, "birth class is not forced edge")
    return {"counts": {str(n):len(p) for n,p in paths.items()},
            "degree_two_endpoint_hostile": hostile,
            "Q25_path": paths[25][0], "Q25_endpoints": [[p[0],p[-1]] for p in paths[25]],
            "Q25_tree_block_quotient": {"vertices":next_id,"edges_with_multiplicity":len(qedges),
                 "cut_size":3,"components_after_cut":len(components(core,core_edges))}}


def reconfiguration_components(states, width):
    vertices = list(range(len(states)))
    edges = [(i,j) for i in vertices for j in range(i+1,len(states))
             if sum(a != b for a,b in zip(states[i], states[j])) <= width]
    return sorted(len(c) for c in components(vertices,edges))


def check_graceful():
    fixed = [(5,6),(1,3),(8,11),(3,7),(0,5),(4,11),(1,9),(0,10),(0,12)]
    colors = (6,9,11)
    choices = [[(a,a+d) for a in range(13-d)] for d in colors]
    states = [s for s in product(*choices) if forest(list(range(13)), fixed+list(s))]
    need(len(states)==5, "m12 fiber size")
    need(reconfiguration_components(states,2)==[2,3], "m12 width-two disconnect")
    for s in states:
        need(sorted(b-a for a,b in fixed+list(s))==list(range(1,13)), "rainbow differences")
    comp = components(range(13),fixed)
    owner = {v:i for i,c in enumerate(comp) for v in c}
    contacts = {str(d):sorted({edge(owner[a],owner[b]) for a,b in cs if owner[a]!=owner[b]})
                for d,cs in zip(colors,choices)}
    factorial_rows = []
    for k in range(1,5):
        m = k*k+2*k
        special = {i*(k+1):i for i in range(1,k+1)}
        f = [(special[d],d+special[d]) if d in special else (0,d)
             for d in range(k+1,m+1)]
        cs = components(range(m+1),f)
        owner = {v:i for i,c in enumerate(cs) for v in c}
        free_choices = [[(a,a+d) for a in range(m+1-d) if owner[a]!=owner[a+d]]
                        for d in range(1,k+1)]
        need(all(len(x)==2*k for x in free_choices), "factorial nonloop contacts")
        completions = [s for s in product(*free_choices)
                       if forest(list(range(m+1)),f+list(s))]
        need(len(completions)==factorial(k)*2**k, "factorial completion count")
        one = reconfiguration_components(completions,1)
        two = reconfiguration_components(completions,2)
        need(one==[2**k]*factorial(k), "factorial one-color components")
        need(two==[factorial(k)*2**k], "factorial width-two connectivity")
        # Every request is a coordinate overwrite; no-op replacements are legal.
        # Compute guards on actual arithmetic edges, not predicted spoke types.
        profiles = Counter()
        for state in completions:
            profile = []
            for d, options in enumerate(free_choices):
                for replacement in options:
                    after = list(state)
                    after[d] = replacement
                    profile.append(forest(list(range(m+1)), f+after))
            profiles[tuple(profile)] += 1
        need(len(profiles)==factorial(k), "guard quotient recovers permutation only")
        need(set(profiles.values())=={2**k}, "left/right signs are guard-invisible")
        factorial_rows.append({"k":k,"m":m,"states":len(completions),
             "one_color_component_sizes":one,"two_color_component_sizes":two,
             "one_color_guard_profiles":len(profiles),"guard_fiber_size":2**k})
    return {"m12_fixed":fixed,"m12_components":comp,"m12_free_colors":colors,
            "m12_contact_supports":contacts,"m12_candidate_count":7*4*2,
            "m12_states":states,"m12_width_two_components":[2,3],
            "factorial_rows":factorial_rows}


def petersen_edges():
    return ({edge(i,(i+1)%5) for i in range(5)} |
            {edge(i+5,(i+2)%5+5) for i in range(5)} |
            {(i,i+5) for i in range(5)})


def hamilton_cycles(vertices, edges):
    vertices = tuple(vertices)
    start = min(vertices)
    return [p for p in permutations([v for v in vertices if v!=start])
            if p[0] < p[-1] and all(edge(a,b) in edges
                 for a,b in zip((start,)+p,p+(start,)))]


def check_petersen():
    e = petersen_edges()
    perfect = [set(s) for s in combinations(sorted(e),5)
               if Counter(v for a in s for v in a)==Counter(range(10))]
    need(len(perfect)==6, "Petersen perfect matchings")
    complements = [components(range(10),e-m) for m in perfect]
    need(all(sorted(map(len,c))==[5,5] for c in complements), "all two-factors are 5+5")
    paths = hamilton_paths(range(10),e)
    need(bool(paths), "Petersen traceability")
    k5 = {edge(a%5,b%5) for a,b in e if a%5!=b%5}
    need(k5==set(combinations(range(5),2)), "five-spoke contraction gives K5")
    # Delete outer vertex 0. Its former neighbors 1,4,5 have degree two.
    branch = {2,3,6,7,8,9}
    subdivided = [(2,1,6),(3,4,9),(7,5,8)]
    remaining = {x for x in e if 0 not in x}
    reduced = {x for x in remaining if all(v in branch for v in x)}
    for p in subdivided:
        need(all(edge(a,b) in remaining for a,b in zip(p,p[1:])),"subdivision path")
        reduced.add(edge(p[0],p[-1]))
    bipartitions = [(set(s),branch-set(s)) for s in combinations(sorted(branch),3)
                     if reduced=={edge(a,b) for a in s for b in branch-set(s)}]
    need(len(bipartitions)==2,"Petersen minus vertex is K33 subdivision")
    reduced_cycles = hamilton_cycles(sorted(branch), reduced)
    marks = {edge(p[0],p[-1]) for p in subdivided}
    expanded_marked = []
    for p in reduced_cycles:
        cyc = (min(branch),)+p
        edges = {edge(a,b) for a,b in zip(cyc,cyc[1:]+cyc[:1])}
        if marks <= edges:
            expanded = edges-marks
            for path in subdivided:
                expanded |= {edge(a,b) for a,b in zip(path,path[1:])}
            expanded_marked.append(expanded)
    actual_deleted_cycles = hamilton_cycles(range(1,10),remaining)
    actual_deleted_sets = []
    for p in actual_deleted_cycles:
        cyc = (1,)+p
        actual_deleted_sets.append({edge(a,b) for a,b in zip(cyc,cyc[1:]+cyc[:1])})
    need(len(reduced_cycles)==6 and len(expanded_marked)==2, "K33 marked cycle counts")
    need({frozenset(s) for s in expanded_marked}=={frozenset(s) for s in actual_deleted_sets},
         "marked cycle expansion equals all Petersen-minus-0 cycles")
    # Local two-port rule on each contracted spoke: select equally many outer
    # and inner nonspoke edges, either one of each (+spoke) or two of each (-spoke).
    quotient_edges = sorted(k5)
    lifted = []
    for mask in range(1<<10):
        selected = {quotient_edges[i] for i in range(10) if mask>>i&1}
        out = [sum(edge(i,j) in selected for j in ((i-1)%5,(i+1)%5)) for i in range(5)]
        inn = [sum(edge(i,j) in selected for j in ((i-2)%5,(i+2)%5)) for i in range(5)]
        if not all(out[i]==inn[i] and out[i] in (1,2) for i in range(5)):
            continue
        expanded = {x for x in e if x[0]%5!=x[1]%5 and edge(x[0]%5,x[1]%5) in selected}
        expanded |= {(i,i+5) for i in range(5) if out[i]==1}
        need(all(sum(v in x for x in expanded)==2 for v in range(10)),"local rule gives 2-factor")
        lifted.append(expanded)
    need({frozenset(x) for x in lifted}=={frozenset(e-m) for m in perfect}, "response rule exact")
    k5cycles = hamilton_cycles(range(5),k5)
    k33 = {edge(a,b) for a in range(3) for b in range(3,6)}
    k33cycles = hamilton_cycles(range(6),k33)
    need(len(k5cycles)==12 and len(k33cycles)==6,"nonplanar graphs have cycles")
    labels33 = [0,1,2,3,6,9]
    need(sorted(abs(labels33[a]-labels33[b]) for a,b in k33)==list(range(1,10)),
         "K33 graceful positive control")
    # Any graceful K5 must contain labels0 and10; exhaust the84 possible
    # sorted remaining3 labels, independently of the parity proof.
    graceful5 = []
    for mid in combinations(range(1,10),3):
        labels = (0,)+mid+(10,)
        if sorted(b-a for a,b in combinations(labels,2))==list(range(1,11)):
            graceful5.append(labels)
    need(not graceful5,"K5 graceful negative control")
    # A Hamilton quotient cycle cannot lift: its local degrees must be one
    # outer plus one inner, so the five cycle edges would alternate two colors.
    balanced_cycles = [p for p in k5cycles if all(
          ((b-a)%5 in (1,4)) != ((c-b)%5 in (1,4))
          for a,b,c in zip((0,)+p,p+(0,),p[1:]+(0,p[0])))]
    need(not balanced_cycles,"odd quotient cycle cannot alternate port colors")
    return {"edges":sorted(e),"perfect_matchings":[sorted(m) for m in perfect],
            "two_factor_component_sizes":[sorted(map(len,c)) for c in complements],
            "hamilton_path_count_unoriented":len(paths),"hamilton_path":paths[0],
            "K5_cycles":len(k5cycles),"K33_cycles":len(k33cycles),
            "K33_graceful_labels":labels33,"K5_graceful_label_sets":graceful5,
            "K5_normalized_label_sets_checked":84,
            "K33_branch_parts":[sorted(s) for s in bipartitions[0]],
            "K33_subdivided_paths":subdivided,
            "K33_cycles_using_all_suppression_marks":len(expanded_marked),
            "Petersen_minus_zero_cycles":len(actual_deleted_sets),
            "local_two_port_rule_two_factors":len(lifted),"balanced_K5_cycles":len(balanced_cycles)}


def check_hostile_mutation(result):
    try:
        need(result["graceful"]["m12_width_two_components"]==[5], "reject lost fixed-color sidecar")
    except ValueError:
        return "rejected"
    raise ValueError("hostile mutation was accepted")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hostile", action="store_true", help="exit nonzero on intentionally false receipt claim")
    args = parser.parse_args()
    result = {"schema":"session70.graph-alignment.v1",
              "matching":check_matching_carrier(),"square_sum":check_square_sum(),
              "graceful":check_graceful(),"petersen":check_petersen()}
    result["false_claim_negative_control"] = check_hostile_mutation(result)
    if args.hostile:
        need(result["petersen"]["balanced_K5_cycles"]==12,
             "intentional false claim: bare K5 cycles lift through Petersen ports")
    print(json.dumps(result,sort_keys=True,indent=2))


if __name__ == "__main__":
    main()
