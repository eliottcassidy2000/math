#!/usr/bin/env python3
"""Exact, standard-library controls for relational full-witness repair theory.

All fixtures are declared finite abstract record models, not a native census.
No input is silently authenticated as a tournament record by this checker.
"""
from collections import deque, defaultdict
from hashlib import sha256
from itertools import product, combinations, permutations
from math import factorial
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def powerset(items):
    items = tuple(items)
    for mask in range(1 << len(items)):
        yield frozenset(items[i] for i in range(len(items)) if mask >> i & 1)


def exact_flow(nleft, nright, edges, d, e, middle_capacity=1):
    """Integral residual max-flow, cut and actual selected full pair edges."""
    source, sink = nleft+nright, nleft+nright+1
    n = sink+1
    cap = [[0]*n for _ in range(n)]
    for x in range(nleft):
        cap[source][x] = d[x]
    for x, z in edges:
        cap[x][nleft+z] = middle_capacity
    for z in range(nright):
        cap[nleft+z][sink] = e[z]
    residual = [row[:] for row in cap]
    value = 0
    while True:
        parent = {source: None}
        queue = deque([source])
        while queue and sink not in parent:
            u = queue.popleft()
            for v in range(n):
                if residual[u][v] > 0 and v not in parent:
                    parent[v] = u
                    queue.append(v)
        if sink not in parent:
            break
        increment = sum(d)+sum(e)+middle_capacity+1
        v = sink
        while parent[v] is not None:
            u = parent[v]
            increment = min(increment, residual[u][v])
            v = u
        v = sink
        while parent[v] is not None:
            u = parent[v]
            residual[u][v] -= increment
            residual[v][u] += increment
            v = u
        value += increment
    reached = {source}
    queue = deque([source])
    while queue:
        u = queue.popleft()
        for v in range(n):
            if residual[u][v] > 0 and v not in reached:
                reached.add(v)
                queue.append(v)
    cut = sum(cap[u][v] for u in reached for v in range(n) if v not in reached)
    need(cut == value, 'residual reachable set is exact minimum cut')
    used = {(x,z): cap[x][nleft+z]-residual[x][nleft+z] for x,z in edges}
    return value, frozenset(x for x in range(nleft) if x in reached), frozenset(z for z in range(nright) if nleft+z in reached), used


def source_pins():
    pins = json.loads((HERE/'SOURCEPINS.json').read_text())
    actual = {p.relative_to(HERE).as_posix() for p in (HERE/'inputs').rglob('*') if p.is_file()}
    need(actual == set(pins), 'source pins close exactly over input archives')
    for name, row in pins.items():
        data = (HERE/name).read_bytes()
        need(len(data) == row['bytes'] and sha256(data).hexdigest() == row['sha256'], 'byte-identical source archive '+name)
        need(name.endswith('.txt') and row['source_path'].endswith('.md'), 'historical Markdown stored as archival text')
        need(row['source_commit'] == 'eb085bb428517b88bf6881cdab6d4bd980e305a3', 'fixed exact source commit')
    return {'archives':len(pins), 'commit':'eb085bb428517b88bf6881cdab6d4bd980e305a3'}


def branch_controls():
    # Entire arbitrary predicates on two receiver bits and one repair bit.
    X = list(product((0,1), repeat=2))
    Z = (0,1)
    universe = tuple(product(X,Z))
    cases = 0
    witness_checks = 0
    for modulus in range(2,6):
        Q = set(range(modulus))
        charge_x = lambda x: (x[0]+2*x[1]) % modulus
        charge_z = lambda z: (3*z) % modulus
        for relation in powerset(universe):
            for grade in range(4):
                E = {(x,z) for x,z in relation if sum(x)+z == grade}
                branch = {z:defaultdict(list) for z in Z}
                for x,z in E:
                    branch[z][charge_x(x)].append(x)
                A = {z:set(branch[z]) for z in Z}
                H = {z:Q-A[z] for z in Z}
                actual = {(charge_x(x)+charge_z(z))%modulus for x,z in E}
                union = set().union(*({(charge_z(z)+a)%modulus for a in A[z]} for z in Z))
                intersection = Q.copy()
                for z in Z:
                    intersection &= {(charge_z(z)+h)%modulus for h in H[z]}
                need(actual == union and Q-actual == intersection, 'arbitrary relation-conditioned graded branch identity')
                forbidden = {(h0-h1)%modulus for h0 in H[0] for h1 in H[1]}
                need((len(actual) != modulus) == (((charge_z(1)-charge_z(0))%modulus) in forbidden), 'conditioned two-branch difference criterion')
                for q in actual:
                    z = next(z for z in Z if (q-charge_z(z))%modulus in branch[z])
                    x = branch[z][(q-charge_z(z))%modulus][0]
                    need((x,z) in E and (charge_x(x)+charge_z(z))%modulus == q, 'section reconstructs actual full tuple, branch and grade')
                    witness_checks += 1
                cases += 1
    joint = 0
    branch_deletions = 0
    for modulus in range(1,4):
        Q = frozenset(range(modulus))
        holes = list(powerset(Q))
        for H in product(holes, repeat=3):
            for c1,c2 in product(Q,repeat=2):
                c = (0,c1,c2)
                translated = [{(c[j]+h)%modulus for h in H[j]} for j in range(3)]
                common = set.intersection(*translated)
                D = {((h0-h1)%modulus,(h0-h2)%modulus) for h0,h1,h2 in product(*H)}
                need(bool(common) == ((c1,c2) in D), 'three-branch joint difference obstruction')
                for q in Q:
                    covered = {j for j in range(3) if (q-c[j])%modulus not in H[j]}
                    costs = (1,2,4)
                    least = min((sum(costs[j] for j in banned) for banned in powerset(range(3)) if covered <= banned), default=None)
                    need(least == sum(costs[j] for j in covered), 'weighted whole-branch destructive margin')
                    branch_deletions += 1
                joint += 1
    H = [{0,1},{1,2},{0,2}]
    need(not set.intersection(*H) and all(H[a]&H[b] for a,b in combinations(range(3),2)), 'pairwise hole intersections do not imply a joint hole')
    return {'all_arbitrary_8_pair_predicates':256, 'graded_predicate_cases':cases, 'literal_section_witnesses':witness_checks, 'joint_three_branch_cases':joint, 'weighted_branch_deletion_queries':branch_deletions}


def cover_controls():
    X = tuple(range(3))
    Z = tuple(range(3))
    edges = tuple(product(X,Z))
    weighted = 0
    rectangles = 0
    graphs = 0
    schemes = [((1,1,1),(1,1,1)),((1,2,3),(3,2,1)),((0,2,1),(1,0,3))]
    maximum_cover = 0
    for E in powerset(edges):
        cover_pairs = [(U,V) for U in powerset(X) for V in powerset(Z) if all(x in U or z in V for x,z in E)]
        for a,b in schemes:
            brute = min(sum(a[x] for x in U)+sum(b[z] for z in V) for U,V in cover_pairs)
            value,U,V,_ = exact_flow(3,3,E,a,b,sum(a)+sum(b)+1)
            need(value == brute, 'weighted whole-assignment witness vertex cover')
            cover_left = set(X)-set(U)
            need(all(x in cover_left or z in V for x,z in E), 'cut reconstructs whole full-assignment ban')
            need(sum(a[x] for x in cover_left)+sum(b[z] for z in V) == value, 'actual cover has exact claimed cost')
            weighted += 1
        for a,b in product(range(4),repeat=2):
            direct = any(len(U)<=a and len(V)<=b for U,V in cover_pairs)
            formula = any(len(U)<=a and len({z for x,z in E if x not in U})<=b for U in powerset(X))
            need(direct == formula, 'rectangular separate whole-assignment deletion budgets')
            rectangles += 1
        maximum_cover = max(maximum_cover,min(len(U)+len(V) for U,V in cover_pairs))
        graphs += 1
    # Three alternative full kernel-lift records of a single resource rho.
    E = {(0,0),(1,1),(2,2)}
    value,*_ = exact_flow(3,3,E,(1,1,1),(1,1,1),7)
    atoms = ('rho','x0','x1','x2')
    vulnerable = [{'rho','x'+str(x)} for x,z in E]
    atom_covers = [A for A in powerset(atoms) if all(A & D for D in vulnerable)]
    best = min(atom_covers,key=lambda A:(len(A),sorted(A)))
    need(value == 3 and best == {'rho'}, 'kernel-lift grouped resource ban is not record vertex cover')
    # Exhaust arbitrary witness-atom incidence, including invulnerable witnesses.
    atom_cases = 0
    atoms = (0,1,2)
    for vulnerabilities in product(list(powerset(atoms)),repeat=3):
        for bans in powerset(atoms):
            survivors = [j for j,D in enumerate(vulnerabilities) if not (bans & D)]
            need((not survivors) == all(bans & D for D in vulnerabilities), 'exact witness-atom transversal identity')
            atom_cases += 1
    return {'all_3_by_3_witness_graphs':graphs, 'weighted_cover_cases':weighted, 'separate_budget_cases':rectangles, 'maximum_unweighted_cover':maximum_cover, 'grouped_resource_hostile':{'record_vertex_cover':value,'resource_atom_cover':sorted(best)}, 'exhaustive_atom_incidence_queries':atom_cases}


def degree_signature(nleft,nright,R):
    return (tuple(sum(x==i for x,z in R) for i in range(nleft)),tuple(sum(z==j for x,z in R) for j in range(nright)))


def max_deficit(nleft,nright,E,d,e):
    best = None
    for U in powerset(range(nleft)):
        for V in powerset(range(nright)):
            deficit = sum(d[x] for x in U)-sum(e[z] for z in V)-sum(x in U and z not in V for x,z in E)
            if best is None or deficit > best[0]:
                best = (deficit,tuple(sorted(U)),tuple(sorted(V)))
    return best


def quota_controls():
    universe = tuple(product(range(2),range(3)))
    quota_cases = 0
    robust_cases = 0
    feasible_cases = 0
    avoiding_relations = 0
    for E in powerset(universe):
        tables = defaultdict(list)
        for R in powerset(E):
            tables[degree_signature(2,3,R)].append(R)
        for d in product(range(4),repeat=2):
            for e in product(range(3),repeat=3):
                if sum(d) != sum(e):
                    continue
                M = sum(d)
                value,U,V,used = exact_flow(2,3,E,d,e)
                direct = (d,e) in tables
                deficit = max_deficit(2,3,E,d,e)
                need((value == M) == direct == (deficit[0] <= 0), 'quota feasibility direct subsets, exact flow and all cuts')
                quota_cases += 1
                if not direct:
                    continue
                feasible_cases += 1
                for modulus in (2,3):
                    charge = lambda pair:(pair[0]+pair[1])%modulus
                    for q in range(modulus):
                        E_not = frozenset(pair for pair in E if charge(pair) != q)
                        avoiding = [R for R in tables[(d,e)] if all(charge(pair)!=q for pair in R)]
                        value,U,V,used = exact_flow(2,3,E_not,d,e)
                        deficit = max_deficit(2,3,E_not,d,e)
                        need((not avoiding) == (value<M) == (deficit[0]>0), 'unavoidable charge iff positive Hall deficit')
                        if value == M:
                            R = frozenset(pair for pair,n in used.items() if n)
                            need(all(n in (0,1) for n in used.values()), 'relation table is a subset, not repeated pair capacity')
                            need(degree_signature(2,3,R)==(d,e) and R<=E_not, 'flow reconstructs actual q-avoiding full-pair relation')
                            avoiding_relations += 1
                        else:
                            witness = sum(d[x] for x in U)-sum(e[z] for z in V)-sum(x in U and z not in V for x,z in E_not)
                            need(witness == M-value > 0, 'returned exact unavoidable-charge cut witness')
                        robust_cases += 1
    # Positive full-table fixture, all 24 arbitrary matching predicates.
    left = (0,1,1,1)
    right = (0,0,1,1)
    positive = []
    for p in permutations(range(4)):
        charges = [(left[x]+right[p[x]])%2 for x in range(4)]
        need(set(charges)=={0,1}, 'every four-by-four degree-one relation retains both charges')
        positive.append(charges.count(0))
    E = set(product(range(4),repeat=2))
    cuts = {}
    for q in (0,1):
        E_not = {pair for pair in E if (left[pair[0]]+right[pair[1]])%2 != q}
        deficit,U,V = max_deficit(4,4,E_not,(1,)*4,(1,)*4)
        need(deficit == 1, 'positive full-table charge has exact deficit-one certificate')
        cuts[q] = {'deficit':deficit,'left_rows':U,'right_columns':V}
    hostile = []
    for p in permutations(range(2)):
        hostile.append(sorted({(x+p[x])%2 for x in range(2)}))
    need(hostile == [[0],[1]], 'full unary projections alone do not imply charge coverage')
    # At-least degree quotas have no hidden matching upper bound.
    lower_cases = 0
    for E in powerset(universe):
        for d in product(range(2),repeat=2):
            for e in product(range(2),repeat=3):
                direct = any(all(sum(x==i for x,z in R)>=d[i] for i in range(2)) and all(sum(z==j for x,z in R)>=e[j] for j in range(3)) for R in powerset(E))
                full = all(sum(x==i for x,z in E)>=d[i] for i in range(2)) and all(sum(z==j for x,z in E)>=e[j] for j in range(3))
                need(direct == full, 'mere lower quotas are tested by the whole remaining relation')
                lower_cases += 1
    return {'all_2_by_3_compatibility_graphs':64, 'balanced_exact_quota_cases':quota_cases, 'nonempty_quota_families':feasible_cases, 'charge_robustness_cases':robust_cases, 'constructed_avoiding_relations':avoiding_relations, 'positive_matching_tables':24, 'positive_zero_charge_counts':dict(sorted((i,positive.count(i)) for i in set(positive))), 'positive_hall_cuts':cuts, 'full_projection_hostiles':hostile, 'lower_quota_cases':lower_cases}


def binary_marginal_controls():
    cases = 0
    matchings = 0
    minimum_controls = []
    for n in range(1,8):
        perms = list(permutations(range(n)))
        for a,b in product(range(n+1),repeat=2):
            left = (0,)*a+(1,)*(n-a)
            right = (0,)*b+(1,)*(n-b)
            observed = set()
            counts = defaultdict(int)
            for p in perms:
                zero_pairs = sum(left[x]==right[p[x]] for x in range(n))
                both_zero = sum(left[x]==0 and right[p[x]]==0 for x in range(n))
                need(zero_pairs==n-a-b+2*both_zero, 'actual named pairing has the exact four-type charge count')
                observed.add(zero_pairs)
                counts[zero_pairs] += 1
                matchings += 1
            expected = {n-a-b+2*t for t in range(max(0,a+b-n),min(a,b)+1)}
            expected_counts = {}
            for t in range(max(0,a+b-n),min(a,b)+1):
                numerator = factorial(a)*factorial(n-a)*factorial(b)*factorial(n-b)
                denominator = factorial(t)*factorial(a-t)*factorial(b-t)*factorial(n-a-b+t)
                need(numerator%denominator==0, 'named-table counting coefficient is integral')
                expected_counts[n-a-b+2*t] = numerator//denominator
            need(dict(counts)==expected_counts and sum(counts.values())==factorial(n), 'exact marked-table two-charge polynomial multiplicities')
            need(observed==expected, 'all possible binary marginal charge counts are attained')
            need(min(observed)==abs(n-a-b) and n-max(observed)==abs(a-b), 'sharp two-charge witness minima from all named pairings')
            need((all(0<k<n for k in observed)) == (a!=b and a+b!=n), 'exact all-size binary-marginal anti-collapse condition')
            cases += 1
        minimum_controls.append({'n':n,'pairings_per_marginal':len(perms)})
    return {'orders':list(range(1,8)), 'marginal_cases':cases, 'actual_named_pairing_instances':matchings, 'pairing_counts':minimum_controls}


def components(vertices, adjacency):
    unseen = set(vertices)
    groups = []
    while unseen:
        start = min(unseen)
        group = {start}
        queue = deque([start])
        unseen.remove(start)
        while queue:
            u = queue.popleft()
            for v in adjacency[u]:
                if v in unseen:
                    unseen.remove(v)
                    group.add(v)
                    queue.append(v)
        groups.append(frozenset(group))
    return groups


def path_in(vertices, adjacency, source, target):
    allowed = set(vertices)
    parent = {source:None}
    queue = deque([source])
    while queue:
        u = queue.popleft()
        if u == target:
            route = []
            while u is not None:
                route.append(u)
                u = parent[u]
            return list(reversed(route))
        for v in sorted(adjacency[u]):
            if v in allowed and v not in parent:
                parent[v] = u
                queue.append(v)
    return None


def budget_controls():
    states = list(product((0,1),repeat=3))
    weights = (1,3,6)
    costs = {(a,b):(sum(x!=y for x,y in zip(a,b)),sum(w for x,y,w in zip(a,b,weights) if x!=y)) for a,b in product(states,repeat=2)}
    graph_cases = 0
    endpoint_cases = 0
    lifted_routes = 0
    for F in powerset(states):
        for S,T in product(range(4),(0,1,3,4,6,7,9,10)):
            adj = {a:{b for b in F if b!=a and costs[a,b][0]<=S and costs[a,b][1]<=T} for a in F}
            fullgroups = components(F,adj)
            fullid = {a:i for i,C in enumerate(fullgroups) for a in C}
            fibrecomps = []
            for boundary in (0,1):
                fibre = {a for a in F if a[0]==boundary}
                fibrecomps.extend(components(fibre,{a:adj[a]&fibre for a in fibre}))
            fid = {a:i for i,C in enumerate(fibrecomps) for a in C}
            qadj = {i:set() for i in range(len(fibrecomps))}
            bridges = {}
            for a in F:
                for b in adj[a]:
                    if fid[a] != fid[b]:
                        qadj[fid[a]].add(fid[b])
                        bridges.setdefault((fid[a],fid[b]),(a,b))
            quotientgroups = components(qadj,qadj)
            qid = {i:j for j,C in enumerate(quotientgroups) for i in C}
            for a,b in product(F,repeat=2):
                agrees = fullid[a]==fullid[b]
                need(agrees == (qid[fid[a]]==qid[fid[b]]), 'split full-boundary components preserve simultaneous budget connectivity')
                endpoint_cases += 1
                if not agrees:
                    continue
                qroute = path_in(qadj,qadj,fid[a],fid[b])
                route = [a]
                for i,j in zip(qroute,qroute[1:]):
                    u,v = bridges[i,j]
                    inner = path_in(fibrecomps[i],adj,route[-1],u)
                    route.extend(inner[1:])
                    route.append(v)
                inner = path_in(fibrecomps[fid[b]],adj,route[-1],b)
                route.extend(inner[1:])
                need(route[0]==a and route[-1]==b and all(v in adj[u] for u,v in zip(route,route[1:])), 'literal source/target and every paid event in lifted route')
                lifted_routes += 1
            graph_cases += 1
    # Pareto hostile: two real endpoints, incompatible separate minima.
    vectors = {(1,10),(2,2)}
    need(min(x for x,y in vectors)==1 and min(y for x,y in vectors)==2 and not any(x<=1 and y<=2 for x,y in vectors), 'coordinate minima do not witness a simultaneous budget')
    hostile_states = {(0,0,0),(1,0,0),(1,1,1),(0,1,1)}
    adj = {a:{b for b in hostile_states if sum(x!=y for x,y in zip(a,b))==1} for a in hostile_states}
    need(len(components(hostile_states,adj))==2, 'unsplit intermediate fibre cannot concatenate different crossing witnesses')
    need(any(b in adj[(0,0,0)] for b in {(1,0,0),(1,1,1)}) and any(b in adj[(0,1,1)] for b in {(1,0,0),(1,1,1)}), 'false unsplit quotient has both crossing edges')
    return {'all_3_bit_accepted_sets':256, 'simultaneous_budget_graphs':graph_cases, 'endpoint_connectivity_tests':endpoint_cases, 'literal_lifted_routes':lifted_routes, 'pareto_hostile':sorted(vectors), 'unsplit_quotient_components':2}


def directed_controls():
    vertices = tuple(range(3))
    edges = tuple((a,b) for a in vertices for b in vertices if a!=b)
    cases = 0
    endpoint_cases = 0
    for E in powerset(edges):
        adj = {a:{b for x,b in E if x==a} for a in vertices}
        for labels in product(range(3),repeat=3):
            blocks = []
            for boundary in sorted(set(labels)):
                fibre = {a for a in vertices if labels[a]==boundary}
                unseen = set(fibre)
                while unseen:
                    a = min(unseen)
                    C = {b for b in fibre if path_in(fibre,adj,a,b) is not None and path_in(fibre,adj,b,a) is not None}
                    blocks.append(C)
                    unseen -= C
            bid = {a:i for i,C in enumerate(blocks) for a in C}
            quotient = {i:set() for i in range(len(blocks))}
            for a,b in E:
                if bid[a]!=bid[b]:
                    quotient[bid[a]].add(bid[b])
            for a,b in product(vertices,repeat=2):
                direct = path_in(vertices,adj,a,b) is not None
                lifted = path_in(quotient,quotient,bid[a],bid[b]) is not None
                need(direct==lifted, 'directed fibre SCC quotient preserves directed reachability')
                endpoint_cases += 1
            cases += 1
    # Weakly connected B={1,2} has only 2->1, so arrival at 1 cannot reach exit 2.
    adj = {0:{1},1:set(),2:{1,3},3:set()}
    need(path_in(adj,adj,0,3) is None and adj[0]=={1} and 3 in adj[2], 'weak directed fibre connectivity cannot concatenate crossings')
    return {'all_3_vertex_digraphs':64, 'labelled_fibre_partitions_per_graph':27, 'directed_quotient_cases':cases, 'directed_endpoint_tests':endpoint_cases}


def strict_record_search(F,source,target,budget):
    initial = (source,(0,)*len(source))
    queue = deque([initial])
    parent = {initial:None}
    while queue:
        state = queue.popleft()
        u,departed = state
        if u == target:
            route = []
            cursor = state
            while cursor is not None:
                route.append(cursor[0])
                cursor = parent[cursor]
            return list(reversed(route)),len(parent)
        for v in sorted(F):
            changed = [i for i in range(len(u)) if u[i]!=v[i]]
            if not changed or len(changed)>budget or any(departed[i]>>v[i]&1 for i in changed):
                continue
            history = tuple(departed[i] | ((1<<u[i]) if i in changed else 0) for i in range(len(u)))
            next_state = (v,history)
            if next_state not in parent:
                parent[next_state] = state
                queue.append(next_state)
    return None,len(parent)


def nonreturn_controls():
    a = (1,0,0,1,1)
    b = (0,1,1,0,1)
    F = {tuple((r*x+s*y)%3 for x,y in zip(a,b)) for r,s in product(range(3),repeat=2)}
    source = (0,0,0,0,0)
    target = (1,2,2,1,0)
    adj = {u:{v for v in F if 0<sum(x!=y for x,y in zip(u,v))<=3} for u in F}
    normal = path_in(F,adj,source,target)
    blocked,explored3 = strict_record_search(F,source,target,3)
    strict4,explored4 = strict_record_search(F,source,target,4)
    need(len(F)==9 and len(components(F,adj))==1 and normal is not None, 'nine-state code ordinary budget-three connected')
    need(blocked is None and strict4==[source,target], 'full-record no-return fails at three and direct four works')
    # In this fixed-resource kernel-lift model every slot always owns the same rho_i.
    occupied = tuple(range(5))
    released = set()
    for u,v in zip(normal,normal[1:]):
        before,after = set(occupied),set(occupied)
        need(not ((after-before)&released), 'no reacquisition of a released normalized resource')
        released |= before-after
    need(not released, 'kernel phase changes never release the five fixed distinct resources')
    return {'states':9, 'ordinary_budget3_route':normal, 'strict_record_budget3_reachable':False, 'strict_record_budget3_history_states':explored3, 'strict_record_budget4_route':strict4, 'strict_record_budget4_history_states':explored4, 'released_normalized_units':0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    out = {'scope':'finite abstract full-record controls; no native authentication claimed', 'source_archives':source_pins(), 'branch_holes':branch_controls(), 'witness_covers':cover_controls(), 'quota_relations':quota_controls(), 'binary_marginals':binary_marginal_controls(), 'paid_boundary_lifts':budget_controls(), 'directed_boundary_lifts':directed_controls(), 'temporal_nonreturn':nonreturn_controls(), 'checks':CHECKS}
    out = json.loads(json.dumps(out,sort_keys=True))
    receipt = HERE/'receipt.json'
    if args.write:
        receipt.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    else:
        need(out==json.loads(receipt.read_text()), 'entire frozen exact receipt')
    print('PASS',out['checks'],'exact conditions; abstract relational hole/cut/cost packet')


if __name__ == '__main__':
    main()
