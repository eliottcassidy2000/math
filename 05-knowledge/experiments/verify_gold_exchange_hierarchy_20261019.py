#!/usr/bin/env python3
"""Exact certificates for one fixed Gold exchange hierarchy, not a poset census.

All guards remain active under -O. Standard library only. The receipt stores
matchings and equal-size vertex covers; word IDs refer to lexicographically
sorted legal words of the specified second orientation in the fixed child.
"""
from collections import Counter, deque
from functools import lru_cache
from itertools import combinations, permutations
import json

checks = 0

def require(ok, msg):
    global checks
    if not ok:
        raise RuntimeError(msg)
    checks += 1


def close(pred):
    p = list(pred)
    for k in range(len(p)):
        for v in range(len(p)):
            if p[v] >> k & 1:
                p[v] |= p[k]
    return None if any(p[v] >> v & 1 for v in range(len(p))) else tuple(p)


def add(pred, a, b):
    p = list(pred)
    p[b] |= 1 << a
    return close(p)


@lru_cache(None)
def count(pred):
    if pred is None:
        return 0
    full = (1 << len(pred))-1
    @lru_cache(None)
    def visit(mask):
        if mask == full:
            return 1
        return sum(visit(mask | 1 << v) for v in range(len(pred))
                   if not mask >> v & 1 and not pred[v] & ~mask)
    return visit(0)


def words(pred):
    out = []
    def visit(mask, w):
        if len(w) == len(pred):
            out.append(tuple(w))
        for v in range(len(pred)):
            if not mask >> v & 1 and not pred[v] & ~mask:
                visit(mask | 1 << v, w + [v])
    visit(0, [])
    return out


def upsets(pred):
    return tuple(sum(1 << v for v in range(len(pred)) if pred[v] >> u & 1)
                 for u in range(len(pred)))


def good_candidates(pred):
    U = upsets(pred)
    E = count(pred)
    out = []
    for dual, (D, V) in enumerate(((pred, U), (U, pred))):
        for a, b in permutations(range(len(pred)), 2):
            if (D[a] | V[a]) >> b & 1 or D[a] & ~D[b]:
                continue
            diff = [v for v in range(len(pred)) if (V[b] & ~V[a]) >> v & 1]
            if not all((D[x] >> y | D[y] >> x) & 1 for x, y in combinations(diff, 2)):
                continue
            N = count(add(pred, b, a) if dual else add(pred, a, b))
            out.append((dual, a, b, N))
    return out


def perm(n, cycle):
    sigma = list(range(n))
    for a, b in zip(cycle, cycle[1:] + cycle[:1]):
        sigma[a] = b
    return tuple(sigma)


def matching_cover(adj, nr):
    """Deterministic augmenting paths plus an independently checked cover."""
    nl = len(adj)
    left, right = [-1]*nl, [-1]*nr
    # One augmenting search per left vertex suffices; previous vertices remain
    # matched when a later augmenting path is taken.
    for start in range(nl):
        queue = deque([start])
        seen_left = {start}
        parent_right = {}
        terminal = None
        while queue and terminal is None:
            a = queue.popleft()
            for b in adj[a]:
                if b in parent_right:
                    continue
                parent_right[b] = a
                if right[b] == -1:
                    terminal = b
                    break
                nxt = right[b]
                if nxt not in seen_left:
                    seen_left.add(nxt)
                    queue.append(nxt)
        if terminal is not None:
            b = terminal
            while b != -1:
                a = parent_right[b]
                old = left[a]
                left[a], right[b] = b, a
                b = old
    pairs = [(a, b) for a, b in enumerate(left) if b != -1]
    require(len({b for a, b in pairs}) == len(pairs), 'matching right capacities')
    require(all(b in adj[a] for a, b in pairs), 'matching edge membership')
    # Reachable alternating vertices from unmatched left vertices.
    ZL = {a for a in range(nl) if left[a] == -1}
    ZR = set()
    queue = deque(sorted(ZL))
    while queue:
        a = queue.popleft()
        for b in adj[a]:
            if left[a] == b or b in ZR:
                continue
            ZR.add(b)
            nxt = right[b]
            require(nxt != -1, 'no remaining augmenting path')
            if nxt not in ZL:
                ZL.add(nxt)
                queue.append(nxt)
    CL, CR = set(range(nl))-ZL, ZR
    require(all(a in CL or b in CR for a in range(nl) for b in adj[a]),
            'vertex cover covers every edge')
    require(len(CL)+len(CR) == len(pairs), 'equal-size matching and cover')
    return pairs, sorted(CL), sorted(CR)


def exchange_graph(A, B, cycles, n):
    ids = {w: j for j, w in enumerate(B)}
    edge_cycle = [dict() for w in A]
    nonempty = []
    for cycle in cycles:
        sigma = perm(n, cycle)
        size = 0
        for i, w in enumerate(A):
            z = tuple(sigma[v] for v in w)
            if z in ids:
                j = ids[z]
                edge_cycle[i].setdefault(j, cycle)
                size += 1
        if size:
            nonempty.append({'cycle': cycle, 'edges': size})
    adj = [tuple(sorted(row)) for row in edge_cycle]
    pairs, CL, CR = matching_cover(adj, len(B))
    # Check the stored certificate by direct words and support, independently
    # of the candidate cycle list and the augmenting search.
    certified = []
    for i, j in pairs:
        w, z = A[i], B[j]
        sigma = [0]*n
        for x, y in zip(w, z):
            sigma[x] = y
        support = [v for v in range(n) if sigma[v] != v]
        require(len(support) in (2, 3), 'matched support at most three')
        require(all(sigma[v] in support for v in support), 'support invariant')
        require(w.index(3) < w.index(7) and z.index(3) > z.index(7),
                'matched orientation flips')
        require(tuple(perm(n, edge_cycle[i][j])) == tuple(sigma),
                'stored cycle agrees with independent relative permutation')
        certified.append([i, j, list(edge_cycle[i][j])])
    return {'candidate_cycles': len(cycles), 'nonempty_cycles': nonempty,
            'edge_count': sum(map(len, adj)), 'matching_size': len(pairs),
            'matching': certified, 'vertex_cover_left': CL, 'vertex_cover_right': CR}


def is_module(pred, M):
    U = upsets(pred)
    for z in range(len(pred)):
        if z in M:
            continue
        types = {(bool(pred[v] >> z & 1), bool(U[v] >> z & 1)) for v in M}
        if len(types) != 1:
            return False
    return True


def module_controls():
    # Module {0,1,2} is a two-point chain and an isolated vertex; an outside
    # isolated vertex is interleaved with it, and a global top is added.
    P = (0, 1, 0, 0, 15)
    M = {0, 1, 2}
    require(close(P) == P and is_module(P, M), 'selected nontrivial module')
    W = words(P)
    restrictions = Counter(tuple(v for v in w if v in M) for w in W)
    require(set(restrictions.values()) == {4}, 'uniform module restriction fibers')
    require(len(restrictions) == 3 and len(W) == 12, 'module extension counts')
    # A more asymmetric module also has equal fibers under outside shuffles.
    Q = (0, 0, 1, 0, 0, 31)
    M2 = {0, 1, 2}
    require(is_module(Q, M2), 'second module')
    C = Counter(tuple(v for v in w if v in M2) for w in words(Q))
    require(len(C) == 3 and len(set(C.values())) == 1, 'second uniform fibers')
    return {'first_module_pred': P, 'fibers': [(k, v) for k, v in sorted(restrictions.items())],
            'second_module_pred': Q, 'second_fiber_size': next(iter(C.values()))}


def main():
    Q = (0, 0, 0, 6, 0, 19, 2, 18, 86, 151)
    P = add(Q, 6, 7)
    require(P == (0, 0, 0, 6, 0, 19, 2, 82, 86, 215), 'Qprime masks')
    U = upsets(P)
    degrees = [9-(D|V).bit_count() for D,V in zip(P,U)]
    require(max(degrees) == 7, 'seven-thin')
    eligible = good_candidates(P)
    E = count(P)
    require(E == 7249 and len(eligible) == 42, 'total and eligibility')
    require(min(r[3] for r in eligible) == 3804, 'minimum eligible numerator')
    require(all(2*r[3] > E for r in eligible), 'no good pair in P or dual')
    # Each seed vertex has an incomparable degree-seven neighbor, so a
    # nontrivial chain substitution at any vertex violates seven-thinness.
    saturation = []
    for seed in (Q, P):
        V = upsets(seed)
        d = [9-(D|T).bit_count() for D,T in zip(seed,V)]
        blockers = []
        for v in range(10):
            bs = [u for u in range(10) if u != v and d[u] == 7
                  and not (seed[v]|V[v]) >> u & 1]
            require(bool(bs), 'chain substitution blocked at every seed vertex')
            blockers.append(bs)
        saturation.append(blockers)
    W = words(P)
    require(len(W) == E, 'full enumeration agrees with ideal DP')
    first, second = (2, 6), (3, 7)
    R = add(P, *first)
    require(R == (0, 0, 0, 6, 0, 19, 6, 86, 86, 215), 'child masks')
    RW = words(R)
    b = count(R)
    require(b == 4710 == len(RW), 'heavy count')
    require(3*E < 5*b and 3*b <= 2*E, 'hard probability window')
    A = sorted(w for w in RW if w.index(3) < w.index(7))
    B = sorted(w for w in RW if w.index(3) > w.index(7))
    require((len(A),len(B)) == (2422,2288), 'heavy second cells')
    require(count(add(R,3,7)) == len(A) and count(add(R,7,3)) == len(B),
            'independent cell ideal counts')
    cells = [sum((w.index(2)<w.index(6)) == a and (w.index(3)<w.index(7)) == z
                 for w in W) for a,z in [(True,True),(True,False),(False,True),(False,False)]]
    margins = [E-b-cells[0],E-b-cells[1],b-cells[2],b-cells[3]]
    require(min(margins) > 0, 'actual Gold margins')
    tolerance, need = 2*E-3*b, 2*b-E
    require((tolerance,need) == (368,2171), 'exchange thresholds')
    position = Counter()
    for w in RW:
        i,j = w.index(3),w.index(7)
        position[min(i,j),max(i,j),i<j] += 1
    rows = []
    for i,j in sorted({(i,j) for i,j,o in position}):
        a,z = position[i,j,True],position[i,j,False]
        rows.append([i,j,a,z,a-z])
    dpos = sum(abs(row[4]) for row in rows)
    cap = sum(min(row[2],row[3]) for row in rows)
    require((dpos,cap) == (2034,1338), 'universal pair-swapping positional cap')
    require(cap < need and dpos > tolerance, 'all pair-swapping involution families fail')
    ranks = Counter(w.index(3) for w in RW)
    require((ranks[2],ranks[3],ranks[9]) == (208,476,674), 'compressed rank obstruction')
    require(all(4 <= w.index(7) <= 8 for w in RW), 'four predecessors and successor of seven')
    require(ranks[2]+ranks[3]+ranks[9] == 1358 > tolerance, 'compressed unmatched mass suffices')
    swaps = [c for c in combinations(range(10),2) if set(c)&{3,7}]
    both = [(3,7)] + [cyc for z in range(10) if z not in(3,7)
                      for cyc in [(3,7,z),(3,z,7)]]
    triples = [cyc for c in combinations(range(10),3) if set(c)&{3,7}
               for cyc in [c,(c[0],c[2],c[1])]]
    graphs = {}
    for name,cycles,expected_edges,expected_match in [
            ('all_transpositions',swaps,3210,1360),
            ('pair_swap_and_triples_containing_both',both,6952,2086),
            ('all_triples_meeting_pair',triples,12870,2206)]:
        result = exchange_graph(A,B,cycles,10)
        require(result['edge_count'] == expected_edges, name+': edge count')
        require(result['matching_size'] == expected_match, name+': certified optimum')
        graphs[name] = result
    nu = graphs['all_triples_meeting_pair']['matching_size']
    require(nu >= need, 'three-cycle matching certifies Gold')
    require(b-2*nu == 298 <= tolerance, 'three-cycle defect')
    require(E-2*b+nu == 35, 'certified minimum heavy margin')
    module = module_controls()
    out = {'status':'PASS', 'scope':'One prescribed child of Qprime; no poset census',
           'P_pred':P, 'degrees':degrees, 'E':E, 'eligible_good_orientations':eligible,
           'first':first, 'second':second, 'R_pred':R, 'b':b,
           'cells_11_10_01_00':cells, 'actual_gold_margins':margins,
           'tolerance':tolerance, 'required_matching':need,
           'position_rows_zero_based_i_j_Aij_Bij_difference':rows,
           'positional_defect':dpos, 'positional_matching_cap':cap,
           'rank3_histogram':dict(sorted(ranks.items())),
           'compressed_unmatched_mass':1358,
           'matching_certificate_encoding':'Left/right IDs are zero-based in lexicographically sorted child extension words with 3<7 / 7<3. Each record is [leftID,rightID,directed label cycle].',
           'graphs':graphs, 'three_cycle_certified_defect':b-2*nu,
           'three_cycle_certified_heavy_margin':E-2*b+nu,
           'chain_substitution_blockers_Q_Qprime':saturation, 'module_controls':module,
           'checks':checks}
    print(json.dumps(out,sort_keys=True,indent=2))

if __name__ == '__main__':
    main()
