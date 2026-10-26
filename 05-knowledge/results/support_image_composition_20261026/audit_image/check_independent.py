#!/usr/bin/env python3
"""Independent edge-set/permutation oracle; imports no producer module."""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations, product
import json

GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def make(n, digits, radix):
    es = []
    for (u, v), bit in zip(combinations(range(n), 2), digits):
        if radix == 2:
            es.append((u, v) if bit else (v, u))
        elif bit:
            es.append((u, v) if bit == 1 else (v, u))
    return frozenset(es)


def row_edges(rows):
    return frozenset((u, v) for u, row in enumerate(rows) for v in range(len(rows)) if row >> v & 1)


@lru_cache(None)
def words(n, es):
    return frozenset(p for p in permutations(range(n)) if all((u, v) in es for u, v in zip(p, p[1:])))


def support(n, es):
    return frozenset(edge for p in words(n, es) for edge in zip(p, p[1:]))


def reach(n, es, x):
    done, todo = {x}, [x]
    while todo:
        v = todo.pop()
        for w in range(n):
            if (v, w) in es and w not in done:
                done.add(w)
                todo.append(w)
    return done


def decode(n, es):
    rs = [reach(n, es, x) for x in range(n)]
    components = set(frozenset(y for y in range(n) if y in rs[x] and x in rs[y]) for x in range(n))
    valid_orders = [cs for cs in permutations(components)
                    if all(next(iter(cs[j])) in rs[next(iter(cs[i]))]
                           for i in range(len(cs)) for j in range(i+1, len(cs)))]
    if len(valid_orders) != 1:
        return None
    cs = valid_orders[0]
    out = set(es)
    for i, c in enumerate(cs):
        for j in range(i+1, len(cs)):
            for u in c:
                for v in cs[j]:
                    if ((u, v) in es) != (j == i+1):
                        return None
                    out.add((u, v))
        for u, v in combinations(c, 2):
            if (u, v) in es or (v, u) in es:
                continue
            uv = sum((u, w) in es and (w, v) in es for w in c)
            vu = sum((v, w) in es and (w, u) in es for w in c)
            if (uv >= 2) == (vu >= 2):
                return None
            out.add((u, v) if uv >= 2 else (v, u))
    return frozenset(out)


def defect(n, T, D):
    return Counter(sum(edge not in D for edge in zip(p, p[1:])) for p in words(n, T))


def cover_counts(n, T, D):
    if n == 0:
        return Counter({0:1})
    missing = T-D
    out = Counter()
    def visit(remaining, last, blocks):
        if not remaining:
            out[blocks-1] += 1
            return
        for length in range(1, len(remaining)+1):
            for path in permutations(remaining, length):
                if not all(edge in D for edge in zip(path, path[1:])):
                    continue
                if last is not None and (last, path[0]) not in missing:
                    continue
                visit(remaining-set(path), path[-1], blocks+1)
    visit(set(range(n)), None, 0)
    return out


def relabel(es, perm):
    return frozenset((perm[u], perm[v]) for u, v in es)


def canonical(n, es):
    return min(tuple(sorted(relabel(es, p))) for p in permutations(range(n)))


def partner(es, path):
    kept = set(zip(path, path[1:]))
    return frozenset((u, v) if (u, v) in kept else (v, u) for u, v in es)


def blue(es, path):
    reversal = dict(zip(path, reversed(path)))
    return relabel(es, reversal) == frozenset((v, u) for u, v in es)


def agreement(n, T, S):
    common = T&S
    candidates = [p for p in permutations(range(n)) if frozenset(zip(p, p[1:])) == common]
    need(len(candidates) <= 1, 'at most one agreement path')
    return candidates[0] if candidates else None


def padded(n, es):
    out = {(u+1, v+1) for u, v in es}
    out.update((0, v) for v in range(1, n+2))
    out.update((v, n+1) for v in range(1, n+1))
    return frozenset(out)


def hub(n, h, s, t):
    out = {(u+1, v+1) for u, v in h}
    out.update((0, s+1) if v == s else (v+1, 0) for v in range(n))
    out.update((t+1, n+1) if v == t else (n+1, v+1) for v in range(n))
    out.add((n+1, 0))
    return frozenset(out)


def main():
    census = []
    all_tournaments = {}
    for n in range(5):
        count = n*(n-1)//2
        Ts = [make(n, bits, 2) for bits in product((0, 1), repeat=count)]
        all_tournaments[n] = Ts
        image = {support(n, T):T for T in Ts}
        need(len(image) == len(Ts), 'support injectivity small oracle')
        accepted, decoded = 0, 0
        for digits in product((0, 1, 2), repeat=count):
            D = make(n, digits, 3)
            T = decode(n, D)
            passes = False
            if T is not None:
                decoded += 1
                need(D <= T, 'decoder extension')
                poly = defect(n, T, D)
                need(poly[0] == len(words(n, D)), 'zero coefficient')
                need(cover_counts(n, T, D) == poly, 'independently enumerated ordered covers')
                passes = support(n, D) == D and sum(poly.values()) == poly[0]
            need(passes == (D in image), 'exact image iff')
            if D in image:need(T == image[D], 'decoder returns original tournament')
            accepted += passes
        census.append({'n':n,'partial_graphs':3**count,'decoded':decoded,'image':accepted})
    D = row_edges((10,20,16,6,1));T = row_edges((14,20,16,22,1))
    need(decode(5,D) == T and support(5,D) == D, 'five-vertex false local candidate')
    need(defect(5,T,D) == Counter({0:7,1:2}) and cover_counts(5,T,D) == Counter({0:7,1:2}), 'false image exact polynomial')
    need(T-D == {(0,2),(3,4)}, 'false-image missing arcs')
    D4=row_edges((2,4,8,1));T4=row_edges((6,12,8,1))
    need(defect(4,T4,D4) == Counter({0:4,2:1}) and cover_counts(4,T4,D4) == Counter({0:4,2:1}), 'joint connector interaction')
    need(words(4,T4)-words(4,D4) == {(1,3,0,2)}, 'literal joint word')
    for edge in ((0,2),(1,3)):
        need(words(4,D4|{edge}) == words(4,D4), 'individual connector inert')
    transitive = frozenset(combinations(range(3),2))
    Ds=support(3,transitive)
    need(Ds == {(0,1),(1,2)} and (0,2) not in Ds, 'nonhereditary support')
    need(all(support(2,T) for T in all_tournaments[2]), 'edgeless pair not in image')
    need((0,2) in support(4,T4), 'insertion can reveal old arc')
    line_tests, hub_tests = 0, 0
    for n in (3,4):
        reps = [frozenset(t) for t in sorted({canonical(n,T) for T in all_tournaments[n]})]
        aut = {S:sum(relabel(S,p)==S for p in permutations(range(n))) for S in reps}
        counts={}
        for T in reps:
            expected=Counter()
            joint=defaultdict(list)
            for P in words(n,T):
                S=partner(T,P)
                need(agreement(n,T,S) == P, 'line agreement recovers word')
                need(blue(T,P)==blue(S,P), 'line colour unchanged')
                expected[canonical(n,S),blue(T,P)] += 1
                paddedP=(0,)+tuple(v+1 for v in P)+(n+1,)
                H=partner(padded(n,S),paddedP)
                need(H == hub(n,T,P[0],P[-1]), 'explicit fixed-port hub')
                need(blue(padded(n,S),paddedP) == blue(T,P), 'padding colour')
                need(all(len(reach(n+2,H,x))==n+2 for x in range(n+2)), 'hub strong by literal reachability')
                joint[P[0],P[-1]].append(H)
                hub_tests+=1
            for hs in joint.values():need(len(set(hs))==1,'same ports same hub')
            for S in reps:
                observed=Counter()
                for g in permutations(range(n)):
                    target=relabel(S,g)
                    P=agreement(n,T,target)
                    if P is not None:
                        need(partner(T,P)==target,'agreement sufficiency')
                        observed[blue(T,P)] += 1
                for colour in (False,True):need(observed[colour]==aut[S]*expected[canonical(n,S),colour],'exact transporter division')
                counts[T,S]=expected[canonical(n,S),False]
                line_tests+=1
        for T in reps:
            for S in reps:need(counts[T,S]*aut[S]==counts[S,T]*aut[T],'normalized black incidence symmetric')
    print(json.dumps({'status':'FINITE-EXACT independent PASS','gates':GATES,
                      'all_partial_graphs_n0_to4':census,'false_five_polynomial':[7,2],
                      'joint_four_polynomial':[4,0,1], 'ordinary_type_pairs':line_tests,
                      'padded_hub_paths':hub_tests,
                      'scope':'Independent proof audit I1-I5; finite oracle n<=4 plus explicit five-vertex controls.'},sort_keys=True,indent=2))


if __name__=='__main__':
    main()
