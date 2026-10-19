#!/usr/bin/env python3
"""Exact fixed-instance controls for four-collar balance transport.
No external packages. Run: python3 poset_four_collar_20261019.py
"""
from functools import lru_cache
from fractions import Fraction
from itertools import permutations, product, combinations
import json

CORE = (0, 0, 1, 0, 8, 24, 3, 5, 189)
PORTS = (0, 6)
PAIR = (2, 4)


def validate(p):
    n = len(p)
    for j, ds in enumerate(p):
        assert not (ds >> j & 1)
        for i in range(n):
            if ds >> i & 1:
                assert p[i] & ~ds == 0
    return p


def ecount(p, pair=None):
    p = list(p)
    if pair is not None:
        x, y = pair
        p[y] |= 1 << x
    n = len(p)
    full = (1 << n) - 1
    @lru_cache(None)
    def rec(s):
        if s == full:
            return 1
        return sum(rec(s | 1 << x) for x in range(n)
                   if not (s >> x & 1) and p[x] & s == p[x])
    return rec(0)


def extensions(p):
    n = len(p)
    def rec(word, mask):
        if len(word) == n:
            yield word
        else:
            for x in range(n):
                if not (mask >> x & 1) and p[x] & mask == p[x]:
                    yield from rec(word + (x,), mask | 1 << x)
    return rec((), 0)


def insert_context(A, u, B, a, b, C, v):
    na, nb = len(A), len(B)
    p = list(A)
    p += [(x << na) | ((1 << na) - 1) for x in B]
    p += [(x << (na + nb)) | ((1 << (na + nb)) - 1) for x in C]
    if A:
        assert not any(ds >> u & 1 for ds in A)
        p[na + a] &= ~(1 << u)
    if C:
        assert C[v] == 0
        p[na + nb + v] &= ~(1 << (na + b))
    return validate(tuple(p))


def flags(p, a, b, pair):
    N = [0] * 4
    C = [0] * 4
    x, y = pair
    for w in extensions(p):
        tests = (True, w[0] == a, w[-1] == b, w[0] == a and w[-1] == b)
        event = w.index(x) < w.index(y)
        for k, ok in enumerate(tests):
            N[k] += ok
            C[k] += ok and event
    assert N[0] == ecount(p)
    assert C[0] == ecount(p, pair)
    return N, C


def degree(p):
    n = len(p)
    return [sum(i != j and not (p[i] >> j & 1 or p[j] >> i & 1)
                for j in range(n)) for i in range(n)]


def width(p):
    return max(s.bit_count() for s in range(1 << len(p))
               if all(p[i] & s == 0 for i in range(len(p)) if s >> i & 1))


def cover_edges(p):
    n = len(p)
    return [(i, j) for j in range(n) for i in range(n) if p[j] >> i & 1
            and not any((p[k] >> i & 1) and (p[j] >> k & 1) for k in range(n))]


def automorphisms(p):
    n = len(p)
    sig = [(p[i].bit_count(), sum(ds >> i & 1 for ds in p)) for i in range(n)]
    groups = [[i for i in range(n) if sig[i] == s] for s in sorted(set(sig))]
    ans = []
    for choices in product(*(permutations(g) for g in groups)):
        perm = [0] * n
        for g, c in zip(groups, choices):
            for i, j in zip(g, c):
                perm[i] = j
        if all((p[j] >> i & 1) == (p[perm[j]] >> perm[i] & 1)
               for i in range(n) for j in range(n)):
            ans.append(perm)
    return ans


def bridge_chain(m):
    p = ()
    for _ in range(m):
        if not p:
            p = CORE
        else:
            p = insert_context(p, len(p) - 9 + PORTS[1], CORE, *PORTS, (), 0)
    return p


def mv(v, M):
    return [sum(v[i] * M[i][j] for i in range(2)) for j in range(2)]


def transfer(m, marked=None):
    M = ((504, 160), (714, 230))
    E = ((264, 80), (411, 129))
    v = [1, 0]
    for i in range(m):
        v = mv(v, E if marked == i else M)
    return sum(v)


def frontier_count(p):
    """Sparse ideal flow with the proved displacement-window restriction."""
    n = len(p)
    assert all(ds < (1 << i) for i, ds in enumerate(p)), 'natural reference extension required'
    D = max(degree(p), default=0)
    levels = [{0: 1}]
    max_states = 1
    for k in range(n):
        row = {}
        for s, val in levels[-1].items():
            for x in range(max(0, k - D), min(n, k + D + 1)):
                if not (s >> x & 1) and p[x] & s == p[x]:
                    t = s | 1 << x
                    forced = (1 << max(0, k + 1 - D)) - 1
                    if t & forced != forced:
                        continue
                    row[t] = row.get(t, 0) + val
        levels.append(row)
        max_states = max(max_states, len(row))
    return levels[-1][(1 << n) - 1], D, max_states


def run():
    B = validate(CORE)
    N, C = flags(B, *PORTS, PAIR)
    assert N == [664, 280, 160, 70]
    assert C == [344, 196, 80, 49]
    denominators = [N[0], N[0]+N[1], N[0]+N[2], sum(N)]
    numerators = [C[0], C[0]+C[1], C[0]+C[2], sum(C)]
    probs = [Fraction(c, n) for c, n in zip(numerators, denominators)]
    assert min(probs) == Fraction(53, 103)
    assert max(probs) == Fraction(135, 236)
    assert all(Fraction(1, 3) <= p <= Fraction(2, 3) for p in probs)
    assert width(B) == 3 and max(degree(B)) == 7
    assert len(automorphisms(B)) == 1
    # Independent complete enumeration of the actual four collar posets.
    corner_results = []
    for left, right in product((0, 1), repeat=2):
        P = insert_context((0,) if left else (), 0, B, *PORTS,
                           (0,) if right else (), 0)
        H = ecount(P)
        K = ecount(P, tuple(x+left for x in PAIR))
        expected_H = N[0]+left*N[1]+right*N[2]+left*right*N[3]
        expected_K = C[0]+left*C[1]+right*C[2]+left*right*C[3]
        assert (H, K) == (expected_H, expected_K)
        brute_H = brute_K = 0
        for w in extensions(P):
            brute_H += 1
            brute_K += w.index(PAIR[0]+left) < w.index(PAIR[1]+left)
        assert (H, K) == (brute_H, brute_K)
        corner_results.append({'left':left,'right':right,'extensions':H,'event':K,'p':str(Fraction(K,H))})
    # Nontrivial exact contexts with different endpoint distributions.
    contexts = [(), (0,), (0, 1), (0, 0), (0, 0, 1), (0, 0, 1, 3)]
    context_checks = 0
    for A, Ctx in product(contexts, repeat=2):
        umax = [i for i in range(len(A)) if not any(ds >> i & 1 for ds in A)] or [0]
        vmin = [i for i,ds in enumerate(Ctx) if ds == 0] or [0]
        for u,v in product(umax,vmin):
            P = insert_context(A,u,B,*PORTS,Ctx,v)
            ea,ec = ecount(A),ecount(Ctx)
            ra = Fraction(sum(w[-1]==u for w in extensions(A)),ea) if A else Fraction(0)
            sc = Fraction(sum(w[0]==v for w in extensions(Ctx)),ec) if Ctx else Fraction(0)
            Z = N[0]+ra*N[1]+sc*N[2]+ra*sc*N[3]
            E = C[0]+ra*C[1]+sc*C[2]+ra*sc*C[3]
            assert ecount(P) == ea*ec*Z
            assert ecount(P,tuple(x+len(A) for x in PAIR)) == ea*ec*E
            context_checks += 1
    chain_results = []
    for m in range(1,13):
        P = bridge_chain(m)
        H,D,states = frontier_count(P)
        assert H == transfer(m) == ecount(P)
        assert D == 7
        for j in range(m):
            K=transfer(m,j)
            assert Fraction(53,103) <= Fraction(K,H) <= Fraction(135,236)
            if m <= 3:
                assert K == ecount(P,(9*j+PAIR[0],9*j+PAIR[1]))
        chain_results.append({'blocks':m,'n':len(P),'extensions':H,'frontier_states_max':states,'first_pair_p':str(Fraction(transfer(m,0),H))})
    # Minimal boundary migration control: core Q3 -> two-collar Q5.
    Q3=(0,0,1);Q5=insert_context((0,),0,Q3,0,2,(0,),0)
    assert ecount(Q5)==8
    assert ecount(Q5,(1,2))==6 and ecount(Q5,(2,3))==6
    assert ecount(Q5,(1,0))==3 and ecount(Q5,(4,3))==3
    # Failed displacement-extremum witness selector (exact 7-element example).
    bad=(0,0,0,3,1,21,55)
    words=list(extensions(bad));H=len(words)
    h=[sum(w.index(i)+1 for w in words) for i in range(7)]
    order=sorted(range(7),key=lambda i:h[i]);disp=[h[x]-(i+1)*H for i,x in enumerate(order)]
    assert H==55 and order==[0,1,2,4,5,3,6] and disp==[35,34,-18,-28,19,-27,-15]
    assert ecount(bad,(1,0))==18 and ecount(bad,(2,0))==17
    assert Fraction(18,55)<Fraction(1,3) and Fraction(17,55)<Fraction(1,3)
    assert all(bad[i]&1 for i in (3,4,5,6)) and ecount(bad,(1,2))==28
    out={'status':'exact checks passed; general conjecture remains open',
         'core':{'predecessor_masks':B,'ports':PORTS,'pair':PAIR,'flag_totals':N,'flag_event_counts':C,
                 'corners':corner_results,'width':width(B),'incomparability_degrees':degree(B),'automorphisms':len(automorphisms(B)),
                 'cover_edges':cover_edges(B)},
         'arbitrary_context_checks':context_checks,'bridge_chains':chain_results,
         'boundary_migration':{'core_predecessors':Q3,'collared_predecessors':Q5,'extensions':8,'core_minorities':['1/4','1/4'],'bridge_minorities':['3/8','3/8']},
         'displacement_selector_hostile':{'pred':bad,'extensions':55,'mean_order':order,'displacement_numerators':disp,'all_selected_incomparable_minorities':['18/55','17/55'],'balanced_pair_1_2':'28/55'}}
    print(json.dumps(out,indent=2))

if __name__=='__main__':
    run()
