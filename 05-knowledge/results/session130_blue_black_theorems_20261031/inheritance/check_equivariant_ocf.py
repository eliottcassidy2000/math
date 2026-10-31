#!/usr/bin/env python3
"""Literal-cycle parity quotient and coloured-incidence bridge. Standard library.

No inherited producer/data imports. Assertions remain active under python -O.
The fixed-path convention agrees with the owner's explorer, but the primary
OCF verifier enumerates literal directed cycles and packs their vertex sets.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations
import json


def demand(test, text):
    if not test:
        raise RuntimeError(text)


@lru_cache(None)
def pairs(n):
    return tuple(combinations(range(n), 2))


@lru_cache(None)
def orders(n):
    return tuple(permutations(range(n)))


def canonical(A):
    return min(sum(A[p[i]][p[j]] << (len(pairs(len(A)))-1-k)
                   for k, (i, j) in enumerate(pairs(len(A))))
               for p in orders(len(A)))


def tile_coordinates(n):
    return [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1)]


def board(n, word):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    for k, (x, y) in enumerate(tile_coordinates(n)):
        i, j = n-x, n-y
        if word & (1 << k):
            A[j][i] = 1
        else:
            A[i][j] = 1
    return tuple(map(tuple, A))


def reflected(n, word):
    coords = tile_coordinates(n)
    idx = {p: k for k, p in enumerate(coords)}
    return sum(((word >> k) & 1) << idx[(n+1-y, n+1-x)]
               for k, (x, y) in enumerate(coords))


def rotate_min(c):
    i = c.index(min(c))
    return c[i:] + c[:i]


def odd_cycles(A):
    """Enumerate simple directed cycles modulo rotation, never modulo support.

    DFS fixes the least vertex as root, avoiding factorial generation of words
    that are not directed paths. A cycle and its reversed word are distinct;
    in a tournament at most one orientation of a particular cyclic word exists.
    """
    n = len(A)
    result = []
    for root in range(n):
        def visit(word, used):
            if len(word) >= 3 and len(word) % 2 and A[word[-1]][root]:
                result.append(tuple(word))
            for nxt in range(root+1, n):
                if not (used >> nxt) & 1 and A[word[-1]][nxt]:
                    visit(word+[nxt], used | (1 << nxt))
        visit([root], 1 << root)
    demand(len(result) == len(set(result)), 'Duplicate literal cycle')
    return tuple(result)


def support(c):
    return sum(1 << v for v in c)


def packing_polynomial(n, objects):
    """Objects are (nonempty vertex support, polynomial degree/weight).

    Smallest-vertex recurrence counts every disjoint collection once, retaining
    multiplicity when different literal objects share a support.
    """
    by_vertex = [[] for _ in range(n)]
    for mask, degree in objects:
        for v in range(n):
            if mask & (1 << v):
                by_vertex[v].append((mask, degree))

    @lru_cache(None)
    def recur(available):
        if not available:
            return (1,)
        bit = available & -available
        v = bit.bit_length()-1
        total = list(recur(available ^ bit))
        for mask, degree in by_vertex[v]:
            if available & mask == mask:
                value = recur(available ^ mask)
                while len(total) < len(value)+degree:
                    total.append(0)
                for j, coefficient in enumerate(value):
                    total[j+degree] += coefficient
        return tuple(total)
    ans = list(recur((1 << n)-1))
    while len(ans) > 1 and not ans[-1]:
        ans.pop()
    return ans


def length_refined_polynomial(n, objects):
    """Independent explicit collection recursion, not the scalar vertex DP.

    Each object is (vertex mask, exponent tuple); literal parallel objects
    remain separate choices. This counts actual independent collections.
    """
    dimension = (n-1)//2
    total = Counter({(0,)*dimension:1})
    def extend(first, occupied, exponents):
        for i in range(first, len(objects)):
            mask, powers = objects[i]
            if not mask & occupied:
                new = tuple(a+b for a,b in zip(exponents,powers))
                total[new] += 1
                extend(i+1, occupied | mask, new)
    extend(0, 0, (0,)*dimension)
    return total


def hp_dp(A):
    n = len(A)
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for subset in range(1 << n):
        for end, count in enumerate(dp[subset]):
            if count:
                for nxt in range(n):
                    if not subset & (1 << nxt) and A[end][nxt]:
                        dp[subset | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def verify_board(A, involution=None):
    n = len(A)
    cycles = odd_cycles(A)
    alpha = packing_polynomial(n, [(support(c), 1) for c in cycles])
    def power(c, multiplicity=1):
        return tuple(multiplicity if 2*j+3==len(c) else 0 for j in range((n-1)//2))
    refined = length_refined_polynomial(n, [(support(c),power(c)) for c in cycles])
    independent_alpha = Counter()
    for exponents,count in refined.items():
        independent_alpha[sum(exponents)] += count
    demand(all(independent_alpha[k] == count for k,count in enumerate(alpha)),
           'Independent packing collection recursion differs from vertex DP')
    H = hp_dp(A)
    demand(H == sum(a*(2**k) for k, a in enumerate(alpha)), 'OCF failed')
    result = dict(alpha=alpha, H=H,
                  counts_by_length=dict(sorted(Counter(map(len, cycles)).items())),
                  refined_profile=[dict(multiplicities=p,count=v) for p,v in sorted(refined.items())])
    if involution is None:
        return result
    rho = involution
    demand(all(rho[rho[i]] == i for i in range(n)), 'Not involutory')
    demand(all(A[rho[j]][rho[i]] == A[i][j] for i, j in pairs(n)),
           'Not an anti-automorphism')
    cs = set(cycles)
    starred = {c: rotate_min(tuple(rho[v] for v in reversed(c))) for c in cycles}
    demand(all(d in cs for d in starred.values()), 'Cycle action leaves universe')
    fixed = [c for c in cycles if starred[c] == c]
    fixed_vertices = {i for i in range(n) if rho[i] == i}
    demand(len(fixed_vertices) == n % 2, 'Anti-involution fixed vertex count')
    demand(all(fixed_vertices <= set(c) for c in fixed), 'Fixed cycle misses center')
    blocks = []
    refined_blocks = []
    for c in cycles:
        d = starred[c]
        if c < d and not (support(c) & support(d)):
            blocks.append((support(c) | support(d), 2))
            refined_blocks.append((support(c) | support(d),power(c,2)))
    folded = packing_polynomial(n, blocks+[(support(c), 1) for c in fixed])
    fixed_refined = length_refined_polynomial(n,refined_blocks+[(support(c),power(c)) for c in fixed])
    demand(all((refined[p]-fixed_refined[p]) % 2 == 0
               for p in set(refined)|set(fixed_refined)), 'Multivariate fixed-packing parity')
    for p,count in refined.items():
        odd_multiplicities = sum(v % 2 for v in p)
        if odd_multiplicities > n % 2:
            demand(count % 2 == 0, 'Forbidden multi-odd packing sector')
    padded = max(len(alpha), len(folded))
    a = alpha+[0]*(padded-len(alpha))
    b = folded+[0]*(padded-len(folded))
    demand(all((x-y) % 2 == 0 for x, y in zip(a, b)), 'Packing quotient parity failed')
    demand(len(fixed_vertices) or not fixed, 'Even order has fixed odd cycle')
    if n % 2 == 0:
        demand(all(a[k] % 2 == 0 for k in range(1, len(a), 2)),
               'Even SC odd packing coefficient')
        demand(H % 4 == 1, 'Even SC H mod4')
    demand(H % 4 == (1+2*len(fixed)) % 4, 'Fixed cycle H mod4')
    result.update(fixed_odd_cycles=len(fixed),
                  fixed_by_length=dict(sorted(Counter(map(len, fixed)).items())),
                  admissible_paired_cycle_orbits=len(blocks),
                  fixed_packing_polynomial=folded)
    return result


def metagraph(n):
    m = len(tile_coordinates(n))
    As = [board(n, w) for w in range(1 << m)]
    cc = [canonical(A) for A in As]
    opp = {c: canonical(tuple(zip(*A))) for c, A in zip(cc, As)}
    merged = [min(c, opp[c]) for c in cc]
    codes = sorted(set(merged))
    node_of = [codes.index(c) for c in merged]
    fibres = defaultdict(list)
    for w, u in enumerate(node_of):
        fibres[u].append(w)
    unmerged_sizes = Counter(cc)
    rows = []
    for u in range(len(codes)):
        words = fibres[u]
        blue = [w for w in words if reflected(n, w) == w]
        A = As[min(blue or words)]
        cyc = verify_board(A, tuple(reversed(range(n))) if blue else None)
        H = cyc['H']
        aut = H // unmerged_sizes[cc[words[0]]]
        direct_aut = sum(all(A[p[i]][p[j]] == A[i][j] for i, j in pairs(n))
                         for p in orders(n))
        demand(H % unmerged_sizes[cc[words[0]]] == 0 and aut == direct_aut,
               'Aut and HP fibre law disagree')
        b = len(blue)
        degree = (len(words)-b)//2
        a1 = cyc['alpha'][1] if len(cyc['alpha']) > 1 else 0
        charge_formula = ((a1+(1-aut*b)//2) % 2) if blue else 1
        demand(degree % 2 == charge_formula, 'Cycle/black incidence charge')
        if blue:
            demand(2*aut*degree == H-aut*b, 'SC degree bridge')
        rows.append(dict(node=u,representative_word=min(blue or words),
                         H=H,aut=aut,blue_degree=b,mass=len(words),
                         quotient_black_degree=degree,charge=degree % 2,
                         scores=sorted(map(sum,A)),**{k:v for k,v in cyc.items() if k!='H'}))
    edges = Counter()
    for w in range(1 << (m-1)):
        if reflected(n, w) != w:
            edges[tuple(sorted((node_of[w],node_of[((1 << m)-1)^w])))] += 1
    demand(all(c % 2 == 0 for c in edges.values()), 'Odd black multiplicity')
    quotient = {uv: multiplicity//2 for uv, multiplicity in edges.items()}
    graph_degrees = Counter()
    for (u,v), multiplicity in quotient.items():
        graph_degrees[u] += multiplicity
        graph_degrees[v] += multiplicity
    demand(all(graph_degrees[u] == rows[u]['quotient_black_degree'] for u in fibres),
           'Loop-aware incidence degree disagreement')
    test_subsets = range(1 << len(rows)) if n <= 5 else [0,(1 << len(rows))-1] + [1 << u for u in fibres]
    cut_count = 0
    for subset in test_subsets:
        cut = sum(v for (i,j),v in edges.items() if ((subset >> i) ^ (subset >> j)) & 1)
        rhs = 2*sum(row['charge'] for row in rows if (subset >> row['node']) & 1)
        demand((cut-rhs) % 4 == 0, 'OCF black-cut congruence')
        cut_count += 1
    return dict(n=n,ordinary_classes=len(set(cc)),merged_nodes=len(rows),rows=rows,
                black_line_count=sum(edges.values()),
                quotient_black_edges=[dict(nodes=uv,multiplicity=c) for uv,c in sorted(quotient.items())],
                cut_tests=cut_count)


def symmetric_words(n):
    coords = tile_coordinates(n)
    unseen = set(range(len(coords)))
    orbits = []
    while unseen:
        i = min(unseen)
        x,y = coords[i]
        j = coords.index((n+1-y,n+1-x))
        orb = {i,j}
        unseen -= orb
        orbits.append(sum(1 << k for k in orb))
    for bits in range(1 << len(orbits)):
        yield sum(orb for k,orb in enumerate(orbits) if bits & (1 << k))


def triangle_stack(k):
    n=3*k
    A=[[0]*n for _ in range(n)]
    for i,j in pairs(n):
        if i//3 != j//3:
            A[i][j]=1
        else:
            A[i][j]=int((j-i) % 3 == 1)
            A[j][i]=1-A[i][j]
    # Reverse block order and fix vertex0 in each triangle while swapping1,2.
    rho=tuple(3*(k-1-i//3)+({0:0,1:2,2:1}[i%3]) for i in range(n))
    return tuple(map(tuple,A)),rho


def equivariant_conflict_hostile():
    examples=[]
    rho=tuple(reversed(range(5)))
    for word in (1,26):
        A=board(5,word)
        cycles=odd_cycles(A)
        paths=[p for p in orders(5) if all(A[p[j]][p[j+1]] for j in range(4))]
        blue=[]
        for p in paths:
            rev={p[j]:p[4-j] for j in range(5)}
            if all(A[rev[j]][rev[i]]==A[i][j] for i,j in pairs(5)):
                blue.append(p)
        common=set(range(5))
        for c in cycles:
            common &= set(c)
        examples.append(dict(word=word,adjacency_rows=A,
            out_neighbor_bitmasks=[sum(x << j for j,x in enumerate(row)) for row in A],
            cycles=cycles,hamiltonian_paths=paths,blue_paths=blue,
            common_vertices_of_all_odd_cycles=sorted(common)))
    demand([len(x['hamiltonian_paths']) for x in examples]==[9,9], 'H9 hostile paths')
    demand([len(x['blue_paths']) for x in examples]==[1,3], 'H9 hostile blue counts')
    phi={(0,1,2,3,4):(0,4,1,2,3), (0,2,4):(1,2,3),
         (0,1,4):(0,2,3), (0,3,4):(1,2,4)}
    demand(set(phi)==set(examples[0]['cycles']) and
           set(phi.values())==set(examples[1]['cycles']), 'Conflict map not bijection')
    def star(c):
        return rotate_min(tuple(rho[v] for v in reversed(c)))
    for c,d in combinations(phi,2):
        demand(bool(set(c)&set(d)) == bool(set(phi[c])&set(phi[d])), 'Conflict map adjacency')
    demand(all(len(c)==len(phi[c]) and phi[star(c)]==star(phi[c]) for c in phi),
           'Conflict map length or reflection failure')
    return dict(examples=examples,cycle_bijection=[dict(source=c,target=d) for c,d in phi.items()],
                preserves=['cycle length','pairwise cycle intersection','anti-involution action'],
                loses=['actual vertex-incidence supports','blue Hamiltonian path count'])


def main():
    out={'status':'FINITE-EXACT verification of general elementary theorem',
         'metagraphs':{str(n):metagraph(n) for n in (3,4,5,6)}}
    symmetry=[]
    mixed_length_hostile=None
    for n in (3,4,5,6,7,8):
        tested=0
        # Every reflection-fixed word through n7; deterministic stride sample n8.
        for i,w in enumerate(symmetric_words(n)):
            if n==8 and i % 32:
                continue
            details=verify_board(board(n,w),tuple(reversed(range(n))))
            if n==8 and mixed_length_hostile is None:
                for item in details['refined_profile']:
                    if item['multiplicities']==(1,1,0) and item['count']:
                        mixed_length_hostile=dict(n=n,word=w,**details)
            tested+=1
        symmetry.append(dict(n=n,symmetric_words_tested=tested,
                             scope='exhaustive' if n<=7 else 'every32nd half-cube word'))
    out['reflection_fixed_tests']=symmetry
    demand(mixed_length_hostile is not None,'No nonzero (3,5) packing control')
    out['nonzero_mixed_length_even_control']=mixed_length_hostile
    out['triangle_stacks']=[]
    for k in range(1,6):
        A,rho=triangle_stack(k)
        v=verify_board(A,rho)
        demand(v['alpha']==[__import__('math').comb(k,j) for j in range(k+1)],
               'Triangle stack binomial profile')
        out['triangle_stacks'].append(dict(k=k,**v))
    # The n5 same-score pair is where triangle-only counting first fails.
    rows=out['metagraphs']['5']['rows']
    hostile=[r for r in rows if r['scores']==[1,2,2,2,3]]
    demand(sorted(r['H'] for r in hostile)==[11,13,15], 'n5 same-score branches')
    out['same_score_hostile_n5']=hostile
    out['equivariant_conflict_graph_hostile']=equivariant_conflict_hostile()
    out['checks']={'literal_directed_cycles_not_supports':True,
                   'held_karp_vs_cycle_packings':True,
                   'full_polynomial_mod2_not_only_H':True,
                   'multivariate_packing_and_forbidden_parity_sectors':True,
                   'black_loops_count_twice_in_degree':True,
                   'failure_gates_survive_python_optimization':True}
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
