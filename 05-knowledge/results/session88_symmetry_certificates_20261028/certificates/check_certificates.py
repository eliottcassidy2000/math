#!/usr/bin/env python3
"""Session88 proof-carrying decision controls; stdlib, no assert statements."""
from collections import Counter, deque
from itertools import combinations_with_replacement, permutations, product
from math import isqrt, prod
import json
import random


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def expect_reject(fn, label):
    try:
        fn()
    except (ValueError, RuntimeError, IndexError, KeyError):
        return label
    raise RuntimeError('hostile accepted: '+label)


def maximum_spanning_tree(weights):
    """Dense Prim, O(k^2) comparisons, O(k) extra words, including zero edges."""
    k = len(weights)
    if k <= 1:
        return 0, []
    used = [False] * k
    best = [-1] * k
    parent = [-1] * k
    best[0] = 0
    total = 0
    edges = []
    for _ in range(k):
        v = max((i for i in range(k) if not used[i]), key=lambda i: best[i])
        used[v] = True
        total += best[v]
        if parent[v] >= 0:
            edges.append((parent[v], v))
        for u in range(k):
            if not used[u] and weights[v][u] > best[u]:
                best[u], parent[u] = weights[v][u], v
    return total, edges


def partial_count(gram, ambient, required, forbidden):
    """Count any disjoint required/forbidden membership pattern; interval promise."""
    m = len(gram)
    need(not required & forbidden, "overlapping pattern guards")
    need(0 <= required < (1 << m) and 0 <= forbidden < (1 << m), "pattern range")
    yes = [i for i in range(m) if required >> i & 1]
    no = [i for i in range(m) if forbidden >> i & 1]
    length = min((gram[i][j] for i in yes for j in yes), default=ambient)
    clipped = [min(length, *(gram[a][i] for a in yes)) if yes else gram[i][i]
               for i in no]
    weights = [[min(clipped[a], clipped[b], gram[i][j])
                for b, j in enumerate(no)] for a, i in enumerate(no)]
    weight, tree = maximum_spanning_tree(weights)
    return length - sum(clipped) + weight, [(no[i], no[j]) for i, j in tree]


def pattern_count(gram, ambient, pattern):
    return partial_count(gram, ambient, pattern, ((1 << len(gram))-1) ^ pattern)


def reconstruct_patterns(gram, ambient):
    """Positive-prefix trie: O(s*m) oracle calls, O(s*m^3) arithmetic, s<=2m+1."""
    m = len(gram)
    stack = [(0, 0, 0, ambient)] if ambient else []
    output = {}
    queries = nodes = 0
    while stack:
        depth, required, forbidden, count = stack.pop()
        nodes += 1
        if depth == m:
            output[required] = count
            continue
        bit = 1 << depth
        positive, _ = partial_count(gram, ambient, required | bit, forbidden)
        queries += 1
        negative = count-positive
        need(0 <= positive <= count, 'partial pattern inconsistent with parent')
        if positive:
            stack.append((depth+1,required|bit,forbidden,positive))
        if negative:
            stack.append((depth+1,required,forbidden|bit,negative))
    return output, {'oracle_queries':queries, 'positive_trie_nodes':nodes}


def family_control(family, ambient):
    gram = [[len(a & b) for b in family] for a in family]
    direct = Counter(sum(1 << i for i, interval in enumerate(family) if x in interval)
                     for x in range(1, ambient + 1))
    for a in range(1 << len(family)):
        count, _ = pattern_count(gram, ambient, a)
        need(count == direct[a], ('pattern count', ambient, family, a, count, direct[a]))
    reconstructed, stats = reconstruct_patterns(gram, ambient)
    need(reconstructed == dict(direct), 'entire multiset via positive trie')
    need(len(reconstructed) <= min(ambient, 2*len(family)+1), 'interval pattern support bound')
    need(stats['oracle_queries'] <= len(family)*len(reconstructed), 'output-sensitive query bound')
    return 1 << len(family)


def partial_controls():
    # Every ternary specification (unconstrained/required/forbidden) on every
    # interval family of <=3 rows on lines of length <=4.
    families = queries = 0
    for ambient in range(5):
        intervals = [frozenset()] + [frozenset(range(a,b+1))
                     for a in range(1,ambient+1) for b in range(a,ambient+1)]
        for m in range(4):
            for family in combinations_with_replacement(intervals,m):
                gram = [[len(a&b) for b in family] for a in family]
                for ternary in range(3**m):
                    value, required, forbidden = ternary, 0, 0
                    for i in range(m):
                        digit, value = value%3, value//3
                        if digit == 1:required |= 1<<i
                        if digit == 2:forbidden |= 1<<i
                    actual = 0
                    for x in range(1,ambient+1):
                        pattern = sum(1<<i for i,a in enumerate(family) if x in a)
                        actual += (pattern & required == required and not pattern & forbidden)
                    computed,_ = partial_count(gram,ambient,required,forbidden)
                    need(computed == actual, 'partial pattern count')
                    queries += 1
                families += 1
    # A large ambient line does not force enumeration of its positions.
    n = 10**12
    family = [(1,n//2), (n//3,n), (n//4,3*n//4)]
    gram = [[max(0,min(b,d)-max(a,c)+1) for c,d in family] for a,b in family]
    output, stats = reconstruct_patterns(gram,n)
    need(sum(output.values())==n and len(output)<=7, 'huge ambient symbolic recovery')
    return {'families':families, 'partial_queries':queries,
            'large_ambient':n, 'large_ambient_pattern_count':len(output),
            'large_ambient_oracle_queries':stats['oracle_queries']}


def gram_realization_controls():
    """All 4-tuples of arbitrary subsets of a 4-point universe, not just intervals."""
    n, m = 4, 4
    intervals = [0] + [sum(1<<i for i in range(a,b+1))
                      for a in range(n) for b in range(a,n)]
    def gram(f):
        return tuple((a&b).bit_count() for a in f for b in f)
    def columns(f):
        return tuple(sorted(Counter(sum(1<<i for i,a in enumerate(f) if a>>x&1)
                                    for x in range(n)).items()))
    bank = {}
    for family in product(intervals, repeat=m):
        key, counts = gram(family), columns(family)
        need(key not in bank or bank[key]==counts, 'interval Gram fibre')
        bank[key] = counts
    matched = 0
    for family in product(range(1<<n), repeat=m):
        key = gram(family)
        if key in bank:
            matched += 1
            need(columns(family)==bank[key], 'arbitrary realization of interval Gram')
    return {'ambient':n, 'rows':m, 'arbitrary_families':(1<<n)**m,
            'interval_realizable_gram_matrices':len(bank),
            'arbitrary_families_matching_interval_gram':matched}


def interval_controls():
    families = patterns = 0
    counts = []
    for ambient in range(0, 6):
        intervals = [frozenset()] + [frozenset(range(a, b + 1))
                    for a in range(1, ambient + 1) for b in range(a, ambient + 1)]
        nf = np = 0
        # Multisets suffice: every row ordering is a permutation of one checked family.
        for m in range(0, 5):
            for family in combinations_with_replacement(intervals, m):
                np += family_control(family, ambient)
                nf += 1
        families += nf
        patterns += np
        counts.append({'ambient': ambient, 'families': nf, 'patterns': np})
    rng = random.Random(880031)
    for _ in range(100):
        ambient, m = 30, 10
        family = []
        for _ in range(m):
            a, b = sorted(rng.sample(range(1, ambient+1), 2))
            family.append(frozenset(range(a,b+1)))
        patterns += family_control(family, ambient)
        families += 1
    # Three circular arcs/arbitrary sets: all pair overlaps one, no common point.
    bad = [frozenset((1,2)), frozenset((2,3)), frozenset((1,3))]
    g = [[len(a & b) for b in bad] for a in bad]
    need(pattern_count(g,3,7)[0] == 1, 'three-arc false positive predicted')
    need(not set.intersection(*(set(a) for a in bad)), 'three-arc actual triple empty')
    need(pattern_count(g,3,0)[0] == -1, 'three-arc union formula fails')
    # Reflecting a family preserves Gram and all counts, but loses cut address.
    family = [frozenset((1,2)), frozenset((2,3,4))]
    reflected = [frozenset(6-x for x in a) for a in family]
    g1 = [[len(a&b) for b in family] for a in family]
    g2 = [[len(a&b) for b in reflected] for a in reflected]
    need(g1 == g2 and family != reflected, 'Gram does not recover positions')
    return {'families': families, 'patterns': patterns, 'all_multisets_reconstructed_by_trie':True,
            'partial_patterns':partial_controls(),
            'interval_realizable_gram_fibres':gram_realization_controls(),
            'exhaustive_by_ambient': counts,
            'random_seed': 880031, 'random_families': 100,
            'hostiles': ['three circular arcs / arbitrary sets', 'reflection loses addresses']}


def orbit_certificate(seed, generators, action):
    states, tree, index = [seed], [None], {seed:0}
    for state in states:
        i = index[state]
        for g, gen in enumerate(generators):
            child = action(state, gen)
            if child not in index:
                index[child] = len(states)
                states.append(child)
                tree.append((i,g))
    closure = [[index[action(state,g)] for g in generators] for state in states]
    return states, tree, closure


def verify_orbit(seed, generators, action, cert):
    states, tree, closure = cert
    need(states and states[0] == seed, 'wrong orbit seed')
    need(len(set(states)) == len(states), 'repeated literal state')
    need(len(tree) == len(states) == len(closure), 'orbit dimensions')
    need(tree[0] is None, 'root tree record')
    for i in range(1, len(states)):
        p, g = tree[i]
        need(0 <= p < i and 0 <= g < len(generators), 'invalid tree index')
        need(action(states[p], generators[g]) == states[i], 'invalid transporter edge')
    for i, state in enumerate(states):
        need(len(closure[i]) == len(generators), 'missing closure generator')
        for g, gen in enumerate(generators):
            j = closure[i][g]
            need(0 <= j < len(states), 'closure outside certified set')
            need(action(state,gen) == states[j], 'false closure equality')
    return len(states)


def orbit_controls():
    gens = [(1,2,0,3), (0,1,3,2)]
    action = lambda x,g: sum(((x>>i)&1) << g[i] for i in range(4))
    cert = orbit_certificate(3, gens, action)
    need(verify_orbit(3,gens,action,cert) == 6, 'S4 2-set orbit')
    bad = (cert[0][:-1], cert[1][:-1], cert[2][:-1])
    h1 = expect_reject(lambda:verify_orbit(3,gens,action,bad), 'incomplete orbit')
    bad = ([*cert[0], 0], [*cert[1], (0,0)], [*cert[2], [6,6]])
    h2 = expect_reject(lambda:verify_orbit(3,gens,action,bad), 'closed extra component lacks transport')
    # Quotient equality can be exact in the quotient yet wrong for source orbit size.
    quotient_action = lambda x,g: x
    coarse = orbit_certificate(2, gens, quotient_action)
    need(verify_orbit(2,gens,quotient_action,coarse) == 1, 'coarse orbit valid only in quotient')
    return {'literal_orbit_size':6, 'unfaithful_weight_quotient_orbit_size':1,
            'rejected':[h1,h2]}


def tournament(n, code):
    rows=[0]*n
    bit=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>bit&1: rows[i] |= 1<<j
            else: rows[j] |= 1<<i
            bit += 1
    return tuple(rows)


def encode(rows, order):
    code = bit = 0
    for i,u in enumerate(order):
        for v in order[i+1:]:
            code |= ((rows[u]>>v)&1)<<bit
            bit += 1
    return code


def canonical_pointed(rows, root):
    others = [v for v in range(len(rows)) if v != root]
    best = None
    best_order = None
    for order in permutations(others):
        key = (encode(rows,order), sum(((rows[root]>>v)&1)<<i for i,v in enumerate(order)))
        if best is None or key < best:
            best, best_order = key, order
    # image permutation source vertex -> canonical target vertex
    image=[0]*len(rows)
    for i,v in enumerate(best_order): image[v]=i
    image[root]=len(rows)-1
    return best, tuple(image)


def extend(parent, mask):
    n=len(parent)
    return tuple(r | ((1-((mask>>i)&1))<<n) for i,r in enumerate(parent))+(mask,)


def verify_row(rows, bank, witnesses):
    n=len(rows)
    need(len(witnesses)==n, 'missing root record')
    need({w[0] for w in witnesses}==set(range(n)), 'duplicate or missing root')
    result=Counter()
    for root,key,image in witnesses:
        need(key in bank, 'unknown named bank state')
        target=bank[key]
        need(len(image)==n and sorted(image)==list(range(n)), 'not a bijection')
        need(image[root]==n-1, 'wrong transported root')
        for u in range(n):
            for v in range(n):
                need(((rows[u]>>v)&1)==((target[image[u]]>>image[v])&1), 'arc or mask alignment failure')
        result[key]+=1
    return result


def inverse(p):
    out=[0]*len(p)
    for i,j in enumerate(p): out[j]=i
    return tuple(out)


def compose(p,q):
    return tuple(p[q[i]] for i in range(len(p)))


def group_closure(generators,n):
    seed=tuple(range(n))
    return set(orbit_certificate(seed,generators,lambda x,g:compose(g,x))[0])


def row_controls():
    counts=[]
    for n in range(1,6):
        seen=set()
        rows_checked=0
        for code in range(1<<(n*(n-1)//2)):
            rows=tournament(n,code)
            canonical=min(encode(rows,p) for p in permutations(range(n)))
            if canonical in seen: continue
            seen.add(canonical)
            records=[]
            bank={}
            for v in range(n):
                key,image=canonical_pointed(rows,v)
                bank[key]=extend(tournament(n-1,key[0]),key[1])
                records.append((v,key,image))
            row=verify_row(rows,bank,records)
            r=0
            key,f=records[r][1:]
            # Exact stabilizer independently enumerated; this represents known parent-bank data.
            full={p for p in permutations(range(n)) if encode(rows,p)==encode(rows,range(n))}
            stabilizer={p for p in full if p[r]==r}
            gen=list(stabilizer)
            gen += [compose(inverse(fv),f) for v,k,fv in records if k==key]
            need(group_closure(gen,n)==full, 'transversal/stabilizer reconstructs full child group')
            need(row[key]*len(stabilizer)==len(full), 'row multiplicity times parent mask stabilizer')
            rows_checked += 1
            if n>1:
                expect_reject(lambda:verify_row(rows,bank,records[:-1]), 'missing root')
                duplicate = [records[0]] + records[:-1]
                expect_reject(lambda:verify_row(rows,bank,duplicate), 'duplicate roots')
                badbank = dict(bank)
                k = records[0][1]
                badbank[k] = extend(tournament(n-1,k[0]), k[1]^1)
                expect_reject(lambda:verify_row(rows,badbank,records), 'unaligned mask')
            if n>2:
                root, k, image = records[0]
                swapped = tuple(1 if x==0 else 0 if x==1 else x for x in image)
                corrupted = [(root,k,swapped)] + records[1:]
                expect_reject(lambda:verify_row(rows,bank,corrupted), 'false deletion permutation')
        counts.append({'n':n,'ordinary_classes':rows_checked})
    return {'rows':counts, 'full_group_reconstruction':True, 'incomplete_rows_rejected':True,
            'additional_hostiles_rejected':['duplicate roots','unaligned mask','false deletion permutation'],
            'scope':'all ordinary tournament classes through order 5; direct exhaustive S_n audit'}


def arithmetic_controls():
    primes=[2,5,7,19,29,59,89]
    row_primes=[2,3,5]
    columns=[]
    for p in primes:
        value=p+1
        column=0
        for j,q in enumerate(row_primes):
            while value%q==0:
                value//=q
                column ^= 1<<j
        need(value==1,'complete row factor bank')
        columns.append(column)
    for x in range(1<<len(primes)):
        value=prod(p+1 for i,p in enumerate(primes) if x>>i&1)
        syndrome=0
        for i,c in enumerate(columns):
            if x>>i&1:syndrome ^= c
        need((syndrome==0)==(isqrt(value)**2==value),'square-class certificate')
    path=[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9]
    need(sorted(path)==list(range(1,16)),'square-sum path permutation')
    need(all(isqrt(a+b)**2==a+b for a,b in zip(path,path[1:])),'square-sum path local guard')
    labels=[0,1,2]
    oriented=[labels[0]-labels[1],labels[2]-labels[1]]
    need(len(set(oriented))==2 and len({abs(x) for x in oriented})==1,'oriented graceful hostile')
    return {'square_class_subsets':128,'square_sum_path_order':15,
            'graceful_hostile_oriented':oriented,'graceful_hostile_absolute':[abs(x) for x in oriented]}


def main():
    print(json.dumps({'intervals':interval_controls(),'orbits':orbit_controls(),
                      'root_rows':row_controls(),'arithmetic_guards':arithmetic_controls()},indent=2))

if __name__=='__main__':main()
