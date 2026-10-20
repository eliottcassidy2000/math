"""Independent audit: no producer-code imports; parsed literal certificates.

Run from any directory. Uses only Python's standard library. Checks do not
depend on assertions and remain enabled under python -O.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
import ast
import hashlib
import json
import re
import sys

HERE = Path(__file__).resolve().parent
FROZEN = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE / "frozen"


def check(ok, message):
    if not ok:
        raise ValueError(message)


def matrix(a, b):
    n = a + b + 1
    return tuple(tuple(False if u == v else
                       (u > v if {u, v} in ({0, a}, {a, a+b}) else u < v)
                       for v in range(n)) for u in range(n))


def legal(A, p):
    return set(p) == set(range(len(A))) and len(p) == len(A) and all(
        A[u][v] for u, v in zip(p, p[1:]))


def all_paths(A):
    # List/set recursion independent from producer bitmask DFS.
    result = []
    def extend(prefix, available):
        if not available:
            result.append(tuple(prefix))
        else:
            for v in available:
                if not prefix or A[prefix[-1]][v]:
                    extend(prefix + [v], available - {v})
    extend([], set(range(len(A))))
    return result


def positional_blue(A, p):
    # Test the relabelled matrix directly, with positions as indices.
    n = len(p)
    return all(A[p[i]][p[j]] == A[p[n-1-j]][p[n-1-i]]
               for i in range(n) for j in range(n))


def coordinates(p, a, b):
    middle = p.index(a)
    left = frozenset(v for v in p[:middle] if v < a)
    right = frozenset(v for v in p[middle+1:] if v > a)
    I = None if left == frozenset(range(a)) else left
    Kraw = None if right == frozenset(range(a+1, a+b+1)) else right
    if a == b and Kraw is not None:
        Kraw = frozenset(2*a-v for v in Kraw)
    return I, Kraw


def literal_decoder(p, a, b):
    # Reconstruct from subsets extracted from the actual path, not its index.
    cut = p.index(a)
    I = set(v for v in p[:cut] if v < a)
    J = set(v for v in p[:cut] if v > a)
    check(I == set(range(a)) or 0 not in I, 'left decoder domain')
    check(not J or a+b in J, 'right decoder domain')
    recovered = tuple(sorted(I) + sorted(J) + [a] +
                      sorted(set(range(a))-I) +
                      sorted(set(range(a+1, a+b+1))-J))
    check(recovered == p, 'literal decoder inverse')


def rectangular():
    maps = deletions = 0
    details = []
    for a, b in product(range(2, 6), repeat=2):
        A = matrix(a, b)
        source = all_paths(A)
        check(len(source) == (2**(a-1)+1)*(2**(b-1)+1), 'rectangle count')
        for p in source:
            literal_decoder(p, a, b)
        if a+b+1 <= 7:
            brute = {p for p in permutations(range(len(A))) if legal(A, p)}
            check(set(source) == brute, 'permutation/DFS independent check')
        for r in range(a-1):
            for s in range(b-1):
                for dl in combinations(range(1, a), r):
                    for dr in combinations(range(a+1, a+b), s):
                        gone = set(dl+dr)
                        keep = [v for v in range(len(A)) if v not in gone]
                        reindex = {v:i for i,v in enumerate(keep)}
                        targetA = matrix(a-r, b-s)
                        check(tuple(tuple(A[u][v] for v in keep) for u in keep) == targetA,
                              'induced graph equality')
                        buckets = Counter()
                        for p in source:
                            q = tuple(reindex[v] for v in p if v not in gone)
                            check(legal(targetA, q), 'actual deletion HP')
                            buckets[q] += 1
                        target = all_paths(targetA)
                        check(set(buckets) == set(target), 'deletion surjectivity')
                        for q in target:
                            factor = 2**((0 if q[0] == 0 else r) +
                                         (0 if q[-1] == a+b-r-s else s))
                            check(buckets[q] == factor, 'all individual deletion fibres')
                        maps += 1
                        deletions += len(source)
        details.append([a, b, len(source)])
    check((maps, deletions) == (676, 112896), 'reported finite universe')
    return dict(maps=maps, source_path_deletions=deletions, cases=details)


def predicted_lifts(k, s, t, blue=False):
    # Direct transcription of equations (4)-(5); concrete paths are independent.
    emb = dict((v, v + (0 if v < k else 1 if v == k else 2)) for v in range(2*k+1))
    ports = Counter()
    for x in ([None] if s == 0 else [0, 1]):
        for y in ([None] if t == 2*k else [0, 1]):
            if blue and x != y:
                continue
            start, end = emb[s], emb[t]
            if s >= k:
                start = k if x == 1 else k+2 if x == 0 and y == 0 else start
            if t <= k:
                end = k+2 if y == 1 else k if y == 0 and x == 0 else end
            ports[start, end] += 1
    return ports


def endpoints():
    report = []
    ambiguous = None
    for k in range(2, 7):
        A, AA = matrix(k, k), matrix(k+1, k+1)
        old, new = all_paths(A), all_paths(AA)
        grouped = defaultdict(list)
        for p in new:
            q = tuple(v if v < k else v-1 if v == k+1 else v-2
                      for v in p if v not in {k, k+2})
            grouped[q].append(p)
        check(set(grouped) == set(old), 'endpoint-map surjection')
        for q in old:
            check(Counter((p[0], p[-1]) for p in grouped[q]) ==
                  predicted_lifts(k, q[0], q[-1]), 'endpoint kernel per literal parent')
            if positional_blue(A, q):
                check(Counter((p[0], p[-1]) for p in grouped[q] if positional_blue(AA, p)) ==
                      predicted_lifts(k, q[0], q[-1], blue=True), 'blue lift kernel')
            else:
                check(not any(positional_blue(AA, p) for p in grouped[q]),
                      'black parents have only black lifts')
        for p in old + new:
            kk = (len(p)-1)//2
            x, y = coordinates(p, kk, kk)
            check(positional_blue(matrix(kk,kk), p) == (x == y), 'actual blue diagonal')
        by_ports = defaultdict(dict)
        for p in old:
            by_ports[p[0], p[-1]][positional_blue(A,p)] = p
        if ambiguous is None:
            for port, witnesses in by_ports.items():
                if len(witnesses) == 2:
                    ambiguous = dict(k=k, endpoints=port, blue=witnesses[True], black=witnesses[False])
                    break
        report.append(dict(k=k, parents=len(old), children=len(new)))
    check(ambiguous is not None, 'endpoint alone is insufficient hostile')
    return dict(per_parent=report, endpoint_blue_ambiguity=ambiguous)


def expected_port_matrix(k):
    n = 2*k+1
    A = [[0]*n for _ in range(n)]
    for i, j in product(range(k), repeat=2):
        aa = 1 if i == 0 else 2**(k-1-i)
        bb = 1 if j == 0 else 2**(k-1-j)
        A[i][2*k-j] = aa*bb
    for alpha in range(1,k):
        for beta in range(k):
            c = 1 if beta == alpha-1 else 2**(beta-alpha-1) if beta > alpha else 0
            A[alpha][beta] = c
            A[2*k-beta][2*k-alpha] = c
    A[0][k] = A[k][2*k] = A[k+1][k-1] = 1
    return A


def determinant(A):
    # Rational Gaussian elimination, not producer Bareiss.
    Q = [[Fraction(x) for x in row] for row in A]
    result = Fraction(1)
    for j in range(len(Q)):
        pivot = next((i for i in range(j,len(Q)) if Q[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            Q[j], Q[pivot] = Q[pivot], Q[j]
            result *= -1
        v = Q[j][j]
        result *= v
        Q[j] = [x/v for x in Q[j]]
        for i in range(j+1,len(Q)):
            if Q[i][j]:
                w = Q[i][j]
                Q[i] = [x-w*y for x,y in zip(Q[i],Q[j])]
    check(result.denominator == 1, 'integer-matrix determinant is integral')
    return int(result)


def unique_matching(A):
    # Independently certify uniqueness by repeatedly finding any forced column.
    rows, cols = set(range(len(A))), set(range(len(A)))
    matching = {}
    while cols:
        possibilities = [(j, [i for i in rows if A[i][j]]) for j in cols]
        check(all(ii for _,ii in possibilities), 'matching not singular')
        forced = next(((j,ii[0]) for j,ii in possibilities if len(ii)==1), None)
        check(forced is not None, 'support peeling always has a forced column')
        j,i = forced
        check(A[i][j] == 1, 'forced weight equals one')
        matching[i] = j
        rows.remove(i)
        cols.remove(j)
    inversions = sum(matching[i] > matching[j] for i in matching for j in matching if i<j)
    check(inversions % 2 == 0, 'matching sign positive')
    return [matching[i] for i in range(len(A))]


def full_ports():
    out = []
    for k in range(2, 13):
        A = expected_port_matrix(k)
        if k <= 7:
            literal = Counter((p[0],p[-1]) for p in all_paths(matrix(k,k)))
            check(A == [[literal[i,j] for j in range(2*k+1)] for i in range(2*k+1)],
                  'full matrix coefficients against actual paths')
        # Check c(alpha,beta) separately against actual subsets.
        for alpha in range(1,k):
            for beta in range(k):
                count = sum(bool(bits) and min(bits)==alpha and
                            max(set(range(k))-set(bits))==beta
                            for r in range(k) for bits in combinations(range(1,k),r))
                check(A[alpha][beta] == count, 'coefficient subset enumeration')
        check(sum(x!=0 for row in A for x in row) == 2*k*k-k+3, 'support count')
        match = unique_matching(A)
        check(determinant(A)==1, 'independent exact determinant')
        out.append(dict(k=k, support=2*k*k-k+3, unique_matching=match, determinant=1, permanent=1))
    return out


def literal_text(text=None):
    if text is None:
        text = (FROZEN / 'BALANCED_DELETION_EXTENSIONS.md').read_text()
    section = text.split('## 5. ')[1].split('## 6. ')[0]
    labels = re.search(r'\(u,v,a,b,c\)=\(([^)]+)\)', section)
    check(labels is not None, 'literal label map present')
    label = dict(zip('uvabc', ast.literal_eval('('+labels.group(1)+')')))
    row_match = re.search(r'outgoing bit rows are\s+\(([^)]+)\)',section)
    check(row_match is not None, 'literal bit rows present')
    rows = ast.literal_eval('('+row_match.group(1)+')')
    A = tuple(tuple(bool(rows[u] & (1<<v)) for v in range(5)) for u in range(5))
    check(all(not A[i][i] for i in range(5)) and
          all(A[i][j] != A[j][i] for i,j in combinations(range(5),2)), 'literal rows form tournament')
    source = re.search(r'The Hamiltonian path\s+\(([^)]+)\)',section)
    target = re.search(r'deletes to \(([^)]+)\)',section)
    check(source is not None and target is not None, 'literal paths are present')
    p = tuple(label[x.strip()] for x in source.group(1).split(','))
    q = tuple(label[x.strip()] for x in target.group(1).split(','))
    check(legal(A,p), 'printed source witness is a Hamiltonian path')
    check(q == tuple(v for v in p if v != label['v']), 'printed remainder is actual deletion')
    bad_arcs = [(u,v) for u,v in zip(q,q[1:]) if not A[u][v]]
    check(bad_arcs == [(label['a'],label['b'])], 'printed remainder has stated bad joining arc')
    u,v = label['u'],label['v']
    check(A[u][v] and all(A[u][w]==A[v][w] for w in range(5) if w not in {u,v}),
          'literal numerical rows certify homogeneous twins')
    # Regression: a wrong source path must be rejected even if generic routines pass.
    corrupted = (p[1],p[0])+p[2:]
    check(not legal(A,corrupted), 'negative literal-path parser control')
    return dict(rows=rows, source=p, remainder=q, bad_arcs=bad_arcs,
                intentionally_corrupted_source=corrupted, corrupted_source_rejected=True)


def parsed_text_regression():
    original = (FROZEN / 'BALANCED_DELETION_EXTENSIONS.md').read_text()
    bad = original.replace('    (a,v,b,c,u)', '    (v,a,b,c,u)')
    check(bad != original, 'negative test actually mutates printed example')
    try:
        literal_text(bad)
    except ValueError as error:
        check(str(error) == 'printed source witness is a Hamiltonian path',
              'negative parser regression rejected for correct reason')
        return dict(mutated_printed_path='(v,a,b,c,u)',
                    outcome='REJECTED: printed source witness is not a Hamiltonian path')
    raise ValueError('parser failed to reject a deliberately corrupted printed witness')


def minimal_twins():
    controls = 0
    for n in range(2,5):
        edges = list(combinations(range(n),2))
        for orientations in product((0,1),repeat=len(edges)):
            A = [[False]*n for _ in range(n)]
            for bit,(u,v) in zip(orientations,edges):
                A[u][v],A[v][u] = bool(bit),not bit
            P = [p for p in permutations(range(n)) if legal(A,p)]
            for u,v in edges:
                if all(A[u][w]==A[v][w] for w in range(n) if w not in {u,v}):
                    for removed in (u,v):
                        for p in P:
                            q = tuple(x for x in p if x!=removed)
                            check(all(A[x][y] for x,y in zip(q,q[1:])), 'no smaller twin orphan')
                            controls += 1
    check(controls==412, 'small twin universe count')
    return controls


def join_and_marking():
    A = matrix(3,3)
    P = all_paths(A)
    n = len(A)
    cycle = (0,1,2,4,5,6,3)
    check(set(cycle)==set(range(n)) and all(A[cycle[i]][cycle[(i+1)%n]] for i in range(n)),
          'literal directed Hamiltonian cycle proves source strong')
    autos = [g for g in permutations(range(n))
             if all(A[i][j] == A[g[i]][g[j]] for i,j in combinations(range(n),2))]
    antis = [g for g in permutations(range(n))
             if all(A[i][j] == A[g[j]][g[i]] for i,j in combinations(range(n),2))]
    check(autos == [tuple(range(n))] and antis == [tuple(reversed(range(n)))],
          'literal source rigidity and unique anti-automorphism')
    AA = tuple(tuple(A[u%n][v%n] if u//n==v//n else u<n
                     for v in range(2*n)) for u in range(2*n))
    Q = all_paths(AA)
    expected = {p+tuple(v+n for v in q) for p,q in product(P,repeat=2)}
    check(set(Q)==expected and len(Q)==625, 'literal join Cartesian product')
    out = Counter()
    for p in Q:
        bulk = p[:n] != tuple(range(n)) and p[n:] != tuple(range(n,2*n))
        colour = 'blue' if positional_blue(AA,p) else 'black'
        stratum = 'bulk' if bulk else 'boundary'
        out[colour]+=1
        out[stratum]+=1
        out[stratum+'_'+colour]+=1
    check(out==Counter(blue=25,black=600,bulk=576,boundary=49,bulk_blue=24,
                       bulk_black=552,boundary_blue=1,boundary_black=48), 'all eight join counts')
    # The supplied JSON is read and verified, never regenerated as a substitute.
    marked = json.loads((FROZEN/'marked_f5_paths.json').read_text())
    check(len(marked)==25 and {tuple(x['path']) for x in marked}==set(P), 'marked literal path completeness')
    label = {None:0,frozenset():1,frozenset({1}):2,frozenset({2}):3,frozenset({1,2}):4}
    mapping = {}
    for item in marked:
        p = tuple(item['path'])
        check(legal(A,p), 'every printed marked path legal')
        I,K = coordinates(p,3,3)
        vector = (label[I],label[K])
        check(vector==tuple(item['vector']), 'marked vector agrees with literal path inverse')
        states = [4 if x is None else sum(2**(v-1) for v in x) for x in (I,K)]
        check(states==item['states'], 'marked state masks agree with inverse')
        check(item['blue']==positional_blue(A,p), 'marked colour is actual predicate')
        mapping[p]=vector
    check(set(mapping.values())==set(product(range(5),repeat=2)), 'marked F5 bijection')
    for p,v in mapping.items():
        r = tuple(6-x for x in reversed(p))
        check(mapping[r]==v[::-1], 'actual anti-reversal becomes swap')
        I,K = coordinates(p,3,3)
        kappas = tuple(None if x is None else frozenset({1,2})-x for x in (I,K))
        J = next(q for q in P if coordinates(q,3,3)==kappas)
        check(mapping[J]==tuple(-x%5 for x in v), 'J becomes sign negation')
    check(mapping[tuple(range(7))]==(0,0), 'standard increasing path is zero')
    kappa = (0,4,3,2,1)  # indexed by F5 labels, after this marking
    marking_count = sum(all(f[kappa[x]] == -f[x]%5 for x in range(5))
                        for f in permutations(range(5)))
    check(marking_count==8, 'all equivariant one-coordinate markings')
    sign_orbits = {min(v,tuple(-x%5 for x in v)) for v in mapping.values() if v!=(0,0)}
    check(len(sign_orbits)==12, 'nonzero sign quotient')
    return dict(join=dict(out), marked_paths=25, equivariant_markings=marking_count,
                nonzero_sign_orbits=len(sign_orbits), ordered_sign_pairs=len(sign_orbits)**2,
                source_hamiltonian_cycle=cycle, source_automorphisms=len(autos),
                source_anti_automorphisms=len(antis))


def algebra():
    cases = 0
    for q in range(1,101):
        for D in range(1,17):
            qp = D*(q-1)+1
            total = D*D*q*q-2*D*(D-1)*q+(D-1)**2
            black = D*D*q*(q-1)-D*(D-1)*(q-1)
            strata = ((q-1)*(q-2)*D*D + 2*(q-1)*D + (q-1)*D*(D-1))
            check(total==qp*qp and black==qp*(qp-1)==strata, 'all-depth black strata identity')
            cases+=1
    # General nonuniform fibre law, including singleton/unequal weights.
    fibres = 0
    for length in range(1,6):
        for w in product(range(1,5),repeat=length):
            lhs = sum(w)*(sum(w)-1)
            rhs = sum(x*y for i,x in enumerate(w) for j,y in enumerate(w) if i!=j)
            rhs += sum(x*(x-1) for x in w)
            check(lhs==rhs, 'general product-fibre identity')
            fibres+=1
    return dict(all_depth_integer_controls=cases, nonuniform_fibre_controls=fibres)


def main():
    result = dict(rectangular=rectangular(), endpoints=endpoints(), full_ports=full_ports(),
                  literal_text=literal_text(), small_twin_controls=minimal_twins(),
                  join_and_marking=join_and_marking(), algebra=algebra(),
                  parsed_text_regression=parsed_text_regression())
    result['audited_document_sha256'] = hashlib.sha256((FROZEN/'BALANCED_DELETION_EXTENSIONS.md').read_bytes()).hexdigest()
    result['audit_source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
