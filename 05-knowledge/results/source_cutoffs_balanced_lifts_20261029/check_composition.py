"""Exact Pfaffian kernel and heterogeneous regular-composition controls."""
from functools import lru_cache
from itertools import combinations, product, permutations
from math import gcd
from functools import reduce
import random
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def matrix(n, mask):
    A = [[0]*n for _ in range(n)]
    for bit, (i,j) in enumerate(combinations(range(n), 2)):
        A[i][j] = 1 if mask >> bit & 1 else -1
        A[j][i] = -A[i][j]
    return A


def kernel(A):
    @lru_cache(None)
    def pf(ids):
        if not ids:
            return 1
        return sum((-1)**(j+1)*A[ids[0]][ids[j]]*
                   pf(ids[1:j]+ids[j+1:]) for j in range(1,len(ids)))
    n = len(A)
    return [(-1)**i*pf(tuple(j for j in range(n) if j != i)) for i in range(n)]


def prime(A):
    n = len(A)
    for i,j in combinations(range(n), 2):
        C = {i,j}
        while True:
            added = {x for x in range(n) if x not in C and
                     len({A[x][y] for y in C}) > 1}
            if not added:
                break
            C |= added
        if len(C) != n:
            return False
    return True


def normalize(v):
    if all(x < 0 for x in v):
        v = [-x for x in v]
    if not all(x > 0 for x in v):
        return None
    d = reduce(gcd,v)
    return [x//d for x in v]


def inflate(A, weights):
    vertices = [(i,x) for i,w in enumerate(weights) for x in range(w)]
    return [[(0 if x == y else
              (1 if (y-x) % weights[i] <= (weights[i]-1)//2 else -1))
             if i == j else A[i][j] for j,y in vertices] for i,x in vertices]


def automorphisms(A):
    n=len(A)
    degrees=[sum(x==1 for x in row) for row in A]
    return [pi for pi in permutations(range(n))
            if all(degrees[i]==degrees[pi[i]] for i in range(n)) and
            all(A[i][j]==A[pi[i]][pi[j]] for i,j in combinations(range(n),2))]


def main():
    rng = random.Random(43)
    tested = 0
    positive5 = 0
    witness = None
    for n, masks in [(5, range(1<<10)),
                     (7, [rng.randrange(1<<21) for _ in range(2000)]),
                     (9, [rng.randrange(1<<36) for _ in range(300)])]:
        for mask in masks:
            A = matrix(n,mask)
            v = kernel(A)
            require(all(x % 2 for x in v), 'every cofactor is odd and nonzero')
            require(all(sum(a*x for a,x in zip(row,v)) == 0 for row in A),
                    'literal integer kernel equation')
            u = normalize(v)
            if n == 5 and u:
                positive5 += 1
            # Switching is an existence probe, not a native metagraph operation.
            d = [1 if x > 0 else -1 for x in v]
            B = [[d[i]*d[j]*A[i][j] for j in range(n)] for i in range(n)]
            weights = normalize([abs(x) for x in v])
            require(all(sum(a*x for a,x in zip(row,weights)) == 0 for row in B),
                    'positive switched kernel')
            if witness is None and len(set(weights)) > 1 and prime(B):
                expanded = inflate(B,weights)
                require(all(sum(row) == 0 for row in expanded), 'literal regular inflation')
                require(not prime(expanded), 'inflation has genuine child modules')
                # Strong positive control: multiplying all child sizes by any odd t works.
                expanded3 = inflate(B,[3*w for w in weights])
                require(all(sum(row) == 0 for row in expanded3), 'common odd multiplier')
                altered = weights[:];altered[0] += 2
                require(any(sum(row) != 0 for row in inflate(B,altered)),
                        'independent child-size change destroys regularity')
                aut=automorphisms(B)
                require(all(weights[i]==weights[pi[i]] for pi in aut for i in range(n)),
                        'primitive kernel sizes retained by automorphisms')
                require(len(aut)==1 and weights.count(3)==2 and weights.count(1)==5,
                        'explicit rigid skeleton / two triangles / five singleton ports')
                # Exact action on the actual eleven labelled vertices, not a cycle-index guess.
                starts=[];at=0
                for w in weights:starts.append(at);at+=w
                groups=[]
                for turns in product(range(3),repeat=2):
                    pi=list(range(at));t=iter(turns)
                    for start,w in zip(starts,weights):
                        if w==3:
                            h=next(t)
                            for j in range(3):pi[start+j]=start+(j+h)%3
                    require(all(expanded[i][j]==expanded[pi[i]][pi[j]]
                                for i,j in combinations(range(at),2)), 'native child action')
                    groups.append(pi)
                def orbit_count(items,act):
                    unseen=set(items);answer=0
                    while unseen:
                        x=min(unseen);unseen-={act(pi,x) for pi in groups};answer+=1
                    return answer
                subsets=orbit_count(range(1<<at),lambda pi,m:sum(1<<pi[i] for i in range(at) if m>>i&1))
                points=orbit_count(range(at),lambda pi,x:pi[x])
                require((subsets,points)==(512,7),'actual subset and vertex orbit counts')
                witness = {'prime_order':n, 'sign_matrix':B, 'primitive_sizes':weights,
                           'regular_composite_order':sum(weights),
                           'regular_rows_verified':len(expanded),
                           'triple_inflation_order':len(expanded3),
                           'prime_automorphism_order':len(aut),
                           'composite_automorphism_order':len(groups),
                           'binary_subset_orbits':subsets,'vertex_orbits':points,
                           'unmarked_one_vertex_children':subsets-points}
            tested += 1
    require(witness is not None, 'nonuniform positive prime-kernel witness')
    return {'exact_kernel_cases':tested,'positive_kernel_labelled_order5':positive5,
            'nonuniform_prime_witness':witness,
            'hostile_control':'Increasing only the first child size by two destroys regularity.'}


if __name__ == '__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
