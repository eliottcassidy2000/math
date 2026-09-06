#!/usr/bin/env python3
"""Native parity and Peck-source hostile; complete small/exact type controls."""
from collections import Counter, deque
from itertools import product, permutations
from math import factorial, gcd
import hashlib
import json
import sympy as sp
import eulerian_uniform_orbit_sampler_overnight_hexagon_sep05 as fixed
import even_graph_triangle_quotient_spectrum_thm4078 as primal


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def parity_trivial(parts):
    ones = sum(a==1 for a in parts)
    if any(a>1 and a%2 for a in parts):
        return False
    return ones==0 or (ones<=2 and all(a==1 or a%4==0 for a in parts))


def main():
    gates = 0
    for n in range(1,17):
        for parts in fixed.partitions(n):
            basis,_,_ = fixed.fixed_kernel(n,parts)
            literal = all(v.bit_count()%2==0 for v in basis)
            require(literal==parity_trivial(parts), ('literal parity kernel',parts,literal))
            gates += 1
    require(not parity_trivial((2,1)), 'fixed triangle hostile')
    gates += 1
    rows = []
    for n in range(1,31):
        total = twisted = 0
        for parts in fixed.partitions(n):
            b = sum(a//2 for a in parts)+sum(gcd(a,c) for i,a in enumerate(parts) for c in parts[i+1:])
            f = b-len(parts)+int(any(a%2 for a in parts))
            w = fixed.class_size(n,parts)*2**f
            total += w
            if parity_trivial(parts):
                twisted += w
        require(total%factorial(n)==twisted%factorial(n)==0, 'integer Burnside projections')
        q, difference = total//factorial(n), twisted//factorial(n)
        require((q+difference)%2==0 and difference>=0, 'parity class counts')
        if n>=3:
            require((difference==0)==(n%4==3), 'exact vanishing locus')
        rows.append((n,q,(q+difference)//2,(q-difference)//2,difference))
        gates += 3
    print('TWISTED BURNSIDE n,q,q_even,q_odd,imbalance')
    for row in rows[:16]:
        print(row)
    for n in range(3,7):
        edges,ei = primal.edge_system(n)
        reps,orbit,sizes,_ = primal.quotient_data(n)
        literal = Counter(rep.bit_count()%2 for rep in reps)
        require((literal[0],literal[1])==rows[n-1][2:4], 'literal orbit parity count')
        triangles = primal.simple_cycles(n,3,ei)
        M = primal.weighted_operator(reps,orbit,triangles)
        B = primal.boolean_support(M)
        require(all(not B[i][j] or (reps[i].bit_count()-reps[j].bit_count())%2 for i in range(len(reps)) for j in range(len(reps))), 'native bipartition')
        nullity = len(reps)-sp.Matrix(B).rank()
        require(nullity>=rows[n-1][-1], 'native adjacency index lower bound')
        print('LITERAL BOOLEAN n',n,'NULLITY',nullity,'IMBALANCE',rows[n-1][-1])
        gates += 3
        if n==5:
            distance = {0:0}
            queue = deque([0])
            while queue:
                i = queue.popleft()
                for j in range(len(reps)):
                    if B[i][j] and j not in distance:
                        distance[j]=distance[i]+1
                        queue.append(j)
            ranks = tuple(Counter(distance.values())[j] for j in range(max(distance.values())+1))
            require(ranks==(1,1,2,2,1), 'empty-rooted rank hostile')
            print('N5 INTRINSIC EMPTY-DISTANCE RANKS',ranks,'NOT RANK SYMMETRIC')
            # This alternative is a positive small quotient control, not
            # an inherited unitary-Peck source or an all-order construction.
            rank = {i:1 if rep.bit_count()%2 else 0 if rep.bit_count()<5 else 2 for i,rep in enumerate(reps)}
            bottom = [i for i in rank if rank[i]==0]
            middle = [i for i in rank if rank[i]==1]
            top = [i for i in rank if rank[i]==2]
            up0 = sp.Matrix([[B[j][i] for i in bottom] for j in middle])
            up1 = sp.Matrix([[B[j][i] for i in middle] for j in top])
            require((len(bottom),len(middle),len(top))==(2,3,2) and abs((up1*up0).det())==1, 'alternative quotient Peck positive control')
            print('N5 ALTERNATIVE QUOTIENT RANKS',(2,3,2),'ABS DET U2',abs((up1*up0).det()))
            gates += 2
        if n==4:
            states = primal.cycle_space(n,ei)
            even = [x for x in states if x.bit_count()%2==0]
            odd = [x for x in states if x.bit_count()%2]
            incidence = sp.Matrix([[int((a^b) in triangles) for a in even] for b in odd])
            require(incidence==sp.ones(4), 'labelled triangle graph is K4,4')
            gradings = [r for r in product(range(3),repeat=3) if min(r)==0 and abs(r[0]-r[1])==abs(r[1]-r[2])==1]
            require(set(gradings)=={(0,1,0),(1,0,1),(0,1,2),(2,1,0)}, 'complete invariant grading universe')
            for r in gradings:
                profile = tuple(sum(size for size,rank in zip((1,4,3),r) if rank==j) for j in range(max(r)+1))
                require(profile!=profile[::-1] or incidence.rank()==1<4, 'unitary Peck source obstruction')
            print('N4 LABELLED K4,4; FOUR S4-INVARIANT GRADINGS ALL FAIL UNITARY PECK')
            gates += 6
    print('ALL FIXED-SPACE PARITY CONTROLS',914,'PURE TYPE COUNT RANGE n1..30')
    for k in range(1,4):
        arcs = [(i,j) for i in range(k) for j in range(k) if i!=j]
        index = {arc:i for i,arc in enumerate(arcs)}
        maps = [[index[p[i],p[j]] for i,j in arcs] for p in permutations(range(k))]
        representatives = set()
        for state in range(1 << (2*len(arcs))):
            images = []
            for mapping in maps:
                image = sum(((state >> (2*i)) & 3) << (2*j) for i,j in enumerate(mapping))
                images.append(image)
            representatives.add(min(images))
        count = len(representatives)
        require(count==rows[4*k][-1], 'literal ordered digraph-pair equinumerosity')
        print('ORDERED LOOPLESS DIGRAPH PAIRS k',k,'ORBITS',count,'= q-even-minus-odd at n',4*k+1)
        gates += 1
    print('GATES',gates,'SEMANTIC SHA256',hashlib.sha256(json.dumps(rows).encode()).hexdigest())
    print('PROVED parity index / exact type criterion; CITED Stanley theorem has missing source hypothesis at n4')
    print('OPEN exact Boolean adjacency nullity and Laplacian gap; no general Peck grading claimed')


if __name__=='__main__':
    main()
