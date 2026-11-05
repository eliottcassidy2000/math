#!/usr/bin/env python3
"""Exact integer Heisenberg substitution tests; standard library only."""
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent

def require(condition, message):
    if not condition:
        raise ValueError(message)

def omega(u, v):
    return u[0]*v[1]-u[1]*v[0]+u[2]*v[3]-u[3]*v[2]

def graph(points):
    rows = [0]*len(points)
    for i, j in combinations(range(len(points)), 2):
        z = omega(points[i], points[j])
        require(z != 0, 'tie')
        a, b = (i,j) if z > 0 else (j,i)
        rows[a] |= 1 << b
    return rows

def substitute_graph(outer, inners):
    labels = [(a,i) for a,S in enumerate(inners) for i in range(len(S))]
    rows = [0]*len(labels)
    for u, v in combinations(range(len(labels)), 2):
        a,i = labels[u]
        b,j = labels[v]
        forward = (inners[a][i] >> j & 1) if a == b else (outer[a] >> b & 1)
        x,y = (u,v) if forward else (v,u)
        rows[x] |= 1 << y
    return rows

def compile_substitution(outer, inners):
    require(len(outer) == len(inners) and all(inners), 'nonempty matching blocks')
    graph(outer)
    for block in inners:
        graph(block)
    require(all(v[0] == 1 for v in outer+sum(inners, [])), 'common affine chart')
    D = E = 0
    for a,b in combinations(range(len(outer)), 2):
        va,vb = outer[a],outer[b]
        for u,v in product(inners[a], inners[b]):
            B = (vb[3]-va[3])*(u[2]+v[2])+(va[2]-vb[2])*(u[3]+v[3])
            D = max(D,abs(B))
            E = max(E,abs(omega(u,v)))
    M = D+E+1
    points = []
    for (_,p,q,r), block in zip(outer, inners):
        points.extend((1,M*M*p+M*(r*Q-q*R)+P,M*q+Q,M*r+R)
                      for _,P,Q,R in block)
    require(graph(points) == substitute_graph(graph(outer),list(map(graph,inners))),
            'full labelled substitution equality')
    offsets = [sum(map(len,inners[:a])) for a in range(len(inners))]
    for a,block in enumerate(inners):
        for i,j in combinations(range(len(block)),2):
            require(omega(points[offsets[a]+i],points[offsets[a]+j]) == omega(block[i],block[j]),
                    'internal pairings unchanged exactly')
    return points, {'D':D,'E':E,'M':M}

def subset_path_counts(rows):
    n = len(rows)
    table = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        table[1<<v][v] = 1
    for mask in range(1,1<<n):
        for last,c in enumerate(table[mask]):
            targets = rows[last] & ~mask
            while c and targets:
                bit = targets & -targets
                targets -= bit
                table[mask|bit][bit.bit_length()-1] += c
    return [1]+[sum(row) for row in table[1:]]

def brute_H(rows):
    return sum(all(rows[a] >> b & 1 for a,b in zip(p,p[1:]))
               for p in permutations(range(len(rows))))

def path_cover_profile(rows):
    """Unordered components, choosing the component containing the least name."""
    n = len(rows)
    h = subset_path_counts(rows)
    pc = [[0]*(n+1) for _ in range(1<<n)]
    pc[0][0] = 1
    for mask in range(1,1<<n):
        first = mask & -mask
        rest = mask ^ first
        sub = rest
        while True:
            block = sub|first
            other = mask^block
            for d in range(1, n+1):
                pc[mask][d] += h[block]*pc[other][d-1]
            if sub == 0:
                break
            sub = (sub-1)&rest
    fact=1
    F=[0]
    for d in range(1,n+1):
        fact *= d
        F.append(fact*pc[-1][d])
    return F

def c3_diagonal_H(F):
    F=F+[0]
    return 3*sum(F[r]**3+F[r+1]*F[r]**2+F[r+1]**2*F[r]
                 for r in range(1,len(F)-1))

def transitive(rows, chosen):
    return not any(rows[a]>>b&1 and rows[b]>>c&1 and rows[c]>>a&1
                   for a,b,c in permutations(chosen,3))

def feedback_number(rows):
    for size in range(len(rows),-1,-1):
        for chosen in combinations(range(len(rows)),size):
            if transitive(rows,chosen):
                return len(rows)-size,list(chosen)

def rank(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    row = 0
    for col in range(len(a[0])):
        pivot = next((r for r in range(row,len(a)) if a[r][col]),None)
        if pivot is None:
            continue
        a[row],a[pivot] = a[pivot],a[row]
        z=a[row][col]
        a[row]=[x/z for x in a[row]]
        for r in range(row+1,len(a)):
            z=a[r][col]
            a[r]=[x-z*y for x,y in zip(a[r],a[row])]
        row += 1
    return row

def word_points(word):
    K=3*len(word)+1
    return [(1,0,K,0),(1,0,0,K)]+[
        (1,3*i,-1 if c&2 else 1,1 if c&1 else -1)
        for i,c in enumerate(word,1)]

def left_translate(v,u):
    _,p,q,r=v
    _,P,Q,R=u
    return (1,p+P+r*Q-q*R,q+Q,r+R)

def main():
    singleton=[(1,0,0,0)]
    c3=[(1,0,1,0),(1,0,0,1),(1,0,-1,-1)]
    require(graph(c3)==[2,4,1],'triangle orientation')
    t2,bound2=compile_substitution(c3,[c3,c3,c3])
    r2=graph(t2)
    h2=subset_path_counts(r2)[-1]
    require(h2==3159==brute_H(r2),'independent nine-vertex Hamilton count')
    feedback,witness=feedback_number(r2)
    require(feedback==5,'feedback vertex number exactly five')
    require(rank([[omega(u,v) for v in t2] for u in t2])==4,'exact alternating rank four')
    F2=path_cover_profile(r2)
    t3,bound3=compile_substitution(c3,[t2,t2,t2])
    h3=c3_diagonal_H(F2)
    require(h3==417382500592116023859,'27-vertex inherited target independently recalculated')
    require(h3==3**13*1309949*199849717,'factorization multiplication')
    require(all(p%d for p in (1309949,199849717) for d in range(2,__import__('math').isqrt(p)+1)),
            'two non-139-smooth factors are prime by trial division')
    require(rank([[omega(u,v) for v in t3] for u in t3])==4,'27-vertex exact rank four')

    # Exhaustive all-word geometric controls, not a census of rank-four tournaments.
    compiled=0
    path_controls=0
    for m in range(4):
        for word in product(range(4),repeat=m):
            outer=word_points(word)
            inners=[c3 if a == 0 else singleton for a in range(len(outer))]
            out,_=compile_substitution(outer,inners)
            compiled+=1
            if len(out)<=6:
                r=graph(out)
                require(subset_path_counts(r)[-1]==brute_H(r),'direct count control')
                path_controls+=1
    # Every ordered pair of the four one-letter outer/inner words.
    for u,v in product(range(4),repeat=2):
        inner=word_points((v,))
        out,_=compile_substitution(word_points((u,)),[inner,inner,inner])
        compiled+=1

    # Uniform Euclidean shrinking omits the quadratic central scale and fails.
    v=(1,0,1,0)
    naive=[(1,v[1]+P,v[2]+Q,v[3]+R) for _,P,Q,R in c3]
    naive_pairs=[omega(naive[i],naive[j]) for i,j in combinations(range(3),2)]
    require(graph(naive) != graph(c3) and transitive(graph(naive),range(3)),
            'naive translation hostile changes a triangle to a transitive tournament')
    shifted=[left_translate(v,u) for u in c3]
    require(graph(shifted)==graph(c3),'Heisenberg compensator repairs translation')

    # Same H, different first path-cover jet: transitive orders one and two.
    chain2=[(1,0,0,0),(1,1,0,0)]
    one_wrap,_=compile_substitution(c3,[singleton,singleton,singleton])
    two_wrap,_=compile_substitution(c3,[chain2,singleton,singleton])
    require(subset_path_counts(graph(one_wrap))[-1]==3,'H=1 one-vertex wrapper')
    require(subset_path_counts(graph(two_wrap))[-1]==5,'H=1 two-vertex wrapper differs')

    # Exact H(C3[S,1,1]) = 3H(S)+F_S(2), tested on all short word geometries.
    wrap_controls=0
    for m in range(3):
        for word in product(range(4),repeat=m):
            inner=word_points(word)
            out,_=compile_substitution(c3,[inner,singleton,singleton])
            F=path_cover_profile(graph(inner))
            require(subset_path_counts(graph(out))[-1]==3*F[1]+F[2], 'wrapper law')
            wrap_controls+=1
    receipt={
        'status':'PROVED geometry; FINITE-EXACT controls; no completeness claim',
        'triangle_squared':{'points':t2,'rows':r2,'bound':bound2,'H':h2,
                            'H_methods':['subset DP','literal permutations'],
                            'feedback_vertex_number':feedback,'max_transitive_subset':witness,
                            'ordered_path_cover_profile':F2[1:], 'alternating_rank':4},
        'triangle_cubed':{'points':t3,'rows':graph(t3),'bound':bound3,'H':h3,
                         'H_method':'proved quotient-walk law with independently counted nine-vertex path-cover profile',
                         'not_full_27_vertex_DP':True,'factorization':{'3':13,'1309949':1,'199849717':1},
                         'alternating_rank':4,'feedback_vertex_number_by_proved_recursion':19},
        'controls':{'geometric_substitutions':compiled,'small_direct_counts':path_controls,
                    'c3_wrapper_law':wrap_controls,'naive_translation_pairings':naive_pairs,
                    'same_H_wrapper_outputs':[3,5]},
        'source_hashes':{}
    }
    for path in ['session192-harmonics/PROOF.md','session191-grammar/END_COMMUTATION.md',
                 'session191-grammar/SPARSITY.md',
                 'math/01-canon/theorems/THM-3243-contact-deformation-blowup-equivariance-and-full-orbit-resolution.md',
                 'math/01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md',
                 'math/01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md']:
        receipt['source_hashes'][path]=hashlib.sha256((ROOT.parent/path).read_bytes()).hexdigest()
    output=json.dumps(receipt,sort_keys=True,indent=2)+'\n'
    (ROOT/'receipt.json').write_text(output)
    print(output,end='')

if __name__=='__main__':
    main()
