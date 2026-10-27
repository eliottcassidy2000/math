#!/usr/bin/env python3
"""Standalone exact symmetry-collision certificates (Python standard library).

Run normally or with python3 -O; every mathematical check stays active.
Vertices are F3^4 in lexicographic order, indexed 0,...,80. Orient g -> gs
for s in the explicit connection set. Matrix multiplication constructs the
nonabelian group independently. No input files, packages, or network required.

Root-pinned discrete refinement and all-pair module closures certify full
regular automorphism groups and primality. Index-three subgroup kernels and
coset-union orbit enumeration certify the binary stabilizer discrepancy.
Inflation counts use the intrinsic prime-quotient blocks and rigid children,
not brute-force automorphism enumeration on 108, 135, or 324 vertices.

Scope: these are fixed-quotient families, not complete prime censuses.
"""
from itertools import product, permutations
from collections import Counter
from fractions import Fraction
import json

CONNECTION_SETS = {
    'E81': (
        1, 6, 7, 8, 9, 10, 11, 15, 16, 17,
        24, 25, 26, 27, 36, 40, 41, 43, 44, 48,
        49, 51, 52, 55, 56, 57, 58, 59, 60, 61,
        62, 63, 64, 65, 67, 70, 73, 74, 75, 78,
    ),
    'Heisenberg27_times_C3': (
        1, 5, 6, 8, 9, 10, 13, 14, 15, 16,
        17, 19, 24, 28, 29, 30, 31, 32, 34, 35,
        42, 43, 46, 48, 49, 51, 53, 54, 57, 65,
        67, 69, 70, 72, 73, 74, 75, 76, 77, 79,
    ),
}

def need(ok, message):
    if not ok: raise RuntimeError(message)

def matrix(a,b,c): return (1,a,c,0,1,b,0,0,1)
def mm(A,B): return tuple(sum(A[3*i+k]*B[3*k+j] for k in range(3))%3 for i in range(3) for j in range(3))
def cycles(p):
    seen=set(); out=[]
    for i in range(len(p)):
        if i in seen: continue
        j=i; n=0
        while j not in seen: seen.add(j); n+=1; j=p[j]
        out.append(n)
    return tuple(sorted(out))

def audit(name, connection_set):
    nonabelian=name!='E81'
    elements=list(product(range(3),repeat=4)); idx={e:i for i,e in enumerate(elements)}
    matrices=[matrix(*e[:3]) for e in elements]
    table=[]
    for i,x in enumerate(elements):
        row=[]
        for j,y in enumerate(elements):
            if nonabelian:
                z=mm(matrices[i],matrices[j]); out=(z[1],z[5],z[2],(x[3]+y[3])%3)
            else: out=tuple((u+v)%3 for u,v in zip(x,y))
            row.append(idx[out])
        table.append(row)
    need(all(table[0][i]==table[i][0]==i for i in range(81)), 'identity')
    need(all(table[table[i][j]][k]==table[i][table[j][k]] for i,j,k in product(range(81),repeat=3)), 'associativity')
    inverse=[next(j for j in range(81) if table[i][j]==table[j][i]==0) for i in range(81)]
    need(all(table[table[i][i]][i]==0 for i in range(81)), 'exponent three')
    abelian=all(table[i][j]==table[j][i] for i,j in product(range(81),repeat=2))
    need(abelian != nonabelian,'abelian versus nonabelian')
    S=set(connection_set)
    need(len(S)==40 and 0 not in S and all((i in S)!=(inverse[i] in S) for i in range(1,81)), 'skew connection set')
    # Build arcs as g -> gs rather than testing inverse(g)h membership.
    adjacency=[{table[g][s] for s in S} for g in range(81)]
    rows=[sum(1<<j for j in a) for a in adjacency]
    need(all(i not in adjacency[i] for i in range(81)) and all((j in adjacency[i])!=(i in adjacency[j]) for i in range(81) for j in range(i)), 'tournament')
    for p in table:
        need(all({p[j] for j in adjacency[i]}==adjacency[p[i]] for i in range(81)), 'translation not automorphism')
    partition=[{0},set(range(1,81))]; trace=[2]
    while len(partition)<81:
        refined=[]
        for cell in partition:
            buckets={}
            for v in sorted(cell):
                sig=tuple(len(adjacency[v]&c) for c in partition)
                buckets.setdefault(sig,set()).add(v)
            refined.extend(buckets.values())
        need(len(refined)>len(partition),'refinement stalled before rigidity')
        partition=refined; trace.append(len(partition))
    # Incremental splitter closure with a fixed reference member: any outside
    # vertex seeing that member differently from another must enter the module.
    universe=(1<<81)-1; closure_sizes=Counter()
    for a in range(81):
        for b in range(a):
            inside=(1<<a)|(1<<b); queue=[b]; cursor=0
            while cursor<len(queue):
                x=queue[cursor]; cursor+=1
                split=(rows[a]^rows[x]) & (universe^inside)
                while split:
                    bit=split&-split; split-=bit; v=bit.bit_length()-1
                    inside|=bit; queue.append(v)
            closure_sizes[inside.bit_count()]+=1
    need(closure_sizes=={81:3240}, 'proper module found')
    action_cycles=Counter(cycles(p) for p in table)
    need(action_cycles=={(1,)*81:1,(3,)*27:80},'cycle index')
    # Test every linear functional, rather than preselecting abelianization coordinates.
    kernels=set(); homs=0
    for coeff in elements[1:]:
        f=[sum(a*b for a,b in zip(coeff,x))%3 for x in elements]
        if all(f[table[i][j]]==(f[i]+f[j])%3 for i,j in product(range(81),repeat=2)):
            homs+=1; kernels.add(frozenset(i for i in range(81) if not f[i]))
    need(homs==(26 if nonabelian else 80), 'homomorphism count')
    need(len(kernels)==(13 if nonabelian else 40), 'kernel count')
    # All index-three subgroups are kernels: their coset action is a
    # transitive 3-subgroup of S3, hence C3. The abelianizations are C3^4
    # and C3^3, respectively; in the latter the c-axis is the commutator.
    commutators={table[table[table[inverse[i]][inverse[j]]][i]][j]
                 for i,j in product(range(81),repeat=2)}
    expected_commutators={idx[(0,0,c,0)] for c in range(3)} if nonabelian else {0}
    need(commutators==expected_commutators, 'commutator subgroup')
    # Exact conjugacy and colored-subset orbits for all nonconstant index-3 coset unions.
    subsets=set()
    for K in kernels:
        need(len(K)==27 and all(table[i][j] in K for i in K for j in K),'kernel subgroup')
        cosets={frozenset(table[k][g] for k in K) for g in range(81)}
        need(len(cosets)==3,'coset number')
        for C in cosets:
            subsets.add(C); subsets.add(frozenset(set(range(81))-C))
    orbits=set()
    for C in subsets:
        images={tuple(sorted(table[g][v] for v in C)) for g in range(81)}
        need(len(images)==3,'exact stabilizer not 27')
        orbits.add(min(images))
    need(len(orbits)==(26 if nonabelian else 80),'orbit count')
    weighted=Counter(81+len(C) for C in orbits)
    need(weighted=={108:len(kernels),135:len(kernels)},'unequal-size inflation counts')
    return {'root_refinement_cell_counts':trace,'automorphism_order':81,'pair_module_closures':dict(closure_sizes),'nonzero_C3_homomorphisms':homs,'index3_kernels':len(kernels),'binary_orbits_stabilizer27':len(orbits),'inflated_order108_Aut27':weighted[108],'inflated_order135_Aut27':weighted[135],'inflated_order324_Aut27':len(orbits),'abelian':abelian,'cycle_index':'(p1^81+80*p3^27)/81'}

def ordered(*parts):
    out=[]; offset=0; total=sum(map(len,parts))
    for part in parts:
        later=set(range(offset+len(part),total))
        out.extend({v+offset for v in row}|later for row in part)
        offset+=len(part)
    return out

def tt(n): return [set(range(i+1,n)) for i in range(n)]
def blown_cycle():
    return [{j for j in range(6) if (j//2-i//2)%3==1 or (j//2==i//2 and j>i)} for i in range(6)]
def tournament(A):
    n=len(A)
    need(all(i not in A[i] for i in range(n)), 'loop')
    need(all((j in A[i])!=(i in A[j]) for i in range(n) for j in range(i)), 'tournament arcs')

def component_word(A):
    # Reachability gives strong components; tournament condensation is a chain.
    reach=[set(row)|{i} for i,row in enumerate(A)]
    for k in range(len(A)):
        for i in range(len(A)):
            if k in reach[i]: reach[i]|=reach[k]
    components={frozenset(j for j in range(len(A)) if j in reach[i] and i in reach[j]) for i in range(len(A))}
    ordered_components=sorted(components,key=lambda C:-len(reach[min(C)]))
    return tuple(map(len,ordered_components))

def inventory():
    C=[{1},{2},{0}]; one=tt(1); D=blown_cycle()
    L=[ordered(C,C,one,one),ordered(C,one,C,one),ordered(one,C,C,one),tt(8)]
    R=[ordered(C,tt(5)),ordered(one,C,tt(4)),ordered(D,tt(2)),ordered(one,D,one)]
    words=[component_word(A) for A in L+R]
    need(len(set(words))==8, 'inventory repeats an isomorphism class')
    result=[]
    for family in [L,R]:
        agg=Counter(); orders=[]
        for A in family:
            tournament(A)
            autos=[p for p in permutations(range(8)) if all({p[j] for j in A[i]}==A[p[i]] for i in range(8))]
            orders.append(len(autos))
            for p in autos: agg[cycles(p)]+=Fraction(1,len(autos))
        result.append((agg,orders))
    need(result[0][0]==result[1][0],'inventory cycle index differs')
    need(result[0][0]=={(1,)*8:Fraction(4,3),(1,)*5+(3,):Fraction(4,3),(1,1,3,3):Fraction(4,3)},'aggregate polynomial')
    need(result[0][1]==[9,9,9,1] and result[1][1]==[3]*4,'inventory orders')
    return {'L_automorphism_orders':result[0][1],'R_automorphism_orders':result[1][1],'aggregate':{str(k):str(v) for k,v in result[0][0].items()}}

if __name__=='__main__':
    output={name:audit(name,S) for name,S in CONNECTION_SETS.items()}
    output['order8_inventory']=inventory()
    palette=[tt(4),[{1,2},{2},{3},{0,1}]]
    # Second palette element: blocks {0,1}->{2}->{3}->{0,1}, with 0->1.
    for A in palette:
        tournament(A)
        need(sum(all({p[j] for j in A[i]}==A[p[i]] for i in range(4)) for p in permutations(range(4)))==1,'order4 palette rigidity')
    need(sorted(map(len,palette[0]))!=sorted(map(len,palette[1])),'palette nonisomorphism')
    output['order4_palette_rigid_nonisomorphic']=True
    print(json.dumps(output,indent=2,sort_keys=True))
