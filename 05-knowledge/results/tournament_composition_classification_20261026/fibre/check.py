"""Exact cyclic-triple fibres, decorated symmetry, and the sharp n=12 hostile.

No producer imports, external libraries, or Python assertion gates.
"""
import itertools
import json
import sys
from collections import Counter
from fractions import Fraction
from math import factorial

sys.stdout.reconfigure(newline="\n")


def need(test, message):
    if not test:
        raise RuntimeError(message)


def mask(rows):
    return sum(((rows[i] >> j) & 1) << k
               for k, (i, j) in enumerate(itertools.combinations(range(len(rows)), 2)))


def frommask(n, value):
    rows = [0]*n
    for k, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        a,b=(i,j) if value >> k & 1 else (j,i)
        rows[a] |= 1 << b
    return tuple(rows)


def converse(rows):
    full=(1 << len(rows))-1
    return tuple(full ^ row ^ (1 << i) for i,row in enumerate(rows))


def relabel(rows, perm):
    out=[0]*len(rows)
    for i,row in enumerate(rows):
        for j in range(len(rows)):
            if row >> j & 1:
                out[perm[i]] |= 1 << perm[j]
    return tuple(out)


def cyclic(rows):
    return frozenset((i,j,k) for i,j,k in itertools.combinations(range(len(rows)),3)
                     if ((rows[i]>>j)&1)==((rows[j]>>k)&1)==((rows[k]>>i)&1))


def lift(root_bit, triangle_bits, child_types):
    """C3[A_0,A_1,A_2], A_0=C3->1 and A_1=1->C3 on fixed labels."""
    rows=[0]*12
    for i,j in itertools.combinations(range(12),2):
        bi,ri=divmod(i,4); bj,rj=divmod(j,4)
        if bi != bj:
            forward=((bj-bi)%3==1) ^ bool(root_bit)
        elif rj==3:
            forward=child_types[bi]==0
        else:
            forward=((rj-ri)%3==1) ^ bool(triangle_bits[bi])
        a,b=(i,j) if forward else (j,i)
        rows[a] |= 1 << b
    return tuple(rows)


def automorphisms(rows, fixed=()):
    """All permutations by degree/adjacency constrained backtracking, independently
    of modules and the predicted hypergraph automorphism group."""
    n=len(rows); degrees=[row.bit_count() for row in rows]
    choices=[[j for j in range(n) if degrees[j]==degrees[i] and (i not in fixed or j==i)] for i in range(n)]
    order=sorted(range(n),key=lambda i:len(choices[i]))
    mapping=[-1]*n; total=0
    def walk(k,used):
        nonlocal total
        if k==n:
            total+=1; return
        v=order[k]
        for w in choices[v]:
            if used >> w & 1:
                continue
            if all(((rows[v]>>u)&1)==((rows[w]>>mapping[u])&1) for u in order[:k]):
                mapping[v]=w
                walk(k+1,used|(1<<w))
                mapping[v]=-1
    walk(0,0)
    return total


def hamiltonian(rows):
    n=len(rows); size=1<<n
    dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for subset in range(1,size):
        for v,count in enumerate(dp[subset]):
            if not count: continue
            choices=rows[v] & ~subset
            while choices:
                bit=choices & -choices; choices-=bit
                dp[subset|bit][bit.bit_length()-1]+=count
    return sum(dp[-1])


def sign(perm):
    return (-1)**sum(perm[i]>perm[j] for i in range(len(perm)) for j in range(i+1,len(perm)))


def path_cover_profile(rows):
    n=len(rows); size=1<<n
    dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for subset in range(1,size):
        for v,count in enumerate(dp[subset]):
            if not count: continue
            choices=rows[v] & ~subset
            while choices:
                bit=choices & -choices; choices-=bit
                dp[subset|bit][bit.bit_length()-1]+=count
    hp=[sum(a) for a in dp]
    covers=[[0]*(n+1) for _ in range(size)]; covers[0][0]=1
    for subset in range(1,size):
        first=subset & -subset
        rest=subset ^ first; part=rest
        while True:
            block=part|first
            for k,value in enumerate(covers[subset ^ block][:-1]):
                covers[subset][k+1]+=hp[block]*value
            if not part: break
            part=(part-1)&rest
    return [factorial(k)*value for k,value in enumerate(covers[-1])]


def prime_quotient(rows):
    n=len(rows)
    for subset in range(1,1<<n):
        if subset.bit_count() in (1,n): continue
        if all((rows[z]&subset) in (0,subset) for z in range(n) if not (subset>>z)&1):
            return False
    return True


def substitute_four(q,types):
    n=4*len(q); out=[0]*n
    for i,j in itertools.combinations(range(n),2):
        bi,ri=divmod(i,4); bj,rj=divmod(j,4)
        if bi!=bj: forward=(q[bi]>>bj)&1
        elif rj==3: forward=types[bi]==0
        else: forward=(rj-ri)%3==1
        a,b=(i,j) if forward else (j,i)
        out[a]|=1<<b
    return tuple(out)


def main():
    # Independently exhaustive small lemma: the sole n=4 C3-fibre with
    # multiple ordinary tournament types is a triangle and an isolated vertex.
    four={}
    for value in range(64):
        rows=frommask(4,value)
        key=min(mask(relabel(rows,p)) for p in itertools.permutations(range(4)))
        four.setdefault(cyclic(rows),set()).add(key)
    ambiguous={h:types for h,types in four.items() if len(types)>1}
    need(len(ambiguous)==4 and all(len(h)==1 and len(types)==2 for h,types in ambiguous.items()),'four-vertex ambiguity')
    need(len(set.union(*four.values()))==4,'four tournament classes')

    fibre={mask(lift(root,tri,types)):lift(root,tri,types)
           for root in range(2) for tri in itertools.product(range(2),repeat=3)
           for types in itertools.product(range(2),repeat=3)}
    need(len(fibre)==128,'fibre cardinality')
    h=cyclic(next(iter(fibre.values())))
    need(len(h)==67 and all(cyclic(t)==h for t in fibre.values()),'literal cyclic-triple equality')
    # Pair codegrees recover the three blocks. This proves these generated
    # hypergraph permutations are all of Aut(H), not just a subgroup.
    codegree={(i,j):sum(i in e and j in e for e in h)
              for i,j in itertools.combinations(range(12),2)}
    need(all((d<4)==(i//4==j//4) for (i,j),d in codegree.items()),'block recovery by codegree')
    perms=[]
    for blocks in itertools.permutations(range(3)):
        for inside in itertools.product(list(itertools.permutations(range(3))),repeat=3):
            p=tuple(4*blocks[b]+(inside[b][r] if r<3 else 3) for b in range(3) for r in range(4))
            need(frozenset(tuple(sorted(p[x] for x in e)) for e in h)==h,'H automorphism')
            perms.append(p)
    need(len(set(perms))==1296,'hypergraph group')
    done=set(); classes=[]; profiles=[]
    for number_b in range(4):
        types=(1,)*number_b+(0,)*(3-number_b)
        t=lift(0,(0,0,0),types)
        orbit={mask(relabel(t,p)) for p in perms}
        need(orbit<=fibre.keys() and not orbit&done,'disjoint lift orbits')
        done |= orbit
        aut=automorphisms(t)
        need(aut==(81 if number_b in (0,3) else 27),'literal automorphism count')
        need(len(orbit)*aut==1296,'orbit-stabilizer')
        hp=hamiltonian(t)
        profiles.append(path_cover_profile(t))
        classes.append({'B_children':number_b,'canonical_in_H_orbit':str(min(orbit)),
                        'representative_mask':str(mask(t)),'representative_outrows':list(t),
                        'automorphism_order':aut,'H_orbit_size':len(orbit),
                        'Hamiltonian_paths':hp,'tc':str(Fraction(hp,aut))})
    need(done==fibre.keys(),'four classes exhaust fibre')
    all_paths={hamiltonian(t) for t in fibre.values()}
    need(len(all_paths)==1,'H invariant on all128 lifts')
    need(all(p==profiles[0] for p in profiles),'full path-cover profile invariant')
    need(sum(Fraction(1,c['automorphism_order']) for c in classes)==Fraction(128,1296),'groupoid mass')

    # Independent adjacent correction: even permutation can send T to T*.
    t5=tuple(sum(1<<j for j in range(5) if (j-i)%5 in (1,2)) for i in range(5))
    anti=tuple((-i)%5 for i in range(5))
    need(sign(anti)==1 and relabel(t5,anti)==converse(t5),'A5/C3 chirality distinction')
    prime_cases=[]
    for p in (3,5,7,11,13):
        q=tuple(sum(1<<j for j in range(p) if 1<=(j-i)%p<=(p-1)//2) for i in range(p))
        need(prime_quotient(q),'prime-order circulant is prime')
        a=automorphisms(q); point=automorphisms(q,fixed=(0,))
        need(a==p*point,'prime orbit stabilizer')
        uniform=substitute_four(q,(0,)*p)
        marked=substitute_four(q,(1,)+(0,)*(p-1))
        need(cyclic(uniform)==cyclic(marked),'sharp4p same labelled C3')
        prime_cases.append({'p':p,'vertices':4*p,'quotient_aut':a,'point_stabilizer':point,
                            'uniform_aut':3**p*a,'one_changed_child_aut':3**p*point,'ratio':p})
    return {'status':'FINITE-EXACT','universes':{'all_labelled_order4':64,'labelled_order12_lifts':128,
            'hypergraph_automorphisms':1296,'direct_tournament_automorphism_representatives':4},
            'order4':{'C3_hypergraphs':len(four),'ambiguous_labelled_hypergraphs':len(ambiguous),'ambiguous_structure':'one triple plus one isolated vertex'},
            'order12':{'C3_edges':sorted(map(list,h)),'C3_edge_count':len(h),
                        'classes':classes,'all_lift_Hamiltonian_counts':sorted(all_paths),
                        'AutH':1296,'weighted_mass':'8/81','ordered_path_cover_profile':profiles[0]},
            'sharp_prime_threshold_controls':prime_cases,
            'adjacent_chirality_control':{'T5_rows':t5,'anti_permutation':anti,'permutation_sign':sign(anti)}}


if __name__=='__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
