#!/usr/bin/env python3
"""Complete pair inventory and scoped controls for H4 type (4)(2).

Status: RESERVED pending independent audit.  The proof is uniform in the
ambient degree: a pair has union support at most twelve, while a joint
centralizer and the proved single-cycle theorem pay the four-generator step.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json


gates=0


def need(condition,label):
    global gates
    gates+=1
    if not condition:
        raise RuntimeError(label)


def mul(p,q):
    return tuple(p[x] for x in q)


def inv(p):
    result=[0]*len(p)
    for i,j in enumerate(p):result[j]=i
    return tuple(result)


def cyc(p,block):
    result=list(p)
    for x,y in zip(block,block[1:]+block[:1]):result[x]=y
    return tuple(result)


def cycles(p):
    seen=set();result=[]
    for x in range(len(p)):
        if x in seen or p[x]==x:continue
        row=[];y=x
        while y not in seen:
            row.append(y);seen.add(y);y=p[y]
        result.append(tuple(row))
    return tuple(sorted(result,key=lambda row:(-len(row),row)))


def flags(p,q):
    pq=mul(p,q);qp=mul(q,p)
    return (mul(pq,p)==mul(qp,q),
            mul(mul(pq,pq),p)==mul(mul(qp,qp),q),pq==qp)


def letter_relation(p,q,length):
    # This direct nested letter action does not call mul or flags.
    for x in range(len(p)):
        left=right=x
        for j in range(length-1,-1,-1):
            left=(p if j%2==0 else q)[left]
            right=(q if j%2==0 else p)[right]
        if left!=right:return False
    return True


def joint_orbits(p,q):
    seen=set();result=[]
    for x in range(len(p)):
        if x in seen:continue
        seen.add(x);row=[x]
        for y in row:
            for g in (p,q):
                z=g[y]
                if z not in seen:
                    seen.add(z);row.append(z)
        if len(row)>1:result.append(tuple(sorted(row)))
    return tuple(sorted(result,key=lambda row:(len(row),row)))


def class42(degree):
    identity=tuple(range(degree))
    for four in combinations(range(degree),4):
        rest=tuple(x for x in range(degree) if x not in four)
        # The minimum of the four-cycle is first; all six cyclic orders
        # occur once.  The unordered two-set has only one transposition.
        for order in permutations(four[1:]):
            p=cyc(identity,(four[0],)+order)
            for two in combinations(rest,2):yield cyc(p,two)


D=12;identity=tuple(range(D))
sigma=cyc(cyc(identity,(0,1,2,3)),(4,5))
A4=set(range(4));A2={4,5};A=A4|A2
universe=tuple(class42(D))
need(len(universe)==len(set(universe))==83160,'full unfiltered class of partners')
need(83160==495*6*28==924*90,'two independent combinatorial universe counts')

expected_cells=(
    {(2,0,0,1):480,(2,0,0,2):60,(2,2,2,0):4,(4,0,0,1):60,(4,0,0,2):5},
    {(2,0,0,2):120,(3,0,0,2):24,(3,1,1,1):16,(4,0,0,2):1},
    {(0,0,0,0):90,(0,0,0,2):90,(4,0,0,0):30,(4,0,0,2):2},
)
expected_orbits=(
    {(2,4):5,(2,6):60,(3,4):60,(3,6):480,(6,):4},
    {(2,4):1,(2,5):24,(2,6):120,(6,):16},
    {(2,2,4):30,(2,2,4,4):90,(2,4):2,(2,4,4):90},
)
expected_cell_orbits=(
    {(2,0,0,1):(3,6),(2,0,0,2):(2,6),(2,2,2,0):(6,),
     (4,0,0,1):(3,4),(4,0,0,2):(2,4)},
    {(2,0,0,2):(2,6),(3,0,0,2):(2,5),(3,1,1,1):(6,),
     (4,0,0,2):(2,4)},
    {(0,0,0,0):(2,2,4,4),(0,0,0,2):(2,4,4),
     (4,0,0,0):(2,2,4),(4,0,0,2):(2,4)},
)
cell_hist=[Counter() for _ in range(3)]
orbit_hist=[Counter() for _ in range(3)]
admitted=[[],[],[]];data={};raw=[]
for tau in universe:
    cc=cycles(tau)
    need(tuple(map(len,cc))==(4,2),'literal cycle type')
    ff=flags(sigma,tau)
    for length,flag in zip((3,5),ff):
        need(letter_relation(sigma,tau,length)==flag,'independent odd letter action')
    need(all(sigma[tau[x]]==tau[sigma[x]] for x in range(D))==ff[2],
         'independent commuting letter action')
    B4,B2=map(set,cc)
    cell=(len(A4&B4),len(A4&B2),len(A2&B4),len(A2&B2))
    oo=joint_orbits(sigma,tau)
    sizes=tuple(map(len,oo))
    raw.append((tau,ff,cell,sizes))
    data[tau]=(B4,B2,cell,oo)
    for k,yes in enumerate(ff):
        if yes:
            admitted[k].append(tau)
            cell_hist[k][cell]+=1;orbit_hist[k][sizes]+=1
            need(expected_cell_orbits[k].get(cell)==sizes,
                 'full typed-cell and joint-orbit coupling')
    if ff[0]:
        need(len(A4&B4)>=2,'ordinary four-block overlap')
        need(len(A4&B4)==4 or 6 in sizes,'ordinary unequal four-blocks create joint six-orbit')
        if sizes==(2,6):
            need(A2==B2,'two-six ordinary pair shares transposition')
    if ff[1]:
        need(len(A4&B4)>=2,'fifth four-block overlap')
        need(A2==B2 or cell==(3,1,1,1),'fifth transposition-sharing dichotomy')
        if len(A4&B4)==4:
            need(tau==sigma,'fifth common four-block forces exact equality')
    if ff[2]:
        need(not(A4&B2) and not(A2&B4),'commuting unequal-size blocks do not mix')
        need(A4==B4 or not(A4&B4),'commuting four-block equal or disjoint')
        need(A2==B2 or not(A2&B2),'commuting two-block equal or disjoint')

for k in range(3):
    need(dict(cell_hist[k])==expected_cells[k],'complete literal cycle-intersection table')
    need(dict(orbit_hist[k])==expected_orbits[k],'complete joint-orbit size table')
need(tuple(map(len,admitted))==(609,161,212),'all admitted partner totals')

# A commuting permutation permutes joint orbits of equal size.  The
# ordinary inventory has no repeated nontrivial sizes.  These always-active
# controls retain the full literal joint centralizers inside S_12; they are
# checks of the separate all-degree argument, not a bound on a whole action.
ordinary,fifth,commuting=admitted
centralizer_trials=0;centralizer_edges=0;fixed_six_edges=0
propagation_trials=0
for tau in ordinary:
    B4,B2,cell,oo=data[tau]
    sizes=tuple(map(len,oo))
    need(len(sizes)==len(set(sizes)),'ordinary joint sizes all distinct')
    for d in commuting:
        centralizer_trials+=1
        if mul(tau,d)!=mul(d,tau):continue
        centralizer_edges+=1
        for orbit in oo:
            need({d[x] for x in orbit}==set(orbit),'literal joint orbit preserved')
            if len(orbit) in (3,6):
                need(all(d[x]==x for x in orbit),'literal three-six pointwise fixed')
                if len(orbit)==6:fixed_six_edges+=1
    # This is exactly the source of the shared-transposition propagation.
    # It includes every c commuting with sigma inside the twelve-label
    # universe.  The proof explains why the same implication is uniform.
    if sizes==(2,6):
        for c in commuting:
            if not letter_relation(tau,c,3):continue
            C4,C2,_,_=data[c]
            if C4==A4:continue
            propagation_trials+=1
            need(not(C4&A),'commuting outside four-block avoids both original blocks')
            need(C2==A2,'whole two-six propagation to c transposition')
            need(len(C4&B4)==2,'two-six c four-block overlap exactly two')

# The small joint-centralizer obstruction has no allowed semiregular
# restriction on six or three points for a global type (4)(2) permutation.
allowed_restrictions=[(),(2,),(4,),(4,2)]
for size in (3,6):
    semiregular=[]
    for typ in allowed_restrictions:
        if sum(typ)==size and len(set(typ))==1:semiregular.append(typ)
    need(semiregular==[],'no nonidentity semiregular type on three or six points')
need(6%4!=0 and 6//2==3,'six-orbit order-four/order-two numerical obstruction')

# Minimal ambient-size hostiles: extracting the transposition from an odd
# pair need not preserve that odd relation.  The global proof must pay the
# common-transposition step before deleting it.
hostile3=cyc(cyc(identity,(0,4,2,5)),(1,3))
hostile5=cyc(cyc(identity,(0,2,1,4)),(3,5))
trans_sigma=cyc(identity,(4,5))
hostile_rows=[]
for length,tau in ((3,hostile3),(5,hostile5)):
    need(letter_relation(sigma,tau,length),'full mixed pair satisfies named odd relation')
    trans_tau=cyc(identity,cycles(tau)[1])
    need(not letter_relation(trans_sigma,trans_tau,length),'isolated transpositions lose odd relation')
    need(len(set().union(*map(set,cycles(tau)),A))==6,'hostile reaches minimal ambient support six')
    hostile_rows.append((length,cycles(tau),data[tau][2]))

# The equal tuple is a real positive control, but is always intransitive:
# its two nontrivial cycle orbits already have different sizes.
need(flags(sigma,sigma)==(True,True,True),'equal four-generator positive control')
need(tuple(map(len,joint_orbits(sigma,sigma)))==(2,4),'equal tuple is not transitive')

print('h4_mixed42: PASS')
print('scope: all-degree equality for H4 generators of type(4)(2); no ambient-degree tuple census')
print('full unfiltered pair universe:',len(universe),'on twelve labels')
for name,k in (('ordinary',0),('fifth',1),('commuting',2)):
    print(name,'partners:',len(admitted[k]))
    print(name,'cells (4/2 rows and columns):',sorted(cell_hist[k].items()))
    print(name,'joint nontrivial orbit sizes:',sorted(orbit_hist[k].items()))
print('literal centralizer controls (trials,edges,fixed-six):',
      (centralizer_trials,centralizer_edges,fixed_six_edges))
print('two-six common-transposition propagation controls:',propagation_trials)
print('minimal split-projection hostiles:',hostile_rows)
print('gates:',gates)
print('raw semantic sha256:',sha256(json.dumps(raw,separators=(',',':')).encode()).hexdigest())
