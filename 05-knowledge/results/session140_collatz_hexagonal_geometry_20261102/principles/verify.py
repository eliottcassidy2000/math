#!/usr/bin/env python3
"""Independent exact controls for the Session140 principles note.
No repository source is imported. Assertions are deliberately not used.
"""
from fractions import Fraction as Q
from itertools import product, permutations
import json

def check(test, message):
    if not test:
        raise RuntimeError(message)

def T(n):
    return n//2 if n % 2 == 0 else 3*n+1

def inv_candidates(n):
    es = (2,4,6) if n % 3 == 1 else (3,5,7)
    out = []
    for e in es:
        z = (2**e*n-1)//3
        if z % 3:
            out.append((e,z))
    return out

# Complete residue proof checks: n mod 18 determines all divisibility guards.
residue_controls=[]
for n in (5,7,11,13,17,19):
    out=inv_candidates(n)
    check(len(out)==2, 'exactly two modulo9 choices')
    for e,z in out:
        check(z%2==1 and z%3 and z>n, 'positive increasing odd inverse')
        check(3*z+1==2**e*n, 'inverse identity')
        v=z
        for j in range(e+1):
            v=T(v)
        check(v==n, 'exact standard time')
    residue_controls.append({'n':n,'children':out})

# Literal, distinct-labelled macro tree. This does NOT assert internally
# disjoint route paths: even segments can and do overlap.
level={5}; seen={5}; levels=[]; max_time=0
for depth in range(1,17):
    next_level=set()
    for n in level:
        for e,z in inv_candidates(n):
            check(z not in next_level and z not in seen, 'global child collision')
            check((3*z+1)>>e==n and (3*z+1)%(2**e)==0, 'route arithmetic')
            v=z
            for _ in range(e+1):
                v=T(v)
            check(v==n, 'literal T route')
            next_level.add(z)
            max_time=max(max_time,e+1)
    check(len(next_level)==2**depth, 'binary level count')
    seen.update(next_level); level=next_level
    levels.append({'depth':depth,'leaves':len(level),'total_odd_vertices':len(seen),
                   'max_standard_radius_bound':8*depth,
                   'triangular_ball_at_radius_bound':1+3*(8*depth)*(8*depth+1)})
check(max_time==8, 'eight-step upper bound attained')
check(levels[-1]['leaves']>levels[-1]['triangular_ball_at_radius_bound'],
      'one-edge embedding count contradiction')

# Triangular lattice unit distances, proper 3-coloring, exact small balls.
steps=((1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1))
check(all(a*a+a*b+b*b==1 for a,b in steps), 'unit step norm')
check(all((a-b)%3 for a,b in steps), 'proper 3-coloring')
balls=[]
ball={(0,0)}
for r in range(0,21):
    if r:
        ball |= {(a+da,b+db) for a,b in tuple(ball) for da,db in steps}
    check(len(ball)==1+3*r*(r+1),'triangular exact ball count')
    balls.append(len(ball))
# Honeycomb is the induced graph on colors 1 and 2.
for a,b in product(range(-4,5),repeat=2):
    c=(a-b)%3
    if c:
        nbr=[((a+da-b-db)%3) for da,db in steps if (a+da-b-db)%3]
        check(len(nbr)==3 and all(d!=c for d in nbr),'honeycomb bipartition')
check([T(n) for n in (1,4,2)]==[4,2,1],'standard Collatz directed triangle')

# First-return lifting theorem: all maps on at most4 states and all clocks
# tau in {1,2}. Canonical suspension is independently enumerated.
def orbit(f,start):
    seen={}; seq=[]; v=start
    while v not in seen:
        seen[v]=len(seq);seq.append(v);v=f(v)
    return seq,seen[v],len(seq)-seen[v]
clock_systems=0;clock_root_tests=0
for n in range(1,5):
    for fv in product(range(n),repeat=n):
        for tv in product((1,2),repeat=n):
            def R(v):
                x,j=v
                return (x,j+1) if j+1<tv[x] else (fv[x],0)
            for x in range(n):
                fseq,mu,period=orbit(lambda z:fv[z],x)
                rseq,rmu,rperiod=orbit(R,(x,0))
                landmarks=[v[0] for v in rseq if v[1]==0]
                check(landmarks==fseq,'first-return orbit/target preservation')
                check(rperiod==sum(tv[z] for z in fseq[mu:]),'least cycle clock sum')
                clock_root_tests+=1
            clock_systems+=1
check(clock_systems==4330,'exact clock universe')
# Endpoint-only routes fail: F=id, R=swap, proposed tau=2.
check(all((x^1)^1==x for x in (0,1)),'endpoint route equations')
check((0^1)==1 and 0!=1,'new reachability despite endpoint equations')

# Explicit ring-order site bijection; triangular graph radius is max(|a|,|b|,|a+b|).
ring_points=[(0,0)];ring_dirs=((-1,1),(-1,0),(0,-1),(1,-1),(1,0),(0,1))
for r in range(1,101):
    a,b=r,0
    for da,db in ring_dirs:
        for _ in range(r):
            check(max(abs(a),abs(b),abs(a+b))==r,'site ring radius')
            ring_points.append((a,b));a+=da;b+=db
    check((a,b)==(r,0),'ring closes')
check(len(ring_points)==len(set(ring_points))==30301,'ring-order injection/count')
site_index={v:i+1 for i,v in enumerate(ring_points)}
for n in range(1,10001):
    image=ring_points[T(n)-1]
    check(site_index[image]==T(n),'exact site-conjugate successor')
for n in range(1,100001):
    check((n.bit_length()-1)%3!=(T(n).bit_length()-1)%3,'explicit bitlength3-colouring')

# Proper marginals lose exactly the parity coordinate; a pinned context
# recovers it. This checks loss without assuming any unprescribed HN gadget.
B=list(product((0,1),repeat=3))
even={b:int(sum(b)%2==0) for b in B}
odd={b:int(sum(b)%2==1) for b in B}
for missing in range(3):
    ix=[i for i in range(3) if i!=missing]
    for vals in product((0,1),repeat=2):
        em=sum(v for b,v in even.items() if tuple(b[i] for i in ix)==vals)
        om=sum(v for b,v in odd.items() if tuple(b[i] for i in ix)==vals)
        check(em==om,'parity tables equal proper counting marginals')
check(even[(0,0,0)]!=odd[(0,0,0)],'pin distinguishes parity tables')

# Subdivision preserves path existence but not chromatic obstruction:
# K4 is not 3-colorable; once subdividing every edge is bipartite.
K4=[(i,j) for i in range(4) for j in range(i+1,4)]
check(not any(all(c[i]!=c[j] for i,j in K4) for c in product(range(3),repeat=4)),
      'K4 chromatic lower bound')
sub_edges=[(i,4+k) for k,(i,j) in enumerate(K4)]+[(j,4+k) for k,(i,j) in enumerate(K4)]
sub_colors=[0]*4+[1]*6
check(all(sub_colors[i]!=sub_colors[j] for i,j in sub_edges),'subdivision bipartite')

# Chronological Collatz odd words (1,7),(7,1) have same diagonal/trace,
# distinct carries, and different source guards.
def affine(word):
    a,b,d=1,0,1
    for e in word:
        a,b,d=3*a,3*b+d,2**e*d
    return a,b,d
u=affine((1,7)); v=affine((7,1))
check(u==(9,5,256) and v==(9,131,256),'ordered carry values')
check(u[0]+u[2]==v[0]+v[2] and u[1]!=v[1],'trace blindness')
check(Q(u[0]*483+u[1],u[2])==Q(v[0]*469+v[1],v[2])==17,'carry-controlled merge')

# All-modulus 2-variable jet theorem controls, independently from coefficients.
# P=c00+c10*x+c01*y+c11*x*y, top c11 in a mod m, others in mZ.
jet_checks=0
for m in range(2,8):
    for a in range(1,m):
        delta=min(a,m-a)
        for z00,z10,z01,ztop in product((-1,0,1),repeat=4):
            c00,c10,c01,c11=m*z00,m*z10,m*z01,a+m*ztop
            for x,y in product((Q(-1),Q(-1,2),Q(0),Q(1,2),Q(1)),repeat=2):
                val=c00+c10*x+c01*y+c11*x*y
                px=c10+c11*y; py=c01+c11*x
                energy=val*val+px*px+py*py+c11*c11
                floor=delta*delta*(1+x*x)*(1+y*y)
                check(energy>=floor,'full all-modulus tensor norm floor')
                jet_checks+=1
check((Q(3)*Q(2,3)-2)==0,'single value can vanish')
check(Q(0)**2+Q(1)**2<2,'x-1 violates floor without proper coefficient divisibility')

# Real/integral repair distinction for the triangle Laplacian.
L=((2,-1,-1),(-1,2,-1),(-1,-1,2))
def matvec(A,x): return tuple(sum(a*b for a,b in zip(row,x)) for row in A)
v=(1,-1,0)
check(matvec(L,v)==(3,-3,0),'triangle torsion multiplied by3 is repairable')
check(matvec(L,(Q(1,3),Q(-1,3),Q(0)))==v,'real repair exists')
# If Lz=v with integer z then first minus second equation gives
# 3(z0-z1)=2, impossible. This exact congruence certifies all integer z.
check((v[0]-v[1])%3!=0,'integral repair obstruction')

result={'status':'PASS','source_imports':False,
        'inverse_residue_controls':residue_controls,'macro_tree_levels':levels,
        'max_macro_standard_time':max_time,
        'triangular_ball_counts_r0_to20':balls,
        'tensor_jet_checks':jet_checks,'first_return_systems':clock_systems,
        'first_return_root_tests':clock_root_tests,'ring_sites':len(ring_points),
        'site_conjugacy_checks':10000,'bitlength_colouring_checks':100000,
        'controls':{'parity_marginal_hostile':True,'K4_subdivision_hostile':True,
                    'Collatz_trace_carry_hostile':{'u':u,'v':v,'common_trace':265,'merge':17},
                    'triangle_real_integral_repair_hostile':True},
        'scope':'Finite exact controls supplement the written proofs; no convergence or chromatic-bound advance.'}
print(json.dumps(result,sort_keys=True,indent=2))
