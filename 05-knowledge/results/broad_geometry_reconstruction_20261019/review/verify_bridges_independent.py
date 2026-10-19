#!/usr/bin/env python3
"""Independent finite controls: no producer script is imported.
Uses direct moment signatures, full Padé systems, row-deletion separators,
parity-syndrome dynamic programming, and local Smith certificates.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, prod
from collections import Counter
import json
import sympy as sp


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def rref(a, columns=None):
    a=[[Q(x) for x in row] for row in a]
    n=columns if columns is not None else len(a[0])
    row=0; pivots=[]
    for j in range(n):
        pivot=next((i for i in range(row,len(a)) if a[i][j]),None)
        if pivot is None: continue
        a[row],a[pivot]=a[pivot],a[row]
        v=a[row][j]; a[row]=[x/v for x in a[row]]
        for i in range(len(a)):
            if i!=row and a[i][j]:
                v=a[i][j]; a[i]=[x-v*y for x,y in zip(a[i],a[row])]
        pivots.append(j); row+=1
        if row==len(a): break
    return a,pivots


def rank(a): return len(rref(a)[1])


def endpoint_runs(s):
    s=set(s)
    return sorted((x,next(y for y in range(x+1,max(s)+2) if y not in s))
                  for x in s if x-1 not in s) if s else []


def coefficients_from_moments(m,k):
    d=[0]+[sum(comb(n,j)*m[j] for j in range(n)) for n in range(1,2*k+1)]
    # Exponentiate the series by multiplying exp(D_n z^n/n) factors,
    # rather than using the producer's differential recurrence.
    c=[Q(1)]+[Q(0)]*(2*k)
    for n in range(1,2*k+1):
        exponential=[Q(0)]*(2*k+1); exponential[0]=1
        value=Q(1)
        for j in range(1,2*k//n+1):
            value*=Q(d[n],n*j); exponential[n*j]=value
        c=[sum(c[i]*exponential[h-i] for i in range(h+1)) for h in range(2*k+1)]
    return c


def decode(m,k):
    c=coefficients_from_moments(m,k)
    for d in range(k+1):
        rows=[]
        for n in range(1,2*k+1):
            # Unknowns are p_1..p_d,q_1..q_d simultaneously.
            rows.append([int(n==j) for j in range(1,d+1)]
                        +[-c[n-j] if n>=j else 0 for j in range(1,d+1)]+[c[n]])
        reduced,pivots=rref(rows,2*d)
        if any(not any(row[:2*d]) and row[-1] for row in reduced): continue
        require(len(pivots)==2*d,'minimal Pade pair not unique')
        sol=[Q(0)]*(2*d)
        for i,j in enumerate(pivots): sol[j]=reduced[i][-1]
        return [Q(1)]+sol[:d],[Q(1)]+sol[d:]
    raise RuntimeError('no bounded-degree endpoint ratio')


def polynomial(endpoints):
    p=[Q(1)]
    for endpoint in endpoints:
        nxt=[Q(0)]*(len(p)+1)
        for i,c in enumerate(p): nxt[i]+=c; nxt[i+1]-=endpoint*c
        p=nxt
    return p


def moment_controls():
    counts={}; zero=negative=0
    grid=list(range(-3,5))
    for k in (1,2,3):
        signatures={}; tested=0
        for mask in range(1<<len(grid)):
            s=[x for i,x in enumerate(grid) if mask>>i&1]
            runs=endpoint_runs(s)
            if len(runs)>k: continue
            moments=tuple(sum(x**j for x in s) for j in range(2*k))
            require(moments not in signatures or signatures[moments]==s,'moment signature collision')
            signatures[moments]=s
            P,R=decode(moments,k)
            require(P==polynomial(a for a,b in runs),'entrance polynomial')
            require(R==polynomial(b for a,b in runs),'exit polynomial')
            require(len(P)-1==len(R)-1==len(runs),'homogeneous degree / zero endpoint')
            zero+=any(a==0 or b==0 for a,b in runs)
            negative+=any(a<0 for a,b in runs)
            tested+=1
        counts[k]=tested
    first=[1,5,6]; second=[2,3,7]
    m1=[sum(x**j for x in first) for j in range(4)]
    m2=[sum(x**j for x in second) for j in range(4)]
    require(m1[:3]==m2[:3] and m1[3]!=m2[3],'discrete k2 sharpness')
    vandermonde=[]
    for k in range(1,7):
        endpoints=list(range(-k,k))
        J=sp.Matrix([[(-1 if col%2==0 else 1)*x**j for col,x in enumerate(endpoints)]
                      for j in range(2*k)])
        expected=(-1)**k*prod(endpoints[j]-endpoints[i] for i in range(2*k) for j in range(i+1,2*k))
        require(J.det()==expected and expected!=0,'signed Vandermonde')
        require(J[:-1,:].rank()==2*k-1,'continuous submersion')
        vandermonde.append(k)
    return {'grid':grid,'decodes_by_k':counts,'zero_endpoint_cases':zero,
            'negative_endpoint_cases':negative,'k2_sharpness_moments':[m1,m2],
            'continuous_Jacobian_k':vandermonde,
            'boundary':'No all-k discrete sharpness or stability claim is tested.'}


def evaluation(v,q):
    return [[x**a*y**b*z**(q-a-b) for a in range(q+1) for b in range(q-a+1)] for x,y,z in v]


def conductor_controls():
    records=[]; triple_count=0; separator_tests=0
    for r in range(6,14):
        for off in (False,True):
            vectors=[(1,t,t*t) for t in range(r)]
            if off: vectors[-1]=(1,r-1,(r-1)**2+1)
            for ids in combinations(range(r),3):
                a,b,c=(vectors[i] for i in ids)
                determinant=(a[0]*(b[1]*c[2]-b[2]*c[1])
                             -a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))
                require(determinant!=0,'uniform rank3 matroid')
                triple_count+=1
            hs=[]; separators=[None]*r
            for q in range((r+1)//2+1):
                matrix=evaluation(vectors,q); h=rank(matrix); hs.append(h)
                expected=min(r,2*q+1) if not off or q<2 else min(r,2*q+2)
                require(h==expected,'Hilbert-rank formula')
                for i in range(r):
                    if separators[i] is not None: continue
                    # e_i is in the evaluation image iff evaluation row i
                    # is not in the span of the remaining rows.
                    without=matrix[:i]+matrix[i+1:]
                    if h-rank(without)==1: separators[i]=q
                    separator_tests+=1
            expected=[r//2]*r if not off else [((r-2)+1)//2]*(r-1)+[2]
            require(separators==expected,'branch conductor exponent')
            defect=sum(r-h for h in hs[1:])
            expected_defect=(r-2)**2//4-((r-4)//2 if off else 0)
            require(defect==expected_defect,'seminormal defect formula')
            records.append({'r':r,'off_conic':off,'hilbert_ranks':hs,'defect':defect,'separator_degrees':separators})
    return {'r_range':[6,13],'triple_determinants':triple_count,
            'row_deletion_separator_tests':separator_tests,'records':records}


EDGES=[(2,13),(2,17),(2,19),(2,21),(3,11),(3,16),(7,22),(10,11),(16,17),(18,19),(21,22),
       (2,10),(3,22),(7,18),(10,16),(10,22),(11,13),(13,18),(13,22),(16,21),(17,22),(19,22)]
SMALL={2,11,16,18,22}
VERTICES=sorted({v for e in EDGES for v in e})
INDEX={v:i for i,v in enumerate(VERTICES)}


def syndrome(edges):
    result=0
    for u,v in edges: result^=(1<<INDEX[u])^(1<<INDEX[v])
    return result


def subset_dp(orbits,maxcost=6):
    dp=[Counter() for k in range(maxcost+1)]; dp[0][0]=1
    for orbit in orbits:
        cost=len(orbit); boundary=syndrome(orbit)
        for k in range(maxcost,cost-1,-1):
            for old,count in list(dp[k-cost].items()): dp[k][old^boundary]+=count
    target=syndrome(EDGES)
    return [dp[k][target] for k in range(maxcost+1)]


def laplacian(duplicates):
    L=[[0]*len(VERTICES) for v in VERTICES]
    for u,v in EDGES:
        i,j=INDEX[u],INDEX[v]; w=1+int((u,v) in duplicates)
        L[i][i]+=w;L[j][j]+=w;L[i][j]-=w;L[j][i]-=w
    return L


def mod_rref(a,p):
    a=[[x%p for x in row] for row in a]; piv=[]; row=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(row,len(a)) if a[i][j]),None)
        if pivot is None: continue
        a[row],a[pivot]=a[pivot],a[row]
        inv=pow(a[row][j],-1,p);a[row]=[x*inv%p for x in a[row]]
        for i in range(len(a)):
            if i!=row:
                c=a[i][j];a[i]=[(x-c*y)%p for x,y in zip(a[i],a[row])]
        piv.append(j);row+=1
        if row==len(a):break
    return a,piv


def kernel_mod2(a):
    R,piv=mod_rref(a,2);free=[j for j in range(len(a[0])) if j not in piv]
    basis=[]
    for j in free:
        x=[0]*len(a[0]);x[j]=1
        for i,col in enumerate(piv):x[col]=R[i][j]
        basis.append(x)
    return [[sum(c*v[j] for c,v in zip(coeff,basis))%2 for j in range(len(a[0]))]
            for coeff in product((0,1),repeat=len(basis))]


def local_smith(a):
    determinant=abs(int(sp.Matrix(a).det(method='bareiss')))
    valuations={}
    for p,total in sp.factorint(determinant).items():
        p=int(p);total=int(total)
        nullity=len(a)-len(mod_rref(a,p)[1])
        if nullity==1: vals=[total]
        elif total==nullity:vals=[1]*nullity
        elif p==2:
            K=kernel_mod2(a);rank_mod2=len(a)-nullity;liftable=0
            for x in K:
                Ax=[sum(v*y for v,y in zip(row,x)) for row in a]
                require(all(y%2==0 for y in Ax),'kernel lift divisibility')
                rhs=[-y//2%2 for y in Ax]
                augmented=[row+[b] for row,b in zip(a,rhs)]
                if len(mod_rref(augmented,2)[1])==rank_mod2:liftable+=1
            kernel4=liftable*(2**nullity)
            candidates=[v for v in combinations(range(1,total+1),nullity)]
            # Enumerate nondecreasing positive valuations (repetitions allowed).
            candidates=[v for v in product(range(1,total+1),repeat=nullity)
                        if tuple(sorted(v))==v and sum(v)==total
                        and 2**sum(min(2,e) for e in v)==kernel4]
            require(len(candidates)==1,'local Smith information not decisive')
            vals=list(candidates[0])
        else:raise RuntimeError('unhandled local Smith certificate')
        valuations[p]=vals
    width=max(map(len,valuations.values()))
    invariants=[1]*width
    for p,vals in valuations.items():
        for j,e in enumerate([0]*(width-len(vals))+vals):invariants[j]*=p**e
    return {'determinant':determinant,'prime_valuations':valuations,'invariants':invariants}


def graph_controls():
    swap=lambda v:21 if v==17 else 17 if v==21 else v
    image=lambda edge:tuple(sorted(swap(v) for v in edge))
    require({image(e) for e in EDGES}==set(EDGES),'core involution')
    unseen=set(EDGES);orbits=[]
    while unseen:
        e=min(unseen);orbit={e,image(e)};orbits.append(sorted(orbit));unseen-=orbit
    ordinary=subset_dp([[e] for e in EDGES]);symmetric=subset_dp(orbits)
    require(ordinary[4]==2 and all(n==0 for n in ordinary[:4]),'unconstrained repair DP')
    require(symmetric[6]==27 and all(n==0 for n in symmetric[:6]),'invariant repair DP')
    sets={'original':[], 'four_edges':[(2,17),(3,11),(18,19),(21,22)],
          'six_symmetric_edges':[(2,10),(2,17),(2,21),(3,11),(10,22),(18,19)]}
    expected={'original':[74748],'four_edges':[2,4,49156],'six_symmetric_edges':[4,4,46620]}
    groups={}
    for tag,edges in sets.items():
        L=laplacian(edges);c=[int(v in SMALL) for v in VERTICES]
        harmonic=all(sum(x*y for x,y in zip(row,c))%2==0 for row in L)
        require(harmonic==(tag!='original'),'harmonic character')
        certificate=local_smith([row[:-1] for row in L[:-1]])
        require(certificate['invariants']==expected[tag],'Smith factors')
        groups[tag]=certificate
    return {'parity_syndrome_method':'dynamic programming on edge or automorphism orbits',
            'counts_by_cost_0_to_6':ordinary,'invariant_counts_by_cost_0_to_6':symmetric,
            'involution_orbit_sizes':sorted(map(len,orbits)),
            'smith_method':'Bareiss determinant, finite-field ranks, and mod-4 kernel lifts; no Smith-library call',
            'groups':groups}


def theta_controls():
    records=[]
    for n in range(2,9):
        samples=[([1]*n,True),([3]*n,True),([1]*(n-1)+[2],False),
                 ([(-1)**i*(i+1) for i in range(n)],False)]
        for ratios,trivial in samples:
            # Unknowns: n component constants and two vertex-fibre values.
            constraints=[]
            for i,r in enumerate(ratios):
                p=[0]*(n+2);p[i]=1;p[n]=-1;constraints.append(p)
                q=[0]*(n+2);q[i]=r;q[n+1]=-1;constraints.append(q)
            rk=rank(constraints);h0=n+2-rk;h1=2*n-rk
            require(h0==int(trivial),'theta h0')
            require(h1==(n-1 if trivial else n-2),'theta h1')
            records.append({'n':n,'ratios':ratios,'h0':h0,'h1':h1})
    # Coordinate-free toric character: the last-coordinate functional is
    # fixed by every block shear [[I,B],[0,1]].
    for B in ([[1,0],[0,1]],[[1,0],[0,7]],[[3,-2],[5,4]]):
        for u,v in ((0,1),(1,0),(2,-3)):
            shift=[B[0][0]*u+B[0][1]*v,B[1][0]*u+B[1][1]*v]
            shear=sp.Matrix([[1,0,shift[0]],[0,1,shift[1]],[0,0,1]])
            require(sp.Matrix([[0,0,1]])*shear==sp.Matrix([[0,0,1]]),'toric fibre character')
    return {'endpoint_complex_cases':records,'toric_character_checks':9,
            'boundary':'Descent to T^1(W) requires the stipulated local unramified analytic quotient.'}


if __name__=='__main__':
    print(json.dumps({'moments':moment_controls(),'conductor':conductor_controls(),
                      'critical_graph':graph_controls(),'theta_holonomy':theta_controls(),
                      'result':'PASS'},indent=2,sort_keys=True))
