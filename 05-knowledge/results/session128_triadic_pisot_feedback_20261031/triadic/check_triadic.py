#!/usr/bin/env python3
"""Portable exact checks for the fixed-degree field/companion/recurrence triad.
No external modules, repository imports, or assertion-dependent checks.
The finite checks corroborate the general proof in THEORY.md.
"""
from fractions import Fraction
from itertools import product
import hashlib
import json

checks = 0

def check(value, label):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(label)

def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]

def matmul(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]

def matvec(a, v):
    return [sum(x*y for x, y in zip(row, v)) for row in a]

def plus(a, b):
    return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]

def scale(a, c):
    return [[c*x for x in row] for row in a]

def determinant(a):
    a = [list(map(Fraction,row)) for row in a]
    sign = 1
    ans = Fraction(1)
    for j in range(len(a)):
        p = next((i for i in range(j,len(a)) if a[i][j]),None)
        if p is None: return 0
        if p != j: a[j],a[p]=a[p],a[j]; sign = -sign
        pivot=a[j][j]; ans*=pivot
        for i in range(j+1,len(a)):
            q=a[i][j]/pivot
            for m in range(j,len(a)): a[i][m]-=q*a[j][m]
    ans *= sign
    return int(ans) if ans.denominator == 1 else ans

def inverse(a):
    n=len(a); a=[list(map(Fraction,row))+list(map(Fraction,er)) for row,er in zip(a,eye(n))]
    for j in range(n):
        p=next(i for i in range(j,n) if a[i][j]); a[j],a[p]=a[p],a[j]
        pivot=a[j][j]; a[j]=[x/pivot for x in a[j]]
        for i in range(n):
            if i != j:
                q=a[i][j]; a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    out=[row[n:] for row in a]
    return [[int(x) if x.denominator == 1 else x for x in row] for row in out]

def power(a,n):
    if n<0: return power(inverse(a),-n)
    out=eye(len(a))
    while n:
        if n%2: out=matmul(out,a)
        a=matmul(a,a); n//=2
    return out

def companion(k):
    a=[[0]*k for _ in range(k)]
    for j in range(k-1): a[j+1][j]=1
    for i in range(k): a[i][-1]=1
    return a

def polymul(a,b):
    k=len(a); c=[0]*(2*k-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    for m in range(2*k-2,k-1,-1):
        for j in range(k): c[m-k+j]+=c[m]
    return c[:k]

def coefficient_inverse(k):
    b=[[0]*k for _ in range(k)]
    b[0][0]=1
    for j in range(1,k):
        b[k-j][j]=1
        for i in range(1,j): b[k-j][i]=-1
    return b

def stream_forward(initial,count):
    s=list(initial); k=len(initial)
    while len(s)<count: s.append(sum(s[-k:]))
    return s

def constant_stream(c,count):
    x=list(c); k=len(c); alpha=[0,1]+[0]*(k-2); s=[]
    for _ in range(count): s.append(x[0]); x=polymul(x,alpha)
    return s

def multiplication_matrix(c):
    k=len(c); e=eye(k)
    return [list(row) for row in zip(*(polymul(c,v) for v in e))]

def rank(a):
    a=[list(map(Fraction,row)) for row in a]; row=0
    for col in range(len(a[0])):
        p=next((i for i in range(row,len(a)) if a[i][col]),None)
        if p is None: continue
        a[row],a[p]=a[p],a[row]; q=a[row][col]; a[row]=[x/q for x in a[row]]
        for i in range(len(a)):
            if i != row:
                q=a[i][col]; a[i]=[x-q*y for x,y in zip(a[i],a[row])]
        row+=1
        if row==len(a): break
    return row

results=[]
for k in (3,4,5):
    C=companion(k); I=eye(k)
    O=[power(C,n)[0] for n in range(k)]
    B=coefficient_inverse(k)
    R=[[0]*k for _ in range(k)]
    for i in range(k-1): R[i][i+1]=1
    R[-1]=[1]*k
    J=[list(reversed(row)) for row in I]
    A=matmul(matmul(J,R),J)
    U=[[0]*k for _ in range(k)]
    for j in range(k):
        U[j][j]=1
        for i in range(j): U[i][j]=-1
    Bdescending=matmul(J,O)
    check(matmul(U,A)==matmul(C,U),'substitution-to-power-basis intertwiner')
    check(matmul(U,Bdescending)==inverse(C),'convention difference is one inverse time step')
    check(inverse(U)==matmul(Bdescending,C),'substitution coordinates as forward stream state')
    check(matmul(O,B)==I and matmul(B,O)==I,'observability inverse')
    check(abs(determinant(O))==1,'unimodular observation')
    check(matmul(O,C)==matmul(R,O),'shift naturality')
    check(matmul(matmul(J,O),C)==matmul(A,matmul(J,O)),'reversed-state paper convention')
    sum_powers=[[0]*k for _ in range(k)]
    for n in range(k): sum_powers=plus(sum_powers,power(C,n))
    check(power(C,k)==sum_powers,'annihilating polynomial')
    # Centralizer constraints B*C-C*B=0; nullity k is an independent full-hom check.
    eq=[]
    for i in range(k):
        for j in range(k):
            row=[0]*(k*k)
            for m in range(k): row[i*k+m]+=C[m][j]; row[m*k+j]-=C[i][m]
            eq.append(row)
    nullity=k*k-rank(eq)
    check(nullity==k,'centralizer dimension')
    coeffs=list(product(range(-1,2),repeat=k))
    for c in coeffs:
        v=list(c); w=matvec(O,v)
        check(matvec(B,w)==v,'coefficient reconstruction')
        check((all(x==0 for x in v[2:]))==(all(x==0 for x in w[1:k-1])),'sparse gate')
        if not any(v[2:]):
            check(matvec(inverse(U),v)==[v[0]+v[1],v[1]]+[0]*(k-2),
                  'sparse shear is substitution basis coordinate')
        check(constant_stream(v,3*k)==stream_forward(w,3*k),'independent polynomial recurrence stream')
        M=multiplication_matrix(v)
        check(matmul(M,C)==matmul(C,M),'multiplication intertwiner')
        P=[[0]*k for _ in range(k)]
        for j,x in enumerate(v): P=plus(P,scale(power(C,j),x))
        check(M==P,'full hom inverse first-column')
        H=matmul(matmul(O,M),B)
        check(matmul(H,R)==matmul(R,H),'recurrence hom shift commute')
        check(determinant(H)==determinant(M),'norm survives coordinate conjugacy')
    # Composition/naturality checks on k+1 positive and hostile morphism generators.
    morphs=[multiplication_matrix(v) for v in I]+[multiplication_matrix([-1]*k)]
    for M in morphs:
        for N in morphs:
            HM=matmul(matmul(O,M),B); HN=matmul(matmul(O,N),B)
            check(matmul(HM,HN)==matmul(matmul(O,matmul(M,N)),B),'functor composition')
            check(matmul(HM,O)==matmul(O,M),'naturality square')
    # Rank-two K-linear determinant / rational-realization identity.
    # Distinct from checking only scalar multiplication norms.
    norm_det_tests=0
    entries=[I[0],I[1],[-1]*k,[2]+[0]*(k-1),[0]*k]
    for a,b,c,d in product(entries,repeat=4):
        ad=polymul(a,d); bc=polymul(b,c)
        delta=[x-y for x,y in zip(ad,bc)]
        ma,mb,mc,md=map(multiplication_matrix,(a,b,c,d))
        block=[ra+rb for ra,rb in zip(ma,mb)]+[rc+rd for rc,rd in zip(mc,md)]
        check(determinant(block)==determinant(multiplication_matrix(delta)),
              'rank-two determinant equals field norm of determinant')
        norm_det_tests+=1
    D=[[0]*k for _ in range(k)]; D[0][0]=1
    check(matmul(D,C)!=matmul(C,D),'non-K-linear feedback hostile')
    two=scale(I,2)
    check(determinant(plus(I,scale(two,-1)))!=0,'large feedback exact invertibility')
    e0=I[0]
    source=[]
    for n in range(25):
        v=matvec(power(C,-n),e0); w=matvec(O,v)
        check(matvec(power(C,n),v)==e0,'unit-source orbit inverse')
        if not any(v[2:]):
            source.append({'n':n,'coefficients':v,'stream_initial':w,'paper_state':list(reversed(w)),'substitution_state':matvec(inverse(U),v)})
    hostile=matvec(power(C,-1),e0)
    check(any(hostile[2:]),'inverse unit not sparse hostile')
    norm_hostile=matvec(power(C,-2),e0)
    check(determinant(multiplication_matrix(norm_hostile))==1,'same norm one hostile')
    check(norm_hostile != e0,'norm-one distinct source')
    check(any(norm_hostile[2:]),'norm-one not sparse hostile')
    scalar_collision=I[1]
    check(scalar_collision[0]==0 and scalar_collision!=[0]*k,'one recurrence value not invertible')
    results.append({'k':k,'C':C,'O':O,'O_inverse':B,'R_ascending':R,'A_paper_descending':A,
                    'U_substitution_to_power':U,'det_O':determinant(O),'centralizer_dimension':nullity,
                    'coefficient_cube_size':len(coeffs),'rank_two_norm_det_tests':norm_det_tests,'sparse_source_n_0_to_24':source,
                    'norm_one_n2_hostile_coefficients':norm_hostile})
check(results[0]['A_paper_descending']==[[1,1,1],[1,0,0],[0,1,0]],'paper cubic matrix')
check(next(x for x in results[0]['sparse_source_n_0_to_24'] if x['n']==16)['stream_initial']==[-103,0,56],'cubic exceptional stream')
check(next(x for x in results[0]['sparse_source_n_0_to_24'] if x['n']==16)['substitution_state']==[-47,56,0],'cubic exceptional substitution state')
payload={'status':'FINITE-EXACT corroboration of general algebraic proof','orders':[3,4,5],
         'checks':checks,'results':results,
         'scope':'No symbolic substitution/attractor equivalence, all-unit classification, or cross-degree field embedding is asserted.'}
print(json.dumps(payload,indent=2,sort_keys=True))
