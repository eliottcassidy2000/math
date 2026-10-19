#!/usr/bin/env python3
"""Exact controls, with literal deletion multiplicities; no third-party packages."""
from itertools import combinations, permutations
from functools import lru_cache
from math import comb, factorial, gcd
from collections import Counter
import json

checks = Counter()
def ck(value, tag, detail=None):
    checks[tag] += 1
    if not value:
        raise RuntimeError((tag, detail))

@lru_cache(None)
def pairs(n): return tuple(combinations(range(n), 2))
@lru_cache(None)
def index(n): return {ij:k for k,ij in enumerate(pairs(n))}
def arc(n, mask, i, j):
    return ((mask >> index(n)[(i,j)]) & 1) if i < j else 1-((mask >> index(n)[(j,i)]) & 1)
def induced(n, mask, keep):
    return sum(arc(n,mask,keep[i],keep[j]) << k for k,(i,j) in enumerate(pairs(len(keep))))
def sources(n, mask, count):
    # Sources have successive final labels. Every new source beats all earlier vertices.
    return sum(((mask >> k)&1) << index(n+count)[ij] for k,ij in enumerate(pairs(n)))

classes=[]; canonical=[]; automorphisms=[]
for n in range(6):
    mapping={}; reps=[]; aut=[]
    for mask in range(1 << (n*(n-1)//2)):
        if mask in mapping: continue
        orbit=Counter(induced(n,mask,order) for order in permutations(range(n)))
        rep=min(orbit)
        reps.append(rep); aut.append(orbit[rep])
        for other in orbit: mapping[other]=rep
    classes.append(reps); canonical.append(mapping); automorphisms.append(aut)
class_index=[{T:i for i,T in enumerate(cs)} for cs in classes]

def indicator_row(n,T):
    row=[0]*len(classes[n]);row[class_index[n][canonical[n][T]]]=1;return row

def deletion_row(n, m, T):
    row=[0]*len(classes[n])
    for keep in combinations(range(n+m),n):
        sub=induced(n+m,T,keep)
        row[class_index[n][canonical[n][sub]]]+=1
    return row

@lru_cache(None)
def D(n,m):
    if n < 0: return tuple(() for _ in classes[n+m])
    return tuple(tuple(deletion_row(n,m,T)) for T in classes[n+m])

def multiply(A,B):
    if not A: return []
    width=len(B[0]) if B else 0
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A] if B else [[0]*width for row in A]
def times(a,A,p=None):
    return [[(a*x)%p if p else a*x for x in row] for row in A]
def matvec(A,v,p):return [sum(a*b for a,b in zip(row,v))%p for row in A]

def rank_and_kernel(A,width,p):
    M=[[x%p for x in row] for row in A]; piv=[];r=0
    for col in range(width):
        at=next((i for i in range(r,len(M)) if M[i][col]),None)
        if at is None:continue
        M[r],M[at]=M[at],M[r]
        inv=pow(M[r][col],-1,p);M[r]=[(inv*x)%p for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][col]:
                a=M[i][col];M[i]=[(x-a*y)%p for x,y in zip(M[i],M[r])]
        piv.append(col);r+=1
    basis=[]
    for col in range(width):
        if col in piv:continue
        v=[0]*width;v[col]=1
        for i,pc in enumerate(piv):v[pc]=(-M[i][col])%p
        basis.append(v)
    return r,basis

@lru_cache(None)
def composite_row(n, out_n, T, outer, source_count, inner):
    """Row of D_outer h^source_count D_inner on C_n, evaluated at T of out_n."""
    if out_n-outer < 0:return tuple([0]*len(classes[n]))
    mid_n=out_n-outer+source_count
    ck(mid_n-inner==n,'typed_composition')
    row=[0]*len(classes[n])
    for keep in combinations(range(out_n),out_n-outer):
        U=induced(out_n,T,keep)
        V=sources(out_n-outer,U,source_count)
        for remain in combinations(range(mid_n),n):
            sub=induced(mid_n,V,remain)
            row[class_index[n][canonical[n][sub]]]+=1
    return tuple(row)

# Integral divided powers, using independent matrix multiplication.
for n in range(6):
    for a in range(6-n):
        for b in range(6-n-a):
            left=multiply(D(n+b,a),D(n,b))
            right=times(comb(a+b,a),D(n,a+b))
            ck(left==right,'integral_divided_power',(n,a,b))
# Integral source/subset relation, including degree zero.
for n in range(6):
    for m in range(1,7-n):
        out=n+m-1
        if out>5:continue
        for T in classes[out]:
            lhs=list(composite_row(n,out,T,0,1,m))
            rhs=list(composite_row(n,out,T,m,1,0)) if n else [0]*len(classes[n])
            extra=deletion_row(n,m-1,T)
            ck(lhs==[a+b for a,b in zip(rhs,extra)],'integral_source_relation',(n,m,T))

pairs_to_test=[(2,0),(3,0),(5,0),(7,0),(2,1),(2,2),(3,1)]
contraction_summary=[]
for p,e in pairs_to_test:
    s=p**e;terms=0
    for n in range(6):
        for T in classes[n]:
            row=[0]*len(classes[n])
            for j in range(p):
                a=p-1-j
                r=composite_row(n,n,T,a*s,(p-1)*s,j*s)
                c=-factorial(a)*factorial(j)
                row=[(x+c*y)%p for x,y in zip(row,r)]
                terms+=1
            ck(row==indicator_row(n,T),'full_digit_contraction',(p,e,n,T))
    contraction_summary.append({'prime':p,'digit':e,'step':s,'base_degrees':list(range(6)), 'terms':terms})

rank_summary=[]
for p,e in pairs_to_test:
    s=p**e
    for k in range(1,p):
        for n in range(6-k*s):
            A=times(factorial(k),D(n,k*s),p)
            rank,ker=rank_and_kernel(A,len(classes[n]),p)
            low=n-(p-k)*s
            if low<0:
                irank=0;P=[]
            else:
                B=times(factorial(p-k),D(low,(p-k)*s),p)
                irank,_=rank_and_kernel(B,len(classes[low]),p)
                P=[]
                for T in classes[low]:
                    row=[0]*len(classes[n])
                    for j in range(k):
                        a=k-1-j
                        r=composite_row(n,low,T,a*s,(p-1)*s,j*s)
                        c=-factorial(a)*factorial(j)
                        row=[(x+c*y)%p for x,y in zip(row,r)]
                    P.append(row)
                for f in ker:
                    g=matvec(P,f,p)
                    ck(matvec(B,g,p)==f,'explicit_kernel_primitive',(p,e,k,n,f))
            ck(len(classes[n])-rank==irank,'kernel_image_rank',(p,e,k,n))
            rank_summary.append({'p':p,'e':e,'k':k,'degree':n,'dimension':len(classes[n]),'kernel':len(classes[n])-rank,'image':irank})

# Scalar carry transducer and independent Frobenius coefficient law.
for p in (2,3,5,7):
    for a in range(65):
        for b in range(65):
            x,y=a,b;local=1;carry=0;ever=False
            while x or y or carry:
                aa,bb=x%p,y%p
                local=local*(comb(aa+bb,aa) if aa+bb<p else 0)%p
                carry=(aa+bb+carry)//p;ever=ever or bool(carry)
                x//=p;y//=p
            c=comb(a+b,a)%p
            ck(c==local,'lucas_product',(p,a,b))
            ck((c!=0)==(not ever),'carry_support',(p,a,b))
    for e in range(4):
        s=p**e
        scalar=factorial(p*s)//(factorial(s)**p)
        # Avoid gigantic factorials past s=343? This largest is 2401!, still exact and modest.
        ck(scalar%p==0 and scalar%(p*p)!=0,'digit_nilpotence_valuation',(p,e))

# Explicit smallest odd-prime invisible host and labelled-isolate erasure.
cycle=next(T for T,aut in zip(classes[3],automorphisms[3]) if aut==3)
# Enumerate all oriented graphs of order 3, absent/forward/backward on each pair.
embedding_counts=[]
for ternary in range(27):
    q=ternary;edges=[]
    for i,j in pairs(3):
        digit=q%3;q//=3
        if digit:edges.append((i,j) if digit==1 else (j,i))
    count=sum(all(arc(3,cycle,order[i],order[j]) for i,j in edges) for order in permutations(range(3)))
    embedding_counts.append(count)
    ck(count%3==0,'triangle_embedding_blindness',ternary)
ck(deletion_row(0,3,cycle)==[1],'hasse_escape')
ck(factorial(3)%3==0,'labelled_isolate_erasure')
for n,auts in enumerate(automorphisms):
    ck(all(a%2 for a in auts),'odd_tournament_automorphisms',n)

# The integral spanning-embedding lattice is weighted by host automorphisms.
embedding_deletions=[]
for n in range(5):
    A=[]
    for i,row in enumerate(D(n,1)):
        out=[]
        for j,value in enumerate(row):
            numerator=value*automorphisms[n][j]
            ck(numerator%automorphisms[n+1][i]==0,'embedding_lattice_integrality',(n,i,j))
            out.append(numerator//automorphisms[n+1][i])
        A.append(out)
    embedding_deletions.append(A)
E3=multiply(embedding_deletions[2],multiply(embedding_deletions[1],embedding_deletions[0]))
cyclic_index=class_index[3][cycle]
ck(E3[cyclic_index][0]%3==2,'unsaturated_abstract_reduction_hostile')
ck(any(row[0]%3 for row in E3),'abstract_embedding_not_p_nilpotent')
ck(all(x%3==0 for row in multiply(D(2,1),multiply(D(1,1),D(0,1))) for x in row),
   'pointwise_reduction_p_nilpotent')

# Integral source-free induced-count basis, explicit inverse, and Smith controls.
def has_source(n,T):
    return any(all(v==u or arc(n,T,u,v) for v in range(n)) for u in range(n))
cores=[[(n,T) for T in classes[n] if not has_source(n,T)] for n in range(6)]
core_labels=[];bases=[];inverses=[]
for n in range(6):
    labels=[entry for rr in range(n+1) for entry in cores[rr]]
    core_labels.append(labels)
    B=[]
    for T in classes[n]:
        row=[]
        for r,C in labels:
            row.append(sum(canonical[r][induced(n,T,U)]==C for U in combinations(range(n),r)))
        B.append(row)
    bases.append(B)
    I=[]
    for r,C in labels:
        row=[0]*len(classes[n])
        for size in range(r+1):
            for keep in combinations(range(r),size):
                sub=induced(r,C,keep)
                aug=sources(size,sub,n-size)
                row[class_index[n][canonical[n][aug]]]+=(-1)**(r-size)
        I.append(row)
    inverses.append(I)
    ident=[[int(i==j) for j in range(len(classes[n]))] for i in range(len(classes[n]))]
    ck(len(labels)==len(classes[n]),'source_free_basis_size',n)
    ck(multiply(I,B)==ident,'integral_basis_left_inverse',n)
    ck(multiply(B,I)==ident,'integral_basis_right_inverse',n)

def det(M):
    A=[row[:] for row in M];n=len(A)
    if n==0:return 1
    previous=1;sign=1
    for k in range(n-1):
        at=next((i for i in range(k,n) if A[i][k]),None)
        if at is None:return 0
        if at!=k:A[k],A[at]=A[at],A[k];sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=A[i][j]*pivot-A[i][k]*A[k][j]
                ck(numerator%previous==0,'bareiss_exact_division')
                A[i][j]=numerator//previous
            A[i][k]=0
        previous=pivot
    return sign*A[-1][-1]

for n,B in enumerate(bases):ck(abs(det(B))==1,'unimodular_basis_determinant',n)
smith_summary=[]
for n in range(5):
    for m in range(1,6-n):
        A=[list(row) for row in D(n,m)]
        conjugated=multiply(inverses[n+m],multiply(A,bases[n]))
        diag=[comb(n-r+m,m) for r,C in core_labels[n]]
        expected=[[0]*len(classes[n]) for _ in classes[n+m]]
        for j,label in enumerate(core_labels[n]):
            i=core_labels[n+m].index(label);expected[i][j]=diag[j]
        ck(conjugated==expected,'integral_diagonalization',(n,m))
        # Independent determinantal divisors, without using the basis inverse.
        divs=[1];pred=[1]
        for k in range(1,len(classes[n])+1):
            d=0
            for R in combinations(range(len(A)),k):
                for C in combinations(range(len(A[0])),k):
                    d=gcd(d,abs(det([[A[i][j] for j in C] for i in R])))
            pd=0
            for Q in combinations(diag,k):
                z=1
                for q in Q:z*=q
                pd=gcd(pd,z)
            ck(d==pd,'smith_determinantal_divisor',(n,m,k,d,pd))
            divs.append(d);pred.append(pd)
        factors=[divs[j]//divs[j-1] for j in range(1,len(divs))]
        smith_summary.append({'source_degree':n,'deletion_size':m,'diagonal_entries':diag,
            'smith_invariant_factors':factors,'free_cokernel_rank':len(classes[n+m])-len(classes[n])})

print(json.dumps({'status':'PASS','class_counts':[len(x) for x in classes],
  'checks_by_family':dict(sorted(checks.items())),'total_checks':sum(checks.values()),
  'contractions':contraction_summary,'rank_and_kernel_cases':rank_summary,
  'integral_source_free_counts':[len(x) for x in cores],'integral_smith_cases':smith_summary,
  'abstract_embedding_d_cubed':E3,'triangle_host_mask':cycle,'triangle_embedding_counts':embedding_counts,
  'scope':'All class-basis functions at degrees 0..5 for seven prime/digit pairs; integral matrix identities through degree five; bounded ranks and explicit primitives when target degree <=5; all 27 oriented order-three source graphs at cyclic host; scalar Lucas/carry atlas 0<=a,b<=64 and p=2,3,5,7; integral source-free induced-count basis and two-sided inverse through degree five; all positive subset-deletion Smith matrices with target degree at most five, independently checked through every determinantal divisor. Universal theorems are proved in PADDING_FROBENIUS_PROOF.md.'},indent=2,sort_keys=True))
