#!/usr/bin/env python3
"""Exact rational-resolvent/norm bridge checks. Standard library only.
Finite controls supplement the proof in THEORY.md, not the cited global cubic classification.
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


from math import comb,gcd

def F(k,p,q):
    return p**k-sum(p**j*q**(k-j) for j in range(k))

def affine_transform_poly(k,b,a):
    # Ascending coefficients of F_k(bY,b+aY).
    p=[0]*(k+1); p[k]=b**k
    for j in range(k):
        for m in range(k-j+1): p[j+m]-=comb(k-j,m)*b**(k-m)*a**m
    return p

def evaluate_matrix(p,M):
    k=len(M); out=[[0]*k for _ in range(k)]
    for j,x in enumerate(p): out=plus(out,scale(power(M,j),x))
    return out

def integral_matrix(M):
    return all(Fraction(x).denominator==1 for row in M for x in row)

def clean(x):
    if isinstance(x,Fraction):
        return int(x) if x.denominator==1 else f'{x.numerator}/{x.denominator}'
    if isinstance(x,list): return [clean(y) for y in x]
    if isinstance(x,dict): return {k:clean(v) for k,v in x.items()}
    return x

results=[]
controls=[]
for k in (3,4,5):
    C=companion(k); I=eye(k); e0=I[0]; alpha=I[1]
    gains=[(a,b) for b in range(1,7) for a in range(-6,7) if gcd(a,b)==1]
    if k==3: gains.append((56,103))
    for a,b in gains:
        norm=F(k,b,a)
        check(norm!=0,'rational gain not a pole')
        check(gcd(abs(norm),b)==1,'coprime denominator norm')
        D=plus(scale(I,b),scale(C,-a))
        R=scale(inverse(D),b); T=matmul(R,C)
        check(determinant(D)==norm,'denominator matrix norm')
        check(matmul(D,R)==scale(I,b),'exact resolvent equation')
        check(matmul(plus(I,scale(C,Fraction(-a,b))),T)==C,'implicit update equation')
        check(matmul(T,C)==matmul(C,T),'field-linear resolvent')
        check(determinant(T)==Fraction((-1)**(k+1)*b**k,norm),'transformed eigenvalue norm')
        u=Fraction(b*a**(k-1),norm)
        check(T[0][0]==u,'generating function is constant coefficient of feedback generator')
        check(u.denominator==1 if abs(norm)==1 else u.denominator!=1,
              'GF scalar integrality exactly detects full operator integrality')
        check(integral_matrix(T)==(abs(norm)==1),'Cgamma integral iff norm unit')
        check(integral_matrix(R)==(abs(norm)==1),'Rgamma integral iff norm unit')
        raw=affine_transform_poly(k,b,a)
        check(raw[-1]==norm and raw[0]==-b**k,'minpoly endpoint coefficients')
        monic=[Fraction(x,norm) for x in raw]
        check(evaluate_matrix(monic,T)==[[0]*k for _ in range(k)],'Möbius minimal polynomial annihilates')
        check(all(x.denominator==1 for x in monic)==(abs(norm)==1),'monic polynomial integrality')
        check(matmul(scale(T,b),inverse(plus(scale(I,b),scale(T,a))))==C,'inverse Möbius operator')
        powers=[matvec(power(T,j),e0) for j in range(k)]
        basis=[list(row) for row in zip(*powers)]
        basis_det=determinant(basis)
        check(basis_det==Fraction(b**(k*(k-1)),norm**(k-1)),'exact power-order basis determinant')
        check(basis_det!=0,'unchanged degree cyclicity')
        if k==3 and (a,b) in ((0,1),(1,2),(56,103)):
            coeff=matvec(T,e0)
            expected=alpha if a==0 else matvec(scale(power(C,4),2),e0) if a==1 else matvec(scale(power(C,17),-103),e0)
            check(coeff==expected,'closed power form of three integral cubic gains')
            check(abs(basis_det)==b**6,'cubic order index')
            expected_poly=([-1,-1,-1,1] if a==0 else [-8,-20,-22,1] if a==1
                           else [1092727,2875039,3249959,1])
            check(monic==expected_poly,'published cubic minimal polynomial row')
            check(abs(determinant(T))==b**3,'integral but nonunimodular transformed matrix')
            results.append({'gain':[a,b],'F':norm,'beta_gamma_power':
                            'beta' if a==0 else '2 beta^4' if a==1 else '-103 beta^17',
                            'coefficient_vector':coeff,'C_gamma':T,'R_gamma':R,
                            'minimal_polynomial_ascending':monic,'norm':determinant(T),
                            'order_index':abs(basis_det)})
    controls.append({'k':k,'primitive_gain_count':len(gains),'b_range':[1,6],'a_range':[-6,6],
                     'extra_gains':[[56,103]] if k==3 else []})
# Global conclusion uses the cited classification, independently check its exact witnesses/filter.
cubic_norm_one=[(0,-1),(1,0),(2,1),(-103,-56)]
for b,a in cubic_norm_one: check(F(3,b,a)==1,'cited complete cubic norm-one witness')
pm=[(b,a) for b,a in cubic_norm_one]+[(-b,-a) for b,a in cubic_norm_one]
filtered=sorted((Fraction(a,b) for b,a in pm if b>0))
check(filtered==[Fraction(0),Fraction(1,2),Fraction(56,103)],'global cited list gain filter')
# Exact Pisot magnitude bounds, no floating roots.
check(F(3,9,5)<0 and F(3,2,1)>0,'9/5 < beta < 2 sign certificate')
check(Fraction(50,81)<1,'half-gain internal modulus upper bound')
check(9**17>103**2*5**17,'exceptional internal modulus squared upper bound')
check(matvec(power(companion(3),-16),[1,0,0])==[-103,56,0],'exceptional gain beyond real pole')
# Explicit finite hostile to universal rational-gain integrality.
check(F(3,1,1)==-2,'gamma one norm obstruction')
print(json.dumps(clean({'status':'PASS','checks':checks,'finite_universe':controls,
                        'cubic_integral_gains':results,
                        'cited_classification_gain_filter':filtered,
                        'pisot_bounds':{'beta_lower':[9,5],'exceptional_square_upper':
                                        Fraction(103**2*5**17,9**17),
                                        'half_gain_modulus_upper':Fraction(50,81)},
                        'scope':'Global cubic gain completeness uses cited Delone-Nagell classification; k4,k5 gain completeness not asserted.'}),indent=2,sort_keys=True))
