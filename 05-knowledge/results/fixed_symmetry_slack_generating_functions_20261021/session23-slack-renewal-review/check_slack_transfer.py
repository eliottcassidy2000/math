#!/usr/bin/env python3
"""Exact finite controls, not a proof of asymptotic laws or a new tournament census."""
import hashlib, itertools, json, math
from collections import Counter
from fractions import Fraction
from pathlib import Path

CHECKS=Counter()
def ck(value, label):
    CHECKS[label]+=1
    if not value: raise RuntimeError(label)

def add(*ps):
    out=Counter()
    for p in ps:
        for key,v in p.items(): out[key]+=v
    return {key:v for key,v in out.items() if v}
def scale(p,s): return {key:v*s for key,v in p.items() if v*s}
def mul(p,q,N=10000):
    out=Counter()
    for (n,k),v in p.items():
        for (m,l),w in q.items():
            if n+m<=N:out[n+m,k+l]+=v*w
    return {key:v for key,v in out.items() if v}
def power(p,e,N=10000):
    out={(0,0):1}
    for _ in range(e):out=mul(out,p,N)
    return out
def monomial(a,b):return {(a,b):1}
def subst(p):return {(2*n,2*k):v for (n,k),v in p.items()}
def word_poly(N):
    return {(n,k):math.comb(n-2*k,k) for n in range(N+1) for k in range(n//3+1)}
def palindrome_poly(N):
    out={}
    for h,e in [(0,0),(1,0),(3,1)]:
        if N>=h:out=add(out,mul(monomial(h,e),subst(word_poly((N-h)//2)),N))
    return out
p={(0,0):1,(1,0):-1,(3,1):-1}
q=subst(p)
center={(0,0):1,(1,0):1,(3,1):1}

def rational_coeffs(num,den,N):
    # Independent formal division in vertex degree, with a polynomial in feedback.
    ck(den.get((0,0))==1 and all(n or k==0 for n,k in den),'unit_constant')
    by_n=[{} for _ in range(N+1)]
    for (n,k),v in num.items():
        if n<=N:by_n[n][k]=by_n[n].get(k,0)+v
    for n in range(N+1):
        for (j,l),v in den.items():
            if 0<j<=n:
                for k,w in by_n[n-j].items():by_n[n][k+l]=by_n[n].get(k+l,0)-v*w
        by_n[n]={k:v for k,v in by_n[n].items() if v}
    return {(n,k):v for n,row in enumerate(by_n) for k,v in row.items()}

def rational_moments(num,den,N):
    # Differentiate QF=N twice with the Euler derivative y*d/dy.
    def layers(poly):
        out=[Counter(),Counter(),Counter()]
        for (n,k),v in poly.items():
            for r in range(3):out[r][n]+=v*(k**r)
        return out
    ns,ds=layers(num),layers(den)
    out=[[0]*(N+1) for _ in range(3)]
    deg=max(n for n,k in den)
    for n in range(N+1):
        out[0][n]=ns[0][n]-sum(ds[0][j]*out[0][n-j] for j in range(1,min(n,deg)+1))
        out[1][n]=ns[1][n]-sum(ds[0][j]*out[1][n-j]+ds[1][j]*out[0][n-j] for j in range(1,min(n,deg)+1))
        out[2][n]=ns[2][n]-sum(ds[0][j]*out[2][n-j]+2*ds[1][j]*out[1][n-j]+ds[2][j]*out[0][n-j] for j in range(1,min(n,deg)+1))
    return out

def compiler(configs,d,N,self_converse=False):
    # Arbitrary synthetic finite configurations. These are NOT actual slack inventories.
    M=max(((ell+1)//2 if self_converse else ell) for a,b,ell in configs)
    den=power(q if self_converse else p,M)
    D=word_poly(d)
    tail_num=add({(0,0):1},scale(mul(D,p),-1))
    pair_num=subst(tail_num)
    pal_num=add(center,scale(mul(palindrome_poly(d),q),-1))
    num={};direct={}
    L={key:v for key,v in word_poly(N).items() if key[0]>d}
    LP=subst({key:v for key,v in word_poly(N//2).items() if key[0]>d})
    PP={key:v for key,v in palindrome_poly(N).items() if key[0]>d}
    for a,b,ell in configs:
        if not self_converse:
            np=mul(power(tail_num,ell),power(p,M-ell))
            dp=power(L,ell,N)
        else:
            pairs,mid=divmod(ell,2)
            np=mul(mul(power(pair_num,pairs),power(pal_num,mid)),power(q,M-pairs-mid))
            dp=mul(power(LP,pairs,N),power(PP,mid,N),N)
        num=add(num,mul(monomial(a,b),np))
        direct=add(direct,mul(monomial(a,b),dp,N))
    got=rational_coeffs(num,den,N)
    for n in range(N+1):
        ck({k:v for (m,k),v in got.items() if m==n}=={k:v for (m,k),v in direct.items() if m==n},'synthetic_exact_bivariate_row')
    return num,den

# Literal palindrome enumeration independently tests optional center and finite-head removal.
def literal_words(N):
    out=[()]
    for n in range(1,N+1):
        out.extend((1,)+w for w in list(out) if sum(w)==n-1)
        if n>=3:out.extend((3,)+w for w in list(out) if sum(w)==n-3)
    return out
words=literal_words(16)
for d in [0,1,4,7]:
    got=Counter((sum(w),w.count(3)) for w in words if sum(w)>d and w==w[::-1])
    expected={key:v for key,v in palindrome_poly(16).items() if key[0]>d}
    ck(dict(got)==expected,'literal_palindrome_tail')
    # Explicit parity extraction with different half cutoffs for the three centers.
    reconstructed={}
    for h,e in [(0,0),(1,0),(3,1)]:
        cutoff=(d-h)//2
        half={key:v for key,v in word_poly((16-h)//2).items() if key[0]>cutoff}
        reconstructed=add(reconstructed,mul(monomial(h,e),subst(half),16))
    ck(reconstructed==expected,'center_specific_half_cutoff')

cases={}
cases['synthetic_ordinary']=compiler([(0,0,1),(2,1,2),(5,3,4),(7,0,0),(1,6,4)],4,48)
cases['synthetic_sc_parity_unequal_pole_orders']=compiler([(0,0,2),(1,2,4),(2,1,1),(5,2,0)],4,48,True)
cases['synthetic_sc_even_only_unbounded']=compiler([(0,0,4),(5,2,0)],4,48,True)

# Actual inherited low-slack formulas; the theorem's full compiler is not enumerated here.
actual={
 'delta0':({(0,0):1},p),
 'delta0_sc':(center,q),
 'delta3':(mul(monomial(2,1),add({(0,0):1},scale(mul({(0,0):1,(1,0):1},p),-1))),power(p,3)),
 'delta3_sc':(mul(monomial(2,1),add(center,scale(mul({(0,0):1,(1,0):1},q),-1))),power(q,2)),
 'delta4':(monomial(5,3),power(p,2)),
 'delta4_sc':(monomial(5,3),q),
}
N=48;A=word_poly(N);P=palindrome_poly(N)
actual_direct={
 'delta0':A,
 'delta0_sc':P,
 'delta3':mul(monomial(2,1),mul({key:v for key,v in A.items() if key[0]>=2},power(A,2,N),N),N),
 'delta3_sc':mul(monomial(2,1),mul({key:v for key,v in P.items() if key[0]>=2},subst(word_poly(N//2)),N),N),
 'delta4':mul(monomial(5,3),power(A,2,N),N),
 'delta4_sc':mul(monomial(5,3),subst(word_poly(N//2)),N),
}
low_rows={}
for name,(num,den) in actual.items():
    got=rational_coeffs(num,den,N)
    for n in range(N+1):
        ck({k:v for (m,k),v in got.items() if m==n}=={k:v for (m,k),v in actual_direct[name].items() if m==n},'actual_low_slack_exact_bivariate_row')
    low_rows[name]=[sum(v for (m,k),v in got.items() if m==n) for n in range(8)]
expected_small={
 'delta0':[1,1,1,2,3,4,6,9], 'delta0_sc':[1,1,1,2,1,2,2,3],
 'delta3':[0,0,0,0,1,4,10,22], 'delta3_sc':[0,0,0,0,1,2,2,4],
 'delta4':[0,0,0,0,0,1,2,3], 'delta4_sc':[0,0,0,0,0,1,0,1],
}
for name in actual:ck(low_rows[name]==expected_small[name],'actual_known_census_through_seven')
cases.update(actual)

rho=.6823278038280193
mu=rho*rho/(1+3*rho*rho)
sigma2=rho*rho*(1+rho*rho)/(1+3*rho*rho)**3
moment_rows=[]
for name,(num,den) in cases.items():
    mom=rational_moments(num,den,1081)
    # Small exact moment identities independently match full bivariate division.
    small=rational_coeffs(num,den,48)
    for n in range(49):
        for r in range(3):ck(mom[r][n]==sum(v*k**r for (m,k),v in small.items() if m==n),'exact_raw_moment_recurrence')
    is_sc='_sc' in name
    for n in [120,121,360,361,1080,1081]:
        Z,S,S2=[row[n] for row in mom]
        if not Z:
            ck(is_sc and (name=='delta4_sc' and n%2==0 or name=='synthetic_sc_even_only_unbounded' and n%2==1),'unsupported_parity')
            continue
        E=Fraction(S,Z);V=Fraction(S2,Z)-E*E
        row={'family':name,'n':n,'mean':float(E),'variance':float(V),'mean_residual':float(E)-mu*n,'variance_residual':float(V)-(2 if is_sc else 1)*sigma2*n}
        ck(abs(row['mean_residual'])<12,'finite_mean_residual_sanity')
        ck(abs(row['variance_residual'])<12,'finite_variance_residual_sanity')
        moment_rows.append(row)

# Actual order-six inflation C3[11,11,11]. Its canonical block action is diagonal C3,
# not two independent rotations of two visible triangles.
T=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(6):
        if i!=j:T[i][j]=int((i//2==j//2 and i<j) or (i//2!=j//2 and (j//2-i//2)%3==1))
autos=[g for g in itertools.permutations(range(6)) if all(T[i][j]==T[g[i]][g[j]] for i in range(6) for j in range(6))]
f=min(sum(T[o[j]][o[i]] for i in range(6) for j in range(i+1,6)) for o in itertools.permutations(range(6)))
s=sum(any(g[i]!=i for g in autos) for i in range(6))
ck((len(autos),f,s,3*f-s)==(3,4,6,6),'actual_diagonal_semidirect_group_control')
ck(len(autos)!=3**f and len(autos)!=9,'not_feedback_exponential_exact_not_independent_triangles')

# Exact bounded perturbation hostile: bounded does NOT imply O(1) variance correction.
# X_n=±sqrt(n), B_n=sign(X_n), so Var(X_n+B_n)-Var(X_n)=2sqrt(n)+1.
for m in [1,2,7,20,100]:
    n=m*m
    ck(Fraction((m+1)**2+(-m-1)**2,2)-n==2*m+1,'bounded_perturbation_covariance_hostile')

root=Path(__file__).resolve().parent
source=root.parent/'session23-root'/'SLACK_RENEWAL_TRANSFER.md'
compiler_source=root.parent/'session23-connections'/'SYMMETRY_SLACK_RATIONALITY.md'
result={
 'status':'PASS', 'total_assertions':sum(CHECKS.values()), 'exact_assertions':sum(CHECKS.values())-CHECKS['finite_mean_residual_sanity']-CHECKS['finite_variance_residual_sanity'], 'numerical_sanity_assertions':CHECKS['finite_mean_residual_sanity']+CHECKS['finite_variance_residual_sanity'], 'checks_by_kind':dict(sorted(CHECKS.items())),
 'scope':'Synthetic finite compiler models are not actual slack inventories. Actual low-slack formula checks concern delta 0, 3, 4 only; one actual delta 6 diagonal-group control. Finite tests do not prove asymptotics.',
 'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
 'compiler_sha256':hashlib.sha256(compiler_source.read_bytes()).hexdigest(),
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'rho':rho,'mu':mu,'sigma2':sigma2,'low_slack_counts_order_0_through_7':low_rows,
 'moment_controls':moment_rows,
}
print(json.dumps(result,indent=2,sort_keys=True))
