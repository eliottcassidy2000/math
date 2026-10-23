#!/usr/bin/env python3
"""Independent Session37 controls. Standard-library exact matrices; SymPy only factors/root isolation.
Does not import the producer's code. Exceptions remain active under python -O.
"""
from fractions import Fraction as Q
from itertools import combinations, permutations
from collections import Counter
import json
import sympy as sp

COUNT=0
def need(p, why):
    global COUNT
    COUNT+=1
    if not p: raise RuntimeError(why)
def eye(n): return [[Q(i==j) for j in range(n)] for i in range(n)]
def tr(A): return sum(A[i][i] for i in range(len(A)))
def scale(A,q): return [[q*v for v in row] for row in A]
def add(A,B): return [[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def transpose(A): return list(map(list,zip(*A)))
def mul(A,B): return [[sum(a*b for a,b in zip(ar,bc)) for bc in zip(*B)] for ar in A]
def power(A,k):
    out=eye(len(A))
    for _ in range(k):out=mul(out,A)
    return out
def inverse(A):
    n=len(A); M=[[Q(x) for x in r]+e for r,e in zip(A,eye(n))]
    for j in range(n):
        h=next(k for k in range(j,n) if M[k][j])
        M[h],M[j]=M[j],M[h]; pivot=M[j][j]
        M[j]=[x/pivot for x in M[j]]
        for k in range(n):
            if k!=j:
                c=M[k][j]; M[k]=[a-c*b for a,b in zip(M[k],M[j])]
    return [r[n:] for r in M]
def det(A):
    M=[[Q(x) for x in r] for r in A]; d=Q(1);n=len(A)
    for j in range(n):
        pivots=[k for k in range(j,n) if M[k][j]]
        if not pivots:return Q(0)
        k=pivots[0]
        if k!=j: M[j],M[k]=M[k],M[j];d=-d
        p=M[j][j];d*=p
        for k in range(j+1,n):
            c=M[k][j]/p
            for h in range(j+1,n):M[k][h]-=c*M[j][h]
    return d
def char_coeff(A):
    # Faddeev-LeVerrier, coefficients descending, independent of symbolic determinant.
    n=len(A);B=eye(n);out=[Q(1)]
    for k in range(1,n+1):
        AB=mul(A,B);c=-tr(AB)/Q(k);out.append(c); B=add(AB,scale(eye(n),c))
    need(all(c.denominator==1 for c in out),'monic integer characteristic')
    return tuple(int(c) for c in out)
def tournament(n, mask):
    M=[[0]*n for _ in range(n)]
    for b,(i,j) in enumerate(combinations(range(n),2)):
        M[i][j]=1 if (mask>>b)&1 else -1;M[j][i]=-M[i][j]
    return M
def characteristic_eval(c,q):
    v=Q(0)
    for a in c:v=v*q+a
    return v
def cut(r,m):
    v=(r-1)/(r+1)
    lo=2*sum((v**(2*j+1)/Q(2*j+1) for j in range(m)),Q(0))
    return lo,lo+2*v**(2*m+1)/(Q(2*m+1)*(1-v*v))
def paths(S):
    return sum(all(S[p[j]][p[j+1]]>0 for j in range(len(S)-1)) for p in permutations(range(len(S))))

census={}; signatures={}; positive_factors={}; x=sp.symbols('x'); cases=0
for n in range(1,6):
    hist=Counter(); value_cuts=set()
    for mask in range(1<<(n*(n-1)//2)):
        S=tournament(n,mask);I=eye(n);K=scale(mul(S,S),-1); c=char_coeff(K)
        cases+=1;hist[c]+=1
        need(K==mul(transpose(S),S),'Gram identity')
        need(tr(K)==n*(n-1),'trace normalization')
        if n%2==0:
            need(det(S).numerator%2==1,'even skew tournament determinant odd')
        else:
            need(det(S)==0,'odd skew determinant zero')
            sub=[r[:-1] for r in S[:-1]]
            need(det(sub).numerator%2==1,'odd rank via invertible principal minor')
        U=mul(add(I,S),inverse(add(I,scale(S,-1))))
        need(mul(U,transpose(U))==I,'Cayley rational orthogonality')
        need(mul(add(U,scale(I,-1)),inverse(add(U,I)))==S,'full Cayley inverse retains S')
        C=mul(add(I,scale(K,25)),inverse(add(I,K)))
        need(C==add(scale(I,13),scale(add(U,transpose(U)),-6)),'Cayley charge identity')
        need(mul(add(C,scale(I,-1)),inverse(add(scale(I,25),scale(C,-1))))==K,'full charge inverse retains K')
        z=Q(1,2*max(1,(n-1)**2));F=det(add(I,scale(K,z)));D=det(add(I,scale(K,-z)));B=1/D
        need(F*D==det(add(I,scale(mul(K,K),-z*z))),'exact even determinant identity')
        if n==1:need((F,B)==(1,1),'one-vertex log exception')
        else:
            need(1<F<B,'positive partition values and difference')
            need(0<D<1 and F*D<1,'even log argument positive')
            for r in (F,B,F*B,B/F):value_cuts.add(r)
    for r in value_cuts:
        previous=None
        for m in range(1,10):
            L,U=cut(r,m)
            need(0<L<U,'strict cut endpoints')
            if previous:need(previous[0]<L<U<previous[1],'nested rational enclosures')
            previous=(L,U)
    census[str(n)]={'labelled_tournaments':sum(hist.values()),'distinct_spectra':len(hist),'distinct_log_arguments':len(value_cuts)}
    signatures[str(n)]={','.join(map(str,k)):v for k,v in sorted(hist.items())}
    for c in hist:
        p=sp.Poly.from_list(c,x)
        for f,m in sp.factor_list(p)[1]:
            if f.as_expr()==x:continue
            need(m%2==0,'each nonzero energy has even multiplicity')
            intervals=sp.polys.polytools.intervals(f,eps=sp.Rational(1,1000000))
            need(sum(mult for _,mult in intervals)==f.degree(),'all conjugates real')
            need(all(a>0 and b>0 for (a,b),_ in intervals),'all conjugates strictly positive')
            positive_factors[str(f.as_expr())]=str(intervals)

# Rational energy torsion candidates are obtained from rational integral traces in (-2,2).
torsion_mu=sorted({Q(2-c,2+c) for c in [-1,0,1]})
need(torsion_mu==[Q(1,3),Q(1),Q(3)],'complete rational energy torsion candidates')
need([v for v in torsion_mu if v.denominator==1]==[1,3],'integral torsion candidates')

# Transitive family: its full rational Cayley matrix is an exact signed cycle.
transitive=[]
for n in range(1,13):
    S=tournament(n,(1<<(n*(n-1)//2))-1);I=eye(n)
    U=mul(add(I,S),inverse(add(I,scale(S,-1))))
    need(power(U,n)==scale(I,(-1)**(n-1)),'transitive Cayley signed-cycle power')
    transitive.append({'n':n,'U_to_n':(-1)**(n-1)})
S4=tournament(4,63);I4=eye(4);K4=scale(mul(S4,S4),-1)
need(char_coeff(K4)==(1,-12,38,-12,1),'T4 conjugate irrational energies')
U4=mul(add(I4,S4),inverse(add(I4,scale(S4,-1))))
need(power(U4,4)==scale(I4,-1),'T4 primitive eighth-root phase hostile')
z=Q(1,18);F4=det(add(I4,scale(K4,z)));B4=1/det(add(I4,scale(K4,-z)))
need(F4==Q(433,324)**2 and B4==Q(324,217)**2,'T4 explicit fugacity values')

# The same spectral data can have different graph observables.
ST=tournament(3,7);SC=[[0,1,-1],[-1,0,1],[1,-1,0]]
need(char_coeff(scale(mul(ST,ST),-1))==char_coeff(scale(mul(SC,SC),-1)),'T3 spectral collision')
need(paths(ST)==1 and paths(SC)==3,'T3 paths differ')

# Paley primes are explicit exact finite controls, not a DRT existence assertion.
paley={}
for p in [3,7,11,19]:
    squares={a*a%p for a in range(1,p)}
    S=[[0 if i==j else 1 if (j-i)%p in squares else -1 for j in range(p)] for i in range(p)]
    need(transpose(S)==scale(S,-1),'Paley skew')
    K=scale(mul(S,S),-1);need(K==[[p*(i==j)-1 for j in range(p)] for i in range(p)],'Paley DRT Gram')
    c=Q(25*p+1,p+1);tau=Q(2*(1-p),p+1)
    need(c==13-6*tau,'DRT charge')
    need((tau.denominator==1)==(p==3),'DRT root-of-unity obstruction')
    paley[str(p)]={'phase_trace':str(tau),'charge':str(c)}

# Regular two-sink extension: block identities and strict regularity hostile.
extensions=[]
for n in [1,3,5,7,9,11]:
    # Regular cyclic tournament, not necessarily doubly regular.
    S=[[0 if i==j else 1 if 1<=(j-i)%n<=(n-1)//2 else -1 for j in range(n)] for i in range(n)]
    need(all(sum(r)==0 for r in S),'cyclic base regular')
    T=[row+[1,1] for row in S]+[[-1]*n+[0,1],[-1]*n+[-1,0]]
    c_old=char_coeff(S);c_new=char_coeff(T)
    p_old=sp.Poly.from_list(c_old,x).as_expr();p_new=sp.Poly.from_list(c_new,x).as_expr()
    need(sp.expand(p_new-p_old*(x*x+2*n+1))==0,'two-sink exact characteristic lift')
    block=[[0,1,1],[-n,0,1],[-n,-1,0]]
    need(char_coeff(block)==(1,0,2*n+1,0),'two-sink residual block')
    mu=2*n+1;c=Q(25*mu+1,mu+1)
    need(c==25-Q(12,n+1),'two-sink charge')
    need((mu in [1,3])==(n==1),'two-sink torsion boundary')
    extensions.append({'regular_base_n':n,'added_energy':mu,'added_charge':str(c)})
S5=tournament(5,1021)
need(char_coeff(S5)==(1,0,10,0,21,0),'sharp first nontorsion S5 characteristic')
need(char_coeff(scale(mul(S5,S5),-1))==(1,-20,142,-420,441,0),'sharp first nontorsion K5 characteristic')
need(char_coeff(tournament(5,1023))!=(1,0,10,0,21,0),'nonregular transitive3 extension hostile')

# All eight quotient orientations: independent subset DP path count.
def subset_path_counts(S):
    n=len(S);dp=[[0]*n for _ in range(1<<n)];total=[0]*(1<<n)
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask>>v)&1:continue
            previous=mask^(1<<v)
            if previous:dp[mask][v]=sum(dp[previous][u] for u in range(n) if (previous>>u)&1 and S[u][v]>0)
        total[mask]=sum(dp[mask])
    return total

def strongly_connected(S):
    n=len(S)
    for start in range(n):
        seen={start};stack=[start]
        while stack:
            v=stack.pop()
            for w in range(n):
                if S[v][w]>0 and w not in seen:seen.add(w);stack.append(w)
        if len(seen)!=n:return False
    return True

all_eight=[]
for n in [1,3,5,7]:
    S=[[0 if i==j else 1 if 1<=(j-i)%n<=(n-1)//2 else -1 for j in range(n)] for i in range(n)]
    totals=subset_path_counts(S);full=(1<<n)-1;H=totals[full]
    pc2=sum(totals[m]*totals[full^m] for m in range(1,full))//2
    if n<=5:need(H==paths(S),'base DP versus permutation path census')
    base_poly=sp.Poly.from_list(char_coeff(S),x).as_expr();row=[]
    for alpha in [-1,1]:
        for beta in [-1,1]:
            for gamma in [-1,1]:
                T=[r+[alpha,beta] for r in S]+[[-alpha]*n+[0,gamma],[-beta]*n+[-gamma,0]]
                quotient=[[0,alpha,beta],[-alpha,0,gamma],[-beta,-gamma,0]]
                cyclic=(alpha==gamma and beta==-alpha)
                need(paths(quotient)==(3 if cyclic else 1),'quotient type')
                target=(3*H+2*pc2) if cyclic else H
                actual=subset_path_counts(T)[-1]
                need(actual==target,'all-eight Hamiltonian lift formula')
                need(strongly_connected(T)==cyclic,'cyclic extension strong, transitive extension reducible')
                need(not cyclic or actual>=(n+2)*H,'unbounded-ratio cut-count bound')
                if n<=3:need(actual==paths(T),'extension DP versus permutation census')
                extension_poly=sp.Poly.from_list(char_coeff(T),x).as_expr()
                need(sp.expand(extension_poly-base_poly*(x*x+2*n+1))==0,'all-eight common spectral factor')
                row.append({'alpha':alpha,'beta':beta,'gamma':gamma,'cyclic':cyclic,'Hamiltonian_paths':actual,'strong':strongly_connected(T),'score_multiset':sorted(sum(v>0 for v in r) for r in T)})
    need(sum(r['cyclic'] for r in row)==2,'two cyclic six transitive quotient choices')
    all_eight.append({'regular_base_n':n,'H_base':H,'unordered_two_path_covers':pc2,'extensions':row})

# Outside the Bose radius positive rational B need not be a positive partition series.
K2=[[1,0],[0,1]];need(1/det(add(eye(2),scale(K2,-3)))==Q(1,4),'Bose outside-radius hostile')
# Charge is field-preserving but is not generally an algebraic integer.
need(Q(25*7+1,8)==22,'n7 integer charge benign');need(Q(25*11+1,12)==23,'n11 integer charge benign')
need(Q(25*19+1,20).denominator>1,'n19 nonintegral rational charge hostile')

receipt={'status':'PASS','assertion_style':'explicit exceptions; active under python -O','tournaments':cases,'checks':COUNT,'census':census,'signatures':signatures,'root_isolation':positive_factors,'rational_torsion_energies':[str(v) for v in torsion_mu],'transitive_signed_cycle':transitive,'regular_two_sink_controls':extensions,'all_eight_extension_controls':all_eight,'paley':paley,'T4':{'F':str(F4),'B':str(B4),'mu_polynomial':'x^2-6x+1','phase_order':8},'collision':{'n':3,'paths':[paths(ST),paths(SC)]},'scope':'Exact finite controls only; no numerical test of transcendence; all-order proofs in REPORT.md.'}
print(json.dumps(receipt,indent=2,sort_keys=True))
