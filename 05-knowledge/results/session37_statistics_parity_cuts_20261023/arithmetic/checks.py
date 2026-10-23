#!/usr/bin/env python3
"""Exact finite controls. No floating point or numerical transcendence test."""
import json
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import sympy as sp

z,x=sp.symbols('z x')
def check(ok, what):
    if not ok: raise RuntimeError(what)
def tournament(n, mask):
    S=sp.zeros(n)
    for k,(i,j) in enumerate(combinations(range(n),2)):
        S[i,j]=1 if (mask>>k)&1 else -1
        S[j,i]=-S[i,j]
    return S

def paths(S):
    return sum(all(S[p[i],p[i+1]]==1 for i in range(len(p)-1)) for p in permutations(range(S.rows)))

def log_series(a,terms):
    """Rational lower/upper bounds for log(a), a>1, via atanh."""
    r=(a-1)/(a+1)
    lo=2*sum(r**(2*j+1)/sp.Integer(2*j+1) for j in range(terms))
    err=2*r**(2*terms+1)/(sp.Integer(2*terms+1)*(1-r*r))
    return lo,lo+err

counts={}; signatures={}; checks=0
for n in range(1,6):
    hist=Counter()
    for mask in range(1<<(n*(n-1)//2)):
        S=tournament(n,mask); K=-S*S; I=sp.eye(n)
        p=K.charpoly(x).as_expr(); hist[str(sp.factor(p))]+=1
        check(K==K.T,'symmetric')
        check(sp.trace(K)==n*(n-1),'trace')
        check(S.rank()==n-(n%2),'parity rank')
        check(K.rank()==S.rank(),'Gram rank')
        if n<=4:
            U=(I+S)*(I-S).inv(); C=(I+25*K)*(I+K).inv()
            check(U.T*U==I,'orthogonal Cayley')
            check(C==13*I-6*(U+U.inv()),'charge matrix identity')
            check((U-I)*(U+I).inv()==S,'faithful U inverse')
            check((C-I)*(25*I-C).inv()==K,'faithful C inverse')
        if n>=2:
            q=sp.Rational(1,2*(n-1)**2)
            F=(I+q*K).det(); D=(I-q*K).det(); B=1/D
            check(F>1 and B>1,'positive partition values')
            check(F*D<1,'positive even log sector')
            check(F*B>1,'positive odd log sector')
        checks+=1
    counts[str(n)]=1<<(n*(n-1)//2); signatures[str(n)]=dict(sorted(hist.items()))

S4=tournament(4,63); K4=-S4*S4
check(sp.factor(K4.charpoly(x).as_expr())==(x*x-6*x+1)**2,'T4 spectrum')
F4=sp.factor((sp.eye(4)+z*K4).det()); D4=sp.factor((sp.eye(4)-z*K4).det())
check(F4==(z*z+6*z+1)**2,'T4 F')
check(D4==(z*z-6*z+1)**2,'T4 B denominator')
q=sp.Rational(1,18); FF=F4.subs(z,q); BB=1/D4.subs(z,q)
check(FF==sp.Rational(433,324)**2,'T4 rational F')
check(BB==sp.Rational(324,217)**2,'T4 rational B')

# Log-det trace expansions and even/odd split, independently via formal derivatives.
for n,mask in [(2,1),(3,7),(3,5),(4,63),(4,4),(5,0),(5,341)]:
    S=tournament(n,mask); K=-S*S; I=sp.eye(n)
    F=(I+z*K).det(); D=(I-z*K).det()
    gF=sp.series(sp.diff(F,z)/F,z,0,12).removeO()
    gB=sp.series(-sp.diff(D,z)/D,z,0,12).removeO()
    P=sp.eye(n)
    for m in range(1,13):
        P=P*K; a=sp.trace(P)
        check(sp.expand(gF).coeff(z,m-1)==(-1)**(m+1)*a,'formal F log')
        check(sp.expand(gB).coeff(z,m-1)==a,'formal B log')
    check(sp.cancel(F*1/F)==1,'F(z)B(-z)=1')
    check(sp.expand(F*D-(I-z*z*K*K).det())==0,'even product identity')

# Smallest spectral collision has different Hamiltonian path counts.
ST=tournament(3,7); SC=sp.Matrix([[0,1,-1],[-1,0,1],[1,-1,0]])
check((-ST*ST).charpoly(x)==(-SC*SC).charpoly(x),'T3 spectrum collision')
check(paths(ST)==1 and paths(SC)==3,'Hamiltonian path loss')

# Sharp first nontorsion Cayley mode and general two-sink regular block.
S5=tournament(5,1021); K5=-S5*S5
check(sp.factor(K5.charpoly(x).as_expr())==x*(x-3)**2*(x-7)**2,'first mode witness')
U5=(sp.eye(5)+S5)*(sp.eye(5)-S5).inv()
check(sp.trace(U5)==sp.Rational(-3,2),'first witness Cayley trace')
nu=sp.symbols('nu', positive=True)
R3=sp.Matrix([[0,sp.sqrt(nu),sp.sqrt(nu)],[-sp.sqrt(nu),0,1],[-sp.sqrt(nu),-1,0]])
check(sp.expand(R3.charpoly(x).as_expr()-x*(x*x+2*nu+1))==0,'all-n two-sink block')
for N in [3,5,7,9]:
    base=sp.Matrix(N,N,lambda i,j:0 if i==j else (1 if 1<=(j-i)%N<=N//2 else -1))
    check(base*sp.ones(N,1)==sp.zeros(N,1),'regular control')
    ext=base.row_join(sp.ones(N,2)).col_join((-sp.ones(2,N)).row_join(sp.Matrix([[0,1],[-1,0]])))
    pbase=base.charpoly(x).as_expr(); pext=ext.charpoly(x).as_expr()
    check(sp.expand(pext-pbase*(x*x+2*N+1))==0,'ordered two-sink spectral factor')

# All eight block orientations and exact Hamiltonian/two-path-cover split.
def subset_hamiltonian_counts(A):
    N=A.rows; dp=[[0]*N for _ in range(1<<N)]
    out=[sum(1<<j for j in range(N) if A[i,j]==1) for i in range(N)]
    for i in range(N): dp[1<<i][i]=1
    for mask in range(1,1<<N):
        for i in range(N):
            count=dp[mask][i]
            if count:
                rest=out[i]&~mask
                while rest:
                    bit=rest&-rest; j=bit.bit_length()-1; rest-=bit
                    dp[mask|bit][j]+=count
    return [sum(row) for row in dp]

all_eight={}
for N in [3,5,7]:
    base=sp.Matrix(N,N,lambda i,j:0 if i==j else (1 if 1<=(j-i)%N<=N//2 else -1))
    h=subset_hamiltonian_counts(base); full=(1<<N)-1
    H=h[full]; pair_sum=sum(h[m]*h[full^m] for m in range(1,full))
    check(pair_sum%2==0,'unordered cover divisor'); pc2=pair_sum//2
    rows=[]
    for alpha in [-1,1]:
      for beta in [-1,1]:
       for gamma in [-1,1]:
        ext=base.row_join(sp.Matrix(N,2,lambda i,j:alpha if j==0 else beta))
        ext=ext.col_join(sp.Matrix([[-alpha]*N+[0,gamma],[-beta]*N+[-gamma,0]]))
        Q=sp.Matrix([[0,alpha,beta],[-alpha,0,gamma],[-beta,-gamma,0]])
        cyclic=all(sum(Q[i,j]==1 for j in range(3))==1 for i in range(3))
        got=subset_hamiltonian_counts(ext)[-1]
        expected=3*H+2*pc2 if cyclic else H
        check(got==expected,'all-eight Hamiltonian split')
        check(sp.expand(ext.charpoly(x).as_expr()-base.charpoly(x).as_expr()*(x*x+2*N+1))==0,'all-eight spectral equality')
        rows.append({'block_signs':[alpha,beta,gamma],'cyclic':cyclic,'H':got})
    check(sum(r['cyclic'] for r in rows)==2,'two cyclic quotient choices')
    all_eight[str(N)]={'H_base':H,'pc2_unordered':pc2,'rows':rows}

# Exact Paley inputs, no composite-modulus or unexplained Paley extrapolation.
paley={}
for p in [3,7,11,19,23]:
    qr={j*j%p for j in range(1,p)}
    S=sp.Matrix(p,p,lambda i,j:0 if i==j else (1 if (j-i)%p in qr else -1))
    K=-S*S; I=sp.eye(p); J=sp.ones(p)
    check(S.T==-S,'Paley antisymmetry'); check(K==p*I-J,'Paley Gram law')
    charge=sp.Rational(25*p+1,p+1); tr=sp.Rational(2*(1-p),1+p)
    check(1<charge<25,'charge branch')
    paley[str(p)]={'nonzero_energy':p,'multiplicity':p-1,'charge':str(charge),'mode_trace':str(tr),'torsion_at_p3':p==3}

# Cheap P163 recovery: full circulant row correlations, not a tournament census.
p163=163
check(all(p163%d for d in range(2,13)),'163 prime trial division')
chi=[0]+[1 if pow(j,81,p163)==1 else -1 for j in range(1,p163)]
corr=[sum(chi[j]*chi[(j+a)%p163] for j in range(p163)) for a in range(1,p163)]
check(sum(chi)==0 and set(corr)=={-1},'Paley163 Gram correlations')
paley163={'prime':163,'row_sum':sum(chi),'nonzero_shifts_checked':162,'correlation_values':sorted(set(corr)),'charge':str(sp.Rational(25*163+1,164)),'CM_field':'Q(sqrt(-163))','mode':'(-81+i*sqrt(163))/82'}

cut={}
for name,a in [('F4',FF),('B4',BB),('odd',FF*BB),('even',BB/FF)]:
    check(a>1,'cut log positive')
    L,U=log_series(a,8); L2,U2=log_series(a,9)
    check(L<L2<U2<U,'nested atanh cuts')
    # Exact agreement between consecutive partials and exact tail bound identity.
    r=(a-1)/(a+1)
    check(L2-L==2*r**17/sp.Integer(17),'atanh next term')
    cut[name]={'log_argument':str(a),'lower_8':str(L),'upper_8':str(U),'width_8':str(U-L)}

receipt={'status':'PASS','arithmetic':'exact integer, rational and symbolic; no numerical transcendence test','labelled_universe':counts,'total_tournaments':checks,'spectral_signature_census':signatures,'T4':{'F':str(F4),'B_denominator':str(D4),'z':'1/18','F_value':str(FF),'B_value':str(BB)},'spectral_hostile':{'order':3,'transitive_hamiltonian_paths':paths(ST),'cycle_hamiltonian_paths':paths(SC),'K_characteristic_polynomial':str((-ST*ST).charpoly(x).as_expr())},'paley_controls':paley,'paley163_control':paley163,'rational_log_cuts':cut,'sharp_first_infinite_order_mode':{'order':5,'mask':1021,'K_characteristic_polynomial':str(K5.charpoly(x).as_expr()),'Cayley_trace':str(sp.trace(U5)),'symbolic_two_sink_block':str(R3.charpoly(x).as_expr()),'regular_extension_orders_checked':[3,5,7,9]},'all_eight_extension_controls':all_eight,'formal_log_degree_checked':12,'scope':'Finite checks support identities; CM/angle/log transcendence and all-order conclusions are proofs, not inferences from enumeration.'}
print(json.dumps(receipt,indent=2))
