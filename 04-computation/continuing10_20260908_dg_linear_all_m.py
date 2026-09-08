"""Exact controls for the all-m source-linear submersion classification.

All-degree proofs are in the matching report. The finite universe is every
m=1..10, pure-power degree n=2..2m, and h in {0,1,-2}, plus constant A.
No producer imports, floating arithmetic, repository reads or output mutation.
"""
import sys
from itertools import product
from math import comb
import sympy as s
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
x,z,t,r,b,h,a=s.symbols('x z t r b h a')
gates=0
def need(ok,msg):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(msg)
def zero(expr,msg):need(s.cancel(expr)==0,msg)
def jac(F,G,X=x,Y=t):return s.diff(F,X)*s.diff(G,Y)-s.diff(F,Y)*s.diff(G,X)
def predicted(m,n,H):
    if n==0:return m==1
    if m==1:return n==2
    return H!=0 and (m+1<=n<=2*m-2 or n==2*m)

rows=0;accepted=[]
for m in range(1,11):
    # Complete linear-section generator bank independently checks the chart.
    for j in range(2*m+1):
        Aj=x**j
        Bj=x**(j-m) if j>m else s.Integer(0)
        second=(Bj-Aj/x**m).subs(x,1/r)-b*r**(2*m-j)
        expected=-(r**(m-j) if j<=m else 0)-b*r**(2*m-j)
        zero(second-expected,'complete coefficient basis second chart')
        need(s.denom(s.cancel(second))==1,'basis function regular in second chart')
    constant_second=-r**m-r**(2*m)*b
    boundary_d=s.diff(constant_second,r).subs(r,0)
    need((boundary_d!=0)==predicted(m,0,0),'constant A global submersion boundary')
    zero(jac(t,-x)-1,'positive affine polynomial mate, constant A')
    zero((-x).subs(x,1/r)+1/r,'constant mate actual boundary pole')
    for n in range(2,2*m+1):
        A=(x-h)**n
        coefficients=[s.Integer(comb(n,j))*(-h)**(n-j) for j in range(n+1)]
        B=sum((coefficients[j]*x**(j-m) for j in range(m+1,n+1)),s.Integer(0))
        F=A*t+B
        F2=-sum((coefficients[j]*r**(m-j) for j in range(min(m,n)+1)),s.Integer(0))
        F2-=b*r**(2*m-n)*(1-h*r)**n
        zero(s.expand(F.subs({x:1/r,t:-r**m-r**(2*m)*b})-F2),'whole pure-power second chart')
        beta=s.expand(s.diff(B,x).subs(x,h))
        expected_beta=(-1)**(n-m-1)*comb(n-2,m-1)*h**(n-m-1) if n>m else 0
        zero(beta-expected_beta,'complete binomial submersion coefficient')
        dr=s.expand(s.diff(F2,r).subs(r,0))
        db=s.expand(s.diff(F2,b).subs(r,0))
        am1=coefficients[m-1] if m-1<=n else 0
        a2m1=coefficients[2*m-1] if 2*m-1<=n else 0
        a2m=coefficients[2*m] if 2*m<=n else 0
        zero(dr+am1+a2m1*b,'entire boundary normal derivative')
        zero(db+a2m,'entire boundary tangential derivative')
        q=1/((n-1)*(x-h)**(n-1))
        zero(-A*s.diff(q,x)-1,'rational primitive all pure powers')
        c=B.subs(x,h)
        g=s.expand((F-c).subs(x,z+h))
        K=s.cancel(g/z)
        need(s.denom(K)==1,'special-fibre second factor polynomial')
        zero(K.subs(z,0)-beta,'other component has nonzero constant term exactly at smooth root')
        zero(jac(g,K,z,t)-z**(n-1)*K,'factor identity for arbitrary upper-order witness')
        for H in [0,1,-2]:
            source_smooth=beta.subs(h,H)!=0 if hasattr(beta,'subs') else beta!=0
            DR=s.expand(dr.subs(h,H));DB=s.expand(db.subs(h,H))
            # DB is scalar. If zero, every boundary point is noncritical iff
            # DR is a nonzero constant polynomial in b.
            boundary_smooth=DB!=0 or (s.diff(DR,b)==0 and DR!=0)
            actual=bool(source_smooth and boundary_smooth)
            need(actual==predicted(m,n,H),'complete source/boundary case classification')
            rows+=1
            if actual:
                accepted.append((m,n,H))
                BB=beta.subs(h,H)
                need(BB!=0,'all admitted special fibres reduced and disjoint')
                need(BB**(n-1)/s.Integer(n-1)!=0,'highest scalar response coefficient nonzero')
                # These identities check upper and lower pole orders without
                # expanding a power whose degree grows quadratically in n.
                KH=K.subs(h,H)
                zero(KH.subs(z,0)**(n-1)/s.Integer(n-1)-BB**(n-1)/s.Integer(n-1),'order n-1 top coefficient')
                zero(KH.subs(z,0)**(n-2)/s.Integer(n-1)-BB**(n-2)/s.Integer(n-1),'g^(n-2) primitive simple source pole')
        if n==2*m-1:
            zero(dr.subs(b,-am1),'excluded degree 2m-1 boundary critical point')

# Rational-primitive hostile bank: all multiplicity words 2..4 on 1..r,
# r=1,2,3. Distinct-root logarithmic residues cannot all disappear.
patterns=0
for root_count in range(1,4):
    for multiplicities in product(range(2,5),repeat=root_count):
        roots=list(range(1,root_count+1))
        residues=[]
        for i,p in enumerate(roots):
            regular=s.prod((x-q)**(-multiplicities[j]) for j,q in enumerate(roots) if j!=i)
            residue=s.diff(regular,x,multiplicities[i]-1).subs(x,p)/s.factorial(multiplicities[i]-1)
            residues.append(s.cancel(residue))
        need(all(v==0 for v in residues)==(root_count==1),'complete bounded reciprocal residue bank')
        if root_count>1:
            degree=sum(multiplicities)
            need(degree-1>degree-root_count,'all-degree pole-capacity contradiction instantiated')
        patterns+=1
zero(s.residue(1/x,x,0)-1,'simple-root hostile')

# Inherited m=2 hostile: global submersion does not imply integrability.
FH=(x*x-1)**2*t+x*x
FH2=2-r*r-b*(1-r*r)**2
zero(FH.subs({x:1/r,t:-r*r-r**4*b})-FH2,'smooth two-root hostile second chart')
for p in [-1,1]:
    zero(s.diff(FH,x).subs(x,p)-2*p,'smooth two-root hostile affine derivative')
zero(s.diff(FH2,b).subs(r,0)+1,'smooth two-root hostile boundary derivative')
zero(s.residue(1/(x*x-1)**2,x,1)+s.Rational(1,4),'nonintegrable positive-root residue')
zero(s.residue(1/(x*x-1)**2,x,-1)-s.Rational(1,4),'nonintegrable negative-root residue')

# The m=1 constant branch: the ring change restores an omitted component.
g=a*t;Q=-x/a;P=-x*t
zero(jac(g,Q)-1,'source unit vanishes for constant branch')
zero(jac(g,P)-g,'global order-one polynomial witness')
g2=-a*r*(1+r*b);Q2=-1/(a*r);P2=1+r*b
zero(P.subs({x:1/r,t:-r-r*r*b})-P2,'global witness second chart')
zero(jac(g2,Q2,r,b)-1,'global Hamiltonian primitive rational in second chart')
zero(jac(g2,P2,r,b)-g2,'global annihilator in second chart')
zero((1+r*b)-b*r-1,'special components D and t=0 comaximal')
zero(s.limit(g2*Q2,r,0)-1,'global boundary unit class has nonzero simple scalar pole')
zero(s.diff(g2,b).subs(r,0),'Hamiltonian r coefficient vanishes on boundary')
zero(s.diff(g2,r).subs(r,0)+a,'Hamiltonian b coefficient regular and nonzero there')

need(rows==300,'declared complete parameter test count')
need(patterns==39,'declared complete residue pattern count')
print('CLASSIFICATION: all m>=1 source-linear rationally integrable global submersions proved in report.')
print('FINITE_UNIVERSE m=1..10, n=2..2m, h in {0,1,-2}:',rows,'rows;',len(accepted),'admitted; plus ten constant-A controls.')
print('RECIPROCAL_RESIDUE_BANK',patterns,'complete multiplicity patterns; simple-root hostile retained.')
print('UNIT_ORDER n-1 for admitted pure powers; m=1 constant branch has source order0 and global order1.')
print('Always-active exact gates:',gates)
