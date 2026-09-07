"""Exact native DG linear-carrier and rational-primitive controls.

All-degree conclusions are proved in the report. No source producer imports.
The finite universe is symbolic identities plus separately solved finite
polynomial coefficient systems; it is not a degree-bounded claim of closure.
"""
import hashlib
from pathlib import Path
import sys
import sympy as s
sys.stdout.reconfigure(newline='\n')
x,t,r,b,h,k,a,c,z=s.symbols('x t r b h k a c z')
G=0
def need(value,label):
    global G
    G+=1
    if not bool(value): raise RuntimeError(label)
def zero(value,label): need(s.cancel(value)==0,label)
def jac(F,H): return s.diff(F,x)*s.diff(H,t)-s.diff(F,t)*s.diff(H,x)
def chart(F): return s.cancel(F.subs({x:1/r,t:-r*r-r**4*b},simultaneous=True))

# Full six-dimensional native layer, including arbitrary coefficient symbols.
Aco=s.symbols('a0:5'); B0=s.symbols('B0')
A=sum(Aco[i]*x**i for i in range(5))
B=Aco[4]*x*x+Aco[3]*x+B0
F=A*t+B
Finf=B0-Aco[2]-Aco[1]*r-Aco[0]*r*r-b*sum(Aco[i]*r**(4-i) for i in range(5))
zero(chart(F)-Finf,'complete regular second-chart expression')
u=x*x*t; j=x*t; v=x*(1+x*x*t); bb=-x*x-x**4*t
zero(F-(Aco[0]*t+Aco[1]*j+Aco[2]*u+Aco[3]*v-Aco[4]*bb+B0),
     'exact full native generator span')
zero(Finf.subs(r,0)-(B0-Aco[2]-Aco[4]*b),'actual boundary restriction')
zero(s.diff(Finf,r).subs(r,0)+Aco[1]+Aco[3]*b,'actual first normal derivative')

# Independent coefficient ranks recover all Laurent cancellation conditions.
layer=[]
for degree in range(1,9):
    aa=s.symbols('q0:'+str(degree+1)); cc=s.symbols('w0:'+str(degree+1))
    trial=sum(aa[i]*x**i*t+cc[i]*x**i for i in range(degree+1))
    transformed=s.expand(trial.subs({x:1/r,t:-r*r-r**4*b},simultaneous=True))
    bad={}
    for term in s.Add.make_args(transformed):
        power=int(term.as_powers_dict().get(r,0))
        if power<0:
            bd=int(term.as_powers_dict().get(b,0))
            bad[(power,bd)]=bad.get((power,bd),0)+term/(r**power*b**bd)
    variables=aa+cc
    matrix=s.linear_eq_to_matrix(list(bad.values()),variables)[0]
    expected=min(degree+2,6)  # The six global generators enter successively.
    need(len(variables)-matrix.rank()==expected,'full Laurent coefficient-space rank')
    layer.append((degree,expected))

# The highest-t identity underlying the unbounded descent lemma.
for n in range(1,7):
    gn=s.Function('g')(x)
    leading=s.Poly(s.expand(jac(F,gn*t**n)),t).coeff_monomial(t**n)
    zero(leading-(n*s.diff(A,x)*gn-A*s.diff(gn,x)),'highest-t coefficient identity')
    zero(n*s.diff(A,x)*A**n-A*s.diff(A**n,x),'all kernel leading terms')
    zero(jac(F,c*F**n),'exact removable full power, not just leading monomial')

# Bounded independent affine-plane equation solver, without descent pruning.
solved=[]
cases=[('nonconstant_simple',1-x**4,-x*x,False),
       ('critical_free_quartic',(x-1)**4,x*x-4*x,False),
       ('critical_free_cubic',x**3,x,False),
       ('linear_A',x,1,False),
       ('constant_A',s.Integer(2),x**3+1,True)]
for name,Av,Bv,success in cases:
    first=Av*t+Bv
    for max_t in (1,2,3):
        mons=[x**i*t**j0 for j0 in range(max_t+1) for i in range(7)]
        coeff=s.symbols('k0:'+str(len(mons)))
        trial=sum(c0*m for c0,m in zip(coeff,mons))
        eq=s.Poly(s.expand(jac(first,trial)-1),x,t)
        M,rhs=s.linear_eq_to_matrix(eq.coeffs(),coeff)
        attainable=M.rank()==M.row_join(rhs).rank()
        need(attainable==success,'literal full finite polynomial-mate coefficient system')
        solved.append((name,max_t,attainable))
    for n in range(1,5):
        gpoly=first**n+sum(x**i for i in range(1,5))
        zero(jac(first,gpoly)+Av*sum(i*x**(i-1) for i in range(1,5)),
             'exact image intersection positive control')

# The remaining constant-A affine mates fail global regularity.
g=-x/a+(a*t+c)**4
zero(jac(a*t+c,g)-1,'positive polynomial affine-plane mate')
zero(s.limit(r*chart(g),r,0)+1/a,'unremovable global boundary pole')

# Mixed t+lambda*b: all four critical points, without choosing radicals.
lam=s.symbols('lam',nonzero=True)
mixed=t+lam*bb
zero(s.diff(mixed,x).subs(t,-1/(2*x*x)),'mixed pencil critical x derivative')
zero(s.diff(mixed,t)-(1-lam*x**4),'mixed pencil critical t derivative')
zero(mixed.subs(t,-1/(2*x*x))+lam*x*x-(1-lam*x**4)*(-1/(2*x*x)),
     'mixed pencil exact critical value reduction')

# Whole globally smooth submersion family, retaining repeated-root boundary.
AA=a*(x-h)**2*(x-k)**2
BB=a*(x*x-2*(h+k)*x)+c
sub=AA*t+BB
zero(s.diff(sub,x).subs(x,h)+2*a*k,'first repeated A root derivative')
zero(s.diff(sub,x).subs(x,k)+2*a*h,'second repeated A root derivative')
need(s.Poly(chart(sub),r,b).is_multivariate,'second chart is actually polynomial')
zero(s.diff(chart(sub),b).subs(r,0)+a,'no critical point on actual D')
res_h=s.residue(-1/AA,x,h)
zero(res_h-2/(a*(h-k)**3),'nonzero rational-derivative residue at distinct roots')

# Fourfold case admits a rational conjugate but not a regular one.
four=s.expand(sub.subs(k,h))
g_rat=1/(3*a*(x-h)**3)
zero(jac(four,g_rat)-1,'exact generic rational relative primitive')
level=c-3*a*h*h
factor=a*(x-h)*((x-h)**3*t+x-3*h)
zero(four-level-factor,'complete special fibre with both components')
zero(((x-h)**3*t+x-3*h).subs(x,h)+2*h,'special components disjoint when h nonzero')
zero(chart(four).subs(r,0)-(c-6*a*h*h-a*b),'boundary separator remains nonconstant')
zero(s.limit(chart(g_rat),r,0),'rational primitive itself regular along D')

# The second special component is a full A1 including its D point.
# Parameter z=1/(x-h), with source expression t=2h*z^3-z^2.
T=2*h*z**3-z*z
X=h+1/z
zero(four.subs({x:X,t:T},simultaneous=True)-level,'literal nonvertical special component')
R=z/(1+h*z)
BV=s.cancel((-X*X-X**4*T))
need(s.denom(BV)==1,'boundary coordinate polynomial on special component')
zero(BV.subs(z,0)+3*h*h,'retained D point on second special component')
zero(g_rat.subs(x,X)-z**3/(3*a),'no primitive pole on second component')

# Fixed smallest clean control h=k=a=1,c=0.
hostile=(x-1)**4*t+x*x-4*x
hostile_inf=-6+4*r-r*r-b*(1-r)**4
zero(chart(hostile)-hostile_inf,'literal globally critical-free hostile')
zero(s.diff(hostile,x).subs(x,1)+2,'nonvanishing gradient at its only A root')
zero(s.diff(hostile_inf,b).subs(r,0)+1,'nonvanishing gradient on all D')
zero(jac(hostile,1/(3*(x-1)**3))-1,'literal rational-but-not-polynomial hostile')

print('FULL GLOBAL LAYER: O(W) intersect {deg_t<=1}=span(1,t,j,u,v,b).')
print('UNBOUNDED DESCENT (A!=0): J(F,K[x,t]) intersect K[x]=A(x)K[x]; ker J(F,-)=K[F].')
print('CONSEQUENCE: no nonconstant F in that global layer has a global regular nonzero-constant-Jacobian mate.')
print('HOSTILE: F=(x-1)^4*t+x^2-4x is globally critical-free and separates D.')
print('RATIONAL MATE: 1/(3*(x-1)^3); cancelling its vertical pole forces a pole on the other F=-3 component.')
print('DISTINCT-ROOT SUBMERSIONS: residue 2/(a*(h-k)^3) excludes even rational mates.')
print('LAURENT LAYER DIMENSIONS',layer)
print('LITERAL FINITE COEFFICIENT SYSTEMS',solved)
print('PASS',G,'always-active exact gates; no arbitrary global coordinate or Keller conclusion.')
