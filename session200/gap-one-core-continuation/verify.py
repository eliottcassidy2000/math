"""Author-only controls for the separately unreviewed gap-one continuation."""
from fractions import Fraction
from pathlib import Path
import hashlib,json
T=(1,0,2,1)
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def step(s,c):
    D,E,U,d=s
    if c=='O':return D,0,2*U,d
    if c=='Z':return D-E+U,d,2*U-d,d
    if c=='J':return D-E+U,D-E,D+U,D
    raise ValueError(c)
def state(w,s=T):
    for c in w:s=step(s,c)
    return s
def word(t):
    need(t>=2,'nontrivial unit source')
    if t==2:return ''
    return word(t//2)+'O' if t%2==0 else word((t+1)//2)+'Z'
def g(t):
    if t==1:return 1
    if t%2==0:return g(t//2)
    n=(t+1)//2;return g(n)+n-(n%2)
def F(s):return s[0]-s[1]+s[2]
checks=0;examples=[]
for t in range(3,2002,2):
    v=state(word(t));D=g(t)
    need(v==(D,1,t,1),'all-odd-cover physical compiler')
    need(2*D>=t+3,'strong-unit source lower bound')
    q=step(v,'J')
    need((q[0],q[2],q[3])==(D+t-1,D+t,D),'gap-one physical core')
    for a in range(7):
        R=2**(a+1)
        for k in range(7):
            M=2**(k+1)-1;A=M*(R-1)+2;B=M*R+1
            n=A*D+B*t-1
            need(F(state('J'+'O'*(a+1)+'Z'*k,v))==n,'target recipe')
            if k==0:need(n==(R+1)*(D+t)-1,'offset return simplification')
            checks+=1
    if t in (3,5,7,9,11):examples.append(dict(t=t,word=word(t),core=[q[0],q[2],q[3]],a0k0=3*(D+t)-1))
# Formula is affine in D,t on the formal slice E=d=1.
affine=0
for D,t in ((0,0),(1,0),(0,1)):
    for a in range(9):
        for k in range(9):
            R=2**(a+1);M=2**(k+1)-1
            need(F(state('J'+'O'*(a+1)+'Z'*k,(D,1,t,1)))==(M*(R-1)+2)*D+(M*R+1)*t-1,'affine identity')
            affine+=1
# Coarse geometric tail: denominator >=5*2^(a+k).
L=6
partial=sum((Fraction(1,(2**(k+1)-1)*(3*2**(a+1)-1)+4)for a in range(L)for k in range(L)),Fraction())
tail=Fraction(8,5*2**L)
need(partial+tail<Fraction(1,2),'strict rational density upper bound')
short_bound=sum((Fraction(1,d)for d in (9,15,19,37)),Fraction())+Fraction(10,69)+Fraction(32,385)
need(short_bound==Fraction(27199346,56025585),'short rational bound evaluation')
need(Fraction(17,35)-short_bound==Fraction(13081,56025585)>0,'short bound strict gap')
print(json.dumps(dict(status='AUTHOR_CHECK_PASS',scope='UNREVIEWED continuation, not part of the core audit',odd_cover_max=2001,physical_formula_checks=checks,affine_basis_checks=affine,examples=examples,short_density_bound=str(short_bound),short_density_ceiling='17/35',density_partial_square=L,density_partial=str(partial),coarse_tail_bound=str(tail),density_upper=str(partial+tail),density_upper_decimal=float(partial+tail),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
