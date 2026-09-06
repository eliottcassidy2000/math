"""Independent residue-moment identities and original-phase tail audit.
No producer imports; matrices checked from polynomial products and literal
partial fractions. Complete tail rebuilt from the original Laurent rows.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb,prod
import sympy as S
import sys,json,hashlib,argparse
from pathlib import Path
sys.stdout.reconfigure(newline='\n')
N=0

def check(ok,label):
    global N
    N+=1
    if not ok:raise ArithmeticError(label)

x,y,z,d,w,u,t,v=S.symbols('x y z d w u t v')
X=75+d
Cmom=[S.Integer(1),S.Integer(1),S.Integer(3),9+d/3,27+16*d/3-4*y/7,81+54*d-59*y/7+z]
Dmom=[S.Integer(1),S.Integer(2),S.Integer(7),S.Rational(99,4)+7*d/12,S.Rational(347,4)+115*d/12-6*y/7,S.Rational(583,2)+199*d/2-92*y/7+z]
den=1-13*w+55*w*w-X*w**3+y*w**4-z*w**5
numC=1-12*w+45*w*w-S.Rational(2,3)*X*w**3+S.Rational(3,7)*y*w**4
numD=1-11*w+36*w*w-S.Rational(5,12)*X*w**3+S.Rational(1,7)*y*w**4
for mu,num in ((Cmom,numC),(Dmom,numD)):
    product=S.Poly(S.expand(den*sum(a*w**i for i,a in enumerate(mu))-num),w)
    for i in range(6):check(product.nth(i)==0,'formal numerator/denominator moment identity')
HC=S.Matrix([[Cmom[i+j+1] for j in range(3)] for i in range(3)])
HD=S.Matrix([[Dmom[i+j] for j in range(3)] for i in range(3)])
expectC=-X*d*d/27+S.Rational(15,7)*d*y-S.Rational(16,49)*y*y+d*z/3
expectD=(-333+2334*d-49*d*d)/144-S.Rational(18,7)*y
check(S.expand(HC.det()-expectC)==0,'complete shifted C determinant')
check(S.expand(HD.det()-expectD)==0,'complete ordinary D determinant')
check(S.expand(Cmom[1]*Cmom[3]-Cmom[2]**2-d/3)==0,'shifted C order2 gives nonnegative delta')
check(expectC.subs(d,0)==-S.Rational(16,49)*y*y,'zero delta forces zero fourth coefficient')
check(expectD.subs({d:0,y:0})==-S.Rational(37,16),'zero delta rejected by D')
# The anchors imply5M^2-26M-67<=0; derivative positive forM>=13/5.
check(5*F(71,10)**2-26*F(71,10)-67==F(9,20)>0,'strict support endpoint')
check(10*F(13,5)-26==0,'endpoint quadratic monotone above mean')
K=F(25,9)/(F(15,7)+F(71,150))
check(K==F(8750,8241),'strict C slope floor')
coefficient=F(2334)/K-F(2592,7)
check(coefficient>0 and F(333)/coefficient==F(161875,888583),'combined exact fourth-coefficient floor')
check(F(161875,888583)>F(9,50),'coarse floor used for the tail')

# A strict genuine root shape: reconstruct both rational functions from
# literal positive residues, then all six moments directly from those atoms.
roots=tuple(S.Rational(j,5) for j in (1,3,9,22,30))
B=S.prod(v-a for a in roots)
poly=S.Poly(B,v);x0=-poly.nth(2);y0=poly.nth(1);z0=-poly.nth(0)
C=v**4-12*v**3+45*v*v-S.Rational(2,3)*x0*v+S.Rational(3,7)*y0
D=v**4-11*v**3+36*v*v-S.Rational(5,12)*x0*v+S.Rational(1,7)*y0
for row,mu in ((C,Cmom),(D,Dmom)):
    weights=[row.subs(v,a)/prod(a-b for b in roots if b!=a) for a in roots]
    check(all(a>0 for a in weights) and sum(weights)==1,'literal positive probability residues')
    reconstructed=sum(weights[i]*prod(v-b for j,b in enumerate(roots) if j!=i) for i in range(5))
    check(S.expand(row-reconstructed)==0,'complete literal partial-fraction identity')
    for k in range(6):
        mk=sum(weights[i]*a**k for i,a in enumerate(roots))
        check(mk==mu[k].subs({d:x0-75,y:y0,z:z0}),'atom moments equal formal coefficient moments')
    for k in range(1,5):check(mu[k+1].subs({d:x0-75,y:y0,z:z0})<=S.Rational(71,10)*mu[k].subs({d:x0-75,y:y0,z:z0}),'bounded support moment inequalities')

# The old coefficient hostile passes ordinary3 and shifted2 PSD, but fails
# the bounded-support endpoint test; checking just determinants is insufficient.
host={d:S.Integer(29),y:S.Integer(50),z:S.Rational(37435088,3898125)}
for label,mu,loss in (('C',Cmom,S.Rational(2159,105)),('D',Dmom,S.Rational(1091,42))):
    for n,shift in ((3,0),(2,1)):
        H=S.Matrix([[mu[i+j+shift].subs(host) for j in range(n)] for i in range(n)])
        for size in range(1,n+1):
            for ids in combinations(range(n),size):check(H.extract(ids,ids).det()>=0,'all hostile low-PSD principal minors')
    check(S.simplify((mu[4]-S.Rational(71,10)*mu[3]).subs(host))==loss>0,'hostile bounded-support fourth moment excess')
    check(S.Matrix([[mu[i+j+1].subs(host) for j in range(3)] for i in range(3)]).det()<0,'hostile shifted3 already rejects at order5')

# Build the ACTUAL E/O convolution and all carried rows, then eliminate z
# by the original first equation. Do not import the prior envelope/table.
O=sum(comb(14,2*j+1)*t**j for j in range(7))
E=sum(comb(14,2*j)*t**j for j in range(8))
mult=S.expand(O*O+E*E/t)
GB=1+13*t+55*t*t+x*t**3+y*t**4+z*t**5
GC=1+12*t+45*t*t+S.Rational(2,3)*x*t**3+S.Rational(3,7)*y*t**4
GD=1+11*t+36*t*t+S.Rational(5,12)*x*t**3+S.Rational(1,7)*y*t**4
mix=S.expand(GB*GB/t**2+2*GC*GD/t)
qs={j:S.expand(mult.coeff(t,j)*mix.coeff(t,j)) for j in range(-1,9)}
check(qs[-1]==28,'full original lower carry')
P=182-20020/u+2002*x/u**2-3432*y/u**3+2002*z/u**4
replacement=S.Rational(12,7)*y*u-x*u*u+10*u**3-u**4/11
check(S.expand(P.subs(z,replacement))==0,'exact original-root elimination')
R=S.Poly(S.expand(sum(S.Integer(-1)**j*qs[j].subs(z,replacement)*u**(7-j) for j in qs)),x,y,u)
check(R.degree(u)==8 and R.coeff_monomial(y*y)==-S.Rational(26075790,7),'complete eliminated degree and negative constant')
# Positive-term envelope: nonnegative x<123, y>9/50. Every positive term
# has y degree at most 2, so dividing by y^2 is safely bounded at the floor.
floor=S.Rational(9,50);envelope=-S.Rational(26075790,7);positives=0
for (a,b,j),coefficient in R.terms():
    if coefficient<=0:continue
    check(j>0 and b<=2,'every positive monomial fits floor domination')
    envelope+=coefficient*123**a*floor**(b-2)*u**j;positives+=1
envelope=S.Poly(envelope,u)
check(all(coefficient>0 for (j,),coefficient in envelope.terms() if j>0),'tail envelope monotone for positive u')
value=envelope.eval(S.Rational(1,4100))
check(value<-59000,'uniform strict4100 tail endpoint')
# The original-root condition is load-bearing even at large phase.
raw_at4100=sum(S.Integer(-4100)**j*qs[j].subs({x:x0,y:y0,z:z0}) for j in qs)
first_at4100=P.subs({u:S.Rational(1,4100),x:x0,y:y0,z:z0})
check(raw_at4100>0 and first_at4100!=0,'large-phase positivity away from original root')
# Preserve the already-proved correlation delta/y instead of replacing x
# by its independent coefficient cap. This improves the same envelope to2500.
rho=1/floor
ratio1=75*rho+1/S.Rational(8750,8241)
ratio2=ratio1*rho
check(ratio1==S.Rational(10962223,26250) and ratio2==S.Rational(10962223,4725),'correlated x/y and x/y^2 bounds')
coupled=-S.Rational(26075790,7)
for (a,b,j),coefficient in R.terms():
    if coefficient<=0:continue
    check((a,b) in ((1,1),(1,0),(0,1),(0,0)),'coupled envelope covers every positive monomial')
    bound={(1,1):ratio1,(1,0):ratio2,(0,1):rho,(0,0):rho*rho}[(a,b)]
    coupled+=coefficient*bound*u**j
coupled=S.Poly(coupled,u)
check(all(coefficient>0 for (j,),coefficient in coupled.terms() if j>0),'coupled envelope monotone')
coupled_value=coupled.eval(S.Rational(1,2500))
check(coupled_value==-S.Rational(653528391359305169367452997041401,13323669433593750000000000000)<-49000,'strict2500 coupled tail endpoint')
check(sum(S.Integer(-2500)**j*qs[j].subs({x:x0,y:y0,z:z0}) for j in qs)>0
      and P.subs({u:S.Rational(1,2500),x:x0,y:y0,z:z0})!=0,'2500 root-condition omission remains hostile')

# A repeated-root weak C boundary is a genuine positive residue measure after
# cancellation, so the passage to weak interlacing cannot assume five atoms.
B0=v**2*(v-3)*(v-5)**2;C0=v*(v-2)*(v-5)**2
check(S.cancel(C0/B0-S.Rational(2,3)/v-S.Rational(1,3)/(v-3))==0,'weak repeated-root canceled probability measure')
for k in range(6):
    expected=S.Integer(1) if k==0 else S.Integer(3)**k/3
    check(Cmom[k].subs({d:0,y:0,z:0})==expected,'weak-boundary moments including zero atom')

# Pin and independently compare every saved identity in the producer's full
# JSON, treating it solely as data; all expressions above were derived first.
parser=argparse.ArgumentParser()
parser.add_argument('--certificate',type=Path,default=Path(__file__).with_name('continuing3_20260906_residue_floor_certificate.json'))
certificate=parser.parse_args().certificate.read_bytes()
check(hashlib.sha256(certificate).hexdigest()=='c2a20f53670198790529fa980b7d33913b9e5a1f64c3b206f9fc34aa42e9104e','full producer certificate pin')
bank=json.loads(certificate)
def unpack(packet,variables):
    return S.expand(sum(S.Rational(coefficient)*prod(vv**ee for vv,ee in zip(variables,exponents)) for exponents,coefficient in packet))
for name,mu in (('C_moments',Cmom),('D_moments',Dmom)):
    for packet,moment in zip(bank[name],mu):
        check(S.expand(unpack(packet,(x,y,z))-moment.subs(d,x-75))==0,'saved full moment polynomial matches independent identity')
check(S.expand(unpack(bank['C_shifted3'],(x,y,z))-expectC.subs(d,x-75))==0,'saved C determinant identity')
check(S.expand(unpack(bank['D_ordinary3'],(x,y,z))-expectD.subs(d,x-75))==0,'saved D determinant identity')
check(S.expand(unpack(bank['original_root_response'],(x,y,u))-R.as_expr())==0,'saved entire eliminated carried response')
json_envelope=-S.Rational(26075790,7)+sum(S.Rational(c)*u**k for k,c in bank['tail_envelope'])
check(S.expand(json_envelope-envelope.as_expr())==0,'saved entire positive tail envelope')
json_coupled=-S.Rational(26075790,7)+sum(S.Rational(c)*u**k for k,c in bank['coupled_tail_envelope'])
check(S.expand(json_coupled-coupled.as_expr())==0,'saved entire coupled tail envelope')
check(tuple(map(S.Rational,bank['coupled_ratios']))==(ratio1,ratio2),'saved coupled ratio constants')

print('PROOF_AUDIT: positive residue support -> C shifted3 and D ordinary3 -> e4>161875/888583>9/50.')
print('EXACT_SLOPE',K,'EXACT_FLOOR',F(161875,888583))
print('HOSTILE_SUPPORT_EXCESSES C',F(2159,105),'D',F(1091,42),'while ordinary3/shifted2 PSD pass')
print('TAIL_POSITIVE_MONOMIALS',positives,'ENVELOPE_AT_1/4100',value)
print('COUPLED_RATIOS',ratio1,ratio2,'ENVELOPE_AT_1/2500',coupled_value)
print('SCOPE: Q(-s)<0 for original positive phases s>=2500 in the retained two-anchor/two-interlacer model; remaining phases open.')
print('PASS',N,'always-active independent exact gates')
