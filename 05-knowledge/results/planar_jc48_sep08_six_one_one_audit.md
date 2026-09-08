# Independent audit: the complete (6,1,1) boundary class

**Status: PASS — complete analytic proof, source inspection, independent
symbolic reconstruction and normal/optimized/frozen replays accepted.**
Auditor: root; producer: certificate_audit. Audit performed September8,
2026 in the isolated planar-jacobian-sep06 session worktree. No JC(2)
closure or external-priority claim follows.

## 1. Accepted scope and immutable inputs

The [primary proof](planar_jc48_sep08_six_one_one.md) excludes every
polynomial mate in the full binary(6,1,1) DG square-prefix class. For
all-finite leading roots, every rational mate forces L constant.
Both infinity placements already fail the proved leading differential
gate. Explicit all-p constant-L rational examples pay sharpness.

The prepromotion primary is17,752 bytes with SHA256
`9a613f4e1c488649143cdd2cd38935f8ac1bab182f1a21ed61de96f355a81cb3`.
The source is12,727 bytes with SHA256
`8448a49b7a52a2b3a35c13c8666425144781c76ed2343eccd21add012b96c196`.
The frozen output is772 bytes with SHA256
`ecd81171fe3f0dbcb353eecb240b877a0e00e6488d279971dd63efaafe1a907d`.
Its63 semantic gates have digest
`8faee9ffd19ef9d77353e8d6f3156d7ebd2d1831ea948eea08bef2c7e397b7a9`.
Parent promotion may change the primary's status and audit routing;
the audited formulas, scope, source and output must remain unchanged.

## 2. Quantifiers, complete rows and local entry

The theorem quantifies over every original finite position p and every
lower coefficient in the full global section filtration, not a bounded
coefficient sample or a bound on the unknown mate. The degree-eight
leading exactness supplier gives3s²=4h, with h and the leading scalar
nonzero. Hence s!=0, and the actual scaling x=lambda X,z=lambda²Z
validly sets s=2,h=3. Scalar target changes preserve nonzero constant
Jacobians. Writing u=x-p is then an ordinary source coordinate change;
the full shifted numerator rows explicitly keep p. No projective
translation of the compactification is assumed.

I checked the coefficient count and the direct numerator equations.
After the necessary active jets, the free H coordinates are exactly
P2,P3,P4,Q0. Their affine coefficient map from a,c,b,k has determinant
one; the lower constant is independently free. The apparent p/u terms
in the q carrier cancel in the original polynomial H,L. This proves
completeness as well as global membership. The source checks every
original x row rather than only the new rational carrier.

The two simple-root residues of M/N force C|M. The sixfold M-unit
case retains its balanced j=2 exception: it can have second-kind poles,
but its total primitive pole budget is at most two, by the cited
complete Morse bound including contact four. The actual D zero gives
local primitive degree three, contradicting that total budget on the
same component. This is componentwise and needs no generic geometric
integrality. A shared normal unit is regular; the separate finite
first-jet supplier then forces P=P'=M=M'=0 at the sixfold point.
Together these imply M=kappa*u²*C, with kappa=0 exactly the constant-L
boundary. The proof never copies the false assertion that all sixfold
M-unit forms are regular.

For a=P2!=0, I independently checked the moving-centre calculation.
The numerator zero is an analytic Zbar(u), not the fixed root Z0.
With s=u⁴[Zbar(u)+uY], the leading equation is
`a²Y²+m2 Z0³=0`, with two nonzero simple roots. The numerator of eta
has u-order8, its denominator order9 and leading coefficient2a²Y0.
Thus the actual residue is Z0²/(2a²Y0), nonzero. Arbitrary first-unit
coefficients merely change Zbar. The producer's source retains n7,p3
in an independent local check. This pays a=0 before any algebraic
elimination or parameter division.

## 3. Genuine quadratic field and exact trace

After a=0, the explicit quadratic E2 in u recovers the original q,
so the extension degree is at most two. It is genuinely degree two:
the leading f² coefficient of its discriminant simplifies to
`kappa²(c²+12w-12k)`. This is nonsquare in C(w). A polynomial in f
that is a square in C(w)(f) must have a square leading coefficient.
The source-to-(F,H) map is dominant because its actual Jacobian has
nonzero q² coefficient kappa*u*C*(2u-6). Hence the quadratic equation
is a field equation, not an unproved use of two geometric branches.

The actual volume is du wedge dq/u². I reconstructed the transport
independently in the original (u,t) coordinates and verified
`J(F,H)=-u² E2_u/kappa` after substituting f=F,w=H. Before restriction
the producer correctly retains the multiple C' E2/(kappa C); it is
not silently deleted in an off-fibre identity. Therefore
eta=-kappa dw/(u² E2_u) is the correct relative form.

For Y=2Au+B, roots u=(-B±Y)/(2A), I independently computed the
normalized trace of1/(u²Y) as B/(2C0²). The sign and factor two in
the trace differential agree. Since trace commutes with the derivation
holding f fixed, its rational residue must vanish. C0 is a genuinely
linear polynomial in w with derivative -3kappa². Its double-pole
residue is zero exactly when B'(w0)=0, giving
`(c-12p)(k+12p²-2pc)=kappa`. If B(w0)=0 there can still be a simple
pole, so cancellation does not discard this condition. In particular
c-12p is nonzero before it is divided out.

## 4. All inverse coefficients and every denominator branch

I checked the formal exactness identity directly: on F(x,T)=v^-4,
`d_x G=(v⁵/4)d_v T` follows from J(F,G)=1 and F_t T_v=-4v^-5.
Thus each positive T_j has exact differential. A rational denominator
embeds into the Laurent field because its top t coefficient gives a
unique lowest v power; no mate degree is bounded. The normalized
quadratic trace gives T2's rational exactness, and the odd residues
are ordinary residues at the two unramified points above u=0.

The formal residue formula
`T_j=-(1/j)[y^-1]F^(j/4)` has the correct sign and factor. It follows
by integrating the Laurent residue of d(T*v^-j); centering t changes
no y^-1 residue. I read both universal implementations: complete
multinomial extraction and recursive inverse solving through T7.
They use independent unspecialized coefficients and agree before any
branch condition. The source also reconstructs every actual-family
residue from those formulas, separately from its collected expressions.

My additional independent code below does not import the producer or
copy the five displayed T_j formulas. It takes the formal residue of
the centered whole quartic directly, then computes coefficients of
C^e from its differential recurrence instead of calling sympy.series.
It reproduces T1's conic, T3's exact nonzero factor4, and the T5/T7
numerator factors -d² and d² after the prescribed coefficient ratios.
This is a separate reconstruction of the actual necessary equations.

The target translation making the constant part of Z0 vanish is
legitimate and does not change H or the trace condition. T1 gives
3b²+6bc+c²=0; T3 then gives c³(6b+c)=12kappa(c-12p).
Thus c!=0 is proved. The b=c=0 branch is not silently divided away.
The homogeneities are exact polynomial identities with weights
1 on b,c,p,2 on k,3 on kappa. Passing to b/c,p/c,k/c²,kappa/c³ is
an algebraic change of ratios in necessary identities, not a surface
automorphism or a claim of preserved polynomiality of a rational chart.
The only subsequent parameter denominator is d=1-12p, already nonzero.
No factor6b+1, no leading coefficient of a resultant, and no unproved
coefficient is inverted.

## 5. Coprimality, infinity locations and sharp boundary

I independently recover both exact polynomials E5,E7. Their resultants
are3R5 and -6R57 with the stated orientations. Only the implication
from a common b root to zero resultant is used. This implication
remains valid when leading coefficients specialize to zero; no
converse lifting is claimed. The two quartics R5,R57 have an exact
extended Euclidean identity A(V)R5+B(V)R57=1 with coefficient degrees
at most three, independently checked below. This complements the
producer's standard-library Fraction primitive-remainder sequence
ending at1. The proof requires no numerical root separation.

The infinity placements are genuinely closed by their original
finite leading polynomials. Degree two with two simple roots gives
nonzero infinity residues. Degree seven of type6+1 has primitive pole
capacity four at the finite sixfold point, but its unique infinity
point forces local primitive degree five. These are the relevant
entries of the audited radical classification, not a translation of
a finite boundary lemma to a weighted point.

For constant L the factor2H excludes a polynomial mate directly in
the original polynomial ring. The all-p rational positive is exact:
H=C(u+u³t-2p)²+k and G_H=(u+1)/[12u²(u+u³t-2p)]. I recomputed its
Jacobian in original (u,t) coordinates; the producer also verifies
the entire original second chart is polynomial. The two simple
roots are distinct and disjoint from the sixfold point. This is an
actual example in the asserted class, with essential rational poles.

## 6. Replay and independent reconstruction

I inspected the full source, including both universal paths, every
family residue reconstruction, all homogeneous substitutions and
the Fraction polynomial division. Gates use explicit exceptions, so
optimization cannot erase them. Independent invocations of

```sh
python3 -B 04-computation/planar_jc48_sep08_six_one_one.py
python3 -B -O 04-computation/planar_jc48_sep08_six_one_one.py
```

both reproduce the772 frozen bytes and all63 semantic gates exactly.
My separate script below passes81 additional direct-formal and
original-source checks. Its output is:

```text
T3 nonzero scalar: 4
T5 exact scalar: -d**2
T7 exact scalar: d**2
Independent direct-formal and original-source gates: 81
Semantic SHA256: 4677c189a8aacae87ce4991aa0de2fa70f69f36bc7c25f8750ab9906bdc8139b
```

These81 checks are independent audit controls, not an extra producer
program in the session aggregate. Save this literal script and run
with Python/SymPy to reproduce it. No correction remains pending.

```python
from pathlib import Path
from hashlib import sha256
import sympy as S
from math import factorial
u,W,b,c,k,p,K,d,V=S.symbols('u W b c k p K d V')
gates=[]
def need(label,v):
 if not bool(v):raise RuntimeError(label)
 gates.append(label)
def zero(label,v):need(label,S.cancel(v)==0)
# Direct formal residue of the full centered quartic; no producer import,
# no copied T1/T3/T5/T7 formulas and no sympy series call.
N=u**6*W**2;D=u**6*((b*u+c)**2-4*k*W**2)
M=K*u*u*W**2;Z=-K*(c-12*p)/(2*u)
A=-D/(2*N*N);B=M/(N*N);C=D*D/(16*N**4)+Z/(N*N)
def power_coefficient(exponent,n):
 g=[S.Integer(1)]
 for j in range(n):
  g.append(S.cancel(((2*j-2*exponent)*g[j]+(2*exponent-j+1)*(g[j-1] if j else 0))/(3*(j+1))))
 return g[n]
res={}
for n in (1,3,5,7):
 raw=0
 for i in range((n+1)//2+1):
  for j in range((n+1)//3+1):
   for h in range((n+1)//4+1):
    if 2*i+3*j+4*h!=n+1:continue
    total=i+j+h
    mult=S.binomial(S.Rational(n,4),total)*S.Rational(factorial(total),factorial(i)*factorial(j)*factorial(h))
    raw-=u**(3*n)*W**n*mult*A**i*B**j*C**h/n
 value=0
 for term in S.Add.make_args(S.expand(raw)):
  powers=term.as_powers_dict();eu=int(powers.get(u,0));ew=int(powers.get(W,0));idx=-eu-1
  if idx<0:continue
  need('odd radical parity row'+str(n)+' '+str(len(gates)),ew%2==1)
  value+=S.cancel(term/u**eu/W**ew)*3**S.Rational(ew+1,2)*power_coefficient(S.Rational(ew,2),idx)
 res[n]=S.factor(value)
quad=3*b*b+6*b*c+c*c
zero('direct first residue',res[1]-quad/72)
red3=S.rem(S.together(res[3]).as_numer_denom()[0],quad,b)
ratio=S.factor(red3/(c**3*(6*b+c)-12*K*(c-12*p)))
need('direct third residue only constant factor',ratio.is_Rational and ratio!=0)
print('T3 nonzero scalar:',ratio)
Qb=3*b*b+6*b+1
E5=(96*b+16)*V*V+(1296*b+240)*V+900*b+165
E7=(96*b+16)*V**4+(1296*b+240)*V**3+(180*b+33)*V*V+(14952*b+2744)*V+10098*b+1853
sub={c:1,p:(1-d)/12,K:(6*b+1)/(12*d),k:(6*b+1)/(12*d*d)+(1-d*d)/12}
for n,target in ((5,E5),(7,E7)):
 num=S.together(res[n].subs(sub)).as_numer_denom()[0]
 remainder=S.rem(num,Qb,b)
 ratio=S.factor(remainder/target.subs(V,d*d))
 need('higher independent residue factor no hidden parameter'+str(n),not ratio.has(b))
 nu,de=S.fraction(ratio)
 need('higher factor pure monomial'+str(n),len(S.Poly(nu,d).terms())==len(S.Poly(de,d).terms())==1)
 print('T'+str(n)+' exact scalar:',ratio)
R5=256*V**4+3072*V**3-2208*V*V-2880*V+225
R57=192*V**4-4320*V**3+2656*V*V+3252*V-255
zero('resultant necessary identity1',S.resultant(Qb,E5,b)-3*R5)
zero('resultant necessary identity2',S.resultant(E5,E7,b)+6*R57)
bez_a,bez_b,gcd=S.gcdex(R5,R57,V)
zero('explicit Bezout constant',gcd-1)
zero('independent extended Euclidean identity',bez_a*R5+bez_b*R57-1)
need('Bezout coefficient degrees',S.degree(bez_a,V)<=3 and S.degree(bez_b,V)<=3)
# Original field transport, now using literal original t coordinates.
t,f,w,ell=S.symbols('t f w ell')
Cu=u*u-2*u+3;q=1+u*u*t-2*p/u
H=u*u*Cu*q*q+(b*u*u+c*u)*q+k
L=K*Cu*q+ell+6*K*p/u;F=H*H+L
z=f-ell-w*w
E=(u*z-6*K*p)**2+K*(b*u+c)*(u*z-6*K*p)-K*K*Cu*(w-k)
jac=S.diff(F,u)*S.diff(H,t)-S.diff(F,t)*S.diff(H,u)
zero('original t-coordinate exact fibre transport',
 jac+u*u*S.diff(E,u).subs({f:F,w:H},simultaneous=True)/K)
AA,BB,CC,Y=S.symbols('AA BB CC Y')
x=(-BB+Y)/(2*AA)
trace=S.cancel((1/(x*x*Y)+1/(((-BB-Y)/(2*AA))**2*(-Y)))/2)
zero('quadratic normalized trace',trace.subs(Y*Y,BB*BB-4*AA*CC)-BB/(2*CC*CC))
# Sharp control verified without q-coordinate transport.
v=u+u**3*t-2*p;Hp=Cu*v*v+k;Gp=(u+1)/(12*u*u*v)
zero('literal sharp H Jacobian',S.diff(Hp,u)*S.diff(Gp,t)-S.diff(Hp,t)*S.diff(Gp,u)-1)
print('Independent direct-formal and original-source gates:',len(gates))
print('Semantic SHA256:',sha256('\n'.join(gates).encode()).hexdigest())
```
