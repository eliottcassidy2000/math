#!/usr/bin/env python3
"""Exact combined Euler/pencil model hostile; no actual trinomial refutation.

Certificate universe: q5,h4 geometric-root model and its dyadic and
genuine-linear-coefficient scalings;
one actual q2,h1 control; the previous q6 Euler hostile as a rejection
control. No imported producer, floats, root finder, or disabled assertions.
"""
from fractions import Fraction as F
from math import comb,gcd,lcm
from functools import reduce
import hashlib,json
GATES=0

def need(ok,msg):
 global GATES
 GATES+=1
 if not ok:raise RuntimeError(msg)

def trim(p):
 p=list(map(F,p))
 while len(p)>1 and p[-1]==0:p.pop()
 return p

def cv(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c)

def add(*rows):
 return trim([sum(row[i] if i<len(row) else F(0) for row in rows) for i in range(max(map(len,rows)))])

def mul(p,c):return trim([c*x for x in p])
def coeff(p,j):return p[j] if 0<=j<len(p) else F(0)
def at(p,x):
 a=F(0)
 for v in reversed(p):a=a*x+v
 return a

def rem(a,b):
 a=trim(a)
 while len(a)>=len(b) and any(a):
  e=len(a)-len(b);c=a[-1]/b[-1]
  for j,v in enumerate(b):a[e+j]-=c*v
  a=trim(a)
 return a

def interval(p,I):
 lo=hi=F(0);a,b=I
 for c in reversed(p):
  products=[lo*a,lo*b,hi*a,hi*b]
  lo=min(products)+c;hi=max(products)+c
 return lo,hi

def sign_on(p,I,sgn,label):
 lo,hi=interval(p,I)
 need(lo>0 if sgn>0 else hi<0,label)
 return [str(lo),str(hi)]

def kernels(a):
 q=len(a);b=[F(1)]
 for z in a:b=cv(b,[-z,1])
 c=[F(3*(i+1),q+2+2*i)*b[i+1] for i in range(q)]
 d=[F(2+3*i,q+1+2*i)*c[i] for i in range(q)]
 return b,c,d

def interlace_certificate(a,b,c,d):
 bp=[i*b[i] for i in range(1,len(b))];ratios=[]
 for r in a:
  need(at(b,r)==0,'literal B root')
  rc=at(c,r)/at(bp,r);rd=at(d,r)/at(bp,r)
  need(rc>0 and rd>0,'both strict interlacers')
  ratios.append([str(rc),str(rd)])
 return ratios

def model(b,c,d,h):
 q=len(b)-1;x=q-h;g=q+2*h+1;k=2*h+1
 G=[b[q-j]*(-1)**j for j in range(q+1)]
 C=[c[q-1-j]*(-1)**j for j in range(q)]
 D=[d[q-1-j]*(-1)**j for j in range(q)]
 need(all(v>0 for p in (G,C,D) for v in p),'full nonnegative path-type support')
 O=[F(comb(g,2*j+1)) for j in range((g-1)//2+1)]
 E=[F(comb(g,2*j)) for j in range(g//2+1)]
 OO,EE,GG,CD=cv(O,O),cv(E,E),cv(G,G),cv(C,D)
 p=[b[h-j]*comb(g,k-2*(h-j)) for j in range(h+1)]
 rows={name:[] for name in ['V','G1','G2','G3']}
 for j in range(-1,2*h+1):
  sg=-1 if j%2 else 1
  rows['V'].append(sg*coeff(OO,j)*coeff(GG,j+2*x))
  rows['G1'].append(sg*coeff(EE,j+1)*coeff(GG,j+2*x))
  rows['G2'].append(sg*2*coeff(OO,j)*coeff(CD,j-1+2*x))
  rows['G3'].append(sg*2*coeff(EE,j+1)*coeff(CD,j-1+2*x))
  need(coeff(OO,j)+coeff(EE,j+1)==comb(2*g,2+2*j),'complete alpha parity')
 for name in list(rows):rows[name]=trim(rows[name])
 rows['W']=add(rows['V'],rows['G1'])
 rows['K']=mul(rows['G2'],F(1,2))
 rows['R0']=add(rows['G1'],rows['G3'])
 rows['skip']=add(rows['G2'],rows['G3'])
 rows['Q-V']=add(rows['G1'],rows['G2'],rows['G3'])
 rows['Q']=add(rows['V'],rows['Q-V'])
 return p,rows

def even(p,s):
 n=len(p)-1;out=[F(0)]*(2*n+1)
 for i,c in enumerate(p):out[2*i]=c*s**(n-i)
 return out

def carrier(p,s,g):return cv(even(p,s),[F(comb(g,j)) for j in range(g+1)])
def primitive(p):
 den=lcm(*(v.denominator for v in p));nums=[int(v*den) for v in p]
 div=reduce(gcd,map(abs,nums));return [v//div for v in nums]

q,h,x,g,k=5,4,1,14,9
a=[F(64)**i for i in range(q)]
b,c,d=kernels(a);ratios=interlace_certificate(a,b,c,d)
p,rows=model(b,c,d,h)
lam=F(1,2**25)
scaled_a=[lam*v for v in a]
bs,cs,ds=kernels(scaled_a);interlace_certificate(scaled_a,bs,cs,ds)
ps,scaled=model(bs,cs,ds,h)
# A rational bracket for the smallest scaled phase, corresponding to base s=lam*z.
I=(F(4938,1000),F(4939,1000));J=(lam*I[0],lam*I[1])
need(at(ps,I[0])*at(ps,I[1])<0,'original joint zero bracket')
need(at(p,J[0])*at(p,J[1])<0,'same base joint zero bracket')
# A positive derivative on [0,right] proves that this root is the smallest positive root.
deriv=[i*p[i] for i in range(1,len(p))]
need(interval(deriv,(F(0),J[1]))[0]>0 and at(p,0)<0,'unique smallest root')
base_signs={name:sign_on(rows[name],J,sgn,'base '+name) for name,sgn in
 [('V',-1),('W',-1),('G1',-1),('G3',1),('R0',-1),('K',1),('skip',1),('Q-V',-1),('Q',-1)]}
scaled_signs={name:sign_on(scaled[name],I,sgn,'scaled '+name) for name,sgn in
 [('V',-1),('W',-1),('R0',-1),('K',1),('skip',1),('Q-V',1),('Q',1)]}
# Exact coefficient homogeneity; proof applies to every positive real scale.
need(all(ps[j]/lam**j==lam**x*p[j] for j in range(len(p))),'first row scaling')
for name in ['V','G1','W']:
 need(all(coeff(scaled[name],j)/lam**j==lam**(2*x-1)*coeff(rows[name],j) for j in range(max(len(scaled[name]),len(rows[name])))),name+' hit scaling')
for name in ['G2','G3','skip','K']:
 need(all(coeff(scaled[name],j)/lam**j==lam**(2*x-2)*coeff(rows[name],j) for j in range(max(len(scaled[name]),len(rows[name])))),name+' skip scaling')
# The source retains both coupled u-Euler identities at the coefficient level.
for B,C,D in [(b,c,d),(bs,cs,ds)]:
 for i in range(q):
  need(2*(i+1)*B[i+1]==F(2,3)*(q+2+2*i)*C[i],'first Euler identity')
  need((q+1+2*i)*D[i]==(2+3*i)*C[i],'second Euler identity')
# Independent literal ordinary-carrier construction; the retuned carriers are identical.
for s in [J[0],J[1],F(1,2**23)]:
 H=[carrier(P,s,g) for P in (b,c,d)]
 Hs=[carrier(P,s/lam,g) for P in (bs,cs,ds)]
 need(H==Hs,'same ordinary carriers under scale and retuning')
 need(coeff(H[0],k)==s**x*at(p,s),'literal original first coefficient')
 need(coeff(cv(H[0],H[0]),2*k)==s**(2*x-1)*at(rows['W'],s),'literal hit coefficient')
 need(-2*coeff(cv(H[1],H[2]),2*k-2)==s**(2*x-2)*at(rows['skip'],s),'literal skip coefficient and carry')
 # Same-zero lowering identity, kept without an independent-factor root substitution.
 lower=cv([F(comb(g-1,j)) for j in range(g)],add(even(b,s),mul([0]+even(c,s),F(2,3))))
 need(k*coeff(H[0],k)==g*coeff(lower,k-1),'zero-preserving lowering')
# Algebraic cancellation uses lambda0=-skip/W in Q(s*) exactly.
# These polynomial remainders give explicit representatives in the quartic field.
residues={name:[str(v) for v in rem(rows[name],p)] for name in ['W','skip']}
need(interval(rows['skip'],J)[0]>0 and interval(rows['W'],J)[1]<0,'algebraic scale lambda0 positive')
need(interval(add(rows['skip'],mul(rows['W'],lam)),J)[0]>0,'rational scale below lambda0')
# The genuine F_15 linear anchor fixes the remaining relative scale.
# Direct reconstruction checks this separately from the scale proof.
lam_norm=F(13,sum(a))
bn,cn,dn=kernels([lam_norm*v for v in a])
pn,normalized=model(bn,cn,dn,h)
need([-bn[-2],-cn[-2],-dn[-2]]==[13,12,11],'genuine linear anchors')
need(all(pn[j]/lam_norm**j==lam_norm**x*p[j] for j in range(len(p))),'normalized first row scaling')
normalized_I=(J[0]/lam_norm,J[1]/lam_norm)
need(at(pn,normalized_I[0])*at(pn,normalized_I[1])<0,'normalized joint zero')
normalized_signs={name:sign_on(normalized[name],normalized_I,sgn,'normalized '+name) for name,sgn in
 [('V',-1),('W',-1),('K',1),('skip',1),('Q-V',-1),('Q',-1)]}
need(interval(add(rows['skip'],mul(rows['W'],lam_norm)),J)[1]<0,'linear anchor above lambda0')
# Actual inherited q2/h1 control, unrelated to the formal hostile's beta law.
ap,ar=model([F(1),F(-4),F(1)],[F(-3),F(1)],[F(-2),F(1)],1)
need(at(ap,2)==0 and at(ar['Q'],2)==-12610,'actual (-9,1,6) negative doubled control')
# The earlier Euler-only hostile fails the stronger all-weight pencil requirement.
old_a=list(map(F,[1,2,3,5,9,12]));ob,oc,od=kernels(old_a)
need(at(oc,1)*at(od,1)<0,'old Euler hostile excluded by weighted compatibility')

bank={'q':q,'h':h,'x':x,'g':g,'B_roots':list(map(str,a)),'p_primitive':primitive(p),'scaled_roots':list(map(str,scaled_a)),
 'phase_bracket':list(map(str,I)),'lambda':str(lam),'lambda_norm':str(lam_norm),'normalized_phase_bracket':list(map(str,normalized_I)),
 'interlacer_ratios':ratios,'quartic_residues':residues}
# Compact outward-rounded integer intervals are independently inspectable.
# The underlying acceptance tests above use the unrounded rational endpoints.
outer_intervals={}
for regime,signs in [('base',base_signs),('dyadic',scaled_signs),('linear_anchor',normalized_signs)]:
 outer_intervals[regime]={}
 for name,(lo,hi) in signs.items():
  lo,hi=F(lo),F(hi)
  lower=lo.numerator//lo.denominator
  upper=-((-hi.numerator)//hi.denominator)
  need(lower>0 or upper<0,'outward interval retains strict '+regime+' '+name)
  outer_intervals[regime][name]=[lower,upper]
bank['outward_integer_intervals_for_s_times_response']=outer_intervals
print('STATUS: EXACT combined-model full-response sign hostile and algebraic double cancellation; actual binomial transport OPEN')
print('parameters',q,h,x,g,'base_B_roots',list(map(str,a)))
print('first_positive_phase_polynomial_ascending',primitive(p))
print('scale',lam,'scaled_B_roots',list(map(str,scaled_a)),'scaled_phase_bracket',list(map(str,I)))
print('strict_interlacer_residues_C_D',ratios)
print('BASE SIGNS: V,W,G1,R0,Q-V,Q negative; G3,K,skip positive')
print('SCALED SIGNS: V,W,R0 negative; K,skip,Q-V,Q positive')
print('LINEAR-ANCHORED SCALE',lam_norm,'beta_linear_coefficients',[13,12,11])
print('LINEAR-ANCHORED SIGNS: V,W,Q-V,Q negative; K,skip positive')
print('outward_integer_intervals_for_s_times_response',json.dumps(outer_intervals,sort_keys=True))
print('FIELD CERTIFICATE: s* is unique smallest positive root; lambda0=-skip(s*)/W(s*)>0')
print('AT PHASE -s*/lambda0: P_lambda0=Q_lambda0=0 exactly; all retained predicates survive')
print('quartic_remainders_of_sW_and_sSkip',json.dumps(residues,sort_keys=True))
print('gates',GATES)
print('semantic_sha256',hashlib.sha256(json.dumps(bank,sort_keys=True,separators=(',',':')).encode()).hexdigest())
