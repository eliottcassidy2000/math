#!/usr/bin/env python3
"""Read-only rational certificate; no word census and no floating point."""
from fractions import Fraction as F
import json

checks=0

def check(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def log_interval(z,m):
    # ln((1+z)/(1-z)), with positive geometric-tail majorant.
    check(0<z<1 and m>=1,'valid atanh-series parameters')
    lower=2*sum((z**(2*k+1)/F(2*k+1) for k in range(m)),F(0))
    tail=2*z**(2*m+1)/(F(2*m+1)*(1-z*z))
    return lower,lower+tail

P,Q,C0=3**19,2**30,4912087879
D=P-Q
r1=815464443
cl4=C0*(P**4-Q**4)//D
r4=(-cl4*pow(P**4,-1,Q**4))%(Q**4)
check(r4==588943162011782273059357097228697595,'least positive K4 source')
check((P**4*r4+cl4)%(Q**4)==0 and 0<r4<Q**4,'K4 cylinder')
# V=K4F, center maximum in the fixed expanding bank.
sv=3*P**4;mv=4*Q**4;hv=3*cl4+Q**4
C=F(hv,sv-mv)
check(C==F(83829391903560839836419813590128166431,157489106280555890461927692448373459),'bank center')
check(C>F(C0,D),'maximum expansion center')
E=C+1
epsilon=E/r4
eta=epsilon/(1-epsilon)
check(0<epsilon<1,'valid non-descent gap threshold')
# All small expansion words K,K2,K3 descend after one reset.
cl3=C0*(P**3-Q**3)//D
s3=F(3*P**3,4*Q**3)
e3=F(3*cl3+Q**3,4*Q**3)
check(s3<1 and s3*r1+e3<r1,'K3 followed F descends on entire source lower bound')
check(F(81,512)<F(3,4) and F(101,512)<F(1,4),'H affine map below F for positive inputs')
# Two explicit neighboring fractions: a/b < ln2/ln3 < c/d.
a=253438635377865553122983825848491
b=401690733307859070263635245069941
c=326161216350563845551964374972875
d=516953297105243726377334761586514
check(c*b-a*d==1,'unimodular bracketing vectors')
l2,u2=log_interval(F(1,3),160)
l3,u3=log_interval(F(1,2),160)
alpha_lo,alpha_hi=l2/u3,u2/l3
check(l3>1,'ln3 exceeds one')
check(F(a,b)<alpha_lo<alpha_hi<F(c,d),'strict logarithmic bracketing')
check(b*alpha_lo-a>eta,'lower linear-form gap exceeds counterexample threshold')
check(b+d==918644030413102796640970006656455,'length bound')
# Bounds printed as rational intervals with a shared decimal denominator.
scale=10**36

def enclosure(lo,hi):
    lower=(lo.numerator*scale)//lo.denominator
    upper=-((-hi.numerator*scale)//hi.denominator)
    check(F(lower,scale)<=lo<=hi<=F(upper,scale),'reported rational enclosure')
    return [lower,upper]

report=dict(status='PROVED necessary counterexample length; zero-head conjecture remains OPEN',
 scope='Nonempty expansion-bank word followed by nonempty reset-bank word',
 K4_least_source=r4,lower_fraction=[a,b],upper_fraction=[c,d],
 minimum_counterexample_steps=b+d,log_series_terms=160,
 envelope_constant=str(E),epsilon=str(epsilon),eta=str(eta),
 interval_denominator=scale,
 epsilon_interval=enclosure(epsilon,epsilon),eta_interval=enclosure(eta,eta),
 lower_linear_form_interval=enclosure(b*alpha_lo-a,b*alpha_hi-a))
report['checks']=checks
print(json.dumps(report,indent=2,sort_keys=True))
