#!/usr/bin/env python3
"""Exact controls for the fixed-chart alternating cubic contact-chain theorem.

Run normally and under -O; no check relies on Python's removable assertions.
"""
from fractions import Fraction as Q
from math import gcd, isqrt, lcm
from collections import Counter
from hashlib import sha256
import json
import sys
import sympy as s
sys.set_int_max_str_digits(0)
checks=Counter()
def check(condition,name):
    checks[name]+=1
    if not condition: raise RuntimeError(name)
def same(a,b,name): check(s.cancel(a-b)==0,name)
def spread(t): return 4*t*(t*t-1)/(t*t+1)**2
def F(t): return -t*(t*t-4*t+1)/(2*t**3+t*t-2*t+1)
def gp(t): return 1/F(t)
def gm(t): return -(t*t-4*t+1)/(t**3-2*t*t+t+2)
def ppt(t):
    n,m=t.numerator,t.denominator
    delta=2 if n%2 and m%2 else 1
    a,b,c=(n*n-m*m)//delta,2*n*m//delta,(n*n+m*m)//delta
    check(a>0 and b>0 and a*a+b*b==c*c and gcd(a,b)==1,'primitive_ppt')
    return min(a,b),max(a,b),c
def rational_sqrt(x):
    check(x>0,'positive_contact_square')
    n,d=isqrt(x.numerator),isqrt(x.denominator)
    check(n*n==x.numerator and d*d==x.denominator,'rational_square_exact')
    return Q(n,d)
def positive_on_shift(poly,radius,name):
    p=s.Poly(s.expand(poly.subs(t,1+h)),h)
    check(all(x.is_Integer for x in p.all_coeffs()),'integer_sign_coefficients')
    lower=Q(int(p.nth(0)))-sum((Q(-int(p.nth(i)))*radius**i for i in range(1,p.degree()+1) if p.nth(i)<0),Q(0))
    check(lower>0,name)
    return str(lower)
t,h,n,m=s.symbols('t h n m', real=True)
A=t*t-4*t+1;B=t**3-2*t*t+t+2
Pp=t**6+24*t**5-13*t**4-16*t**3+19*t*t-8*t+1
Pm=t**6-8*t**5+19*t**4-16*t**3-13*t*t+24*t+1
same(1+spread(t)+spread(gp(t)),(Pp/((t*t+1)**2*(5*t*t-4*t+1)))**2,'plus_same_port_identity')
same(1-spread(t)-spread(gm(t)),(Pm/((t*t+1)**2*(t*t-4*t+5)))**2,'minus_same_port_identity')
same(gm(t),F(1/t),'reciprocal_relation')
for G in [gp,gm]:
    num,den=s.fraction(s.cancel(G(t)))
    check(max(s.degree(num,t),s.degree(den,t))==3,'degree_exactly_three')
    check(s.degree(s.gcd(num,den),t)==0,'coprime_cubic_forms')
same(gp(t)-1,-(t-1)*(3*t*t-1)/(t*A),'plus_above_boundary')
same(gm(t)-1,(t-1)*(3-t*t)/B,'minus_above_boundary')
same(gp(t)-t,-(t-1)**2*(t*t+1)/(t*A),'plus_growth_identity')
same(gm(t)-t,-(t-1)**2*(t*t+1)/B,'minus_decrease_identity')
same(s.diff(gp(t),t),(t*t+1)*(9*t*t-8*t+1)/(t*t*A*A),'plus_derivative_identity')
same(s.diff(gm(t),t),(t*t+1)*(t*t-8*t+9)/(B*B),'minus_derivative_identity')
same((2*t-1)-gp(t),(t-1)*(2*t**3-5*t*t+2*t-1)/(t*A),'plus_growth_upper_identity')
chart={}
for name,p in [('negative_A',-A),('minus_positive',3-t*t),('plus_derivative',9*t*t-8*t+1),('minus_derivative',t*t-8*t+9),('plus_growth_under_two',-(2*t**3-5*t*t+2*t-1)),('spread_derivative',-(t**4-6*t*t+1))]:
    chart[name]=positive_on_shift(p,Q(1,16),'uniform_chart_sign')
Hp=s.cancel(gm(gp(t)));Hm=s.cancel(gp(gm(t)))
same(Hp,F(F(t)),'two_step_F_squared')
D9=16*t**9-3*t**8-132*t**7+248*t**6-124*t**5-6*t**4+20*t**3-4*t+1
same(Hp-t,-16*t**3*(t-1)**3*(t*t+1)*(t*t+2*t-1)/D9,'two_step_decrease_identity')
check(s.Poly(D9.subs(t,1+h).expand(),h).all_coeffs()==[16,141,420,584,440,280,272,224,96,16],'positive_two_step_denominator')
for H in [Hp,Hm]:
    check(s.series(H.subs(t,1+h),h,0,4).removeO()==1+h-4*h**3,'parabolic_cubic_coefficient')
check(spread(Q(17,16))==Q(35904,297025)<Q(1,2),'entire_chart_admissible')
N=-n*(n*n-4*n*m+m*m)
D=2*n**3+n*n*m-2*n*m*m+m**3
same(N*N+D*D,(n*n+m*m)**2*(5*n*n-4*n*m+m*m),'homogeneous_norm_identity')
same(N+s.I*D,(n-s.I*m)**2*(-n+s.I*(2*n-m)),'gaussian_factorization')
check(s.resultant(-t*(t*t-4*t+1),2*t**3+t*t-2*t+1,t)==-64,'resultant_only_two')
for p in [(5*n*n-4*n*m+m*m)-(n*n+m*m)-4*n*(n-m),(n*n-4*n*m+5*m*m)-(n*n+m*m)-4*m*(m-n)]:
    same(p,0,'new_factor_old_prime_residue')
Rp=(5*t*t-4*t+1)/(t*t+1);Rm=(t*t-4*t+5)/(t*t+1)
for ratio in [Rp,Rm]:
    for q in [ratio-s.Rational(7,8),s.Rational(9,8)-ratio]:
        num,den=s.fraction(s.cancel(q))
        chart['ratio_sign_'+str(len(chart))]=positive_on_shift(num,Q(1,16),'uniform_cubic_height_ratio')
check(Rp.subs(t,s.Rational(17,16))==s.Rational(613,545),'plus_ratio_endpoint')
check(Rm.subs(t,s.Rational(17,16))==s.Rational(481,545),'minus_ratio_endpoint')
normalization_probes=0
for mm in range(1,101):
    for nn in range(mm+1,101):
        if gcd(nn,mm)!=1: continue
        delta=2 if nn%2 and mm%2 else 1
        aa=-nn*(nn*nn-4*nn*mm+mm*mm)
        bb=2*nn**3+nn*nn*mm-2*nn*mm*mm+mm**3
        check(gcd(aa,bb)==delta,'exact_gcd_parity_probe')
        check((aa//delta)%2 == (nn%2) if delta==2 else (aa%2)!=(bb%2),'output_parity_probe')
        normalization_probes+=1
cycles=[]
for eps in [1,-1]:
    for length in [1,6]:
        ts=[Q(65,64)]
        triples=[ppt(ts[0])]
        contact_roots=[]
        for i in range(length):
            sigma=eps*((-1)**i)
            old=ts[-1];new=(gp if sigma==1 else gm)(old)
            check(1<new<=Q(17,16),'orbit_in_fixed_chart')
            if i%2==1: check(1<new<ts[-2]<=Q(65,64),'two_step_orbit_decrease')
            check(0<spread(new)<Q(1,2),'orbit_admissible')
            contact_roots.append(rational_sqrt(1+sigma*spread(old)+sigma*spread(new)))
            nn,mm=old.numerator,old.denominator
            delta=2 if nn%2 and mm%2 else 1
            R=((5*nn*nn-4*nn*mm+mm*mm) if sigma==1 else (nn*nn-4*nn*mm+5*mm*mm))//delta
            T=ppt(new);C=triples[-1][2];Cp=T[2]
            check(Cp==C*C*R and gcd(C,R)==1,'primitive_exact_transport')
            check(Q(7,8)*C**3<Cp<Q(9,8)*C**3,'primitive_cubic_bounds')
            check(Cp>C*C,'primitive_strict_growth')
            ts.append(new);triples.append(T)
        check(len(set(triples))==len(triples),'petal_distinctness')
        H=triples[-1][2]
        check(H==lcm(*(T[2] for T in triples)),'last_hypotenuse_is_common_scale')
        M=2*H*H;labels=[M]
        for i,(a,b,c) in enumerate(triples):
            g=H//c;ports={-1:2*g*g*(c*c-4*a*b),1:2*g*g*(c*c+4*a*b)}
            outgoing=eps*((-1)**i)
            labels.extend([ports[-outgoing],ports[outgoing]])
        check(len(set(labels))==len(labels),'cycle_vertices_distinct')
        check(0<min(labels) and 2*H*H<max(labels)<4*H*H,'cycle_label_interval')
        roots=[]
        for i,x in enumerate(labels):
            y=x+labels[(i+1)%len(labels)];r=isqrt(y)
            check(r*r==y,'integer_cycle_edge_square');roots.append(r)
        encoded=json.dumps({'labels':labels,'roots':roots},separators=(',',':'))
        rec={'first_same_port':eps,'contact_count':length,'petal_count':len(triples),'integer_vertices':len(labels),
             'hypotenuse_bits':[T[2].bit_length() for T in triples],'largest_label_bits':max(labels).bit_length(),
             'certificate_sha256':sha256(encoded.encode()).hexdigest()}
        if length==1:rec.update(labels=labels,integer_square_roots=roots)
        cycles.append(rec)
print(json.dumps({'status':'FINITE-EXACT controls; universal theorem proved separately',
 'checks':dict(checks),'normalization_probes':normalization_probes,'uniform_chart_gaps':chart,'cycles':cycles},indent=2))
