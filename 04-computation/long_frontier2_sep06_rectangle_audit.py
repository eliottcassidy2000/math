#!/usr/bin/env python3
"""Independent ordinary-carrier and Bernstein-tensor audit of the rectangle."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib
import itertools
import json
import subprocess
import sys
import tempfile
import sympy as S

BASE=Path(__file__).resolve().parents[1]
STEM='long_frontier2_sep06_rectangle'
PINS={'source':'4fb5fe880fbc030b3461889d021f7fa24b6376b23cdc1251afd99b0327591ea5',
      'output':'572d5f8192b0a567efd95b6db25d23f8dcd7277477dacfd27d75a98134668bad',
      'certificate':'cd1af1f580b51fc1b556a7cdf44f28203e464a3ba41e3d03e55363278816d412'}
GATES=0


def check(value,label):
    global GATES
    GATES+=1
    if not value: raise RuntimeError(label)


def choose(n,k): return comb(n,k) if 0<=k<=n else 0


def finite_basis(power,degree,lo,hi,index):
    return sum(F(choose(power,j)*choose(degree-power,index-j))*lo**(power-j)*hi**j
               for j in range(power+1))


def chart_coeff(terms,boxes,powers):
    result=F(0)
    for exps,coefficient in terms.items():
        value=coefficient
        for p,(degree,lo,hi),k in zip(exps,boxes,powers):
            value*= (F(choose(p,k))*lo**(p-k) if hi is None and k<=p
                     else F(0) if hi is None
                     else finite_basis(p,degree,lo,hi,k))
        result+=value
    return result


def listed(rows): return {tuple(row['powers']):F(row['coefficient']) for row in rows}


def main():
    src=BASE/'04-computation'/f'{STEM}.py'
    out=BASE/'05-knowledge/results'/f'{STEM}.out'
    certpath=BASE/'05-knowledge/results'/f'{STEM}_certificate.json'
    for name,path in [('source',src),('output',out),('certificate',certpath)]:
        check(hashlib.sha256(path.read_bytes()).hexdigest()==PINS[name], 'frozen '+name+' pin')
    cert=json.loads(certpath.read_text())
    a,b,f,t,s=S.symbols('a b f t s')
    beta=(1+13*t+55*t*t+a*t**3+b*t**4+f*t**5)/t
    c=(1+12*t+45*t*t+2*a*t**3/3+3*b*t**4/7)/t
    d=(1+11*t+36*t*t+5*a*t**3/12+b*t**4/7)/t
    O=sum(S.Integer(comb(14,2*j+1))*t**j for j in range(7))
    E=sum(S.Integer(comb(14,2*j))*t**j for j in range(8))
    # Ordinary expanded products, not the producer's integer28 shortcut or array convolution.
    left=S.Poly(S.expand(t*O**2+E**2),t)
    right=S.Poly(S.expand(t*t*(beta**2+2*t*c*d)),t)
    coefficients={j:S.expand(left.nth(j+1)*right.nth(j+2)) for j in range(-1,9)}
    check(coefficients[-1]==28,'independent negative carry')
    # Explicit parity also keeps the negative carry's sign in exact integers.
    qbar=S.expand(sum(coefficients[j]*(-1 if j%2 else 1)*s**(j+1) for j in range(-1,9)))
    bp=S.Poly(S.expand(t*beta),t)
    pp=S.Poly(O,t)
    Poriginal=S.expand(sum(pp.nth(j)*bp.nth(j+1)*(-s)**j for j in range(5)))
    p=f*s**4-S.Rational(12,7)*b*s**3+a*s*s-10*s+S.Rational(1,11)
    check(S.expand(Poriginal-2002*p)==0,'ordinary original first row')
    fzero=12*b/(7*s)-a/s**2+10/s**3-S.Rational(1,11)/s**4
    h=S.Poly(S.cancel(-S.Rational(11,14)*qbar.subs(f,fzero)),a,b,s)
    terms={powers:F(coeff) for powers,coeff in h.terms()}
    check(terms==listed(cert['H']),'all eliminated coefficients independently reconstructed')
    check(h.degree_list()==(2,2,8) and len(terms)==19,'exact eliminated support')
    aa=(F(167,2),F(169,2));bb=(F(69,2),F(71,2))
    charts=[(F(99,10000),F(1,100)),(F(1,9),F(13,100)),(F(1),F(8,5)),(F(10),None)]
    for index,(lo,hi) in enumerate(charts):
        box=[(2,*aa),(2,*bb),(8,lo,hi)]
        saved=listed(cert['positive_charts'][index]['coefficients'])
        check(len(saved)==81,'complete saved chart')
        for powers in itertools.product(range(3),range(3),range(9)):
            actual=chart_coeff(terms,box,powers)
            check(actual==saved[powers] and actual>0,'independent tensor coefficient '+str((index,powers)))
    def pv(A,B,C,x): return C*x**4-F(12,7)*B*x**3+A*x*x-10*x+F(1,11)
    signs=[(1,-1),(-1,1),(1,-1),(-1,)]
    for rec in cert['phase_endpoint_corners']:
        index=rec['chart'];A,B,C=F(rec['a']),F(rec['b']),F(rec['f'])
        ends=[x for x in charts[index] if x is not None]
        vals=[pv(A,B,C,x) for x in ends]
        check(vals==list(map(F,rec['values'])),'independent phase-corner values')
        for value,sign in zip(vals,signs[index]):check(value*sign>0,'uniform phase endpoint sign')
    upper=[F(0),F(71,2),F(-167,2),F(55),F(-13),F(1)]
    derivative=sum(j*upper[j]*F(2,5)**(j-1) for j in range(1,6))
    check(derivative==F(-81,10),'critical derivative upper bound')
    barrier={(0,):F(5)}
    for j in range(1,6):barrier[j,]=-upper[j]
    saved=listed(cert['upper_product_barrier'])
    for k in range(6):
        actual=chart_coeff(barrier,[(5,F(0),F(2,5))],(k,))
        check(actual==saved[k,] and actual>0,'independent product barrier coefficient')
    nodes=[F(0),F(1,10),F(1),F(3),F(5),F(7)]
    for rec in cert['genuine_B_prism_corners']:
        A,B,C=F(rec['a']),F(rec['b']),F(rec['f'])
        vals=[x**5-13*x**4+55*x**3-A*x*x+B*x-C for x in nodes]
        check(vals==list(map(F,rec['values'])),'independent B-prism values')
        for i,value in enumerate(vals):check(value*(-1 if i%2==0 else 1)>0,'genuine B-prism sign')
    hlo,hhi=map(F,cert['f6_hostile_bracket'])
    check(pv(F(84),F(35),F(6),hlo)*pv(F(84),F(35),F(6),hhi)<0,'f6 original root')
    hs={(k,):F(h.as_expr().subs({a:84,b:35}).expand().coeff(s,k)) for k in range(9)}
    saved=listed(cert['f6_negative_H_transform'])
    for k in range(9):
        actual=chart_coeff(hs,[(8,hlo,hhi)],(k,))
        check(actual==saved[k,] and actual<0,'independent hostile chart coefficient')
    check(pv(F(84),F(35),F(1),F(4))!=0,'centre control is off root')
    check(qbar.subs({a:84,b:35,f:1,s:4})==350398552675052,'independent off-root full response')
    with tempfile.TemporaryDirectory(prefix='rectangle-independent-') as tmp:
        generated=[];outputs=[]
        for optimization in [[],['-O']]:
            target=Path(tmp)/('optimized.json' if optimization else 'normal.json')
            outputs.append(subprocess.check_output([sys.executable,'-B']+optimization+[str(src),'--certificate',str(target)]))
            generated.append(target.read_bytes())
        check(outputs[0]==outputs[1]==out.read_bytes(),'normal optimized frozen producer output')
        check(generated[0]==generated[1]==certpath.read_bytes(),'normal optimized frozen certificate')
    print('INDEPENDENT RECTANGLE ANALYTIC CARRIERS AND COEFFICIENTS PASS')
    print('UNIVERSE complete ordinary P,Q; all19 H monomials; all324 tensor coefficients; all corners, barriers and hostiles')
    print('PRODUCER_REPLAY normal optimized frozen output and certificate match;457 producer gates')
    for name,pin in PINS.items():print(name.upper()+'_SHA256',pin)
    print('INDEPENDENT_EXPLICIT_GATES',GATES)


if __name__=='__main__':main()
