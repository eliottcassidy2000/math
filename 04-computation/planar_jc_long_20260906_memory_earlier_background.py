#!/usr/bin/env python3
"""Recover the actual sixth background row and xi10 on the boundary G_m.

Uses the hash-pinned audited THM4426 prefix reconstruction. Specialization
changes six primitive integer normalizations; these gates are replaced by
nonzero rational pivot gates, rather than discarded. Every raw row9..12
compatibility is checked after the resulting source graphs are substituted.
"""
from pathlib import Path
import hashlib
ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/'04-computation/jc2_source_normal_row14_weight18_two_memory_rational_conic_companion_thm4426.py'
EXPECTED='c647bdf2c894d60413b041197cff446884191bf1469372e6b4a6a55ae96fd52a'

def replace_once(text,old,new):
 if text.count(old)!=1:raise RuntimeError('unique inherited anchor '+old[:60])
 return text.replace(old,new)

def main():
 raw=PATH.read_bytes()
 if hashlib.sha256(raw).hexdigest()!=EXPECTED:raise RuntimeError('canonical dependency drift')
 text=raw.decode()
 text=replace_once(text,'x, t = R8.x, R8.t','''x, t = R8.x, R8.t
_bd={R8.Phi:sp.Integer(0),R8.eta:sp.Integer(0),R8.alpha11:sp.Integer(1)}
for _name in ('H','A3','C3'):
    setattr(R8,_name,sp.expand(getattr(R8,_name).subs(_bd)))
R8.Phi=sp.Integer(0);R8.eta=sp.Integer(0);R8.alpha11=sp.Integer(1)
''')
 for expression,value in [('sp.diff(equation9, c14)','216513000'),
                         ('sp.diff(xi_equation, R8.xi10)','57672000'),
                         ('sp.diff(beta_equation, R8.beta11)','15'),
                         ('sp.diff(c42_equation, c42)','89131914000'),
                         ('sp.diff(equation_c23, c23)','675'),
                         ('sp.diff(equation_c70, c70)','368130186720000')]:
  text=replace_once(text,expression+' == '+value,
                    '('+expression+').is_Rational and '+expression+' != 0')
 anchor='    base5 = {R8.Phi, R8.eta, R8.alpha11, c51, z, h18}'
 injection='''    final_graphs={**bracket_subs,**gate,**resolved12,c51:sp.Rational(1087,135)}
    all_raw=raw9+raw_depth9+raw10+raw_depth10+raw11+raw_depth11+raw12+raw_depth12
    for index,value in enumerate(all_raw):
        check(resolve(value,final_graphs)==0,'specialized prefix raw coefficient '+str(index))
    Xi=sp.Symbol('Xi')
    a5=resolve(a10[5],{**bracket_subs,**gate}).subs(R8.xi10,Xi)
    c5=resolve(c10[5],{**bracket_subs,**gate}).subs(R8.xi10,Xi)
    # a10 has row10 graphs already substituted, so recover the raw universal
    # row5 directly from the original rows4..8 reconstruction instead.
    a5=resolve(arows[5],gate).subs(R8.xi10,Xi)
    c5=resolve(crows[5],gate).subs(R8.xi10,Xi)
    expected_a5=128*x**6/sp.Integer(9)-8576*x**4/sp.Integer(75)-6*x*x*Xi/sp.Integer(11)-203776*x*x/sp.Integer(7425)-x/2+sp.Rational(3528704,6075)
    expected_c5=-32*x**7/sp.Integer(3)-128*x**5/sp.Integer(5)+9*x**3*Xi/sp.Integer(22)+130816*x**3/sp.Integer(275)+3*x*x/sp.Integer(8)+9*x*Xi/sp.Integer(11)-4024832*x/sp.Integer(4455)+sp.Rational(3,4)
    check(exact(a5-expected_a5)==0 and exact(c5-expected_c5)==0,'actual sixth row')
    actual_xi=exact(resolved12[R8.xi10].subs(c51,sp.Rational(1087,135)))
    expected_xi=(108391820625*h18+3765431711250*z+sp.Integer(324974300165767168))/sp.Integer(24542012448000)
    check(actual_xi==expected_xi,'actual boundary xi graph')
    u=801*h18+27826*z
    shift=sp.Rational(85855050266495746048,37533020625)
    slope=sp.Rational(1485,269322496);intercept=sp.Rational(34969998848,55604475)
    check(exact(actual_xi-slope*(u+shift)-intercept)==0,'xi affine in actual Gm parameter')
    print('BACKGROUND_DEPENDENCY_SHA256','c647bdf2c894d60413b041197cff446884191bf1469372e6b4a6a55ae96fd52a')
    print('A5',sp.sstr(a5));print('C5',sp.sstr(c5))
    print('XI_Z_H',sp.sstr(actual_xi))
    print('XI_GM',str(slope)+'*s+'+str(intercept))
    print('RAW_PREFIX_COEFFICIENTS',len(all_raw),'ALL_ZERO')
    print('CHECKS',CHECKS,'PASS')
    return
'''+anchor
 text=replace_once(text,anchor,injection)
 exec(compile(text,str(PATH),'exec'),{'__name__':'__main__','__file__':str(PATH)})
if __name__=='__main__':main()
