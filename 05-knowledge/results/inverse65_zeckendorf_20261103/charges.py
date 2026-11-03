"""Complete signed linear-charge classification and a quadratic 7p family."""
from pathlib import Path
from math import gcd,isqrt,lcm
from functools import reduce
import json
import sympy as sp
from bank import *
from lane_inverse import invert
HERE=Path(__file__).resolve().parent

def primitive(v):
    den=lcm(*(x.q for x in v));a=[int(x*den) for x in v];g=reduce(gcd,a)
    a=[x//g for x in a]
    if next(x for x in a if x)!=abs(next(x for x in a if x)):a=[-x for x in a]
    return tuple(a)

def dot(a,b):return sum(x*y for x,y in zip(a,b))

def modular_clock(p0,p1,r,m):
    states=[];seen={};v=(p0%m,p1%m)
    while v not in seen:
        seen[v]=len(states);states.append(v);v=(v[1],(40*v[1]-13*v[0]-10*r)%m)
    start=seen[v]
    return {'modulus':m,'tail':start,'period':len(states)-start,
            'zero_indices_in_stored_tail_and_cycle':[i for i,v in enumerate(states) if v[0]==0],
            'states':states}

def main():
    x=sp.symbols('X');charges=[]
    need(MATS['r21']==MATS['r38'],'same observer, distinct literal charts')
    for label,raw in MATS.items():
        m=sp.Matrix(raw)
        for sign in (1,-1):
            ns=(m-sign*sp.eye(3)).nullspace()
            if not ns:continue
            need(len(ns)==1,('charge dimension',label))
            ell=primitive(ns[0]);poly=m.charpoly(x).as_expr()
            quad=sp.Poly(sp.cancel(poly/(x-sign)),x)
            t=-int(quad.nth(1));d=int(quad.nth(0));f=sp.Matrix([1,0,0])
            forcing=(m*m-t*m+d*sp.eye(3))*f
            index=next(i for i,v in enumerate(ell) if v)
            kappa=sp.Rational(forcing[index],ell[index])
            need(forcing==kappa*sp.Matrix(ell),'phase forcing proportional to charge')
            charges.append({'letter':label,'sign':sign,'primitive_charge':ell,
                            'quadratic_trace':t,'quadratic_norm':d,'forcing_multiplier':str(kappa),
                            'characteristic_polynomial':str(sp.factor(poly))})
    need(len(charges)==7,'complete signed-charge list')
    families={};mods=(2,3,5,7,11,13,17,19,23,27)
    for seed in ('r29','r43'):
        initial=vector((seed,));J=2*(initial[1]+initial[2])-3*initial[0];need(abs(J)==7,'7-divisible odd charge')
        r=J//7;v=initial;terms=[]
        for k in range(6):
            word=(seed,)+('r21',)*k
            need(vector(word)==v,'literal family word')
            need(2*(v[1]+v[2])-3*v[0]==J,'conserved charge')
            p=v[0]+2*r;child=step(v,MATS['b'])
            need(child[0]==7*p,'exact sevenfold attachment')
            terms.append({'k':k,'H':v[0],'p':p,'seven_p':child[0],'state':v,
                          'processing_word':word+('b',),'vertex_clock':5+5*k+3})
            v=step(v,MATS['r21'])
        for term in terms:
            need(invert(term['p'],seed)==tuple(term['processing_word']),'quadratic lane inverse')
        need(invert(17,seed) is None,'lane nonmember control')
        for k in range(4):
            need(terms[k+2]['H']==40*terms[k+1]['H']-13*terms[k]['H']+6*J,'H recurrence')
            need(terms[k+2]['p']==40*terms[k+1]['p']-13*terms[k]['p']-10*r,'quotient recurrence')
        U=[13*t['H']+3*J for t in terms]
        norms=[U[k+1]**2-40*U[k]*U[k+1]+13*U[k]**2 for k in range(5)]
        need(all(norms[k+1]==13*norms[k] for k in range(4)),'quadratic norm clock')
        clocks=[modular_clock(terms[0]['p'],terms[1]['p'],r,m) for m in mods]
        for c in clocks:
            if c['modulus']!=13:need(c['tail']==0,('invertible clock',seed,c['modulus']))
            else:need(c['tail']==1 and c['period']==13,'norm13 lost coordinate')
        powers=[]
        for e in (1,2,3):
            c=modular_clock(terms[0]['p'],terms[1]['p'],r,13**e)
            need(c['tail']==e and c['period']==13**e,('all-power clock control',seed,e))
            cyc=c['states'][c['tail']:]
            need(len({v[0] for v in cyc})==13**e,'clock visits each residue exactly once')
            powers.append({k:v for k,v in c.items() if k!='states'})
        crt=[]
        for e,n in ((1,11),(2,3),(2,19)):
            c=modular_clock(terms[0]['p'],terms[1]['p'],r,(13**e)*n)
            nc=next(x for x in clocks if x['modulus']==n)
            need(c['tail']==e and c['period']==lcm(13**e,nc['period']),'CRT clock')
            crt.append({k:v for k,v in c.items() if k not in ('states','zero_indices_in_stored_tail_and_cycle')})
        families[seed]={'charge_J':J,'r':r,'terms':terms,'norms':norms,'modular_clocks':clocks,
                        'power13_checks':powers,'CRT_checks':crt}
    prime_certificates={}
    for p in (13,499,1009,19861,21001):
        need(all(p%d for d in range(2,isqrt(p)+1)),('trial-division prime',p))
        prime_certificates[str(p)]={'trial_division_through':isqrt(p),'prime':True}
    need(families['r29']['terms'][1]['p']==23**2,'square hostile at first clock')
    need(all(t['p']%3==0 for t in families['r43']['terms'] if t['k']%2),'odd-clock composite control')
    result={'status':'PASS','signed_charges':charges,'families':families,'prime_certificates':prime_certificates,
            'scope':'All-order exact constructions and recurrences; no assertion of infinitely many prime terms or complete prime coverage.'}
    (HERE/'charges.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','charges':charges,'family_first_terms':{s:[t['p'] for t in f['terms']] for s,f in families.items()},
                      'clock_lengths':{s:[(c['modulus'],c['tail'],c['period']) for c in f['modular_clocks']] for s,f in families.items()}}))
if __name__=='__main__':main()
