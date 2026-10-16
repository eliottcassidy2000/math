"""Marked golden Legendre/Paley interfaces; exact finite controls, not coverage."""
from fractions import Fraction
from itertools import product
from pathlib import Path
from math import isqrt
import json


SIGNS=tuple(product((-1,1),repeat=3))


def prime(p):
    return type(p) is int and p>=3 and p%2==1 and all(p%d for d in range(3,isqrt(p)+1,2))


def chi(a,p):
    if type(a) is not int or not prime(p):raise ValueError('integer and odd prime required')
    r=pow(a%p,(p-1)//2,p)
    return 0 if r==0 else (1 if r==1 else -1)


def signature(p,anchors):
    if not prime(p):raise ValueError('odd prime required')
    if type(anchors) is not tuple or len(anchors)!=3 or any(type(a) is not int or not 0<=a<p for a in anchors) or len(set(anchors))!=3:
        raise ValueError('three distinct canonical integer anchors required')
    a,b,c=anchors
    trace=-sum(chi((x-a)*(x-b)*(x-c),p) for x in range(p))
    cells=[]
    for signs in SIGNS:
        boundary=sum(product3(1+signs[i]*chi(x-anchors[i],p) for i in range(3)) for x in anchors)
        numerator=p-sum(signs[i]*signs[j] for i in range(3) for j in range(i+1,3))-product3(signs)*trace-boundary
        if numerator%8:raise ArithmeticError('nonintegral port count')
        cells.append({'signs':list(signs),'count':numerator//8,'boundary':boundary})
    return {'p':p,'anchors':list(anchors),'trace':trace,'cells':cells,
            'relation':'tournament' if p%4==3 else 'graph'}


def product3(values):
    out=1
    for v in values:out*=v
    return out


def same_typed(a,b):
    if type(a) is not type(b):return False
    if type(a) is dict:return a.keys()==b.keys() and all(same_typed(a[k],b[k]) for k in a)
    if type(a) in (list,tuple):return len(a)==len(b) and all(same_typed(x,y) for x,y in zip(a,b))
    return a==b


def audit(record):
    if type(record) is not dict or type(record.get('anchors')) is not list:raise ValueError('source record required')
    expected=signature(record.get('p'),tuple(record['anchors']))
    if not same_typed(record,expected):raise ValueError('altered marked signature')
    return True


def direct(p,anchors):
    out={s:[] for s in SIGNS}
    for x in range(p):
        if x not in anchors:out[tuple(chi(x-a,p) for a in anchors)].append(x)
    return out


def f2_mul(a,b,p):
    # (a0+a1*s)(b0+b1*s), s^2=5, only used when 5 is nonsquare.
    return ((a[0]*b[0]+5*a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0])%p)


def f2_pow(a,n,p):
    v=(1,0)
    while n:
        if n&1:v=f2_mul(v,a,p)
        a=f2_mul(a,a,p);n//=2
    return v


def f2_sub(a,b,p):return ((a[0]-b[0])%p,(a[1]-b[1])%p)


def f2_chi(a,p):
    if a==(0,0):return 0
    return 1 if f2_pow(a,(p*p-1)//2,p)==(1,0) else -1


def main():
    checks=0
    def mulpoly(a,b):
        c=[0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):c[i+j]+=x*y
        return c
    # Exact coefficient comparison: all six Legendre roots of j=2048.
    lhs=mulpoly(mulpoly([1,-1,1],[1,-1,1]),[1,-1,1])
    for k,v in enumerate([0,0,8,-16,8]):lhs[k]-=v
    rhs=mulpoly(mulpoly([-1,-1,1],[1,-3,1]),[-1,1,1])
    if lhs!=rhs:raise ArithmeticError('golden moduli factorization')
    checks+=1
    # Polynomial identities over Q, proved here by elementary coefficient arithmetic.
    a2,a4,a6=1,-3,-2
    b2,b4,b6,b8=4*a2,2*a4,4*a6,4*a2*a6-a4*a4
    c4=b2*b2-24*b4
    delta=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
    if (c4,delta,Fraction(c4**3,delta))!=(160,2000,Fraction(2048)):raise ArithmeticError('golden rational model')
    checks+=1
    rows=[];tested=0;golden=[];above25=0
    primes=[p for p in range(3,102,2) if prime(p)]
    for p in primes:
        e_trace=-sum(chi((x+2)*(x*x-x-1),p) for x in range(p))
        for lam in range(2,p):
            rec=signature(p,(0,1,lam));audit(rec);literal=direct(p,(0,1,lam))
            checks+=1;tested+=1
            if [len(literal[s]) for s in SIGNS]!=[r['count'] for r in rec['cells']]:raise ArithmeticError('literal signature')
            if sum(len(v) for v in literal.values())!=p-3:raise ArithmeticError('lost field vertex')
            if rec['trace']**2>4*p:raise ArithmeticError('Hasse control')
            checks+=2
            if p>25:
                above25+=1;checks+=1
                if min(r['count'] for r in rec['cells'])<1:raise ArithmeticError('extension bound')
            if p<=19:
                point_count=1+sum((y*y-(x*(x-1)*(x-lam)))%p==0 for x in range(p) for y in range(p))
                checks+=1
                if p+1-point_count!=rec['trace']:raise ArithmeticError('independent point enumeration')
            # Affine change preserves labelled signs only after its character gauge.
            for scale,shift in ((1,2),(2,1)):
                anchors=tuple((shift+scale*x)%p for x in (0,1,lam))
                other=signature(p,anchors);ch=chi(scale,p)
                by_sign={tuple(r['signs']):r['count'] for r in other['cells']}
                checks+=1
                if other['trace']!=ch*rec['trace'] or any(by_sign[tuple(ch*x for x in s)]!=len(literal[s]) for s in SIGNS):
                    raise ArithmeticError('affine gauge lost')
            if (lam*lam-lam-1)%p==0 and p!=5:
                s=(2*lam-1)%p;d=(-s)%p
                j=256*pow(1-lam+lam*lam,3,p)*pow(lam*lam*(1-lam)**2,-1,p)%p
                checks+=1
                if j!=2048%p or e_trace!=chi(d,p)*rec['trace']:raise ArithmeticError('golden twist')
                conjugate=signature(p,(0,1,(1-lam)%p))
                checks+=1
                if conjugate['trace']!=chi(-1,p)*rec['trace']:raise ArithmeticError('conjugate twist')
                if chi(lam,p)!=chi(lam-1,p):raise ArithmeticError('golden anchor type')
                for u in range(p):
                    x=(lam-s*u)%p
                    if ((x+2)*(x*x-x-1)-d**3*u*(u-1)*(u-lam))%p:raise ArithmeticError('root coordinate map')
                checks+=p
                golden.append({'p':p,'lambda':lam,'trace':rec['trace'],'rational_model_trace':e_trace,
                               'twist':chi(d,p),'counts':[r['count'] for r in rec['cells']],
                               'relation':rec['relation']})
        if p in (11,19,31,59):rows.append({'p':p,'golden':[r for r in golden if r['p']==p]})
    inert=[]
    for p in (3,7,13,17,23):
        if chi(5,p)!=-1:raise ArithmeticError('inert test type')
        lam=(pow(2,-1,p),pow(2,-1,p));one=(1,0);total=0
        for x in product(range(p),repeat=2):
            val=f2_mul(f2_mul(x,f2_sub(x,one,p),p),f2_sub(x,lam,p),p)
            total+=f2_chi(val,p)
        trace=-total
        ap=-sum(chi((x+2)*(x*x-x-1),p) for x in range(p))
        expected=chi(-5,p)*(ap*ap-2*p)
        checks+=3
        if trace!=expected or f2_chi((p-1,0),p)!=1 or f2_pow(lam,2,p)!=((lam[0]+1)%p,lam[1]):
            raise ArithmeticError('inert descent / relation type')
        inert.append({'p':p,'q':p*p,'trace':trace,'rational_model_trace':ap,'relation':'graph',
                      'twist':chi(-5,p)})
    # Same j and same abstract transitive anchor tournament, different existence result.
    hostile=[]
    for anchors in ((0,1,4),(8,0,1)):
        rec=signature(11,anchors);pts=direct(11,anchors)[(1,1,1)]
        checks+=1
        if not all(chi(anchors[j]-anchors[i],11)==1 for i in range(3) for j in range(i+1,3)):
            raise ArithmeticError('transitive hostile anchors')
        lam=(anchors[2]-anchors[0])*pow(anchors[1]-anchors[0],-1,11)%11
        j=256*pow(1-lam+lam*lam,3,11)*pow(lam*lam*(1-lam)**2,-1,11)%11
        checks+=1
        if j!=2048%11:raise ArithmeticError('hostile must share j')
        checks+=1
        if rec['cells'][0]['boundary']!=4 or rec['cells'][-1]['boundary']!=4 or 8*len(pts)!=11-7-rec['trace']:
            raise ArithmeticError('transitive extreme-cell formula')
        hostile.append({'record':rec,'common_outneighbours':pts,'normalized_lambda':lam,'j':j})
    if [len(h['common_outneighbours']) for h in hostile]!=[1,0]:raise ArithmeticError('hostile not decisive')
    checks+=1
    good=signature(11,(0,1,4))
    import copy
    bads=[]
    bad=copy.deepcopy(good);bad['cells'][0]['count']=False;bads.append(bad)
    bad=copy.deepcopy(good);bad['trace']=-bad['trace'];bads.append(bad)
    bad=copy.deepcopy(good);bad['anchors']=[1,0,4];bads.append(bad)
    bad=copy.deepcopy(good);bad['cells'][0]['boundary']+=4;bads.append(bad)
    for bad in bads:
        checks+=1
        try:audit(bad)
        except ValueError:pass
        else:raise ArithmeticError('forged source accepted')
    for p,anchors in ((True,(0,1,2)),(9,(0,1,2)),(11,(0,1,1)),(11,(False,1,4)),(11,(0,1,4.0))):
        checks+=1
        try:signature(p,anchors)
        except ValueError:pass
        else:raise ArithmeticError('malformed source accepted')
    out={'status':'PROVED scoped identities with FINITE-EXACT controls; global targets OPEN',
         'universe':{'odd_primes_through':101,'all_normalized_triples':tested,'above25':above25,
                     'inert_prime_squares':[x['q'] for x in inert]},
         'checks':checks,'rational_curve':{'coefficients':[0,1,0,-3,-2],'discriminant':delta,'j':2048},
         'legendre_j2048_factors_ascending':[[-1,-1,1],[1,-3,1],[-1,1,1]],
         'golden_specializations':golden,'selected_tables':rows,'inert_specializations':inert,
         'same_j_same_tournament_hostile':hostile}
    path=Path(__file__).resolve().parents[2]/'05-knowledge/results/golden_legendre_ports_20261016c.json'
    path.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'checks':checks,'triples':tested,'golden':len(golden),'hostile_counts':[1,0]}))


if __name__=='__main__':main()
