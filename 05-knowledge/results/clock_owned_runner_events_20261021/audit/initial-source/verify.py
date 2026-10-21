#!/usr/bin/env python3
"""Exact replay for session19 clock-owned runner certificates. Standard library only."""
from collections import defaultdict
from fractions import Fraction as F
from math import gcd, prod
import json


def dist(x):
    z=x%1
    return min(z,1-z)


def safe_fraction(V,t):
    return all(dist(v*t)>=F(1,len(V)+1) for v in V)


def safe_integer(V,t):
    a,b=t.numerator,t.denominator
    n=len(V)+1
    return all(n*min((v*a)%b,(-v*a)%b)>=b for v in V)


def event_tape(V):
    """Ties are simultaneous, and endpoints themselves belong to no owner tooth."""
    n=len(V)+1
    events=defaultdict(lambda:{'exit':[],'entry':[]})
    for v in V:
        for m in range(v):
            events[F(n*m+1,n*v)%1]['exit'].append(v)
            events[F(n*m-1,n*v)%1]['entry'].append(v)
    times=sorted(events)
    debt=len(V)
    records=[]
    for p,e in enumerate(times):
        row=events[e]
        point=debt-len(row['exit'])
        after=point+len(row['entry'])
        assert point==sum(dist(v*e)<F(1,n) for v in V)
        nxt=times[p+1] if p+1<len(times) else times[0]+1
        mid=(e+nxt)/2
        assert after==sum(dist(v*mid)<F(1,n) for v in V)
        assert (point==0)==safe_integer(V,e)==safe_fraction(V,e)
        records.append({'time':str(e),'exit':row['exit'],'entry':row['entry'],
                        'before':debt,'point':point,'after':after})
        debt=after
    assert debt==len(V)
    return records


def critical(V):
    k=len(V);n=k+1;q=2*k+1
    assert len(set(V))==k and min(V)>0
    assert all(gcd(v,q)==1 for v in V)
    grid=[a for a in range(1,q) if safe_integer(V,F(a,q))]
    classes={min(v%q,q-v%q) for v in V}
    if len(classes)<k:
        assert grid
        return {'q':q,'grid':grid,'repairs':[]}
    assert not grid
    rows=[]
    for v in V:
        a=pow(v,-1,q);m=(a*v-1)//q
        for sign,D in [(1,k),(-1,3*k+2)]:
            t=F(n*m+sign,n*v)
            actual=[w for w in V if dist(w*t)<F(1,n)]
            corridor=[];first_tooth=[];timed=[]
            for w in V:
                if w==v:continue
                r=(a*w)%q
                assert 2<=r<=q-2
                c=n*(q-r)-q if sign==1 else n*r-q
                entry=F(a,q)+sign*F(c,n*q*w)
                if D*w>c*v:
                    corridor.append(w)
                    timed.append({'target':w,'first_entry':str(entry)})
                if c*v<D*w<(c+2*q)*v:
                    first_tooth.append(w)
            assert set(first_tooth)<=set(actual)<=set(corridor)
            assert safe_integer(V,t)==safe_fraction(V,t)==(not actual)
            rows.append({'owner':v,'clock':str(F(a,q)),'sign':sign,
                         'endpoint':str(t),'actual_blockers':actual,
                         'corridor_blockers':corridor,'first_tooth_blockers':first_tooth,
                         'timed_entries':timed})
    forward={r['owner']:set(r['corridor_blockers']) for r in rows if r['sign']==1}
    assert all(v not in forward[w] for v in V for w in forward[v])
    return {'q':q,'grid':grid,'repairs':rows}


def prime(q):
    return q>=2 and all(q%d for d in range(2,int(q**0.5)+1))


def main():
    result={}
    pair_count=0
    for q in range(5,102,2):
        k=(q-1)//2;n=k+1
        for a in range(2,q-1):
            if gcd(a,q)>1:continue
            b=pow(a,-1,q)
            assert 2<=b<=q-2
            product=(n*a-q)*(n*b-q)
            assert product>=k*k
            assert (product==k*k)==({a,b}=={2,n})
            pair_count+=1
    result['inverse_pair_checks']=pair_count
    ap=[]
    for k in range(2,14):
        V=list(range(1,k+1));tape=event_tape(V)
        assert all(r['after']>0 for r in tape)
        at=next(r for r in tape if F(r['time'])==F(1,k+1))
        assert at['point']==0 and at['exit']==[1] and at['entry']==[k]
        ap.append({'k':k,'tied_witness':at})
    result['tight_AP_controls']=ap
    result['composite_unit_positive']=critical([1,2,4,5])
    triangle=[11,19,50]
    tr=critical(triangle)
    assert all(r['actual_blockers'] for r in tr['repairs'] if r['sign']==1)
    assert any(not r['actual_blockers'] for r in tr['repairs'] if r['sign']==-1)
    assert safe_integer(triangle,F(27,200))
    result['forward_triangle']=tr
    V=[4,5,9,12,14]
    hostile=critical(V)
    assert all(r['actual_blockers'] for r in hostile['repairs'])
    assert safe_integer(V,F(1,24))
    tape=event_tape(V)
    assert next(F(r['time']) for r in tape if r['point']==0)==F(1,24)
    result['two_sided_hostile']={'speeds':V,'witness':'1/24','certificate':hostile}
    families=[]
    for q in [11,19,43,59,67,83]:
        assert prime(q) and q%8==3
        k=(q-1)//2
        H=sorted({a*a%q for a in range(1,q)})
        assert len(H)==k and set(H).isdisjoint({q-h for h in H})
        assert {(-2*h)%q for h in H}==set(H)
        assert {(4*h)%q for h in H}==set(H)
        for M in [10,100,1000]:
            V=[q*M+h for h in H]
            assert F(max(V),min(V))<F(3*k+2,2*k+3)
            cert=critical(V)
            assert all(r['actual_blockers'] for r in cert['repairs'])
            by_residue={v%q:v for v in V}
            for row in cert['repairs']:
                target=by_residue[((-2 if row['sign']==1 else 4)*row['owner'])%q]
                assert target in row['first_tooth_blockers']
            t=F(1,2*max(V));margin=F(min(V),2*max(V))
            assert safe_integer(V,t) and all(dist(v*t)>=margin for v in V)
            families.append({'q':q,'k':k,'M':M,'min':min(V),'max':max(V),
                             'ordinary_witness':str(t),'margin':str(margin),
                             'all_repair_endpoints_fail':True})
    result['quadratic_residue_families']=families
    # A fixed-row corridor is NOT its endpoint graph: 50 crosses and later re-exits.
    row=next(r for r in tr['repairs'] if r['owner']==11 and r['sign']==-1)
    assert row['corridor_blockers'] and not row['actual_blockers']
    result['scope']='Exact certificate tests, not an LRC census or general proof.'
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
