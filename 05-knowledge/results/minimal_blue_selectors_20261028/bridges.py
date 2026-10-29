"""Exact support crossings, odd fibers, and retained symmetry controls."""
from collections import Counter
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import json
import blue_minimum as m
r=m.r;core=m.core;need=m.need
HERE=Path(__file__).resolve().parent


def atan_interval(q,N=60):
    value=sum((Fraction((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N)),Fraction(0))
    next_term=Fraction((-1)**N,(2*N+1)*q**(2*N+1))
    return min(value,value+next_term),max(value,value+next_term)


def decimal_bound(x,rounding):
    with localcontext() as ctx:
        ctx.prec=60;ctx.rounding=rounding
        return Decimal(x.numerator)/Decimal(x.denominator)


def support_crossings():
    a,b=atan_interval(5),atan_interval(239)
    pi=(16*a[0]-4*b[1],16*a[1]-4*b[0])
    intervals={}
    for rounding in (ROUND_FLOOR,ROUND_CEILING):
        with localcontext() as ctx:
            ctx.prec=60;ctx.rounding=rounding;total=Decimal(0)
            for j in range(40250):
                total+=Decimal(1)/Decimal(2*j+1)
                if j in (6,7,55,56,40248,40249):intervals[j,rounding]=total
    out=[]
    for label,j,target in [('2',7,(Fraction(2),)*2),('3',56,(Fraction(3),)*2),('2pi',40249,(2*pi[0],2*pi[1]))]:
        lo=decimal_bound(target[0],ROUND_FLOOR);hi=decimal_bound(target[1],ROUND_CEILING)
        need(intervals[j-1,ROUND_CEILING]<lo<=hi<intervals[j,ROUND_FLOOR],'strict certified support crossing')
        out.append(dict(threshold=label,m=j,denominator=2*j+1,
                        before=[str(intervals[j-1,v]) for v in (ROUND_FLOOR,ROUND_CEILING)],
                        after=[str(intervals[j,v]) for v in (ROUND_FLOOR,ROUND_CEILING)],
                        target=[str(lo),str(hi)]))
    return out


def odd_controls():
    out=[]
    for n in range(1,8):
        q,_=r.old.native_classes(n);fibers=Counter(q)
        need(all(v%2 for v in fibers.values()),'odd Hamilton presentation fibers')
        out.append(dict(order=n,classes=len(fibers),smallest=min(fibers.values()),largest=max(fibers.values())))
    t=tuple(sum(1<<((v+d)%7) for d in (1,2,4)) for v in range(7))
    G=r.oracle.isomorphisms(t,t,True)
    antis=[tuple((2*a-i)%7 for i in range(7)) for a in range(7)]
    for p in antis:
        need(all(core.edge(t,p[i],p[j])==core.edge(t,j,i) for i,j in combinations(range(7),2)),'anti-involution')
    centralizers=[sum(all(g[p[i]]==p[g[i]] for i in range(7)) for g in G) for p in antis]
    need(len(G)==21 and centralizers==[3]*7,'Paley anti-involution centralizers')
    frobenius=tuple(2*i%7 for i in range(7))
    need(frobenius in G and all(frobenius[antis[0][i]]==antis[0][frobenius[i]] for i in range(7)),'retained Frobenius phase')
    need(r.oracle.hamilton(t)==189,'Paley Hamilton count')
    return dict(hamilton_fibers=out,paley=dict(automorphisms=21,anti_involutions=7,
                centralizers=centralizers,frobenius=frobenius,Hamilton_paths=189))


def hostile_control():
    controls=r.old.source('minimum_hostile','05-knowledge/results/root_ports_closed_carriers_20261028/check_carriers.py')
    groups=controls.strong_controls()
    a,b=(controls.build_hostile(v) for v in (False,True))
    cyc=lambda t:{v for v in combinations(range(len(t)),3) if all(sum(core.edge(t,i,j) for j in v)==1 for i in v)}
    need(cyc(a)==cyc(b) and len(cyc(a))==67,'hostile labelled cyclic triples')
    need(r.oracle.hamilton(a)==r.oracle.hamilton(b)==462915,'hostile Hamilton counts')
    need([x['transporter_0_to_4'] for x in groups[-2:]]==[27,0],'hostile actual transporter')
    return dict(groups=groups[-2:],common_cyclic_triples=67,common_Hamilton_paths=462915)


def main():
    result=dict(status='FINITE-EXACT',support_crossings=support_crossings(),oddness=odd_controls(),hostile=hostile_control())
    (HERE/'bridge-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(denominators=[v['denominator'] for v in result['support_crossings']],
                         groups=[v['automorphisms'] for v in result['hostile']['groups']],
                         odd_classes=[v['classes'] for v in result['oddness']['hamilton_fibers']])));return result


if __name__=='__main__':main()
