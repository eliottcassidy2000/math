"""Exact 121-support decoder for actual 11+2 LRC entry.

The declared domain includes the actual atlas graph and physical Q^2 box.
No mathematical producer imports. Full mixed-support certificates are retained.
Independent controls enumerate literal bounded relations on small physical rows.
"""
import argparse
import hashlib
import json
import sys
from functools import reduce
from itertools import combinations, product
from math import gcd, isqrt
from pathlib import Path

sys.stdout.reconfigure(newline="\n")
Q=91**6
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise ArithmeticError(label)


def ceildiv(x,y):
    return -((-x)//y)


def allowed_sums():
    out={1}
    for p in range(2,357):
        if p%3!=2 or any(p%d==0 for d in range(2,isqrt(p)+1)):
            continue
        out|={m*p**e for m in tuple(out) for e in (1,2) if m*p**e<=356}
    return out


ATLAS_SUMS=allowed_sums()


def graph_components(speeds):
    adjacency=[set() for _ in speeds]
    edges=[]
    for i,j in combinations(range(len(speeds)),2):
        d=gcd(speeds[i],speeds[j])
        a,b=sorted((speeds[i]//d,speeds[j]//d))
        if a+b in ATLAS_SUMS:
            adjacency[i].add(j)
            adjacency[j].add(i)
            edges.append([i,j,a,b])
    components=[]
    unseen=set(range(len(speeds)))
    while unseen:
        component={min(unseen)}
        frontier=list(component)
        for i in frontier:
            for j in sorted(adjacency[i]-component):
                component.add(j)
                frontier.append(j)
        unseen-=component
        components.append(sorted(component))
    return sorted(components,key=lambda p:(len(p),p)),edges


def one_support(speeds,pair,distinguished,B,orientation):
    i,j=sorted(pair,key=lambda i:speeds[i])
    k=distinguished
    A,Bphys,Y=speeds[i],speeds[j],speeds[k]
    D=gcd(A,Bphys)
    a,b=A//D,Bphys//D
    need(1<=a<b<=B,"internal primitive height admits the minimal-coefficient theorem")
    delta=gcd(D,Y)
    c,x=D//delta,Y//delta
    r0=(pow(a,-1,b)*x)%b
    s0=(x-a*r0)//b
    lo=max(ceildiv(-B-r0,b),ceildiv(s0-B,a))
    hi=min((B-r0)//b,(s0+B)//a)
    found=c<=B and lo<=hi
    certificate=dict(orientation=orientation,pair_indices=[i,j],distinguished_index=k,
                     pair_gcd=D,a=a,b=b,delta=delta,c=c,x=x,r0=r0,s0=s0,
                     lower=lo,upper=hi,crossing=found)
    if found:
        r,s=r0+lo*b,s0-lo*a
        vector=[0]*len(speeds)
        vector[i],vector[j],vector[k]=-r,-s,c
        need(sum(v*w for v,w in zip(vector,speeds))==0,"returned crossing is a literal zero relation")
        need(0<max(map(abs,vector))<=B and sum(v!=0 for v in vector)<=3,"returned coefficient height and support")
        certificate['witness']=vector
    return certificate


def split_test(speeds,core,outer,B):
    heights=[]
    for indices in (core,outer):
        for i,j in combinations(indices,2):
            height=max(speeds[i],speeds[j])//gcd(speeds[i],speeds[j])
            heights.append(dict(indices=[i,j],height=height))
    failures=[h for h in heights if h['height']>B]
    if failures:
        return dict(entry=False,reason='internal_height_failure',internal_heights=heights,
                    failed_heights=failures,supports=[])
    supports=[]
    for pair in combinations(core,2):
        for k in outer:
            supports.append(one_support(speeds,pair,k,B,'two_core_one_pair'))
    for pair in combinations(outer,2):
        for k in core:
            supports.append(one_support(speeds,pair,k,B,'one_core_two_pair'))
    return dict(entry=not any(s['crossing'] for s in supports),reason='complete_support_test',
                internal_heights=heights,supports=supports)


def entry_decoder(speeds):
    # Out-of-domain is distinct from a proved negative entry verdict.
    base=dict(speeds=list(speeds),Q=Q)
    if len(speeds)!=13 or any(not isinstance(v,int) or v<=0 for v in speeds):
        return dict(base,domain=False,reason='requires_thirteen_positive_integers')
    if len(set(speeds))!=13 or reduce(gcd,speeds)!=1:
        return dict(base,domain=False,reason='requires_distinct_primitive_row')
    if sum(speeds)>Q**2:
        return dict(base,domain=False,reason='outside_physical_box')
    components,edges=graph_components(speeds)
    if list(map(len,components))!=[2,11]:
        return dict(base,domain=False,reason='outside_actual_11_plus_2_graph',components=components)
    outer,core=components
    need(all(b<=355<Q for i,j,a,b in edges),"actual decoder edge rows belong to the bounded relation span")
    answer=split_test(speeds,core,outer,Q)
    if answer['reason']=='complete_support_test':
        need(len(answer['supports'])==121,"both crossing orientations retained")
        need(sum(s['orientation']=='two_core_one_pair' for s in answer['supports'])==110,"110 first-orientation supports")
        need(sum(s['orientation']=='one_core_two_pair' for s in answer['supports'])==11,"eleven opposite-orientation supports")
    return dict(base,domain=True,components=components,decoder_edges=edges,**answer)


def literal_relations(speeds,B):
    # Independent complete enumeration: solve for the third coefficient after
    # enumerating the first two, without gcd normalization or interval tests.
    out=set()
    n=len(speeds)
    for i,j,k in combinations(range(n),3):
        for ri,rj in product(range(-B,B+1),repeat=2):
            numerator=-ri*speeds[i]-rj*speeds[j]
            if numerator%speeds[k]:
                continue
            rk=numerator//speeds[k]
            if abs(rk)>B or not (ri or rj or rk):
                continue
            row=[0]*n
            row[i],row[j],row[k]=ri,rj,rk
            out.add(tuple(row))
    return out


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-certificates')
    args=parser.parse_args()
    need(sum(1 for s in ATLAS_SUMS for p in range(1,(s+1)//2) if gcd(p,s)==1)==5855,
         "exact inherited all-scale atlas cardinality")
    universe=partitions=relation_count=0
    for B in range(4,8):
        for speeds in combinations(range(1,10),5):
            if sum(speeds)>B**2 or reduce(gcd,speeds)!=1:
                continue
            relations=literal_relations(speeds,B)
            relation_count+=len(relations)
            for core in combinations(range(5),3):
                outer=tuple(i for i in range(5) if i not in core)
                literal=all(sum(row[i]*speeds[i] for i in core)==0 for row in relations)
                decoded=split_test(speeds,core,outer,B)
                need(decoded['entry']==literal,"full literal bounded span versus exact two-component criterion")
                if decoded['supports']:
                    need(len(decoded['supports'])==9,"all mixed supports in the five-label control")
                    for c in decoded['supports']:
                        if c['crossing']:
                            need(tuple(c['witness']) in relations,"returned crossing occurs in independent full relation bank")
                partitions+=1
            universe+=1
    # Explicit positive small-height component-kernel controls. All bounded
    # relation rows are still literally enumerated; no geometric LRC input.
    for B,speeds,core in [(30,(1,3,181,543),(0,1)),(40,(1,3,4,321,963),(0,1,2))]:
        outer=tuple(i for i in range(len(speeds)) if i not in core)
        need(sum(speeds)<=B**2,"positive toy control remains in its physical box")
        relations=literal_relations(speeds,B)
        need(all(sum(row[i]*speeds[i] for i in core)==0 for row in relations),"literal positive component-kernel equality")
        need(split_test(speeds,core,outer,B)['entry'],"exact positive toy decoder")

    U=(1,4,6,8,10,12,14,15,16,18,22)
    V=(2,3,4,5,6,10,12,14,15,20,30)
    cases=[]
    for name,shape,t,g in [('canonical_unit',U,1,2**45),
                           ('unitless',V,1,60*Q+1),
                           ('opposite_orientation_hostile',U,3*Q+1,1)]:
        speeds=tuple(t*v for v in shape)+(g,3*g)
        answer=entry_decoder(speeds)
        need(answer['domain'],"named control has actual graph and physical box")
        selected=[s for s in answer['supports'] if s['orientation']=='two_core_one_pair']
        opposite=[s for s in answer['supports'] if s['orientation']=='one_core_two_pair']
        if name!='opposite_orientation_hostile':
            need(answer['entry'] and g>2*Q*max(shape),"positive entry agrees with independent all-orientation dominance")
        else:
            need(not answer['entry'] and not any(s['crossing'] for s in selected),"110 tests alone falsely accept the hostile")
            positives=[s for s in opposite if s['crossing']]
            need(len(positives)==1 and positives[0]['distinguished_index']==0,"only the physical unit-core column supplies the omitted crossing")
            expected=[0]*13
            expected[0],expected[11],expected[12]=1,-1,-Q
            need(positives[0]['witness']==expected,"exact lost relation t-1-3Q=0")
            need(sum(speeds)<=Q**2 and t%3==1,"hostile finite box and uniform first-orientation coefficient gate")
        answer['name']=name
        cases.append(answer)

    # Actual connected decoder plus W=V_dec outside the physical box: this
    # invalidates extending the necessary internal-height gate beyond scope.
    powers=tuple(355**i for i in range(11))
    g=2*Q*max(powers)+1
    speeds=powers+(g,3*g)
    components,edges=graph_components(speeds)
    need(components==[[11,12],list(range(11))],"outside-box hostile retains the actual decoder graph")
    need(max(powers)>Q and g>2*Q*max(powers),"large internal height coexists with crossing exclusion")
    need(all(any({e[0],e[1]}=={i,i+1} for e in edges) for i in range(10)),"bounded atlas edges still span the whole core kernel")
    need(sum(speeds)>Q**2 and reduce(gcd,speeds)==1,"only the physical box is deliberately dropped")
    outside=entry_decoder(speeds)
    need(not outside['domain'] and 'entry' not in outside,"out-of-domain does not become a false negative entry verdict")
    outside['name']='physical_box_firewall'
    cases.append(outside)

    certificate=dict(Q=Q,scope='actual primitive thirteen-speed Q^2 box and actual decoder graph11+2',cases=cases)
    serialized=(json.dumps(certificate,indent=2,sort_keys=True)+'\n').encode('utf-8')
    if args.write_certificates:
        Path(args.write_certificates).write_bytes(serialized)
    print('STATUS: PASS; exact actual-entry decoder, both mixed orientations')
    print('LITERAL UNIVERSE:',universe,'physical rows Qtoy4..7;',partitions,'3+2 splits;',relation_count,'complete bounded relation rows')
    print('POSITIVE TOY CONTROLS: Q30 four labels and Q40 five labels, literal full relation spans')
    for c in cases:
        supports=c.get('supports',[])
        print('CASE:',c['name'],'domain=',c['domain'],'entry=',c.get('entry'),'supports=',len(supports),
              'crossings=',sum(s['crossing'] for s in supports))
    print('HOSTILE: all110 first-orientation tests negative, omitted relation t-1-3Q=0')
    print('FIREWALL: outside Q^2 box, actual W=V_dec can have an internal primitive pair height>Q')
    print('CERTIFICATE SHA256:',hashlib.sha256(serialized).hexdigest())
    print('SCOPE:121 fixed-support arithmetic tests; bit cost depends on integer height; no LRC14 closure')
    print('ACTIVE GATES:',GATES)


if __name__=='__main__':
    main()
