#!/usr/bin/env python3
"""Independent proof controls and complete endpoint-fixed W census."""
from core import *
from math import lcm

def order(g):
    seen=set();value=1
    for v in range(len(g)):
        if v in seen:continue
        u=v;length=0
        while u not in seen:
            seen.add(u);length+=1;u=g[u]
        need(u==v,'permutation cycle closure')
        value=lcm(value,length)
    return value

def gate(rows,W,all_replay=False):
    n=len(rows);z=n-1;P=list(range(n));C=list(range(1,n))+[0]
    B=move(rows,P);S=move(rows,C)
    k=rows[z].bit_count();d1=rows[1].bit_count()
    I={v for v in range(1,z) if arc(rows,v,0)};J=I|{1}
    need(all(arc(B,a,b)==arc(S,a,b) for a,b in combinations(range(1,n),2)), 'B-0 equals S-0')
    need(rows[0].bit_count()==n-1-k,'source anti endpoint degrees')
    need(S[1].bit_count()==n-d1 and S[0].bit_count()==k-1,'actual intermediate degrees')
    need(sorted(W)==list(range(1,n)) and W[0]==1 and W[-1]==z,'endpoint-fixed W')
    need(all(arc(B,a,b) for a,b in zip(W,W[1:])),'W directed in B-0')
    Q2=list(W)+[0];Q3=[0]+list(W)
    r=dict(zip(Q3,reversed(Q3)));g=[r[z-v] for v in range(n)]
    need(g[0]==0 and g[z]==z,'reflection product fixes endpoints')
    auto=all(arc(B,a,b)==arc(B,g[a],g[b]) for a,b in combinations(range(n),2))
    q3blue=mismatch(B,Q3) is None
    need(auto==q3blue,'exact blue/automorphism equivalence')
    even=order(g)%2==0;port_mismatch={g[v] for v in J}!=J
    need(not(even and auto),'tournament automorphism has odd order')
    need(not(port_mismatch and auto),'automorphism preserves exact endpoint neighborhood')
    if d1!=k:
        need(mismatch(rows,C) is not None and mismatch(S,Q2) is not None,'degree test for both early words')
    sufficient=d1!=k and (even or port_mismatch)
    if all_replay or sufficient:
        flags,trace,target=literal_ledger(rows,[C,Q2,Q3],sufficient)
        need(flags[2]==q3blue and target==B,'actual last source and target')
    else:
        flags=[mismatch(rows,C) is None,mismatch(S,Q2) is None,q3blue]
    return {'degree_first_neighbor':d1,'k':k,'degree_guard':d1!=k,
            'permutation_order':order(g),'port_set_J':sorted(J),'port_mismatch':port_mismatch,
            'sufficient_gate':sufficient,'actual_blue_flags':flags}

def all_W(rows):
    n=len(rows);B=move(rows,list(range(n)))
    tested=0;accepted=[]
    for middle in permutations(range(2,n-1)):
        tested+=1;W=(1,)+middle+(n-1,)
        if all(arc(B,a,b) for a,b in zip(W,W[1:])):
            accepted.append(W)
    return tested,accepted

def main():
    counts=Counter();small=[]
    for n in range(3,8):
        s=Linear(n);s.anti(list(range(n)))
        for v in range(n-1):s.force(v,v+1,1)
        s.force(n-1,0,1)
        need(s.ok,'nonempty generic closed blue source')
        nc=Counter()
        for rows in s.models():
            nc['closed_blue_sources']+=1
            tried,Ws=all_W(rows);nc['endpoint_fixed_permutations']+=tried
            for W in Ws:
                r=gate(rows,W,True);nc['directed_endpoint_fixed_words']+=1
                nc['sufficient_gate_hits']+=r['sufficient_gate']
                nc['actual_three_black']+=not any(r['actual_blue_flags'])
        small.append({'n':n,'source_dimension':s.dimension(),**dict(nc)})
    data_bytes=Path(__file__).with_name('certificates.json').read_bytes()
    data=json.loads(data_bytes);selected=[]
    for item in data:
        n=item['n'];q1,q2,q3=item['words']
        if q1==list(range(1,n))+[0] and q2[-1]==0 and q3==[0]+q2[:-1]:
            rows=tuple(item['source_rows']);r=gate(rows,tuple(q2[:-1]),True)
            need(not any(r['actual_blue_flags']),'supplied cyclic repair literally black')
            selected.append({'n':n,'ports':item['ports'],'source_rows':rows,**r})
    equality=(126,156,440,1008,1248,1474,1926,1793,1555,1079,15)
    need(equality[1].bit_count()==equality[-1].bit_count()==4,'named degree equality hostile')
    tested,Ws=all_W(equality);hist=Counter();examples={}
    for W in Ws:
        r=gate(equality,W,True);key=','.join(str(int(x)) for x in r['actual_blue_flags'])
        hist[key]+=1;examples.setdefault(key,W)
        need(not r['degree_guard'] and not r['sufficient_gate'],'equality is outside sufficient degree guard')
    need(tested==40320 and len(Ws)==132 and hist=={'0,0,0':127,'0,0,1':5},'complete equality W census')
    noncyclic=next(item for item in data if item['n']==11 and item['words'][0][0]==0)
    R=tuple(noncyclic['source_rows']);tried,Vs=all_W(R);hist2=Counter()
    for W in Vs:
        r=gate(R,W,True);hist2[','.join(str(int(x)) for x in r['actual_blue_flags'])]+=1
    print(json.dumps({'status':'PASS','scope':'Audited all-size sufficient bridge criterion, not a necessary criterion or a W-existence theorem; finite controls explicitly bounded',
       'certificates_sha256':sha256(data_bytes).hexdigest(),'small_exhaustive_controls':small,
       'supplied_cyclic_certificates':selected,'equality_hostile':{'source':equality,'endpoint_fixed_permutations':tested,'directed_W':len(Ws),'actual_blue_flag_histogram':dict(hist),'examples':examples},
       'noncyclic_source_diagnostic':{'source':R,'endpoint_fixed_permutations':tried,'directed_W':len(Vs),'actual_blue_flag_histogram':dict(hist2)}},indent=2,sort_keys=True))

if __name__=='__main__':main()
