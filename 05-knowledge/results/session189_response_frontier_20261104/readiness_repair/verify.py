#!/usr/bin/env python3
"""Targeted independent controls for incoming commit a92e0b49; no census claim."""
from pathlib import Path
from itertools import permutations
import json,hashlib,subprocess
P=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def valid(rows):return all(not(row>>i&1)and all(((row>>j)&1)+((rows[j]>>i)&1)==1 for j in range(len(rows))if i!=j)for i,row in enumerate(rows))
def strong(rows):
    n=len(rows);reach=1
    while True:
        nxt=reach
        for i in range(n):
            if reach>>i&1:nxt|=rows[i]
        if nxt==reach:break
        reach=nxt
    if reach!=(1<<n)-1:return False
    rev=[sum(1<<j for j in range(n)if rows[j]>>i&1)for i in range(n)];reach=1
    while True:
        nxt=reach
        for i in range(n):
            if reach>>i&1:nxt|=rev[i]
        if nxt==reach:break
        reach=nxt
    return reach==(1<<n)-1
def endpoints(rows):
    n=len(rows);d=[[0]*n for _ in range(1<<n)]
    for mask in range(1,1<<n):
        for j in range(n):
            if mask>>j&1:
                prev=mask^(1<<j);d[mask][j]=1 if not prev else sum(d[prev][i]for i in range(n)if prev>>i&1 and rows[i]>>j&1)
    return d[-1]
def delete(rows,a):
    ids=[i for i in range(len(rows))if i!=a];return tuple(sum(1<<j for j,v in enumerate(ids)if rows[u]>>v&1)for u in ids),ids
def C2(rows,a):
    n=len(rows);out=list(rows)+[1<<a,1<<n]
    for i in range(n):
        if i!=a:out[i]|=1<<n
        out[i]|=1<<(n+1)
    return tuple(out)
def cycles(rows):
    n=len(rows);return sum(all(rows[p[i]]>>p[(i+1)%n]&1 for i in range(n))for tail in permutations(range(1,n))for p in [(0,)+tail])
def phi(rows,a):
    root=rows[a].bit_length()-1;need(rows[a]==1<<root,'sole exit');out=list(rows);ids=[i for i in range(len(rows))if i!=a]
    for i in ids:
        for j in ids:
            if i!=j:out[i]^=1<<j
    return tuple(out),root

def main():
    witness=(512,17,11,19,1797,31,63,127,239,1518,2543,1023);a=0;root=9;need(valid(witness)and strong(witness),'strong incoming witness');need(witness[a]==1<<root,'exact sole exit');U,ids=delete(witness,a);ends=endpoints(U);digit=ends[ids.index(root)];h=sum(endpoints(witness));need((h,digit)==(12877,2),'incoming H/digit')
    child=C2(witness,a);need(valid(child)and strong(child),'strong C2 child');child_H=sum(endpoints(child));need(child_H==38635==3*h+2*digit,'actual +4 count')
    # Three-name transitive prefix gives the right order and count, but is not strong.
    m=3;n=len(witness);prefix=tuple(sum(1<<j for j in range(i+1,m+n))for i in range(m));padded=prefix+tuple(row<<m for row in witness);pa=a+m;pr=root+m;need(valid(padded)and not strong(padded),'weak padding boundary');pad_H=sum(endpoints(padded));pu,pids=delete(padded,pa);pad_digit=endpoints(pu)[pids.index(pr)];need((pad_H,pad_digit)==(h,digit),'prefix preserves H and digit')
    # Explicit marked counterexample to the claimed onto-bijection in Corollary2.1.
    Y=(2,4,9,3);mark=1;T,r=phi(Y,mark);need(valid(Y)and strong(Y)and Y[mark]==1<<r and cycles(Y)==1,'counterexample target hypotheses');need(T==(14,4,0,6)and not strong(T),'counterexample inverse not strong');Yu,_=delete(Y,mark);need(not strong(Yu),'missing deleted-core hypothesis')
    eligible=[v for v,row in enumerate(Y)if row.bit_count()==1 and strong(delete(Y,v)[0])];need(eligible==[0],'canonical admissible exit repairs this unmarked class');ready,rr=phi(Y,eligible[0]);ru,rids=delete(ready,eligible[0]);need(strong(ready)and endpoints(ru)[rids.index(rr)]==1 and phi(ready,eligible[0])[0]==Y,'canonical inverse is strong ready and involutive')
    tri=(2,4,1);need(not any(strong(delete(tri,v)[0])for v in range(3)),'order-three boundary')
    Y3=(14,20,24,18,1);need(valid(Y3)and strong(Y3)and cycles(Y3)==3,'digit-three boundary graph');need(not any(row.bit_count()==1 and strong(delete(Y3,v)[0])for v,row in enumerate(Y3)),'not an arbitrary-cycle-count unmarked repair')
    record={'incoming_commit':'a92e0b49befb6ba1ac3fd4f42cca5130db76033b','verified_witness':{'rows':witness,'order':len(witness),'H':h,'a':a,'root':root,'digit':digit,'strong':True},'C2_child':{'rows':child,'order':len(child),'H':child_H,'strong':True},'three_prefix_padding':{'rows':padded,'order':len(padded),'H':pad_H,'a':pa,'root':pr,'digit':pad_digit,'strong':False},'marked_bijection_counterexample':{'Y_rows':Y,'Y_strong':True,'sole_exit':mark,'root':r,'Hamilton_cycles':cycles(Y),'Phi_inverse_rows':T,'inverse_strong':False,'deleted_core_strong':False,'canonical_admissible_exit':eligible[0],'canonical_ready_inverse_rows':ready},'boundary_cycle_three':{'rows':Y3,'Hamilton_cycles':3,'admissible_exit_count':0},'scope':'No ready-census completeness, density limit, all-order heredity, or order15 strong digit2 target is certified.'}
    (P/'controls.json').write_text(json.dumps(record,sort_keys=True,indent=2)+'\n');receipt={'status':'PASS','verified_claims':['order12 strong H12877 digit2 witness','order14 strong C2 H38635 child','order15 weak H12877 digit2 padded witness'],'counterexample':'Corollary2.1 is not onto all marked strong sole-exit tournaments with the specified cycle count; deletion-core strongness is missing.','sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['verify.py','INCOMING.json']},'incoming_source_pins':json.loads((P/'INCOMING.json').read_text()),'artifacts':{'controls.json':hashlib.sha256((P/'controls.json').read_bytes()).hexdigest()}}
    (P/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
