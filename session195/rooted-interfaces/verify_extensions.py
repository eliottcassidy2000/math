#!/usr/bin/env python3
"""Additional positive extension and candidate equal-response-barrier controls.
The all-order claims require the accompanying proofs; finite controls do not.
"""
import json
from collections import deque
from pathlib import Path
from verify import dp,literal,require,seam,signature,transfer,triple,census,tournaments

HERE=Path(__file__).resolve().parent
GADGETS=[
    (67,(52,17,2,7,12,30),3,[( (0,),0)]),
    (163,(44,1,2,6,15,30),3,[((0,1),0)]),
    (173,(56,1,3,6,14,30),1,[((0,),0)]),
    (247,(48,1,3,7,14,30),2,[((0,1),0)]),
    (511,(42,16,3,6,13,30),3,[((0,1,3),2)]),
]
# H67 uses the singleton seed directly, not a singleton preliminary move.
GADGETS[0]=(67,GADGETS[0][1],3,[])

def strong(a):
    n=len(a)
    for v in range(n):
        reach=1<<v
        while True:
            new=reach
            for u in range(n):
                if reach>>u&1:new|=a[u]
            if new==reach:break
            reach=new
        if reach!=(1<<n)-1:return False
    return True

def rooted_closure(seed,bound,bank):
    seen={seed};queue=deque([seed]);values={seed[0]}
    while queue:
        x=queue.popleft()
        for y in bank:
            z=transfer(x,y)
            if z[0]>bound:continue
            values.add(z[0])
            if z[2]<=bound and z not in seen:seen.add(z);queue.append(z)
    return seen,values

def main():
    bank={}
    for z in json.loads((HERE/'inventory.json').read_text()):
        a=tuple(z['rows']);b=z['root'];p=dp(a)[1][b]
        require(p==literal(a)[b],'base representative mismatch')
        require(p[:7]==tuple(z['signature']) and p[7]==z['H_deleted_root'],'base metadata mismatch')
        bank[p[:7]]=(a,b,p)
    records=[]
    for N,a,b,prefix in GADGETS:
        p=dp(a)[1][b];require(p==literal(a)[b],'extension literal mismatch')
        bank[p[:7]]=(a,b,p)
        T=(0,);r=0;steps=[]
        for g,c in prefix+[(a,b)]:
            T,r=seam(T,r,g,c);steps.append({'rows':list(g),'root':c})
        q=dp(T)[1][r]
        require(q[0]==N,'extension output count fails')
        records.append({'H':N,'gadget_rows':list(a),'gadget_root':b,'gadget_profile':list(p),'output_rows':list(T),'output_root':r,'output_state':[q[0],q[3],q[4],q[7]],'steps':steps})
    require(len(bank)==70,'extended bank size mismatch')
    states,vals,*_=census(1000,bank)
    missing=sorted(set(range(1,1001,2))-vals)
    require(missing==[7,21],'extended small census mismatch')
    # Finite controls of the equal-response barrier on every rooted labelled
    # strong graph of orders <=5 whose profile has R=h>=5, and the primitive49.
    seeds=set();rooted_count=0
    for n in range(1,6):
        for a in tournaments(n):
            if not strong(a):continue
            for p in dp(a)[1]:
                x=triple(p)
                if x[0]>=5 and x[0]==x[2]:seeds.add(x);rooted_count+=1
    A=(26,16,3,6,12);r=2;T,r=seam(A,r,(0,1),0)
    p=dp(T)[1][r];primitive=triple(p)
    require(primitive==(49,13,49) and strong(T),'primitive equal-response graph fails')
    seeds.add(primitive)
    total=0
    for seed in seeds:
        ss,vs=rooted_closure(seed,3*seed[0]+4,bank)
        require(vs=={seed[0],2*seed[0]-seed[1],3*seed[0]},'equal-response low-spectrum classification counterexample')
        total+=len(ss)
    result={'status':'PASS','extended_bank_size':len(bank),'bound':1000,'retained_states':len(states),'missing_odds':missing,'five_old_holes_filled':[z['H'] for z in records],'equal_response_finite_inputs':rooted_count,'equal_response_distinct_states':len(seeds),'equal_response_search_states':total,'primitive_hostile':{'rows':list(T),'root':r,'triple':list(primitive),'excluded_targets':[149,151]},'scope':'Five literal gadgets extend the base65-matrix3-state grammar to70 matrices; equal-response finite checks are controls, not an all-order proof.'}
    (HERE/'five_gadget_witnesses.json').write_text(json.dumps(records,indent=2)+'\n')
    (HERE/'extension_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
