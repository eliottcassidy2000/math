#!/usr/bin/env python3
"""Independent structural and short-word controls; not a target census."""
from pathlib import Path
import argparse, hashlib, importlib.util, json
HERE=Path(__file__).resolve().parent
HELPER=HERE.parent/'session195-four-state-audit'/'verify_independent.py'
if not HELPER.exists():HELPER=HERE.parent/'combined-growth-audit'/'verify_independent.py'
spec=importlib.util.spec_from_file_location('independent_profiles',HELPER)
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)

def ready(q):
    h,e,R,d=q;return e==d and abs(R-h-d)==1

def gap(q):
    h,e,R,d=q;c=R-h;f=c+e-d
    return h>=9 and d>=3 and all(x==0 or x>=3 for x in (e,c,f))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,required=True);args=parser.parse_args()
    raw=json.loads((args.source/'inventory67.json').read_text());bank=[]
    for rec in raw:
        a=rec['rows'];r=rec['root'];sig=tuple(rec['signature'])+(rec['H_deleted_root'],)
        v.check(v.signature(a,r)==sig,'physical gadget profile')
        g,s,L,E,Q,B,delta,dg=sig;U=Q+E-g-dg;V=B-L
        labels=[i for i in range(len(a)) if i!=r]
        deleted=[sum(1<<j for j,y in enumerate(labels) if a[x]>>y&1) for x in labels]
        missing=sum(v.signature(deleted,j)[2] for j,x in enumerate(labels) if a[r]>>x&1)
        v.check(V==missing and U>=0,'missing-first-bin identity')
        v.check(2*Q+2*B-g-L-2*dg==(g-E)+(L-E)+2*U+2*V>=0,'DK nonnegative coefficient')
        if a[r]:v.check(V>=dg,'non-sink V bound')
        else:v.check(U==V==0 and g==E==Q==dg and s==delta,'sink identities')
        bank.append(sig)
    states=set();endpoint_checks=endpoint_gap_checks=0
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            a=v.graph(n,code);tab=v.tables(a)
            for r in range(n):
                q=v.state(a,r,tab);states.add(q);h,e,R,d=q
                if e==d:
                    v.check((a[r]==0 and q==(d,d,d,d)) or h>=3*d,'endpoint-complete root ratio')
                    endpoint_checks+=1
                if h>=9:
                    v.check(h-e==0 or h-e>=3,'physical endpoint gap');endpoint_gap_checks+=1
    a=[4,1,2];r=0;q=v.state(a,r)
    for op in ('A','D','K'):a,r=v.unary(a,r,op);q=v.step(q,op)
    v.check(q==(17,3,21,3),'positive d1 reset');states.add(q)
    mature=[x for x in states if ready(x) and x[3]>=3]
    v.check((13,3,15,3) in mature and (17,3,21,3) in mature,'both ready boundary controls')
    gap_inputs=[q for q in states if gap(q)];invariant_checks=entry_checks=0
    for q in gap_inputs:
        for op in ['K','D']+bank:
            child=v.step(q,op);v.check(gap(child),'gap class invariant');invariant_checks+=1
    for q in mature:
        v.check(gap(v.step(q,'D')),'first D enters gap');entry_checks+=1
        for op in bank:
            child=v.step(q,op);v.check(gap(child),'first gadget enters gap');entry_checks+=1
    suffix_checks=0
    for q in mature:
        for op in bank:
            layer={v.step(q,op)}
            for length in range(9):
                for t in layer:v.check(not ready(t),'one gadget plus padding returned');suffix_checks+=1
                layer={v.step(t,letter) for t in layer for letter in ('K','D')}
    # Whole-word controls, retaining whether a non-K operation has occurred.
    layer={(q,False) for q in mature};word_checks=0
    for length in range(4):
        for q,disrupted in layer:
            if disrupted:v.check(not ready(q),'multi-gadget or D return detected')
            word_checks+=1
        if length<3:
            layer={(v.step(q,op),disrupted or op!='K') for q,disrupted in layer for op in ['K','D']+bank}
    receipt=dict(status='PASS',physical_gadgets=len(bank),distinct_small_physical_states=len(states),
                 endpoint_complete_checks=endpoint_checks,endpoint_gap_checks=endpoint_gap_checks,
                 gap_inputs=len(gap_inputs),gap_invariant_checks=invariant_checks,first_entry_checks=entry_checks,
                 mature_ready_inputs=sorted(mature),padding_suffix_max_length=8,padding_checks=suffix_checks,
                 full_word_max_length=3,full_word_checks=word_checks,
                 scope='Finite structural controls for separately proved arbitrary-gadget and all-word no-return theorems; no target coverage census.',
                 inventory_sha256=hashlib.sha256((args.source/'inventory67.json').read_bytes()).hexdigest(),
                 helper_sha256=hashlib.sha256(HELPER.read_bytes()).hexdigest(),
                 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    rendered=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(HERE/'receipt.json').write_text(rendered);print(rendered,end='')
if __name__=='__main__':main()
