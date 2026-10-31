#!/usr/bin/env python3
"""All paired/terminal A,C interfaces, symbolic shifts and concrete signed choices."""
import json,sys
from pathlib import Path
from phase_trees import Trees,arity,code,digest,phase
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031c'))
from independent_ray_channels import hostile

def need(x,msg):
    if not x:raise RuntimeError(msg)

def signature(atoms,m,q):
    at=0;heads={};tails={};fixed={}
    for a in atoms:
        k,c=a[:2]
        if k=='single':fixed[c]=at%2;at+=1
        elif k=='head':heads[c]=at%2;at+=1
        else:
            _,c,i,o=a;length=m+i;start=int(k=='tail')
            shift=(at-start if o==1 else at+length-1)%2
            need(c not in tails,'unique signed tail')
            tails[c]=shift
            if not start:heads[c]=shift
            at+=length-start
    half=(q-1)//2
    need(set(heads)==set(tails)==set(range(-half,half+1)),'all signed channels')
    need(set(fixed)==set(range(1,half+1)),'low singleton coverage')
    return heads,tails,fixed,at

def main():
    t=Trees();records=[];channels=0
    for symbol,q in [('A',25),('C',49)]:
        half=(q-1)//2
        for parity in (0,1):
            m=112+parity
            for d in range(half,half+q):
                n=q*m+d;r=t.c.certificate(symbol,parity,d)
                sig=[signature(w,m,q) for w in r['words']]
                need(sig[0][:3]==sig[1][:3],'per-label output phase shifts')
                need([s[3] for s in sig]==[n,n+1],'integer member clock')
                bits=code(symbol,f'{parity}:{d}')
                pair=t.pair(n,symbol,bits)
                for i,p in enumerate(pair):
                    need(t.decode(p,n,symbol,i)==bits,'every signed member decoder')
                uniform=t.pair(n,symbol,'0'*q)
                need([phase(p) for p in pair]==[phase(p) for p in uniform],'common full family phase')
                terminal=t.c.certificate(symbol,parity,d,True)
                w=terminal['word'] if symbol=='C' else terminal['words'][0]
                need(signature(w,m,q)[3]==n,'terminal integer clock')
                cyc=t.cycle(n,symbol,bits)
                records.append([symbol,parity,d,digest(cyc)])
                channels+=3*q
    import hashlib
    print(json.dumps({'status':'passed','pair_states':len(records),'terminal_states':len(records),
        'signed_readable_channels':channels,'cycles_receipt_sha256':hashlib.sha256(json.dumps(records).encode()).hexdigest(),
        'hostile_different_phase':hostile(),
        'scope':'Complete A/C interface states; two members plus terminal per state; phase-fixed choices independently assigned to every signed offset.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
