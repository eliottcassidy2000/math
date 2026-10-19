#!/usr/bin/env python3
"""Bounded exact tests of an OPEN reachable-state refinement.
Eventually periodic sequences are represented exactly, not by finite previews.
No assert statements, floating arithmetic, or third-party dependencies.
"""
from itertools import product
import json

def require(p,msg):
    if not p:raise RuntimeError(msg)

def canonical(pre,cyc):
    pre=list(pre);cyc=list(cyc)
    for d in range(1,len(cyc)+1):
        if len(cyc)%d==0 and all(cyc[i]==cyc[i%d]for i in range(len(cyc))):cyc=cyc[:d];break
    while pre and pre[-1]==cyc[-1]:
        pre.pop();cyc=cyc[-1:]+cyc[:-1]
    return tuple(pre),tuple(cyc)

def T(seq,b):
    pre,cyc=seq;out=[];q=b
    for x in pre:q=abs(q-x);out.append(q)
    seen={};tail=[];i=0
    while (i,q)not in seen:
        seen[i,q]=len(tail);q=abs(q-cyc[i]);tail.append(q);i=(i+1)%len(cyc)
    cut=seen[i,q]
    return canonical(out+tail[:cut],tail[cut:])

def tau(seq):
    pre,cyc=seq
    return next((i for i,x in enumerate(pre+cyc)if x>=2),None)

def phases(seq):
    pre,cyc=seq;h=0;result=set()
    # Two cycles cover both possible parity phases at every periodic position.
    for x in pre+cyc+cyc:
        if x==2:result.add(h)
        h^=x%2
    return result

def C_periodic(cyc):
    require(len(cyc)%2==0,'pair period must be even')
    S={(0,0),(0,1),(1,0),(1,1),(2,0),(2,1),(2,2)}
    q=0;pairs=list(zip(cyc[::2],cyc[1::2]))
    for x in pairs*3:
        if x not in S:return False
        if x==(0,0):continue
        if x==(2,0):
            if q==2:return False
            q=1
        elif x==(2,2):
            if q==1:return False
            q=2
        else:q=0
    return True

def parity_complexity(cyc):
    v=tuple(x%2 for x in cyc);d=0
    while any(v):
        v=tuple(v[i]^v[(i+1)%len(v)]for i in range(len(v)));d+=1
        require(d<=len(v),'expected a power-of-two period')
    return d

def main():
    # Exact general phase criterion on eventually periodic finite controls.
    controls=0
    for n in range(1,9):
        for v in product(range(3),repeat=n):
            seq=((),v)
            killing=[b for b in(0,1)if tau(T(seq,b))is None]
            require(bool(killing)==(len(phases(seq))<=1),'phase equivalence')
            if 2 in v:require(len(killing)==(1 if len(phases(seq))==1 else 0),'unique killing bit')
            controls+=1
    # Exact infinite rotated sequences, through front length 14.
    states=[(((),(3,)),False)];rows=[]
    for L in range(15):
        phase_checks=0;positive=0;maxperiod=0;maxtransient=0
        for seq,hasone in states:
            pre,cyc=seq;maxperiod=max(maxperiod,len(cyc));maxtransient=max(maxtransient,len(pre))
            require(not pre and len(cyc)==1<<L.bit_length(),f'pure dyadic age period failure at L={L}')
            t=tau(seq);positive+=int(t is not None and t>0)
            if L%2==0 and hasone:
                require(max(pre+cyc)<=2,'alphabet after first one')
                require(len(phases(seq))<=1,f'even front phase failure at L={L}')
                if t is not None:
                    child=[tau(T(seq,b))for b in(0,1)]
                    require(child.count(t)==1 and child.count(None)==1,'branch annihilation')
                phase_checks+=1
        rows.append({'front_length':L,'words':len(states),'even_phase_checks':phase_checks,'positive_finite_words':positive,'maximum_exact_period':maxperiod,'maximum_exact_transient':maxtransient})
        if L<14:states=[(T(s,b),h or bool(b))for s,h in states for b in(0,1)]
    # Minimal six-entry finite hostile, with an all-zero continuation.
    finite_hostile=((0,0,2,1,2,0),(0,))
    require(phases(finite_hostile)=={0,1},'finite hostile phases')
    require([tau(T(finite_hostile,b))for b in(0,1)]==[2,4],'finite hostile branches')
    # A stronger period-eight hostile has correct parity complexity L+1=7.
    cyc=(0,0,1,0,2,1,2,0);seq=((),cyc)
    require(C_periodic(cyc)and parity_complexity(cyc)==7,'periodic hostile sidecars')
    require([tau(T(seq,b))for b in(0,1)]==[6,4],'periodic hostile branches')
    alpha=(1,0,1,1,0,0);actual=((),(3,))
    for b in alpha:actual=T(actual,b)
    require(actual==((),(0,0,1,0,0,1,0,0)),'same-parity reachable control')
    require(tuple(x%2 for x in cyc)==actual[1],'same parity control')
    p=actual[1];spike=(p+(2,),p[1:]+p[:1])
    require(tau(spike)==8,'late spike hostile')
    sample=spike[0]+spike[1]*4
    require(all(x%2==p[i%8]for i,x in enumerate(sample)),'late spike parity')
    pairs=list(zip(sample[::2],sample[1::2]))
    require(all(x in {(0,0),(0,1),(1,0),(1,1),(2,0)}for x in pairs),'late spike pair invariant')
    # Minimality only within the explicitly stated dyadic periods 2,4,8.
    minimal=[]
    for P in(2,4,8):
        witness=None
        for v in product(range(3),repeat=P):
            if max(v[:2])<2 and 2 in v and C_periodic(v)and parity_complexity(v)%2==1 and len(phases(((),v)))>1:
                witness=v;break
        minimal.append({'period':P,'first_hostile':witness})
    print(json.dumps({'status':'FINITE-EXACT evidence for an OPEN all-length conjecture','periodic_phase_controls':controls,'reachable_front_checks':rows,'legal_finite_hostile':'002120 followed by zeros','legal_periodic_hostile':'(00102120)^infinity','periodic_hostile_parity_complexity':7,'same_parity_actual_front_word':'101100','same_parity_actual_age_cycle':'00100100','same_parity_late_spike_first_defect':8,'restricted_minimality':minimal},indent=2,sort_keys=True))
if __name__=='__main__':main()
