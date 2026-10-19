#!/usr/bin/env python3
"""Independent finite controls: exact cycle return, phase criterion, and GF indexing."""
from collections import Counter
from itertools import product
import json
checks=0

def check(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def cycle_output(cyc,b,Q):
    # Two whole proposed periods suffice here because input period divides Q.
    # Matching blocks imply identical running state at Q and 2Q, certifying
    # every future block by determinism, rather than guessing from a preview.
    check(Q%len(cyc)==0,('input period divisibility',len(cyc),Q))
    q=b;out=[]
    for j in range(2*Q):q=abs(q-cyc[j%len(cyc)]);out.append(q)
    check(out[:Q]==out[Q:],('exact periodic return',cyc,b,Q,out))
    out=tuple(out[:Q]);least=next(d for d in range(1,Q+1)if Q%d==0 and all(out[i]==out[i%d]for i in range(Q)))
    return out[:least]

def first_defect(cyc):return next((j for j,x in enumerate(cyc)if x>=2),None)
def phase_set(cyc):
    h=0;ph=set()
    for x in cyc+cyc:
        if x==2:ph.add(h)
        h^=x&1
    return ph

def gf_parity(word,Q):
    L=len(word);num=1;power=3
    for b in word:
        if b:num^=power
        power^=power<<1
    return tuple(sum((num>>k&1) for k in range(j+1)if ((j-k)&L)==0)%2 for j in range(Q))

def complexity(cyc):
    a=tuple(x&1 for x in cyc);k=0
    while any(a):a=tuple(x^y for x,y in zip(a,a[1:]+a[:1]));k+=1;check(k<=len(a),'dyadic nilpotence')
    return k

rows=[];states=[((),(3,))];hist=[]
for L in range(15):
    Q=1<<L.bit_length();H=Counter();annihilators=Counter()
    for word,cyc in states:
        check(len(cyc)==Q,('least integer period',word))
        check(tuple(x&1 for x in cyc)==gf_parity(word,Q),('GF coefficient indexing',word))
        check(complexity(cyc)==L+1,('parity complexity',word))
        tau=first_defect(cyc);H[tau]+=1
        if tau is not None:check(tau%2==0 and tau<=2*((Q-1)//2),('conditional horizon finite control',word,tau))
        if L%2==0 and any(word):
            check(max(cyc)<=2,('post-one alphabet',word))
            ph=phase_set(cyc);check(len(ph)<=1,('even reachable phase',word))
            children=[cycle_output(cyc,b,1<<(L+1).bit_length())for b in(0,1)]
            childtau=[first_defect(c)for c in children]
            if tau is None:check(childtau==[None,None],('binary parents',word))
            else:check(sorted(childtau,key=lambda x:-1 if x is None else x)==[None,tau],('including age-zero parents',word,childtau))
            annihilators['age_zero'if tau==0 else'positive'if tau is not None else'binary']+=1
    rows.append({'L':L,'words':len(states),'exact_period':Q,'histogram':{str(k):v for k,v in H.items()},'even_annihilator_parent_types':dict(annihilators)})
    hist.append(H)
    if L<14:
        Qnext=1<<(L+1).bit_length();states=[(word+(b,),cycle_output(cyc,b,Qnext))for word,cyc in states for b in(0,1)]
for L in range(0,14,2):
    check({r:c for r,c in hist[L].items()if r not in(0,None)}=={r:c for r,c in hist[L+1].items()if r not in(0,None)},('paired positive histograms',L))
# A direct finite phase check, independent of eventual-cycle machinery.
finite_controls=0
for n in range(9):
    for v in product(range(3),repeat=n):
        h=0;ph=set()
        for x in v:
            if x==2:ph.add(h)
            h^=x&1
        good=[]
        for b in(0,1):
            q=b;out=[]
            for x in v:q=abs(q-x);out.append(q)
            if all(x<=1 for x in out):good.append(b)
        expected=[0,1]if not ph else[1-next(iter(ph))]if len(ph)==1 else[]
        check(good==expected,('finite phase exact bit',v,good,expected));finite_controls+=1
# Explicit tails for the non-invariance hostile, not zero tails at every stage.
v=(2,1,0,1,2)+(0,)*20
q=0;first=[]
for x in v:q=abs(q-x);first.append(q)
q=0;second=[]
for x in first:q=abs(q-x);second.append(q)
check(first[:5]==[2,1,1,0,2]and first[5:]==[2]*20,'first hostile tail')
check(second[:5]==[2,1,0,0,2]and second[5:]==[0,2]*10,'second hostile tail')
check(phase_set((2,1,0,1,2,0))=={0},'initial hostile phase on zero tail')
check(phase_set((2,1,0,0,2,0))=={0,1},'second hostile phase')
print(json.dumps({'status':'PASS finite controls; no all-length hypothesis proved','checks':checks,'finite_phase_words':finite_controls,'reachable_rows':rows,'noninvariance_tails':{'input':'21012 then 0 forever','first_T0':'21102 then 2 forever','second_T0':'21002 then (02) forever'}},indent=2,sort_keys=True))
