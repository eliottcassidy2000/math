#!/usr/bin/env python3
"""Independent semantic verification of the two 17-state DFA certificates."""
from pathlib import Path
from itertools import combinations, product
from math import isqrt
import json
import sys

GATES=0

def need(ok,msg):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(msg)


def fibonacci(n):
    a,b=0,1
    for _ in range(n):a,b=b,a+b
    return a


def value(msb,kind,shift=0):
    return sum(d*fibonacci(2*i+2+2*shift if kind=='cg' else i+2+shift) for i,d in enumerate(reversed(msb)))


def valid(word,kind):
    if kind=='z':return all(not(a==b==1) for a,b in zip(word,word[1:]))
    # Direct forbidden-factor language; no automaton recurrence.
    for i,d in enumerate(word):
        if d!=2:continue
        j=i+1
        while j<len(word) and word[j]==1:j+=1
        if j<len(word) and word[j]==2:return False
    return True


def locked(word):
    for d in reversed(word):
        if d!=1:return int(d==2)
    return 0


def semantic(word,kind):
    a=[x for x,y in word];b=[y for x,y in word]
    if not valid(a,kind) or not valid(b,kind):return (-1,-1,-1,-1)
    if kind=='cg':l,r=locked(a),locked(b)
    else:l,r=(a[-1] if a else 0),(b[-1] if b else 0)
    return (l,r,(value(a,kind)-value(b,kind))%2,
            (value(a,kind,1)-value(b,kind,1))%2)


def accepted(word,kind):
    s=semantic(word,kind)
    return s!=(-1,-1,-1,-1) and s[2]==0


def check_certificate(cert,kind):
    states={tuple(row['state']):row for row in cert['states']}
    need(len(states)==17,'17 distinct states')
    allstates={(-1,-1,-1,-1)}|set(product(range(2),repeat=4))
    need(set(states)==allstates,'full state universe')
    symbols=cert['alphabet']
    for state,row in states.items():
        reach=row['reach']
        need(semantic(reach,kind)==state,'semantic access word')
        need(accepted(reach,kind)==row['accept'],'semantic acceptance')
        need(len(row['next'])==len(symbols),'transition completeness')
        for symbol,dest in zip(symbols,row['next']):
            need(semantic(reach+[symbol],kind)==tuple(dest),'semantic transition')
    wanted={tuple(sorted(pair)) for pair in combinations(states,2)}
    actual=set()
    for witness in cert['distinguish']:
        left,right=tuple(witness['left']),tuple(witness['right'])
        suffix=witness['suffix']
        need(left!=right and left in states and right in states,'distinction endpoints')
        need(accepted(states[left]['reach']+suffix,kind)!=accepted(states[right]['reach']+suffix,kind),'literal distinguishing words')
        actual.add(tuple(sorted((left,right))))
    need(actual==wanted,'all 136 state pairs distinguished')
    tested=0
    max_length=4 if kind=='cg' else 7
    for size in range(max_length+1):
        for word in product(symbols,repeat=size):
            state=tuple(cert['start'])
            for symbol in word:
                row=states[state]
                state=tuple(row['next'][symbols.index(symbol)])
            need((state!=(-1,-1,-1,-1) and state[2]==0)==accepted(list(word),kind),'full language semantics')
            tested+=1
    return {'states':len(states),'distinguished_pairs':len(actual),'explicit_word_comparisons':tested}


def main():
    if '--test-failing-gate' in sys.argv:need(False,'intentional failure control')
    source=Path(__file__).resolve().parent.parent/'zeckendorf'/'normal.json'
    receipt=json.loads(source.read_text())
    results={kind:check_certificate(receipt['minimality_certificates'][kind],kind) for kind in ('cg','z')}
    cgnumbers={}
    for word in product(range(3),repeat=8):
        if not valid(word,'cg'):continue
        n=value(word,'cg');need(n not in cgnumbers,'CG collision');cgnumbers[n]=word
        expected=0 if n==0 else (3*n+isqrt(5*n*n))//2+1
        need(value(word,'cg',1)==expected,'canonical shift vs exact Beatty')
    need(set(cgnumbers)==set(range(fibonacci(18))),'eight-digit CG exact interval')
    interwoven=0
    znumbers={}
    for word in product(range(2),repeat=16):
        if not valid(word,'z'):continue
        znumbers[value(word,'z')]=word
        z=list(reversed(word));e=z[::2]+[0];o=z[1::2]+[0]
        c=[e[i]+(o[i-1] if i else 0)-o[i] for i in range(9)]
        need(all(d in (-1,0,1) for d in c),'signed even digit bound')
        need(sum(d*fibonacci(2*i+2) for i,d in enumerate(c))==value(word,'z'),'interwoven value')
        recovered=[];old=0
        for d in c:
            candidates=[(x,y) for x,y in product(range(2),repeat=2) if not x*y and not old*x and x+old-y==d]
            need(len(candidates)==1,'interwoven inverse uniqueness')
            x,y=candidates[0];recovered.extend([x,y]);old=y
        need(old==0 and recovered==z+[0,0],'interwoven exact inverse')
        interwoven+=1
    return_checks={}
    for kind,numerals in [('cg',cgnumbers),('z',znumbers)]:
        def shift(n):
            return (1 if n>=0 else -1)*value(numerals[abs(n)],kind,1)
        for n in range(-100,101):
            a=shift(n);b=shift(a);c=shift(b)
            need(b==(3*a-n if kind=='cg' else a+n),'two-shift recurrence')
            need(c==(8*a-3*n if kind=='cg' else 2*a+n),'three-shift recurrence')
            need((c-n)%2==0,'three-shift preserves value parity')
        pairs=0
        for s in range(-40,41):
            for d in range(-40,41):
                if (s-d)%2:continue
                initial=(shift(s)-shift(d))%2
                ss,dd=s,d
                for k in range(4):
                    need(((ss-dd)%2==0)==(not initial or k%3==0),'exact gluing return classification')
                    if k<3:ss,dd=shift(ss),shift(dd)
                pairs+=1
        ss,dd=3,1;path=[]
        for k in range(4):
            path.append([ss,dd])
            if k<3:ss,dd=shift(ss),shift(dd)
        need([(a-b)%2 for a,b in path]==[0,1,1,0],'minimal three-shift block witness')
        return_checks[kind]={'signed_inputs':201,'equal_parity_pairs':pairs,'minimality_witness':path}
    print(json.dumps({'status':'PASS' ,'gates':GATES,'dfa':results,'CG_eight_digit_interval':len(cgnumbers),
        'interwoven_16_digit_words':interwoven,'shift_returns':return_checks},sort_keys=True,indent=2))

if __name__=='__main__':main()
