"""Proof-derived positive Hamiltonian-path transfer for two marked apices."""
from itertools import product, permutations
from collections import defaultdict
import json
from pathlib import Path

# Bit i of a letter says apex i -> this spine vertex. Apex 0 -> apex 1.
ORDERS = ((0,1),(1,0))
STATES = list(product(range(3),range(3),range(2)))
INDEX = {s:i for i,s in enumerate(STATES)}
ZERO = INDEX[(0,0,0)]

def moves(order, state, letter):
    for run in range(3):
        if state[run] == 0 and run and not (letter >> order[run-1]) & 1:
            continue  # first vertex must follow its left apex
        new = list(state)
        new[run] = 1 if run == 2 or not ((letter >> order[run]) & 1) else 2
        yield INDEX[tuple(new)]

def accepting(order, state):
    return state[0] != 2 and state[1] != 2 and (state[1] != 0 or order == (0,1))

def initial():
    return tuple(int(i == ZERO) for order in ORDERS for i in range(18))

def step(v, letter, cap=None):
    out = [0]*36
    for block,order in enumerate(ORDERS):
        for i,state in enumerate(STATES):
            for j in moves(order,state,letter):
                out[18*block+j] += v[18*block+i]
    if cap is not None: out = [min(cap,a) for a in out]
    return tuple(out)

def value(v):
    return sum(v[18*b+i] for b,o in enumerate(ORDERS) for i,s in enumerate(STATES) if accepting(o,s))

def h(word):
    v=initial()
    for letter in word: v=step(v,letter)
    return value(v)

def brute(word):
    n=len(word)+2
    def edge(a,b):
        if a<2 and b<2: return a<b
        if a<2: return bool(word[b-2]>>a & 1)
        if b<2: return not bool(word[a-2]>>b & 1)
        return a<b
    dp={(1<<i,i):1 for i in range(n)}
    for mask in range(1,1<<n):
        for last in range(n):
            c=dp.get((mask,last),0)
            if not c: continue
            for nxt in range(n):
                if not mask>>nxt&1 and edge(last,nxt):
                    key=(mask|1<<nxt,nxt)
                    dp[key]=dp.get(key,0)+c
    return sum(dp.get(((1<<n)-1,i),0) for i in range(n))

def main():
    witnesses={}; layers=[]
    for n in range(6):
        vals=[]
        for word in product(range(4),repeat=n):
            x=h(word)
            if n<=4: assert x==brute(word)
            # Saturation commutes coordinatewise with every transition.
            for cap in (2,8,14):
                v=initial(); w=initial()
                for a in word:
                    v=step(v,a);w=step(w,a,cap)
                    assert w==tuple(min(cap,t) for t in v)
            vals.append(x)
            if x not in witnesses: witnesses[x]=list(word)
        assert 8*sum(vals)==int(n==0)+3*4**n+3*8**n+12**n
        layers.append({'spine_length':n,'words':4**n,'sum_H':sum(vals),
                       'sum_H_squared':sum(x*x for x in vals),
                       'distinct_H':sorted(set(vals)),
                       'cumulative_missing_odd_through_99':[x for x in range(1,100,2) if x not in witnesses]})
    assert 13 in witnesses
    # 16 length-two words, or 16 four-vertex marked tournaments, are NOT
    # automatically the original session's five-cell decomposition.
    out={'status':'FINITE-EXACT','alphabet':'bit i = apex i points to spine vertex',
         'state_dimension':36,'brute_checked_words':sum(4**n for n in range(5)),
         'layers':layers,'first_13_word':witnesses[13],
         'selected_witnesses':{str(x):witnesses[x] for x in range(1,100,2) if x in witnesses}}
    Path(__file__).with_name('receipts.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='layers'},indent=2))
    print('moment rows:', [(r['spine_length'],r['sum_H'],r['sum_H_squared']) for r in layers])

if __name__=='__main__': main()
