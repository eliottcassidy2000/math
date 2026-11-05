"""Independent Session 196 transfer algebra for two-apex tournaments.
Bit 0/1 of a letter says a/b points to the corresponding base vertex.
Base vertices are in transitive order and a -> b. No repository imports.
"""
from itertools import product
from collections import deque, Counter
from functools import lru_cache

LETTERS = tuple(range(4))
STATES2 = tuple(product(range(3), repeat=2))
STATES3 = tuple(product(range(3), repeat=3))

def parameters(kind):
    # Left/right boundary root: None means an unbounded endpoint.
    if kind == 'cycle': return ((0,1),(1,0)), (True,False)
    if kind == 'path_ab': return ((None,0),(0,1),(1,None)), (True,True,True)
    if kind == 'path_ba': return ((None,1),(1,0),(0,None)), (True,False,True)
    raise ValueError(kind)

@lru_cache(None)
def transfers(kind, letter):
    boundaries, allowed_empty = parameters(kind)
    states = tuple(product(range(3), repeat=len(boundaries)))
    index = {s:i for i,s in enumerate(states)}
    ans = [[] for s in states]
    for i,s in enumerate(states):
        for j,(left,right) in enumerate(boundaries):
            if s[j] == 0 and left is not None and not(letter >> left & 1): continue
            good = right is None or not(letter >> right & 1)
            target = list(s); target[j] = 1 if good else 2
            ans[i].append(index[tuple(target)])
    return tuple(tuple(a) for a in ans)

def step(vector, kind, letter, cap=None):
    ans=[0]*len(vector)
    for i,n in enumerate(vector):
        if n:
            for j in transfers(kind,letter)[i]:
                ans[j]+=n
                if cap is not None: ans[j]=min(cap,ans[j])
    return tuple(ans)

def finish(vector,kind,cap=None):
    boundaries,allowed_empty=parameters(kind)
    states=tuple(product(range(3),repeat=len(boundaries)))
    answer=sum(n for s,n in zip(states,vector)
               if all(x==1 or (x==0 and empty) for x,empty in zip(s,allowed_empty)))
    return min(answer,cap) if cap is not None else answer

def count_word(word,kind):
    vector=(1,)+(0,)*(3**len(parameters(kind)[0])-1)
    for letter in word: vector=step(vector,kind,letter)
    return finish(vector,kind)

def counts(word):
    return count_word(word,'path_ab')+count_word(word,'path_ba'), count_word(word,'cycle')

def adjacency(word):
    m=len(word);n=m+2;adj=[0]*n
    for i in range(m):
        for j in range(i+1,m):adj[i]|=1<<j
    for i,letter in enumerate(word):
        for root in range(2):
            if letter>>root&1:adj[m+root]|=1<<i
            else:adj[i]|=1<<(m+root)
    adj[m]|=1<<(m+1)
    return adj

def build_cycle_dfa():
    initial=(1,)+(0,)*8
    states=[initial];index={initial:0};transitions=[]
    for vector in states:
        row=[]
        for letter in LETTERS:
            nxt=step(vector,'cycle',letter,3)
            if nxt not in index:index[nxt]=len(states);states.append(nxt)
            row.append(index[nxt])
        transitions.append(row)
    outputs=[finish(v,'cycle',3) for v in states]
    return states,transitions,outputs

def minimize_dfa(transitions,outputs):
    partition=list(outputs)
    while True:
        signatures=[(outputs[i],tuple(partition[j] for j in transitions[i])) for i in range(len(outputs))]
        unique={};new=[]
        for sig in signatures:
            if sig not in unique:unique[sig]=len(unique)
            new.append(unique[sig])
        if new==partition:break
        partition=new
    reps={}
    for i,p in enumerate(partition):reps.setdefault(p,i)
    minimized_transitions=[[partition[j] for j in transitions[reps[p]]] for p in sorted(reps)]
    minimized_outputs=[outputs[reps[p]] for p in sorted(reps)]
    return partition,minimized_transitions,minimized_outputs

if __name__=='__main__':
    import json
    from pathlib import Path
    from census import dpstats
    tested=0
    for m in range(6):
        for word in product(LETTERS,repeat=m):
            if m == 0: continue # Cycles in a two-vertex tournament do not exist; same output 0.
            got=counts(word);expected=dpstats(adjacency(word))
            assert got==expected,(word,got,expected)
            tested+=1
        print('verified base length',m,flush=True)
    states,transitions,outputs=build_cycle_dfa()
    p,t,o=minimize_dfa(transitions,outputs)
    print('raw DFA',len(states),'minimal Moore DFA',len(t),'outputs',dict(Counter(o)))
    print('independent graph subset-DP cross-checks',tested)
    packet={'alphabet':'letter bit 0: a -> x; bit 1: b -> x; a -> b; base increasing',
            'initial_state':p[0], 'transitions':t,'outputs':o,
            'raw_states':len(states),'states':len(t),'cross_checks':tested}
    Path(__file__).with_name('cycle_dfa.json').write_text(json.dumps(packet,indent=2)+'\n')
