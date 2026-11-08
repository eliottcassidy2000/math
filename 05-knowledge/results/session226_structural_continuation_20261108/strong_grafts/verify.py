#!/usr/bin/env python3
"""Portable literal controls for exact graft strongness; no input files."""
import itertools
import json

def require(p, message):
    if not p: raise RuntimeError(message)

def tournaments(n):
    pairs=list(itertools.combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        out=[0]*n
        for k,(u,v) in enumerate(pairs):
            if bits>>k&1: out[u]|=1<<v
            else: out[v]|=1<<u
        yield out

def reaches(rows, start):
    seen=todo=1<<start
    while todo:
        bit=todo&-todo; todo-=bit
        new=rows[bit.bit_length()-1]&~seen
        seen|=new; todo|=new
    return seen==(1<<len(rows))-1

def strong(rows):
    return all(reaches(rows,v) for v in range(len(rows)))

def graft(a, root, b):
    n=len(a); B=((1<<len(b))-1)<<n
    return [row|(B if u!=root else 0) for u,row in enumerate(a)]+[(row<<n)|(1<<root) for row in b]

def paths(rows):
    n=len(rows); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for end in range(n):
            if dp[mask][end]:
                todo=rows[end]&~mask
                while todo:
                    bit=todo&-todo; todo-=bit
                    dp[mask|bit][bit.bit_length()-1]+=dp[mask][end]
    return sum(dp[-1])

def main():
    bank=[r for n in range(1,5) for r in tournaments(n)]
    counters={'pairs':0,'strong_grafts':0,'strong_with_nonstrong_gadget':0,'singleton_sources':0,'nonreaching_sources':0}
    for a in bank:
        for b in bank:
            actual=strong(graft(a,0,b)); expected=len(a)>=2 and reaches(a,0)
            require(actual==expected,'strongness iff failed')
            counters['pairs']+=1
            counters['strong_grafts']+=actual
            counters['strong_with_nonstrong_gadget']+=actual and not strong(b)
            counters['singleton_sources']+=len(a)==1
            counters['nonreaching_sources']+=not reaches(a,0)
    require(len(bank)==75 and counters['pairs']==5625,'universe')
    examples={}
    for name,a,expected in [('A',[14,28,24,16,1],[423,605,357,185,261]),('B',[10,28,25,16,1],[315,541,551,217,207])]:
        require(paths(a)==9 and strong(a),'source witness')
        out=[]
        for r in range(5):
            g=graft(a,r,[8,1,3,6]); require(strong(g),'strong H9 extension')
            out.append(paths(g))
        require(out==expected,'H9 graph count')
        examples[name]=out
    require(set(examples['A']).isdisjoint(examples['B']),'separator')
    # Positive control: a nonstrong directed-edge gadget gives a strong graft.
    require(not strong([2,0]) and strong(graft([2,0],0,[2,0])),'nonstrong gadget control')
    # Hostile: no gadget repairs an old sink root.
    require(not strong(graft([2,0],1,[2,4,1])),'sink-root hostile')
    print(json.dumps({'status':'PASS','universe':counters,'H9_graft_counts':examples},sort_keys=True,indent=2))

if __name__=='__main__': main()
