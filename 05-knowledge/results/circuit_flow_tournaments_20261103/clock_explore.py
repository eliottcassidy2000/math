"""Bounded discovery only; do not promote incomplete saturations."""
from collections import deque
import json


def step(f,letter,M):
    n=len(f);g=[]
    for r in range(n):
        if letter=='L':
            z=-r*(r-1)*(2*r-1)*f[(r-1)%n]+r*r*(r*r-2)*f[r]+r*(r+1)*(2*r+1)*f[(r+1)%n]
        else:
            z=-2*r*(r*r-r+1)*f[(r-1)%n]+r**4*f[r]+2*r*(r*r+r+1)*f[(r+1)%n]
        g.append(z%M)
    return tuple(g)


def census(M,limit=20000):
    seed=tuple(range(M));states=[seed];index={seed:0};queue=deque([0]);edges=[]
    while queue:
        i=queue.popleft();row=[]
        for letter in 'LU':
            q=step(states[i],letter,M)
            if q not in index:
                if len(states)>=limit:return {'M':M,'complete':False,'states':len(states)}
                index[q]=len(states);states.append(q);queue.append(len(states)-1)
            row.append(index[q])
        edges.append(row)
    colors=[s[1] for s in states]
    while True:
        keys=[(colors[i],*(colors[j] for j in edges[i])) for i in range(len(states))]
        bank={};nxt=[bank.setdefault(k,len(bank)) for k in keys]
        if len(set(nxt))==len(set(colors)):break
        colors=nxt
    commuting=all(nxt[edges[edges[i][0]][1]]==nxt[edges[edges[i][1]][0]] for i in range(len(states)))
    return {'M':M,'complete':True,'states':len(states),'observable_classes':len(set(nxt)),
            'observable_commuting':commuting,'H_residues':sorted(set(s[1] for s in states))}

if __name__=='__main__':
    for M in (8,9,12,24,72,5,7,11,13):print(json.dumps(census(M)),flush=True)
