from math import isqrt
from replay_mersenne import A, join, need, strong, table

def inflate(block):
    labels = [(0,k) for k in range(len(block))]+[(i,0) for i in range(1,len(A))]
    def arc(i,j):
        x,u = labels[i]
        y,v = labels[j]
        return (block[u] >> v & 1) if x==y==0 else (A[x] >> y & 1)
    return [sum(arc(i,j)<<j for j in range(len(labels))) for i in range(len(labels))], len(block)+4

for block, expect in [
    ((22,28,24,17,0),(5,64,4091)),
    ((58,116,65,70,12,28,49),(133,1628,104059)),
]:
    _, f = table(block)
    h=f[-1]
    p=sum(x*f[(len(f)-1)^s] for s,x in enumerate(f))
    q=64*p-h
    need((h,p,q)==expect and q>1 and all(q%d for d in range(2,isqrt(q)+1)))
    a,v=inflate(block)
    d,F=table(a)
    b=sum(row[v]*F[(len(F)-1)^s] for s,row in enumerate(d) if s>>v&1)
    R=p-2*h
    need((F[-1],b)==(129*h+64*R,190*h+96*R))
    child=join(a,v)
    need(strong(child) and table(child)[1][-1]==7*q)
    print(h,p,q,7*q)
