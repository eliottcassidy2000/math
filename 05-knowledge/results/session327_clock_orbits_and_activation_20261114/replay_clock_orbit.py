"""Literal equal-clock switch certificates and modular lattice controls."""
from heapq import heappush, heappop
from math import gcd
import sys
from pathlib import Path
from itertools import combinations
import importlib.util
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root))
from profile_tools import profile, strong
def need(ok):
    if not ok:
        raise RuntimeError("clock certificate failed")

P=(
 (73,85,1088,(2,12,49,69,11,27,55),(16,5,65,7,46,15,59),0),
 (5,19,256,(2,4,9,3),(2,20,9,3,45,15),3),
 (13,41,416,(2,12,9,17,7),(6,12,24,49,35,7),2),
 (23,79,1184,(8,5,9,50,7,23),(6,12,24,97,11,23,55),2),
)
def orbit(M):
    H={1}; orders=[]
    for a,b,q,*_ in P:
        x=b*pow(a,-1,M)%M
        y=x; d=1
        while y not in H: y=y*x%M; d+=1
        H={z*pow(x,i,M)%M for i in range(d) for z in H}
        orders.append(d)
    return H,orders

def least_costs(M):
    R=[(b*pow(a,-1,M)%M,q) for a,b,q,*_ in P]
    heap=[(1,1)]; costs={}
    while heap:
        q,x=heappop(heap)
        if x in costs:continue
        costs[x]=q
        for r,e in R:
            for y in (r,pow(r,-1,M)):
                z=x*y%M
                if z not in costs: heappush(heap,(q*e,z))
    return costs

def main():
    for a,b,q,A,B,t in P:
        need(strong(A) and strong(B) and profile(A)==(a,q>>t) and profile(B)==(b,q))
    for M,orders in ((8712,[330,2,2,2]),(19008,[720,2,2,2])):
        H,D=orbit(M); need(D==orders and H=={i for i in range(M) if gcd(i,M)==1})
    Q=least_costs(396)
    need(Q[199]==303104 and Q[395]==278528)
    need((1501-115)==1386 and (1615+365)==1980)
    for h in range(1,396,2):
        if gcd(h,396)>1: continue
        S={h*c%396 for c in range(1,198,2)}
        Z={i for i in range(1,396,2)}
        need({(h+198)*c%396 for c in range(1,198,2)}==Z-S)
        need({(-h)*c%396 for c in range(1,198,2)}==Z-S)
        u=pow(h,-1,198);s=min(u,198-u)
        nmax=max(next(i for i in range(199) if (x+2*i)%396 not in S) for x in S)
        need(nmax==1+98//s)
    from json import loads
    from collections import Counter
    D={int(c):r for c,r in loads((root/'phase_table.json').read_text()).items()}
    need([Counter(D.values())[r] for r in range(5)]==[33,34,18,11,3])
    for H in (5,19,23,79,91,151):
        E=set()
        for c in range(1,198,2):
            step=396 if c<=99 else 792
            E.update((H*c+j)%792 for j in range(0,792,step))
        need(len(E)==149)
        E=set();J=set()
        for c,r in D.items():
            step=99*(1<<(r+2))
            E.update((H*c+j)%6336 for j in range(0,6336,step))
            J.update((H*c+j*step)%19008 for j in range(19008//step) if j%3)
        need(len(E)==897 and len(J)==1794)
    print('four literal equal-clock switches, full unit orbits, both minimum phase costs, all rotation caps')
if __name__=='__main__': main()
