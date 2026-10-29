#!/usr/bin/env python3
from itertools import product
import json

def need(p,m):
    if not p:raise RuntimeError(m)
def step(b,n):
    t=b*n+1;a=(t&-t).bit_length()-1
    return t>>a,a

def state(b,w):
    A=B=1;C=0;cap=None
    for a in w:
        A,B,C=b*A,B*2**a,b*C+B
        need(C>0,'positive carry')
        if B>A:cap=C//(B-A) if cap is None else min(cap,C//(B-A))
    M=2*B;r=((B-C)*pow(A,-1,M))%M
    return r,M,cap

def main():
    micro=[7];a=b=1;c=0
    for _ in range(8):
        n=micro[-1]
        if n%2:a,b,c=3*a,2*b,3*c+b;micro.append((3*n+1)//2)
        else:b*=2;micro.append(n//2)
    need(micro==[7,11,17,26,13,20,10,5,8],'micro-clock hostile')
    need((a,b,c)==(243,256,347) and min(micro)<7 and micro[-1]>7,'carry vs slope hostile')
    tested=0;hostile=[];replays=0
    for k in range(1,9):
        for w in product(range(1,5),repeat=k):
            r,M,cap=state(3,w);tested+=1
            n=r if r>1 else r+M
            if cap is not None and n<=cap:hostile.append({'word':w,'source':n,'cap':cap})
            if k<=5:
                for q in (0,1,3):
                    n=r+M*q;x=n;survives=True
                    for a in w:
                        x,actual=step(3,x);need(actual==a,'native cylinder')
                        survives=survives and x>=n
                    need(survives==(cap is None or n<=cap),'prefix cutoff equivalence')
                    replays+=1
    r,M,cap=state(3,(2,));need(r==cap==1,'terminal-only finite cell')
    for w,n in [((1,1,5),13),((1,3,3),17)]:
        r,M,cap=state(5,w);need(r==cap==n,'five finite singleton')
    crossings=0;uncrossed=0;scan_hostile=[]
    for n in range(3,100002,2):
        A=B=1;x=n;w=[]
        for k in range(50):
            x,a=step(3,x);A*=3;B*=2**a;w.append(a)
            if B>A:
                crossings+=1
                if x>=n:scan_hostile.append({'source':n,'endpoint':x,'word':w})
                break
        else:uncrossed+=1
    print(json.dumps({'status':'PASS','word_search':{'multiplier':3,'max_length':8,'alphabet':[1,2,3,4],
          'words':tested,'nonterminal_finite_survivors':hostile},'native_cutoff_replays':replays,'micro_clock_hostile':micro,
          'source_search':{'odd_sources':50000,'maximum_source':100001,'horizon':50,
          'crossings':crossings,'not_crossed_within_horizon':uncrossed,'unpaid_crossings':scan_hostile}},sort_keys=True,indent=2))
if __name__=='__main__':main()
