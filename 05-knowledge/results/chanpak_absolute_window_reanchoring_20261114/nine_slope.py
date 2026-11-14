#!/usr/bin/env python3
"""Certified scan of the entire h=9 absorbed left-wing chart at slope D."""
from itertools import combinations
from collections import defaultdict
from two_axis_replay import dp
from complement_gaps import g
BASIS=((14,28,24,16,1),(10,28,25,16,1))
TRI=(2,4,1)

def ordinal(word):
    a=()
    for r in word:
        off=len(a);ext=((1<<len(r))-1)<<off
        a=tuple(x|ext for x in a)+tuple(x<<off for x in r)
    return a

def catalog(n):
    if n>=5:
        for R in BASIS:
            for i in range(n-4):
                yield ordinal([(0,)]*i+[R]+[(0,)]*(n-5-i))
    if n>=6:
        for i,j in combinations(range(n-4),2):
            yield ordinal([TRI if t in (i,j) else (0,) for t in range(n-4)])

def scans(D,sweep=False):
    result=defaultdict(dict); choices=0
    for n in range(5,D.bit_length()):
        for r in catalog(n):
            F=(1<<n)-1
            e=dp(r); s=dp(tuple(F^x^(1<<i) for i,x in enumerate(r)))
            H=[1]+list(map(sum,e[1:]))
            assert H[F]==9
            f=sum(v*H[F^m] for m,v in enumerate(H))
            if f>D:continue
            P=[sum(v[i]*H[F^m] for m,v in enumerate(s)) for i in range(n)]
            Q=[sum(v[i]*H[F^m] for m,v in enumerate(e)) for i in range(n)]
            E=[[sum(e[m][i]*s[F^m][j]+e[m][j]*s[F^m][i] for m in range(F+1)) for j in range(n)] for i in range(n)]
            for U in range(F+1):
                a=[i for i in range(n) if U>>i&1]
                q=f-9-sum(Q[i] for i in a)
                if sweep:
                    B=range(1,(D-q)//f+1)
                else:
                    b,rem=divmod(D-q,f)
                    B=(b,) if not rem and b>=1 else ()
                if not B:continue
                p=sum(P[i] for i in a)
                S=sum(e[F][i] for i in a)
                T=9-sum(s[F][i] for i in a)
                K=9+p-S-sum(E[i][j] for i,j in combinations(a,2))
                for b in B:
                    choices+=1
                    d=p*b+9*(g(b)-2)+K+S*(b%2)
                    key=(f*b+q,T) if sweep else T
                    result[key].setdefault(d,(r,U,b,f,p,q,S,K))
    return result,choices

def runs(values):
    for d in sorted(values):
        if d-2 in values:continue
        t=d
        while t in values:t+=2
        yield list(range(d,t,2))

if __name__=='__main__':
    import sys
    if '--sweep' in sys.argv:
        banks,choices=scans(4095,True)
        wide=[(k,z) for k,values in sorted(banks.items()) for z in runs(values) if len(z)>=4]
        assert choices==86523 and len(banks)==6286
        assert [(k,z[0],len(z)) for k,z in wide]==[
           ((908,0),752,4),((1816,0),1660,4),((1825,9),752,4),
           ((3632,0),3476,4),((3641,9),1660,4),((3650,0),2586,4),((3659,9),752,4)]
        print('PASS: all 86523 marked choices, seven four-windows, no five-window')
        for k,z in wide:print(k,z)
    else:
        D=int(sys.argv[1]) if len(sys.argv)>1 else 908
        banks,choices=scans(D); print('D',D,'marked constructions',choices)
        if D==908:
            assert choices==17 and {T:sorted(x) for T,x in banks.items()}=={0:[624,706,752,754,756,758,772,784]}
        for T,s in sorted(banks.items()):
            best=max(runs(s),key=len,default=[])
            print(T,'offsets',sorted(s),'longest',best)
            if len(best)>=4:
                for d in best: print(d,s[d])
