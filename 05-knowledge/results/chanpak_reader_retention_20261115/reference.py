#!/usr/bin/env python3
from itertools import combinations

def profiles(out, O):
    n=len(out); Z=(1<<n)-1
    end=[[0]*n for _ in range(Z+1)]
    start=[[0]*n for _ in range(Z+1)]
    incoming=[sum(1<<i for i in range(n) if out[i]>>j&1) for j in range(n)]
    for i in range(n):
        end[1<<i][i]=start[1<<i][i]=1
    for s in range(1,Z+1):
        for i in range(n):
            for d,adj in ((end,out),(start,incoming)):
                a=d[s][i]
                if not a:
                    continue
                avail=adj[i]&~s
                while avail:
                    b=avail&-avail
                    avail-=b
                    d[s|b][b.bit_length()-1]+=a
    h=[1]+list(map(sum,end[1:]))
    marks=[i for i in range(n) if O>>i&1]
    E=[sum(x[i] for i in marks) for x in end]
    S=[sum(x[i] for i in marks) for x in start]
    U=sum(S[s]*h[Z^s] for s in range(1,Z+1))
    W=sum(E[s]*S[Z^s] for s in range(1,Z))
    q=sum(h[s]*h[Z^s] for s in range(Z+1))
    return (h[Z]-W,U,h[Z],E[Z]),q

def tt(n):
    return [(1<<n)-(1<<(i+1)) for i in range(n)]

def join(A,B):
    return [x | (((1<<len(B))-1)<<len(A)) for x in A]+[x<<len(A) for x in B]

def close(A,O,X):
    B=join(A,X)
    apex=1<<len(B)
    return [x|(apex if not O>>i&1 else 0) for i,x in enumerate(B)]+[O]

if __name__=='__main__':
    A,O=[18,0,3,7,14],17
    B,P=[10,0,3,6,15],16
    assert profiles(A,O)[0]==(1,24,5,0)
    assert profiles(B,P)[0]==(3,24,3,0)
    assert profiles([20,1,2,7,14],17)[0]==(3,32,9,2)
    c4=[6,12,8,1]
    assert profiles(c4,0)==((5,0,5,0),32)
    assert {profiles(c4,o)[0][1] for o in range(1,16) if profiles(c4,o)[0][0]==3}=={13,14}
    for T in ([0],[2,4,1],c4):
        (_,_,c,_),q=profiles(T,0)
        X=join(T,tt(3))
        A1,B1=close(A,O,X),close(B,P,X)
        o=(1<<(len(T)+5))|(1<<(len(T)+6))
        v=(36*c+720*q,192*q,72*q)
        for G,r0 in ((A1,1),(B1,3)):
            want=(-c*r0,v[0],v[1]+c*r0,v[2])
            assert profiles(G,o)[0]==want
            close_count=profiles(close(G,o,tt(1)),0)[0][2]
            assert close_count==want[0]+2*want[1]
        print('literal synchronized pair and two closings passed',c,q,v[0])
    D=join([20,1,2,7,14],c4)
    assert profiles(D,17)[0]==(15,1024,45,0)
    print('paid cone and dyadic suffix control passed')

