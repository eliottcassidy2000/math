#!/usr/bin/env python3
"""Exact controls only for the new one-query mixed-radix continuation."""
import itertools,json
from pathlib import Path
from controls import S0,step,decode,close,recount,count

def size(B):
    m=0
    while 3**m<B*2**m:m+=1
    return m

def encode(s,m):
    r,M,h,E=s;q=h+r
    return -r+2*(10*3**m*M+3*2**m*q-E)

def nearest(x,y):
    # floor(x/y+1/2), valid also for negative x; proved no ties here.
    assert 2*(x%y)!=y
    return (2*x+y)//(2*y)

def recover(N,r,m):
    assert (N+r)%2==0
    S=(N+r)//2;A=10*3**m;D=3*2**m
    M=nearest(S,A);q=nearest(S-A*M,D);E=A*M+D*q-S
    return r,M,q-r,E

def main():
    depths=[];total=0
    for k in range(11):
        B=3*29**k;m=size(B);values=set()
        assert m<=4+10*k
        for w in itertools.product((3,5),repeat=k):
            s=S0
            for a in w:s=step(s,a)
            assert 0<=s[2]+s[0]<=B and 0<=s[3]<=B
            N=encode(s,m)
            assert N<=1+26*3**m*B
            assert recover(N,s[0],m)==s
            assert decode(recover(N,s[0],m))==w
            assert N not in values
            values.add(N);total+=1
            S=(N+s[0])//2
            assert S%2==1
        depths.append({'depth':k,'m':m,'distinct_counts':len(values),
                       'max_count':max(values),'order_upper_bound':4*k+m+7})
    examples=[([{1},set()],{0,1},0),([{1},set()],{0,1},3),
              ([set()],{0},2),([{1,2},{2},set()],{0,1},1)]
    physical=[]
    for g,O,m in examples:
        s=recount(g,O);q=s[2]+s[0];B=max(1,q,s[3])
        assert 3**m>=B*2**m
        n=len(g);G=close(g,O,m+2);P={n,n+1}
        T=close(G,P,1)
        N=count(T)
        assert N==encode(s,m)
        assert recover(N,s[0],m)==s
        assert len(T)==len(g)+m+5
        physical.append({'initial_state':s,'m':m,'order':len(T),'count':N})
    s=S0;m=0;S=(encode(s,m)+s[0])//2
    assert S//(10*3**m)==s[1]-1
    assert nearest(S,10*3**m)==s[1]
    receipt={'status':'PASS','states':total,'depths':depths,'physical_controls':physical,
             'hostile_floor_division':'fails as predicted at q=0,E=1; nearest rounding passes',
             'scope':'Only new one-query encoding, no inherited census audit'}
    Path(__file__).with_name('one_query_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','states':total,'physical_controls':len(physical)}))
if __name__=='__main__':main()
