#!/usr/bin/env python3
"""Independent subset-DP controls for the unrestricted first ear."""
from offset_check import H, endings, profile, strong, labelled, left, require, add

def restrict(a,verts):
    return tuple(sum(1<<j for j,v in enumerate(verts) if a[u]>>v&1) for u in verts)
def K(a,q):
    n=len(a);full=(1<<n)-1;total=0
    for mask in range(1,1<<n):
        if not(mask>>q&1):continue
        verts=[i for i in range(n) if mask>>i&1]
        rest=[i for i in range(n) if not(mask>>i&1)]
        total+=endings(restrict(a,verts))[verts.index(q)]*H(restrict(a,rest))
    return total

def F(a,q):return left(left(a,q),len(a))
def p2(a):
    n=len(a);s=0
    for mask in range(1,(1<<n)-1):
        v=[i for i in range(n) if mask>>i&1]
        w=[i for i in range(n) if not(mask>>i&1)]
        s+=H(restrict(a,v))*H(restrict(a,w))
    require(s%2==0,'cover parity')
    return s//2

def main():
    graphs=ports=strong_ports=0
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            a=labelled(n,code);h=H(a);js=[];graphs+=1
            for q in range(n):
                k=K(a,q);j=k-h;js.append(j)
                require(j>=0,'negative broken-join slack')
                b=F(a,q);newq=n
                require(H(b)==h+2*k,'unconditional ear count')
                require(profile(b,newq)==(h+2*k,2*k,0),'unconditional output guard profile')
                c=F(b,newq)
                require(profile(c,n+2)==(3*(h+2*k),2*(h+2*k),0),'canonical normalization')
                if strong(a):
                    require(strong(b) and strong(c),'strong normalization')
                    strong_ports+=1
                if sum(a[i]>>q&1 for i in range(n))==1:
                    require(j==profile(a,q)[2],'guard-specialized slack')
                ports+=1
            cover2=p2(a)
            require(sum(js)==2*cover2-(n-1)*h,'sum of rooted slacks')
            wrapper=add(add(a,0),(1<<n)-1)
            require(H(wrapper)==3*h+2*cover2,'triangle wrapper comparison')
            if n>=2:require(min(js)<cover2,'strict rooted improvement')
    print('PASS:',graphs,'labelled tournaments,',ports,'arbitrary ports,',strong_ports,'strong input ports; unrestricted ear, guard specialization, canonical normalization, full-profile slack identity, and triangle-wrapper comparison.')
if __name__=='__main__':main()
