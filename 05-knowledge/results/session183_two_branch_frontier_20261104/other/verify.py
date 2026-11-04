#!/usr/bin/env python3
"""Independent finite controls for the proved all-length chart classification."""
import json
from pathlib import Path

def require(x,msg):
    if not x: raise RuntimeError(msg)

def j(f):
    # Chronological high-end births, independently of proof's low-bit recursion.
    g=e=z=0
    for i in range(f.bit_length()):
        b=(f>>i)&1
        z+=b*(g+e)
        g=2*g+b
        e=b
    return z

def j_pair(f):
    return sum((1<<(k-i-1))+(k==i+1) for k in range(f.bit_length()) for i in range(k) if (f>>i)&1 and (f>>k)&1)

def j_cover(f):
    m=f.bit_length(); full=(1<<m)-1; total=0
    # In TT_m each nonempty subset has exactly one directed path, its sorted list.
    for a in range(1,full):
        b=full^a
        last=a.bit_length()-1
        first=(b&-b).bit_length()-1
        total+= bool((f>>last)&1 and (f>>first)&1)
    return total

def rows_for(f,t):
    m=f.bit_length(); n=1+t+m; rows=[0]*n
    def edge(u,v): rows[u]|=1<<v
    for a in range(1,t+1):
        edge(0,a)
        for b in range(a+1,n): edge(a,b)
    for i in range(m):
        v=1+t+i
        edge(v,0) if (f>>i)&1 else edge(0,v)
        for k in range(i+1,m): edge(v,1+t+k)
    return rows

def hamilton(rows):
    n=len(rows); N=1<<n; dp=[[0]*n for _ in range(N)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,N):
        for v in range(n):
            z=dp[mask][v]
            if z:
                avail=rows[v]&~mask
                while avail:
                    bit=avail&-avail; avail-=bit
                    dp[mask|bit][bit.bit_length()-1]+=z
    return sum(dp[-1])

def strong(rows):
    n=len(rows); full=(1<<n)-1
    for root in range(n):
        reached=1<<root; last=0
        while reached!=last:
            last=reached
            for v in range(n):
                if (reached>>v)&1: reached|=rows[v]
        if reached!=full: return False
    return True

def main():
    covers=0
    for f in range(1,1<<11):
        require(j(f)==j_pair(f)==j_cover(f),('joint',f))
        covers+=1
    bound=1<<18
    js=[j(n) for n in range(3*bound+3)]
    def d(c,n):return 3*js[n]-js[3*n+c]
    hits=[]
    for n in range(bound):
        p=n&1
        require(js[2*n]==js[n] and js[2*n+1]==js[n]+n+p,('low recurrence',n))
        # Check all recurrences where the d-array index is available.
        if 2*n+1<bound:
            actual=[d(0,2*n),d(0,2*n+1),d(1,2*n),d(1,2*n+1),d(2,2*n),d(2,2*n+1)]
            expected=[d(0,n),d(1,n)+4*p-2,d(0,n)-3*n-p,d(2,n)+3*n+3*p,d(1,n),d(2,n)+2*p-2]
            require(actual==expected,('carry',n,actual,expected))
        require(-2*n<=d(0,n)<=2*n and -3*n<=d(1,n)<=4*n and -2*n<=d(2,n)<=2*n,('bounds',n))
        if n:
            odd=n//(n&-n)
            require(d(0,n)!=4 and (d(0,n)==6)==(odd==7),('classification',n))
            if d(0,n) in (4,6):hits.append(n)
    attachments=0
    for f in range(1,1<<7):
        for t in (1,2):
            rows=rows_for(f,t)
            require(strong(rows),('strong',f,t))
            require(hamilton(rows)==1-j(f)+(1<<t)*f,('attachment',f,t))
            attachments+=1
    examples=[]
    for a,b,t in [(0,0,1),(0,0,2),(1,0,1),(0,1,2),(2,1,1)]:
        k=a-b; f=7<<a; g=21<<b
        x=hamilton(rows_for(f,t));y=hamilton(rows_for(g,t+k))
        require(y==3*x+4,('example',a,b,t,x,y))
        examples.append(dict(source_mask=f,target_mask=g,source_clock=t,target_clock=t+k,source_H=x,target_H=y))
    out=dict(status='PASS',cover_checks=covers,carry_and_bounds_range=[0,bound-1],classification_range=[1,bound-1],plus4_masks=hits,plus2_masks=[],actual_attachment_checks=attachments,actual_shift_examples=examples,scope='Finite independent controls; all-length classification proved in PROOF.md.')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
