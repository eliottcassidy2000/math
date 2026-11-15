from random import Random
from itertools import product

def hp(out, start=False):
    n=len(out)
    inc=[sum(1<<j for j in range(n) if out[j]>>i&1) for i in range(n)]
    edges=inc if start else out
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for s in range(1,1<<n):
        for i,x in enumerate(d[s]):
            if not x: continue
            z=edges[i]&~s
            while z:
                v=z&-z
                z-=v
                d[s|v][v.bit_length()-1]+=x
    H=[1]+[sum(row) for row in d[1:]]
    return H,d

def transitive(n): return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def converse(out): return [sum(1<<j for j,a in enumerate(out) if a>>i&1) for i in range(len(out))]
def join(blocks):
    n=sum(map(len,blocks)); out=[]; r=0; all_mask=(1<<n)-1
    for B in blocks:
        out.extend((b<<r)|(all_mask^((1<<(r+len(B)))-1)) for b in B)
        r+=len(B)
    return out
def apex(out,marks):
    n=len(out); return [r|((1<<n) if not (marks>>i&1) else 0) for i,r in enumerate(out)]+[marks]
def state(out,marks):
    H,D=hp(out,True); all_mask=len(H)-1
    z=[sum(x for i,x in enumerate(D[s]) if marks>>i&1)*H[all_mask^s]
       for s in range(1,len(H))]
    return H[-1],hp(apex(out,marks))[0][-1],sum(z),sum(
        x if (i+1).bit_count()%2 else -x for i,x in enumerate(z)),sum(
        h*H[all_mask^s] for s,h in enumerate(H))
def core(u,v):
    return [u,1-v,0,v,u]
def litmirror(out,marks,u,v):
    d=len(out); m=(1<<d)-1
    c=sum(x<<i for i,x in enumerate(core(u,v)))
    return apex(join([out,transitive(5),converse(out)]),marks|(c<<d)|((m^marks)<<(d+5)))
def pred(g,n,s,d,u,v):
    return 2*g*n+(7+2*u-6*v)*g*g+(37-12*v)*g*s+g*d+32*s*s
def p(x):
    z=0
    while x:
        u,v=divmod(x,2)
        if not v:z+=u+(u&1)
        x=u
    return z
def half_inv(w):
    return sum(1<<max(0,j-i-2) for i in range(len(w)) for j in range(i+1,len(w)) if w[i] and not w[j])
def run():
    rng=Random(44073)
    for m in range(1,5):
        for _ in range(16):
            s=[0]*m
            for i in range(m):
                for j in range(i+1,m):
                    if rng.randrange(2):s[i]|=1<<j
                    else:s[j]|=1<<i
            for mk in [rng.randrange(1<<m),(1<<m)-1]:
                g,n,S,d,Q=state(s,mk)
                for u,v in product((0,1),repeat=2):
                    N=hp(litmirror(s,mk,u,v))[0][-1]
                    assert N==pred(g,n,S,d,u,v)
                    if mk==(1<<m)-1:
                        B=27+12*v; C=5+6*v+2*u
                        D=B*B-128*C
                        assert N==32*Q*Q-B*g*Q+C*g*g
                        assert (64*Q-B*g)**2-D*g*g==128*N
                        if v==1 and u==0:
                            x=32*Q-15*g
                            assert x*x-9*x*g-8*g*g==32*N
    for X in range(1,1025):
        w=list(map(int,bin(X)[2:])); wr=[1-x for x in w[::-1]]
        prev=None
        for u,v in product((0,1),repeat=2):
            N=1+2*half_inv(w+core(u,v)+wr)
            assert N==32*X*X+(37-12*v)*X+9+2*u-6*v+(X&1)+2*p(X)
        assert 32*(X+1)**2+25*(X+1)+3+((X+1)&1)+2*p(X+1)>32*X*X+37*X+11+(X&1)+2*p(X)
    assert all(pow(D%37,18,37)==36 for D in [89,-167,113,-143])
    print("512 literal full path counts, 4096 transitive mirror counts, four norm sheets and prime obstruction passed.")

    menus=[
      (1,False,(("0011",-39),("011001",-23),("001100",897))),
      (2,False,(("01101001",-1007),("0",-15),("0000000",15105))),
      (1,True,(("01100001",-663),("0110001",-391),("001100",897))),
      (2,True,(("01001101",-983),("0100110",-87),("001110100",85521))),
    ]
    from math import isqrt,prod
    seen=set()
    for target,times7,words in menus:
        p0=prod(d for _,d in words); assert isqrt(p0)**2==p0
        for w,d0 in words:
            m=len(w); R=1<<m; s=int(w,2); t=R-1-int(w[::-1],2)
            a=s+t; e=(s&1)+1-int(w[0])
            B=a-2*R; C=1+p(s)+e-a+R
            D=B*B-4*R*C
            assert D*d0>0 and D%d0==0 and isqrt(D//d0)**2==D//d0
            assert target in ((R+C+B)%3,(R+C-B)%3)
            if times7: assert pow(D%7,3,7)==1
            # Independently count the literal frame on a single-vertex exterior.
            c=sum(int(x)<<j for j,x in enumerate(w))
            assert hp(apex(transitive(m+2),1|(c<<1)))[0][-1]==4*R+2*B+C
    print("All twelve center choices, joint mod-three gates and literal forms passed.")

if __name__=="__main__":run()
