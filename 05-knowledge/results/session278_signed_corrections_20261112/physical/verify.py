"""Exact signed correction controls. Python stdlib; run from any directory."""
from itertools import product
from functools import lru_cache
import json
from math import comb

@lru_cache(None)
def g(p):
    if not p: return 1
    v=p//2
    return g(v) if p%2 else g(v)+v+v%2

def data(w):
    assert w
    p=int(w,2); q=sum(1<<j for j,b in enumerate(w) if b=='0')
    return ((g(p),q+1),(p+1,1<<len(w))),int(q%2==0),int(p%2==0)
def mv(B,x):return tuple(sum(a*b for a,b in zip(row,x)) for row in B)
def exact(words):
    a,w,d=1,0,0
    for word in words:
        B,e,h=data(word)
        a,w,d=B[0][0]*a+B[0][1]*w-e*d,B[1][0]*a+B[1][1]*w,h*a
    return a,w,d

def layers(words):
    # Coefficients x_m(z) in x_j(z)=B_j x_(j-1)(z)-z D_j x_(j-2)(z).
    prevprev=[]; prev=[(1,0)]; oldeta=0
    for word in words:
        B,eps,eta=data(word)
        cur=[mv(B,v) for v in prev]+[(0,0)]
        for r,v in enumerate(prevprev):
            a,b=cur[r+1];cur[r+1]=(a+eps*oldeta*v[0],b)
        while cur[-1]==(0,0):cur.pop()
        prevprev,prev,oldeta=prev,cur,eta
    return [v[0] for v in prev]

def direct(words):
    starts=[]; rows=[]
    for word in words:
        z=len(rows);starts.append(z);k=len(word)+1
        rows.extend([0]*k)
        for a in range(z):
            for b in range(z,z+k):rows[a]|=1<<b
        for i in range(1,k):
            if word[i-1]=='1':rows[z]|=1<<(z+i)
            else:rows[z+i]|=1<<z
            for j in range(i+1,k):rows[z+i]|=1<<(z+j)
    for a,b in zip(starts,starts[1:]):rows[a]^=1<<b;rows[b]|=1<<a
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for j in range(n):dp[1<<j][j]=1
    for s in range(1,1<<n):
        for j,v in enumerate(dp[s]):
            q=rows[j]&~s
            while q:
                bit=q&-q;q-=bit;dp[s|bit][bit.bit_length()-1]+=v
    return sum(dp[-1])

def check(words,brute=False):
    c=layers(words);h=exact(words)[0];s=0
    assert sum((-1)**r*v for r,v in enumerate(c))==h
    for r,v in enumerate(c):
        s+=(-1)**r*v
        assert (s-h)*(-1)**r>=0,(words,c,h,r)
    if brute:assert direct(words)==h
    return c,h

def main():
    small=['0','1','00','01','10','11'];nchecks=0
    for m in range(1,5):
        for ws in product(small,repeat=m):check(ws);nchecks+=1
    brute=0
    for m in range(1,4):
        for ws in product(small,repeat=m):check(ws,True);brute+=1
    # Repeated-10 characteristic recurrence, independently in integer arithmetic.
    hs=[exact(['10']*m)[0] for m in range(31)]
    assert all(hs[m+3]==7*hs[m+2]-4*hs[m+1]+4*hs[m] for m in range(28))
    # Independently expand each rational error-layer GF (8).
    for r in range(7):
        den=[1]
        for _ in range(r+1):
            nxt=[0]*(len(den)+2)
            for i,a in enumerate(den):
                for j,b in enumerate([1,-7,3]):nxt[i+j]+=a*b
            den=nxt
        coeff=[]
        for m in range(31):
            j=m-2*r
            num=comb(r+1,j)*(-4)**j if 0<=j<=r+1 else 0
            coeff.append(num-sum(den[k]*coeff[m-k] for k in range(1,min(m,len(den)-1)+1)))
            layer=layers(['10']*m) if m else [1]
            assert coeff[-1]==(layer[r] if r<len(layer) else 0)
    out={'algebraic_words_checked':nchecks,'independent_subset_DP':brute,'repeated_10':[]}
    for m in [2,3,4,5,10,20,50,100]:
        c,h=check(['10']*m);s=0;lo=1;hi=c[0];rows=[]
        for r,v in enumerate(c):
            s+=(-1)**r*v
            if r%2:lo=max(lo,s)
            else:hi=min(hi,s)
            if r in sorted(set([0,1,2,3,len(c)//2,len(c)-1])):
                width=hi-lo;b=width.bit_length()
                # Same interval/residue architecture as Session277, physical backend.
                residue=h%(1<<b);candidate=lo+((residue-lo)%(1<<b))
                assert candidate==h and candidate+(1<<b)>hi
                rows.append({'order':r,'lower':str(lo),'upper':str(hi),'residue_bits_sufficient':b})
        out['repeated_10'].append({'blocks':m,'vertices':3*m,'H':str(h),'blind':str(c[0]),'signed_error_blind_minus_H':str(c[0]-h),'correction_degree':len(c)-1,'certificates':rows})
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
