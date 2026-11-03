#!/usr/bin/env python3
from itertools import product

def require(p, msg):
    if not p: raise ValueError(msg)

def bank(limit):
    a=[0]*(limit+1); a[1]=2
    for n in range(2,limit+1):
        m=n//2
        a[n]=a[m]+(m if n%2==0 else int(m%2==0))
    return a

def paths(word):
    # 0 is apex; transitive base positions 1,...,len(word).
    n=len(word)+1; out=[0]*n
    for i,c in enumerate(word,1):
        if c=='1': out[0]|=1<<i
        else: out[i]|=1
        for j in range(i+1,n): out[i]|=1<<j
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            if not dp[mask][i]: continue
            bits=out[i]&~mask
            while bits:
                bit=bits&-bits; bits-=bit; j=bit.bit_length()-1
                dp[mask|bit][j]+=dp[mask][i]
    return sum(dp[-1])

a=bank(100000)
ncheck=0
for L in range(2,9):
    for middle in product('01', repeat=L-2):
        word='1'+''.join(middle)+'0'; m=int(word[:-1],2)
        require(paths(word)==m+a[m], ('DP mismatch',word))
        ncheck+=1
for m in range(1,len(a)):
    require(2<=a[m]<=m+1, ('bound failure',m))
    T=m+a[m]
    require((T-1+1)//2<=m<=T-2, ('search bound',m,T))
for T in [13,29,35,41]:
    lo=T//2; hi=T-2
    candidates=[m for m in range(lo,hi+1) if m+a[m]==T]
    require(not candidates, ('unexpected witness',T,candidates))
    print('restricted absence', T, 'interval',lo,hi,'values', [m+a[m] for m in range(lo,hi+1)])
print('verified binary words via independent subset DP:',ncheck)
print('verified bank states/bounds:',len(a)-1)
# Test the combinatorial dual involution independently on all words <=16 bits.
def dual_word(w):
    return ''.join('1' if c=='0' else '0' for c in w[::-1])
for L in range(2,17):
    for middle in product('01', repeat=L-2):
        w='1'+''.join(middle)+'0'; wd=dual_word(w)
        m=int(w[:-1],2); md=int(wd[:-1],2)
        require(dual_word(wd)==w, ('dual involution',w))
        require(m+a[m]==md+a[md], ('dual count',w,wd))
print('verified dual involution and bank equality for all strong words through length 16')
