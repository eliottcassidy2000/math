from math import isqrt

def require(x):
    if not x: raise ValueError('7p context failed')
A=[2,28,25,81,97,79,7]
B=[2,80,3,7,173,143,189,15]
At=[a|128 for a in A]+[0]
def tournament(mod):
    w='10100010'
    n=len(mod)
    require(n==8)
    m=n+len(w)
    out=[0]*m
    for i in range(m-1):
        for j in range(i+1,m-1):
            if j<n and not mod[i]>>j&1:out[j]|=1<<i
            else:out[i]|=1<<j
        bit=w[0] if i<n else w[i-n+1]
        if bit=='1':out[m-1]|=1<<i
        else:out[i]|=1<<(m-1)
    return out
def count(out):
    n=len(out);full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,full):
        for i,x in enumerate(dp[s]):
            if x:
                t=out[i]&~s
                while t:
                    b=t&-t;t^=b;dp[s|b][b.bit_length()-1]+=x
    return sum(dp[-1])
require(count(tournament(B))==count(tournament(At))==7*33023)
require(all(33023%d for d in range(2,isqrt(33023)+1)))
print('PASS: both literal tournaments realize 231161 = 7*33023')
