from itertools import product
from math import isqrt
from verify import count

def rc(s): return s[::-1].translate(str.maketrans('01','10'))

def exact_rows(w):
    n=1+len(w)
    r=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    for i,b in enumerate(w,1):
        if b=='0':
            r[0]^=1<<i
            r[i]^=1
    return r

def exact_dp(w):
    r=exact_rows(w); n=len(r); N=1<<n
    dp=[[0]*n for _ in range(N)]
    for i in range(n): dp[1<<i][i]=1
    for s,a in enumerate(dp):
        for i,h in enumerate(a):
            if not h: continue
            rest=r[i]&~s
            while rest:
                b=rest&-rest
                rest^=b
                dp[s|b][b.bit_length()-1]+=h
    return sum(dp[-1])

def main():
    for n in range(1,14):
        for s in product('01',repeat=n-1):
            u='1'+''.join(s)
            X=int(u,2); p=int(u[-1])
            l=64*X*X+63*X+2*count(u)+13+p
            assert count(u+'010010'+rc(u))==l
            assert count(u+'100001'+rc(u))==l+2
    # The simultaneous old arc reversals are at central edges 1-2 and 5-6.
    u='1100'; a=u+'010010'+rc(u); b=u+'100001'+rc(u)
    r=exact_rows(a)
    for j in (len(u)+1,len(u)+5):
        r[j]^=1<<(j+1)
        r[j+1]^=1<<j
    perm=list(range(len(r)))
    for j in (len(u)+1,len(u)+5): perm[j],perm[j+1]=perm[j+1],perm[j]
    re=[sum(((r[v]>>z)&1)<<i for i,z in enumerate(perm)) for v in perm]
    assert re==exact_rows(b)
    for w,h in ((a,10007),(b,10009)):
        assert exact_dp(w)==h
        assert all(h%d for d in range(2,isqrt(h)+1))
    print('Formula checks through exterior length 13; exact 15-vertex twin-prime pair 10007/10009.')

if __name__=='__main__':main()
