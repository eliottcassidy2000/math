"""Independent small controls: explicit path covers, subset path DP, differences."""
from itertools import permutations
from math import comb, factorial, isqrt
from random import Random


def ck(condition, tag):
    if not condition:
        raise RuntimeError(tag)


def tt(n):
    return [sum(1 << j for j in range(i + 1, n)) for i in range(n)]


def tournament(n, mask):
    ans = [0]*n
    for j in range(n):
        for i in range(j):
            a, b = (i,j) if mask & 1 else (j,i)
            ans[a] |= 1 << b
            mask >>= 1
    return ans


def converse(a):
    return [sum(1 << i for i in range(len(a)) if a[i] >> j & 1)
            for j in range(len(a))]


def blow(a, blocks):
    s = [0]
    for b in blocks:
        s.append(s[-1]+len(b))
    return [(row << s[i]) | sum(((1 << len(c))-1) << s[j]
             for j,c in enumerate(blocks) if a[i] >> j & 1)
            for i,b in enumerate(blocks) for row in b]


def h(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for S, row in enumerate(dp):
        for i, value in enumerate(row):
            if not value:
                continue
            for j in range(n):
                if a[i] >> j & 1 and not S >> j & 1:
                    dp[S | (1 << j)][j] += value
    return sum(dp[-1])


def q(a):
    n = len(a)
    out = [0]*(n+3)
    for w in permutations(range(n)):
        mandatory = sum(not (a[i] >> j & 1) for i,j in zip(w,w[1:]))
        for cuts in range(mandatory,n):
            out[cuts+1] += comb(n-1-mandatory,cuts-mandatory)
    return out


def run_coeff(q):
    n = len(q)-3
    return [0]+[q[j-1]**2+2*q[j-1]*q[j]+3*q[j]**2+
                2*q[j]*q[j+1]+q[j+1]**2 for j in range(1,n+2)]


def powers(A):
    return [0]+[sum((-1)**(j-i)*comb(j,i)*A[j] for j in range(i,len(A)))
                for i in range(1,len(A))]


def decode(n,m,N):
    b = [0]*(n+2)
    for i in range(n+1,0,-1):
        base = i**m
        b[i] = (2*N+base)//(2*base)
        N -= b[i]*base
    A = [0]*(n+2)
    u = [0]*(n+3)
    for i in range(n+1,0,-1):
        A[i] = b[i]-sum((-1)**(j-i)*comb(j,i)*A[j] for j in range(i+1,n+2))
        x = A[i]-2*u[i]**2-2*u[i]*u[i+1]-u[i+1]**2
        r = isqrt(x)
        ck(r*r == x, 'square')
        u[i-1] = r-u[i]
    ck(N == 0, 'remainder')
    return u


def main():
    cyc = [2,4,1]
    rng = Random(121)
    profiles, direct, differences, modular_switches = 0,0,0,0
    for n in range(1,8):
        masks = range(1 << comb(n,2)) if n <= 4 else [
            rng.randrange(1 << comb(n,2)) for _ in range(15)]
        K = factorial(n)*2**(n-1)
        M = 9*(n+1)*K*K*2**(n+1)
        m = (n+1)*(1+(2*M-1).bit_length())
        for mask in masks:
            x = tournament(n,mask)
            u = q(x)
            b = powers(run_coeff(u))
            ck(sum(map(abs,b)) <= M,'bound')
            ck(decode(n,m,sum(v*i**m for i,v in enumerate(b))) == u,'decode')
            for j in range(1,n+1):
                ck((-1)**n*u[j] == sum((-1)**k*comb(k-1,j-1)*u[k]
                     for k in range(j,n+1)),'reciprocity')
            profiles += 1
            if n <= 5:
                for s in range(1,min(5,13-2*n)):
                    ck(h(blow(cyc,[x,converse(x),tt(s)])) ==
                       sum(v*i**s for i,v in enumerate(b)),'strong count')
                    direct += 1
    for n in range(1,5):
        for mask in range(1 << comb(n,2)):
            a = tournament(n,mask)
            for i in range(n):
                local = a[i]
                active = any(local >> j & 1 and a[j] >> k & 1 and
                             a[k] >> i & 1 for j in range(n) for k in range(n))
                values = []
                for k in range(1,7):
                    values.append(h(blow(a,[tt(k if j==i else 1) for j in range(n)])))
                for order in range(1,4):
                    for j in range(len(values)-order):
                        diff = sum((-1)**(order-t)*comb(order,t)*values[j+t]
                                   for t in range(order+1))
                        ck(diff >= 2**(j+1) if active else diff == 0,'finite difference')
                        differences += 1
                left = blow(a,[tt(3 if j==i else 1) for j in range(n)])
                right = blow(a,[cyc if j==i else tt(1) for j in range(n)])
                f,g,u = q(left),q(right),q(a)
                ck(all(g[j]-f[j] == 2*(u[j] if j<len(u) else 0)
                       for j in range(len(f))), 'all slot switching')
                modular_switches += 1
    print(dict(profiles=profiles,direct=direct,forward_differences=differences,
               full_profile_switches=modular_switches))


if __name__ == '__main__':
    main()
