# Brute-force / DFS check of Z_k = #{residues n mod 4*5^(k-1): last k digits of 2^n (n>=k) all nonzero}
# and the equivalent count #{k-digit zeroless strings x with 2^k | x}.
import sys
def Z_strings(K):
    # DFS over q-process: q_{i+1} = (q_i + a*5^i)/2, a = q_i mod 2
    res=[0]*(K+1); odd=[0]*(K+1)
    stack=[(0,0)]  # (i, q_i)
    while stack:
        i,q=stack.pop()
        res[i]+=1
        if q&1: odd[i]+=1
        if i==K: continue
        p5=5**i
        for a in range(1,10):
            if (a-q)&1: continue
            stack.append((i+1,(q+a*p5)//2))
    return res,odd
def Z_residues(k):
    T=4*5**(k-1); M=10**k
    cnt=0
    x=pow(2,k+T,M)  # start at n = k+T (n>=k), iterate over a full period
    n0=k+T
    for j in range(T):
        s=str(x).zfill(k)
        if '0' not in s: cnt+=1
        x=(2*x)%M
    return cnt
K=int(sys.argv[1]) if len(sys.argv)>1 else 9
res,odd=Z_strings(K)
for k in range(K+1):
    zr = Z_residues(k) if 1<=k<=7 else None
    E=res[k]-odd[k]
    print(k, res[k], 'resid-check', zr, 'O',odd[k],'E',E,'Delta=O-E',odd[k]-E, 'Z*(2/9)^k=%.12f'%(res[k]*(2/9)**k))
