from itertools import product

def formula(w):
    L=len(w); h,p=1,2**(L+1)
    for i,x in enumerate(w):
        for j in range(i+1,L):
            if x and not w[j]:
                d=j-i
                h+=2**(d-1)+(d==1)
                p+=2**(L-d)*(3**(d-1)+1)
    return h,p

def literal(w):
    L=len(w); n=L+1
    rows=[sum(1<<j for j in range(i+1,L))+((not w[i])<<L) for i in range(L)]
    rows.append(sum(x<<i for i,x in enumerate(w)))
    M=1<<n
    dp=[[0]*n for _ in range(M)]
    h=[0]*M; h[0]=1
    for s in range(1,M):
        for j in range(n):
            if (s>>j)&1:
                a=s^(1<<j)
                dp[s][j]=sum(dp[a][i] for i in range(n) if rows[i]>>j&1) if a else 1
        h[s]=sum(dp[s])
    return h[-1],sum(h[s]*h[M-1-s] for s in range(M))

count=0
for L in range(8):
    for w in product(range(2), repeat=L):
        h,p=formula(w)
        if literal(w)!=(h,p): raise ArithmeticError(w)
        idx=[i for i in range(L) if w[i]]
        jdx=[j for j in range(L) if not w[j]]
        if idx and jdx and min(idx)<max(jdx):
            i,j=min(idx),max(jdx); m=j-i+1
            hc,pc=formula(w[i:j+1])
            if (h,p)!=(hc,2**(L-m)*pc) or pc<2*(3**(m-2)+1): raise ArithmeticError(w)
        elif (h,p)!=(1,2**(L+1)): raise ArithmeticError(w)
        count+=1
print(f'PASS: {count} independent literal profile checks')
