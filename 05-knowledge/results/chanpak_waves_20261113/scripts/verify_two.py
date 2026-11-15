import sys
def parse(code,m):
    out=[0]*m; k=0
    for i in range(m):
        for j in range(i+1,m):
            if code[k]=='1': out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
def ham(out,n):
    N=1<<n; dp=[[0]*n for _ in range(N)]
    for v in range(n): dp[1<<v][v]=1
    for S in range(1,N):
        row=dp[S]
        for v in range(n):
            x=row[v]
            if not x: continue
            nb=out[v]&~S
            while nb:
                u=(nb&-nb).bit_length()-1; nb&=nb-1
                dp[S|1<<u][u]+=x
    return sum(dp[N-1])
bad=0; cnt=0
for line in open(sys.argv[1]):
    p=line.split()
    if p[1]=='NONE': print("no witness",p[0]); bad+=1; continue
    N=int(p[0]); mL=int(p[1]); cL=p[2]; VL=int(p[3]); mB=int(p[4]); cB=p[5]; VB=int(p[6])
    oL=parse(cL if cL!='-' else '',mL); oB=parse(cB if cB!='-' else '',mB)
    n=mL+mB; out=[0]*(n+1)
    for i in range(mL): out[i]=oL[i] | (((1<<mB)-1)<<mL)      # L => B
    for i in range(mB):
        out[mL+i]=oB[i]<<mL
    marks=VL | (VB<<mL)
    for v in range(n):
        if not (marks>>v)&1: out[v]|=1<<n
    out[n]=marks
    h=ham(out,n+1); cnt+=1
    if h!=N: print("MISMATCH",line.strip(),h); bad+=1
print(f"verified {cnt} witnesses by direct Hamiltonian-path DP, mismatches/missing {bad}")
