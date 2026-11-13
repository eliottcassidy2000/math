# Monsky nullity nu(T) = dim_F2 ker(M^2 + M + I), M = adjacency mod 2 with diagonal making every row sum 0 (mod 2)
import collections
def decode(code,n):
    out=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1: out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
def nullity(out):
    n=len(out)
    M=[out[i] | ((bin(out[i]).count('1')&1)<<i) for i in range(n)]      # row bitmasks
    # M^2: row i = XOR of rows M[j] for j in row i
    M2=[]
    for i in range(n):
        r=0; x=M[i]
        while x:
            j=(x&-x).bit_length()-1; r^=M[j]; x&=x-1
        M2.append(r)
    A=[M2[i]^M[i]^(1<<i) for i in range(n)]
    rank=0; rows=A[:]
    for bit in range(n):
        piv=None
        for k in range(rank,n):
            if rows[k]>>bit&1: piv=k; break
        if piv is None: continue
        rows[rank],rows[piv]=rows[piv],rows[rank]
        for k in range(n):
            if k!=rank and rows[k]>>bit&1: rows[k]^=rows[rank]
        rank+=1
    return n-rank
strong0=set(); strongall=collections.Counter(); dist={}
for n in range(1,10):
    try: f=open('cls_%d.txt'%n)
    except: continue
    d=collections.Counter()
    for line in f:
        p=line.split('|')[0].split(); code=int(p[0]); H=int(p[1]); st=int(p[3])
        nu=nullity(decode(code,n)); d[nu]+=1
        if st==1 and nu==0: strong0.add(H)
    dist[n]=dict(sorted(d.items()))
    print('n=%d nu distribution'%n,dist[n],flush=True)
B=125   # Busch f(10): strong tournaments with H < 125 have order <= 9
S=set(h for h in strong0 if h<B)|{1}
grew=True
while grew:
    grew=False
    for a in list(S):
        for b in list(S):
            if a>1 and b>1 and a*b<B and a*b not in S: S.add(a*b); grew=True
gaps=[v for v in range(1,B,2) if v not in S]
print('strong nu=0 values below 125:',sorted(h for h in strong0 if h<B))
print('odd values below 125 never H of a nu=0 tournament (any order):',gaps)
