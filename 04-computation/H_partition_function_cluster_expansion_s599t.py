"""PARTITION FUNCTIONS EVERYWHERE. The total H = a partition function: Z_n = sum_{T on n} H(T) =
#(tournament, Ham-path) pairs = n! * 2^(C(n,2)-(n-1)) (each of n! orderings is a Ham path in
2^(C(n,2)-(n-1)) tournaments). H multiplicative over STRONG components (= disjunctive game sum)
=> the SEQUENCE/cluster expansion T_H = 1/(1 - S_H), recursion b_n = sum_k C(n,k) a_k b_{n-k},
a_k = strong (connected) H-sum. Verify. Connect: game-value mex-recursion spectrum (Hamkins
infinite go) <-> H-value multiplicative-recursion spectrum (gaps 7,21). opus-2026-06-03-S599t."""
from itertools import combinations, product
from math import comb, factorial
def Hcount(n,adj):
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        row=dp[mask]
        for v in range(n):
            c=row[v]
            if not c: continue
            av=adj[v]
            for w in range(n):
                if not(mask>>w&1) and av>>w&1: dp[mask|1<<w][w]+=c
    return sum(dp[size-1][v] for v in range(n))
def is_strong(n,adj):
    def reach(nb):
        seen={0}; st=[0]
        while st:
            x=st.pop()
            for y in range(n):
                if y not in seen and nb(x,y): seen.add(y); st.append(y)
        return len(seen)==n
    return reach(lambda x,y:adj[x]>>y&1) and reach(lambda x,y:adj[y]>>x&1)
def main():
    print("(1) Z_n = sum_{T} H(T)  vs closed form  n! * 2^(C(n,2)-(n-1))  [the partition function]")
    a={1:0}; b={}  # a_n strong H-sum, b_n all H-sum
    for n in range(1,6):
        allH=0; strongH=0
        for bits in product((0,1),repeat=n*(n-1)//2):
            adj=[0]*n
            for (i,j),bt in zip(combinations(range(n),2),bits):
                if bt: adj[i]|=1<<j
                else: adj[j]|=1<<i
            h=Hcount(n,adj); allH+=h
            if is_strong(n,adj): strongH+=h
        a[n]=strongH; b[n]=allH
        closed=factorial(n)*2**(comb(n,2)-(n-1))
        print(f"  n={n}: Z_n=sum H={allH}; closed form={closed}; match={allH==closed}; strong-sum a_n={strongH}")
    print("\n(2) CLUSTER EXPANSION recursion b_n = sum_{k=1..n} C(n,k) a_k b_{n-k}  (b_0=1)")
    b0={0:1}; b0.update(b)
    for n in range(1,6):
        rec=sum(comb(n,k)*a[k]*b0[n-k] for k in range(1,n+1))
        print(f"  n={n}: recursion gives {rec}; direct Z_n={b[n]}; match={rec==b[n]}  (a=strong/connected pieces)")
    print("\n(3) FREE ENERGY log Z_n is EXTENSIVE (≈ linear in #vertices => intensive free energy/vertex):")
    import math
    for n in range(2,6):
        print(f"  n={n}: log Z_n={math.log(b[n]):.3f}; per-vertex={math.log(b[n])/n:.4f}; (Z_n=n!2^(C(n,2)-(n-1)))")
    print("\n(4) Rédei parity = partition function at the SIGN fugacity z=-1 (fermionic/Lee-Yang):")
    print("    H(T) is odd for every T (Rédei) => Z evaluated mod 2 / the GF(2) 'permanent->determinant'")
    print("    is the z=-1 slice; the H-spectrum gaps {7,21} are forbidden 'energy levels'.")
if __name__=='__main__': main()
