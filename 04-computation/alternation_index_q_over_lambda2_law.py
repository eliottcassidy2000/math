"""Why the DLO sequence out-alternates the golden ratio: the index tracks Q/lambda^2."""
from math import comb, log, sqrt
def circuit(h):
    R=[h[k]*h[k]/(h[k-1]*h[k+1]) for k in range(1,len(h)-1)]
    return [log(R[i]/R[i-1]) for i in range(1,len(R))]
def alt_index(P):
    s=max(abs(x) for x in P) or 1
    return max(abs(sum((-1)**(j-i)*comb(j,i)*P[i] for i in range(j+1)))/(2**j*s)
               for j in range(1,len(P)))
def seq(P,Q,N=16):
    a=[1,1]
    while len(a)<N: a.append(P*a[-1]+Q*a[-2])
    return a
print("For a_k = P a_{k-1} + Q a_{k-2}, Simson gives a_{k-1}a_{k+1}-a_k^2 = (-Q)^{k-1} c,")
print("so R_k - 1 ~ +-(Q/lambda^2)^k : the alternation MAGNITUDE decays at rate Q/lambda^2.")
print()
print("  (P,Q)     lambda    Q/lambda^2    alternation index")
rows=[(1,1,"Fibonacci/golden"),(2,1,"Pell/silver"),(3,1,"bronze"),(2,5,"DLO count"),(1,2,""),(1,3,""),(2,7,"")]
data=[]
for P,Q,name in rows:
    lam=(P+sqrt(P*P+4*Q))/2
    r=Q/lam**2
    idx=alt_index(circuit(seq(P,Q)))
    data.append((r,idx,P,Q,name))
    print(f"  ({P},{Q})   {lam:7.4f}   {r:9.6f}     {idx:.6f}   {name}")
data.sort()
mono=all(data[i][1]<=data[i+1][1]+1e-9 for i in range(len(data)-1))
print(f"\nindex monotone increasing in Q/lambda^2 across all tested: {mono}")
