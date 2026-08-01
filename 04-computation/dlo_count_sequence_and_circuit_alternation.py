"""The DLO-count sequence a_k = 2a_{k-1} + 5a_{k-2}: arithmetic check, and its
circuit/alternation profile against klein's THM-3010 metallic family."""
from math import comb, log
from fractions import Fraction

# 1. verify the transfer-matrix values and the enumeration
def a(k):
    v=(1,3); M=((1,6),(1,1))
    for _ in range(k-1):
        v=(v[0]*M[0][0]+v[1]*M[1][0], v[0]*M[0][1]+v[1]*M[1][1])
    return v[0]*3+v[1]*1
A=[a(k) for k in range(1,7)]
print("a_k from (1,3)[[1,6],[1,1]]^{k-1}(3,1)^T :", A)
print("recurrence a_k = 2a_{k-1}+5a_{k-2} holds :",
      all(A[i]==2*A[i-1]+5*A[i-2] for i in range(2,len(A))))
print("char poly x^2-2x-5 -> growth 1+sqrt6 =", 1+6**0.5)
print()
# 2. re-derive the realizable list < 1000 (multinomials over component types)
def partitions(n):
    if n==0: yield []; return
    for first in range(n,0,-1):
        for rest in partitions(n-first):
            if not rest or first>=rest[0]: yield [first]+rest
from math import factorial
vals={1}
for k in range(1,7):
    ak=a(k)
    for p in partitions(k):
        mult=factorial(k)
        for part in p: mult//=factorial(part)
        if ak*mult<1000: vals.add(ak*mult)
print("realizable values < 1000 :", sorted(vals))
print("count", len(vals), " sum =", sum(sorted(vals)))
print()
# 3. circuit / alternation index, vs klein's metallic recurrences
def circuit(h):
    R=[h[k]*h[k]/(h[k-1]*h[k+1]) for k in range(1,len(h)-1)]
    return [log(R[i]/R[i-1]) for i in range(1,len(R))]
def alt_index(P):
    s=max(abs(x) for x in P) or 1
    return max(abs(sum((-1)**(j-i)*comb(j,i)*P[i] for i in range(j+1)))/(2**j*s)
               for j in range(1,len(P)))
def metallic(n,N):
    x=[0,1]
    while len(x)<N: x.append(n*x[-1]+x[-2])
    return x[1:]
seq=[1,1]
while len(seq)<16: seq.append(2*seq[-1]+5*seq[-2])
c=circuit(seq)
sg=''.join('+' if v>0 else '-' for v in c)
print(f"DLO sequence (a_k=2a+5a) circuit {sg}")
print(f"   sign changes = {sum(1 for i in range(1,len(sg)) if sg[i]!=sg[i-1])} of {len(sg)-1}"
      f"   alternation index = {alt_index(c):.6f}")
for n in (1,2,3):
    cc=circuit(metallic(n,16))
    ss=''.join('+' if v>0 else '-' for v in cc)
    print(f"metallic n={n}  circuit {ss}   index = {alt_index(cc):.6f}")
