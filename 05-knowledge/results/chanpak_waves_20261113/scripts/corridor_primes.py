# After an O (E = 0) a run of k root changes gives H = A_k H0 + B_k lam0, (A_k, B_k) = (H, lam) entries of the
# first row of J^k.  A prime p is forced (p | H for every source) iff the row e1 J^k is c*(0,1,0) mod p.
# The row e1 J^k mod p is periodic; we scan one full period, so "never" is exact.
from sympy import primerange
J=[[1,-1,1],[1,-1,0],[1,0,1]]
def rowmul(r,p): return tuple(sum(r[t]*J[t][j] for t in range(3))%p for j in range(3))
never=[]; forced=[]
for p in primerange(2,400):
    start=(1,0,0); r=start; k=0; ks=[]
    while True:
        r=rowmul(r,p); k+=1
        if r[0]==0 and r[2]==0: ks.append(k)
        if r==start: period=k; break
    if ks: forced.append((p,ks,period))
    else: never.append(p)
for p,ks,per in forced[:14]:
    print("p=%d: forcing lengths k = %s mod %d"%(p,ks,per))
print("forced primes < 400:",[p for p,_,_ in forced])
print("never forced primes < 400:",never)
