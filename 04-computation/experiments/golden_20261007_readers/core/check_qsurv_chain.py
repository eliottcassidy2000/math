# cross-check: python chain (coal_check.step logic) vs C bignum chain on fixed bit sequences, via 'single' mode is
# not traceable; instead compare survival of trans 1 for a small run with a pure-python simulation with the same law.
import random, sys
def step(k, N, beta):
    if k >= 0:
        if N % 2 == 0:
            if beta == 0: return k, N // 2
            return k, (3*N + 1 - 3**k) // 2
        if beta == 0: return k + 1, (3*N + 1) // 2
        if k >= 1: return k - 1, (N - 3**(k-1)) // 2
        return -1, (3*N - 1) // 2
    a = -k
    if N % 2 == 0:
        if beta == 0: return k, N // 2
        return k, (3*N + 3**a - 1) // 2
    if beta == 0: return k + 1, (N + 3**(a-1)) // 2
    return k - 1, (3*N - 1) // 2
rng = random.Random(1); P=20000; TM=2000
grid=[10,100,1000,2000]; surv=[0]*4
for _ in range(P):
    k,N=0,1; tau=None
    for t in range(TM):
        k,N=step(k,N,rng.getrandbits(1))
        if k==0 and N==0: tau=t+1;break
    for g,T in enumerate(grid):
        if tau is None or tau>T: surv[g]+=1
print([ (T, s/P, round((T**0.5)*s/P,3)) for T,s in zip(grid,surv)])
