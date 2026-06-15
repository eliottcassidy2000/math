# Direct check of the firing logic: an INEXACT relation m != 0 fires at q=p^e iff p^e | m.
# Build a set with a controlled single small relation that is NOT exact, verify D/p^e at
# p^e | m is depressed and at p^e not dividing m is at main, then -> main as e grows.
import math
def deficit(q, S):
    rad = q//14; c=0
    for a in range(q):
        ok=True
        for v in S:
            r=(v*a)%q
            if r<=rad or r>=q-rad: ok=False; break
        if ok: c+=1
    return c
MAIN=(6/7)**13
# A set where speeds have an EXACT relation 1+1 = 2 (i.e. 1*1+1*1-1*2=0 trivially) -- skip.
# Instead: pick S with NO small exact relation but a near-relation m=Sum. Use S so that
# v1+v2-v3 = m for some target m; check firing at q | m only.
# S: 14, and 10,11 with 10+11-21 = 0 EXACT vs 10+11-20 = 1 (m=1, fires nowhere) etc.
# Demonstrate exact (m=0): {10,11,21} relation 10+11-21=0 fires at EVERY q.
# vs inexact: {10,11,22} relation 10+11-22=-1 (m=-1) fires at NO q>1.
S_exact = sorted(set([14,10,11,21] + [3,5,17,33,65,129,257,513,1025]))[:13]
S_inexact = sorted(set([14,10,11,22] + [3,5,17,33,65,129,257,513,1025]))[:13]
print("S_exact (has 10+11-21=0):  ", S_exact)
print("S_inexact (10+11-22=-1):   ", S_inexact)
for label,S in [("exact",S_exact),("inexact",S_inexact)]:
    print(f"  {label}: large-q L =", round(sum(deficit(q,S)/q for q in (14000,15013,16007))/3,4),
          " main=",round(MAIN,4))
    for p in (2,3,5):
        seq=[]
        e=1
        while p**e<=50000:
            if p**e>=28: seq.append(deficit(p**e,S)/p**e)
            e+=1
        print(f"     p={p}: tail {[round(x,4) for x in seq[-3:]]}")
