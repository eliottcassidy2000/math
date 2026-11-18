from itertools import product, permutations
D=[x for x in product(range(3),range(2),range(2)) if x not in [(0,0,0),(2,1,1)]]
C=lambda t:(2-t[0],1-t[1],1-t[2])
Q=lambda t:(t[0],t[2],t[1])
f=lambda t:(0,0) if t[1]==t[2] else (t[1]-t[2],2*abs(t[0]-1)-1)
assert len({f(t) for t in D})==5
for t in D:
 x,y=f(t)
 assert f(C(t))==f(Q(t))==(-x,y)
 assert f(C(Q(t)))==f(t)
counts={}
for bits in product(range(2),repeat=3):
 edges={(0,1),(1,2),(2,0)}|{(3,i) if bits[i] else (i,3) for i in range(3)}
 h=sum(all((a,b) in edges for a,b in zip(p,p[1:])) for p in permutations(range(4)))
 counts[bits]=h
 assert h==(3 if len(set(bits))==1 else 5)
print('PASS five-state quotient and all eight triangle extensions:',counts)
# Pure compatibility counting identity, independent of any physical-gluing claim.
primes=[3,5,7,11,13,17,19,23,29,31]
for E in range(6,63,2):
 pairs=[(p,q) for p in primes for q in primes if p<=q and p+q==E]
 good=[x for x in pairs if x[0]!=7 and x[1]!=7 and x!=(11,11)]
 bad=[x for x in pairs if x not in good]
 assert len(pairs)-len(bad)==len(good)
print('PASS abstract compatibility identity; no physical operation assumed')
