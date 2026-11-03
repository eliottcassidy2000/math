from collections import Counter
from math import gcd,isqrt
P=lambda n:n>1 and all(n%d for d in range(2,isqrt(n)+1))
V=[x for x in range(290) if gcd(x,290)==1]
if not (len(V)==112): raise RuntimeError("independent check failed: line 5")
E=[(x,(x+2)%290) for x in V if gcd((x+2)%290,290)==1]
if not (len(E)==81): raise RuntimeError("independent check failed: line 7")
# Generate components by peeling endpoints, independent of producer's runs logic.
adj={x:set() for x in V}
for x,y in E: adj[x].add(y);adj[y].add(x)
comp=[]; unseen=set(V)
while unseen:
    stack=[unseen.pop()]; C=[]
    while stack:
        x=stack.pop();C.append(x)
        for y in adj[x]&unseen: unseen.remove(y);stack.append(y)
    comp.append(C)
if not (Counter(map(len,comp))=={4:25,3:2,2:2,1:2}): raise RuntimeError("independent check failed: line 18")
if not (sorted(c[0] for c in comp if len(c)==1)==[27,263]): raise RuntimeError("independent check failed: line 19")
for h in range(290):
    if not (sum(gcd(c-h,290)==gcd(c+h,290)==1 for c in range(290))==(4 if h%5==0 else 3)*(28 if h%29==0 else 27)): raise RuntimeError("independent check failed: line 21")
plus=[n for n in range(2,6092) if n%290==1 and P(n)];minus=[n for n in range(2,6092) if n%290==289 and P(n)]
if not (plus[0]==1451 and minus[0]==2029): raise RuntimeError("independent check failed: line 23")
cross=[(n,n+2) for n in range(289,6090,290) if P(n) and P(n+2)]
if not (cross[0]==(6089,6091)): raise RuntimeError("independent check failed: line 25")
if not (sum(P(n) and gcd(n,290)==1 for n in range(2,1452))==227): raise RuntimeError("independent check failed: line 26")
# CRT direct 5 and 29 orbit products.
D={a:((a*a-29)//4)%290 for a in range(1,290,2)}
R={}
for a in D:
    path=[]; x=a
    while x not in path: path.append(x); x=D[x]
    mu=path.index(x); c=path[mu:]; canon=tuple(c[c.index(min(c)):]+c[:c.index(min(c))])
    R.setdefault(canon,Counter())[mu]+=1
if not (len(R)==8): raise RuntimeError("independent check failed: line 35")
if not (Counter(sum(h.values()) for h in R.values())=={5:1,20:7}): raise RuntimeError("independent check failed: line 36")
if not (all(len(c)==3 and h==({0:3,1:2} if sum(h.values())==5 else {0:3,1:7,2:10}) for c,h in R.items())): raise RuntimeError("independent check failed: line 37")
for a in range(-7,8,2):
    x=a
    for _ in range(8): x=(x*x-29)//4
    if not (x in (-7,5,-1)): raise RuntimeError("independent check failed: line 41")
print('PASS: phi=112; edge/paths, all 290 centers, least primes, crossing twin, rank and 145-state dynamics independently reproduced')
print('cycles and tail counts:', R)
