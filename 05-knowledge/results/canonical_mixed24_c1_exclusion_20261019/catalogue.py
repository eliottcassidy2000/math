from itertools import product,combinations

def small_posets(n):
 out=set();edges=list(combinations(range(n),2))
 for bits in range(1<<len(edges)):
  rel={e for i,e in enumerate(edges) if bits>>i&1}
  while True:
   nxt=rel|{(a,d) for a,b in rel for t,d in rel if b==t}
   if nxt==rel:break
   rel=nxt
  out.add(tuple(sum(1<<a for a,b in rel if b==j) for j in range(n)))
 return sorted(out)

def catalogue():
 out=set();l,z,a=2,1,4
 for lp,zp,ap in product(small_posets(l),small_posets(z),small_posets(a)):
  li=[m for m in range(1<<l) if all(not(m>>i&1) or not(lp[i]&~m) for i in range(l))]
  zi=[m for m in range(1<<z) if all(not(m>>i&1) or not(zp[i]&~m) for i in range(z))]
  for zd in product(li,repeat=z):
   if any(zd[i]&~zd[j] for j in range(z) for i in range(z) if zp[j]>>i&1):continue
   for ad in product(li,repeat=a):
    if any(ad[i]&~ad[j] for j in range(a) for i in range(a) if ap[j]>>i&1):continue
    for az in product(zi,repeat=a):
     if any(az[i]&~az[j] for j in range(a) for i in range(a) if ap[j]>>i&1):continue
     if any(zd[k]&~ad[j] for j in range(a) for k in range(z) if az[j]>>k&1):continue
     p=(0,12)+tuple(x<<2 for x in lp)+tuple((zp[i]<<4)|(zd[i]<<2) for i in range(z))+tuple(1|(ap[i]<<5)|(az[i]<<4)|(ad[i]<<2) for i in range(a));out.add(p)
 return out
