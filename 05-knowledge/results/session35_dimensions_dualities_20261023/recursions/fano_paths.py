import itertools,collections,json
V=range(7);QR={1,2,4};Lplus={frozenset((i+j)%7 for j in (0,1,3)) for i in V};Lminus={frozenset((i-j)%7 for j in (0,1,3)) for i in V}
def arc(i,j):return (j-i)%7 in QR
def unused(P,L):return tuple(sorted(tuple(sorted(l)) for l in L if not any({a,b}<=l for a,b in zip(P,P[1:]))))
def blue(P):return all(arc(P[i],P[j])==arc(P[6-j],P[6-i]) for i in V for j in V if i!=j)
def aff(P,a,b):return tuple((a*x+b)%7 for x in P)
paths=[P for P in itertools.permutations(V) if all(arc(a,b) for a,b in zip(P,P[1:]))]
autos=[p for p in itertools.permutations(V) if all(arc(i,j)==arc(p[i],p[j]) for i in V for j in V if i!=j)]
expected={tuple((a*x+b)%7 for x in V) for a in QR for b in V}
if set(autos)!=expected:raise RuntimeError('auto group differs')
for p in autos:
 for L in (Lplus,Lminus):
  if {frozenset(p[x] for x in l) for l in L}!=L:raise RuntimeError('design moved')
counts=collections.Counter();cert=[]
for P in paths:
 u=len(unused(P,Lplus));v=len(unused(P,Lminus));b=blue(P)
 rev=all(P[6-i]==(-P[i]+P[0]+P[-1])%7 for i in V)
 if b!=rev or b!=(u==v) or (not b and abs(u-v)!=1):raise RuntimeError('balance theorem failed')
 counts[u,v,b]+=1
 if P[:2]==(0,1):
  orb={aff(P,a,b) for a in QR for b in V}
  if len(orb)!=21:raise RuntimeError('nonfree orbit')
  cert.append({'path':P,'u+':u,'u-':v,'blue':b,'reversal_b':P[-1] if b else None,'unused+':unused(P,Lplus),'unused-':unused(P,Lminus)})
print('paths',len(paths),'autos',len(autos),'counts',sorted(counts.items()))
print(json.dumps(cert,indent=2))
# exact path-preserving cube and degree preservation for either Fano system
cube_checks=0
for P in paths:
 for L in (Lplus,Lminus):
  U=unused(P,L)
  for bits in range(1<<len(U)):
   flips={frozenset((a,b)) for k,l in enumerate(U) if bits>>k&1 for a,b in itertools.combinations(l,2)}
   def newarc(a,b):return arc(a,b)^(frozenset((a,b)) in flips)
   if not all(newarc(a,b) for a,b in zip(P,P[1:])):raise RuntimeError('path lost')
   if any(sum(newarc(a,b) for b in V if b!=a)!=3 for a in V):raise RuntimeError('regularity lost')
   cube_checks+=1
print('regular path-preserving cube states checked',cube_checks)
def ham_adj(out):
 dp=[[0]*7 for _ in range(128)]
 for i in V:dp[1<<i][i]=1
 for mask in range(1,128):
  for i in V:
   if not dp[mask][i]:continue
   for j in V:
    if not(mask>>j&1) and out[i]>>j&1:dp[mask|1<<j][j]+=dp[mask][i]
 return sum(dp[-1])
base=[sum(1<<j for j in V if i!=j and arc(i,j)) for i in V];cache={};inc=collections.Counter();cube_H=collections.Counter()
for P in paths:
 deg=sum(len(unused(P,L)) for L in (Lplus,Lminus))
 if blue(P)!=(deg in (2,4)) or (not blue(P) and deg!=5):raise RuntimeError('local degree failure')
 for L in (Lplus,Lminus):
  U=unused(P,L)
  for l in U:inc[l]+=1
  for bits in range(1<<len(U)):
   out=base[:]
   for k,l in enumerate(U):
    if bits>>k&1:
     for a,b in itertools.combinations(l,2):out[a]^=1<<b;out[b]^=1<<a
   key=tuple(out)
   if key not in cache:cache[key]=ham_adj(out)
   h=cache[key];cube_H[bool(bits),h]+=1
   if h!=(171 if bits else 189):raise RuntimeError('cube nonorigin returned Paley')
print('triangle edge-avoidance incidences',sorted(inc.items()),'total',sum(inc.values()))
print('cube H census',sorted(cube_H.items()),'distinct presentations',len(cache))
host=[98,84,41,19,37,74,28]
def ha(i,j):return bool(host[i]>>j&1)
tri=[l for l in itertools.combinations(V,3) if all(sum(ha(i,j) for j in l if i!=j)==1 for i in l)]
for P,expect in [((0,1,2,3,4,5,6),True),((0,1,4,2,5,6,3),False)]:
 if not all(ha(a,b) for a,b in zip(P,P[1:])):raise RuntimeError('hostile path invalid')
 b=all(ha(P[i],P[j])==ha(P[6-j],P[6-i]) for i in V for j in V if i!=j)
 deg=sum(not any(set(e)<=set(l) for e in zip(P,P[1:])) for l in tri)
 if b!=expect or deg!=5:raise RuntimeError('hostile control failed')
 print('H171 same-degree hostile',P,'blue',b,'degree',deg,'H',ham_adj(host))
for L in (Lplus,Lminus):
 lines=sorted(tuple(sorted(l)) for l in L);hc=collections.Counter();paleywords=[]
 for bits in range(128):
  out=base[:]
  for k,l in enumerate(lines):
   if bits>>k&1:
    for a,b in itertools.combinations(l,2):out[a]^=1<<b;out[b]^=1<<a
  h=ham_adj(out);hc[h]+=1
  if h==189:paleywords.append(bits)
 if hc!={189:16,171:112}:raise RuntimeError('fixedplane census')
 if any(a^b not in paleywords for a in paleywords for b in paleywords):raise RuntimeError('Paley relative code not linear')
 if min(b.bit_count() for b in paleywords if b)!=3:raise RuntimeError('Paley distance')
 for b in paleywords:
  if b.bit_count()==3 and len(set.intersection(*(set(lines[k]) for k in range(7) if b>>k&1)))!=1:raise RuntimeError('weightthree not pencil')
 print('fixedplane code certified',dict(sorted(hc.items())),'weights',sorted(collections.Counter(b.bit_count() for b in paleywords).items()))
p=11;qr={i*i%p for i in range(1,p)};P=(6,4,5,3,7,0,1,10,8,9,2)
def a11(i,j):return (j-i)%p in qr
if not all(a11(a,b) for a,b in zip(P,P[1:])):raise RuntimeError('P11 invalid path')
bl=all(a11(P[i],P[j])==a11(P[10-j],P[10-i]) for i in range(p) for j in range(p) if i!=j)
cs=[c for c in itertools.combinations(range(p),3) if all(sum(a11(i,j) for j in c if i!=j)==1 for i in c)]
deg=sum(not any(set(e)<=set(c) for e in zip(P,P[1:])) for c in cs)
if bl or deg!=32:raise RuntimeError('P11 parity converse hostile failed')
print('Paley11 even-degree nonblue hostile',P,'degree',deg,'blue',bl)
