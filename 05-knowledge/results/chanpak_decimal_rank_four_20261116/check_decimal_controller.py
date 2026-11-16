import itertools,json,sys
from pathlib import Path
from collections import defaultdict
import numpy as np

HERE=Path(__file__).resolve().parent
if len(sys.argv)>1:
 paths=[Path(sys.argv[1])]
else:
 paths=[HERE/'return_automaton.json',HERE.parent/'eight_braid_root'/'return_automaton.json',HERE.parent/'chanpak_eight_braid_20261116'/'return_automaton.json']
filename=next((p for p in paths if p.exists()),None)
if filename is None:raise SystemExit('Generate return_automaton.json with lattice.py in the previous packet, then give its path as argument.')
D=json.loads(filename.read_text())
f=np.array(D['readout'],dtype=np.int64)
v=np.array(D['start'],dtype=np.int64)
maps={tuple(w):np.array(m,dtype=np.int64).reshape(4,4) for w,m in zip(D['words'],D['matrices'])}
I=np.eye(4,dtype=np.int64)
def need(x):
 if not x:raise AssertionError('exact certificate failure')
def power(A,e,q):
 X=I.copy()
 while e:
  if e&1:X=A@X%q
  A=A@A%q;e>>=1
 return X
def rank(A,p):
 a=np.array(A,dtype=np.int64).copy()%p;k=0
 for i in range(4):
  idx=np.flatnonzero(a[k:,i])
  if not len(idx):continue
  j=k+idx[0];a[[k,j]]=a[[j,k]]
  a[k]*=pow(int(a[k,i]),-1,p);a[k]%=p
  for j in range(k+1,len(a)):
   a[j]=(a[j]-a[j,i]*a[k])%p
  k+=1
 return k
L={tuple(M.ravel()):M for M in maps.values()}
G=np.array([M%800 for M in L.values() if rank(M,5)==4])
def points(q):
 return np.array([t for t in itertools.product(range(q),repeat=4) if next((x for x in t if x%5),None)==1],dtype=np.int64)
need(len(maps)==168 and len(L)==136 and len(G)==111)
def endpoint_strong(w):
 # Contract the known strong old prefix to 0. Every new vertex has
 # an incoming old arc (there is an old vertex more than 3 sites away).
 w=tuple(w)+(0,3,5,2,4)
 A=[set(range(1,9))]+[set() for _ in range(8)]
 for j,b in enumerate(w,1):
  if b>>(j-1):A[j].add(0)
  for i in range(1,j):
   if j-i<=3 and b>>(j-i-1)&1:A[j].add(i)
   else:A[i].add(j)
 seen={0}
 for _ in range(8):seen|={j for j in range(1,9) if A[j]&seen}
 return len(seen)==9
need(all(endpoint_strong(w) for w in maps))
print('PASS: every return block preserves strongness at its endpoint.')

need(list(f)==[32,14,10,27] and list(v)==[-43,-5,4,53])
A=maps[4,3,5]%800
l=np.array([28,24,24,19]);r=np.array([24,10,14,27])
need(np.array_equal(power(A,620,25),I))
P=power(A,4,32)
need(np.array_equal(P,np.outer(l,r)%32))
need(f@l%32==r@l%32==1 and P[3,3]%2==1 and f@P@v%2==1)
need(np.all((f@P@G@l)%2==1))
p5=points(5);need(len(p5)==156)
for u in p5:need(set(f@G@u%5)==set(range(5)))
entries=[((2,1,6),24,(4,2,4,2)),((2,1,7),24,(3,1,1,0)),((3,1,7),20,(0,0,4,4)),((2,3,3),24,(0,0,4,2))]
rows=[]
for w,o,t in entries:
 C=power(maps[w]%25,o,25)
 need(np.array_equal(C%5,I)); row=(f@(((C-I)%25)//5))%5
 need(tuple(row)==t);rows.append(row)
need(rank(rows,5)==4)
for w,o,a in [((3,3,5),15,5),((4,1,6),31,7)]:
 B=maps[w]%800
 need(np.array_equal(power(B,o,5),I) and o%2)
 need((f@P@B@l)%8==a)
print('PASS: arbitrary k,d: rank-one reset, all 156 projective mod-5 starts, four lifting rows, two dyadic unit generators.')

P=I.copy()
for w in [(2,1,7),(4,3,7),(4,5,5)]:P=maps[w]%800@P%800
l0=np.array([7,2,21,9]);r0=np.array([16,30,24,27])
need(rank(P,5)==4 and np.array_equal(P%32,np.outer(l0,r0)%32))
need(f@l0%32==1 and f@P@v%2==1)
signatures=defaultdict(set)
for A in G:
 for B in G:
  q=f@P@B@A%800;g=int(q@l0%32)
  need(g%2==1)
  signatures[g].add(tuple(q%25))
need(len(signatures)==16 and sum(map(len,signatures.values()))==12308)
p25=points(25);need(len(p25)==19500)
for g in range(1,32,2):
 check=np.zeros(len(p25),dtype=np.int64)
 for row in signatures[g]:check|=1<<(p25@row%25)
 need(np.all(check==(1<<25)-1))
print('PASS: mod-800 control in 40 letters, all 16 × 19500 projective (multiplier, start) cases.')
