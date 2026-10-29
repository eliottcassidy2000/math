#!/usr/bin/env python3
"""Positive four-port reconstruction and Construction-A theta/root certificates."""
import contextlib,io,runpy
from pathlib import Path
from itertools import combinations
from collections import Counter
with contextlib.redirect_stdout(io.StringIO()):
 d=runpy.run_path(str(Path(__file__).with_name('check_cyclic_frobenius.py')))
CHECKS=0
def check(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(msg)
def code(g,n):
 words=[0]
 for v in d['code_basis'](g,n):words += [w^v for w in words]
 return words+[w^((1<<(n+1))-1) for w in words]
SIX=[73235,83903,92135,102641,118733,129989]
affine=[]
for g in SIX:
 words=code(g,31);octads=[frozenset(i for i in range(32) if w>>i&1) for w in words if w.bit_count()==8]
 counts=Counter(s for block in octads for s in combinations(sorted(block),4))
 planes=[frozenset(s) for s,k in counts.items() if k==7]
 check(len(planes)==1240,'1240 affine candidate planes')
 add=[[None]*32 for _ in range(32)];zero=31
 for a in range(32):add[a][a]=zero;add[a][zero]=add[zero][a]=a
 for P in planes:
  if zero in P:
   a,b,c=sorted(P-{zero})
   for x,y,z in [(a,b,c),(a,c,b),(b,c,a)]:
    check(add[x][y] is None,'unique fourth point')
    add[x][y]=add[y][x]=z
 for a in range(32):
  for b in range(32):
   check(add[a][b] is not None,'addition defined')
   for c in range(32):check(add[add[a][b]][c]==add[a][add[b][c]],'associative four-port reconstruction')
 # Obtain an explicit five-bit chart from this reconstructed group.
 chart={zero:0};basis=[]
 while len(chart)<32:
  a=min(set(range(32))-chart.keys());bit=1<<len(basis);basis.append(a)
  chart.update({add[x][a]:v^bit for x,v in list(chart.items())})
 check(len(basis)==5 and set(chart.values())==set(range(32)),'F2^5 chart')
 for a in range(32):
  for b in range(32):check(chart[add[a][b]]==chart[a]^chart[b],'chart homomorphism')
 for B in octads:
  a=min(B);W={add[x][a] for x in B}
  check(zero in W and len(W)==8,'eight-element affine span')
  for x in W:
   for y in W:check(add[x][y] in W,'octad affine3flat')
 # Rotation is linear and has one 31-cycle on nonzero vectors.
 rot=list(range(1,31))+[0,31]
 for a in range(32):
  for b in range(32):check(rot[add[a][b]]==add[rot[a]][rot[b]],'Singer rotation linear')
 # Its minimal polynomial is recovered independently from the addition table.
 v=0;seq=[]
 for i in range(6):seq.append(chart[v]);v=rot[v]
 masks=[]
 for mask in range(32):
  z=0
  for i in range(5):
   if mask>>i&1:z^=seq[i]
  if z==seq[5]:masks.append(mask|32)
 check(len(masks)==1,'unique degree5 rotation polynomial')
 affine.append((g,tuple(basis),masks[0]))
check(len({p for g,b,p in affine})==6,'six distinct primitive degree5 polynomials')
# Formal Construction-A series use z^(sum integer coordinate squares).
LIMIT=32
def poly_mul(a,b):
 out=[0]*(LIMIT+1)
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b[:LIMIT+1-i]):
    if y:out[i+j]+=x*y
 return out
def poly_power(a,n):
 z=[1]+[0]*LIMIT
 while n:
  if n&1:z=poly_mul(z,a)
  a=poly_mul(a,a);n//=2
 return z
E=[0]*(LIMIT+1);O=[0]*(LIMIT+1)
for a in range(-5,6):
 if a*a<=LIMIT:(E if a%2==0 else O)[a*a]+=1
def theta(n,g):
 weights=Counter(w.bit_count() for w in code(g,n));out=[0]*(LIMIT+1)
 for w,m in weights.items():
  term=poly_mul(poly_power(E,n+1-w),poly_power(O,w))
  out=[x+m*y for x,y in zip(out,term)]
 return out,weights
thetas=[]
for n,g in [(7,29),(23,7973),(31,75035),(31,92135)]:
 series,weights=theta(n,g);dim=n+1
 check(all(w%4==0 for w in weights),'doubly even')
 check(sum(weights.values())==2**(dim//2),'selfdual size')
 roots=2*dim+16*weights.get(4,0)
 check(series[4]==roots,'root formula')
 check(all(series[j]==0 for j in range(LIMIT+1) if j%4),'even lattice norms')
 thetas.append((n,g,roots,tuple(series[4*i] for i in range(9))))
check(thetas[0][2]==240,'Hamming E8 roots')
check(thetas[1][2]==48,'Golay A1^24 roots not Leech')
check(thetas[2][2]==thetas[3][2]==64,'rank32 A1^32 roots')
check(thetas[2][3]==thetas[3][3],'same theta through norm16')
check(thetas[0][3][2]==2160,'E8 norm4 calibration')
print('CHECKS',CHECKS)
print('Affine reconstruction: generator, chosen basis labels, rotation minimal-polynomial bitmask')
for row in affine:print(row)
print('Construction-A: n, generator, root count, theta coefficients at lattice norms0,2,...,16')
for row in thetas:print(row)
print('Six-code sector: four-port data reconstructs F2^5; octads are exactly all affine3-flats.')
print('Rank32 pair: equal full theta by common weight enumerator; nonisometric by canonical root axes and inequivalent codes.')
