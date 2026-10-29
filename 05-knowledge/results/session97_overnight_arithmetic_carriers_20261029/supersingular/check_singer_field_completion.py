#!/usr/bin/env python3
"""Reconstruct six F32 fields from four-port addition plus a named Singer clock."""
from collections import Counter
from itertools import combinations
SIX={73235:41,83903:55,92135:47,102641:37,118733:61,129989:59}

def require(ok,label):
 if not ok:raise RuntimeError(label)

def codewords(g):
 words=[0]
 for i in range(15):words += [w^(g<<i) for w in words]
 return words+[w^((1<<32)-1) for w in words]

def reconstruct_addition(g):
 octads=[tuple(i for i in range(32) if w>>i&1) for w in codewords(g) if w.bit_count()==8]
 counts=Counter(s for B in octads for s in combinations(B,4))
 add=[[None]*32 for _ in range(32)];zero=31
 for a in range(32):add[a][a]=zero;add[a][zero]=add[zero][a]=a
 for S,k in counts.items():
  if k==7 and zero in S:
   a,b,c=[x for x in S if x!=zero]
   for x,y,z in [(a,b,c),(a,c,b),(b,c,a)]:
    require(add[x][y] is None,'unique four-port addition')
    add[x][y]=add[y][x]=z
 require(all(x is not None for row in add for x in row),'complete addition')
 for a in range(32):
  for b in range(32):
   for c in range(32):require(add[add[a][b]][c]==add[a][add[b][c]],'associative addition')
 return add

def multiplication(clock,one,zero):
 orbit=[];x=one
 while x not in orbit:
  orbit.append(x);x=clock[x]
 require(x==one and len(orbit)==31 and zero not in orbit,'Singer orbit guard')
 log={x:i for i,x in enumerate(orbit)}
 return [[zero if zero in (a,b) else orbit[(log[a]+log[b])%31] for b in range(32)] for a in range(32)]

def binary_product(a,b,f):
 z=0
 while b:
  if b&1:z^=a
  a<<=1;b>>=1
  if a&32:a^=f
 return z

clock=list(range(1,31))+[0,31];zero=31;one=0
rows=[];first_add=None;addition_tables=set()
for g,expected in SIX.items():
 add=reconstruct_addition(g)
 addition_tables.add(tuple(map(tuple,add)))
 if first_add is None:first_add=add
 for a in range(32):
  for b in range(32):require(clock[add[a][b]]==add[clock[a]][clock[b]],'linear clock guard')
 mult=multiplication(clock,one,zero)
 for a in range(32):
  require(mult[a][one]==a and mult[a][zero]==zero,'identity and zero')
  for b in range(32):
   require(mult[a][b]==mult[b][a],'commutative multiplication')
   for c in range(32):
    require(mult[a][add[b][c]]==add[mult[a][b]][mult[a][c]],'distributivity')
    require(mult[mult[a][b]][c]==mult[a][mult[b][c]],'associative multiplication')
  if a!=zero:require(sum(mult[a][b]==one for b in range(32))==1,'unique inverse')
 # The named phase e and clock determine the power basis e,alpha,...,alpha^4.
 chart={zero:0}
 for i in range(5):chart.update({add[x][i]:bits^(1<<i) for x,bits in list(chart.items())})
 require(len(chart)==32 and set(chart.values())==set(range(32)),'primitive power basis')
 minimal=32|chart[5]
 require(minimal==expected,'retained primitive polynomial')
 for a in range(32):
  for b in range(32):
   require(chart[add[a][b]]==chart[a]^chart[b],'additive chart')
   require(chart[mult[a][b]]==binary_product(chart[a],chart[b],minimal),'polynomial field multiplication')
 # Field Frobenius is exponent doubling, not the source's one-step clock.
 frobenius=[mult[a][a] for a in range(32)]
 require(frobenius==list((2*a)%31 for a in range(31))+[zero],'F32 Frobenius clock distinction')
 for a in range(32):
  require(frobenius[clock[a]]==clock[clock[frobenius[a]]],'Frobenius conjugates clock to its square')
  for b in range(32):require(frobenius[add[a][b]]==add[frobenius[a]][frobenius[b]],'additive Frobenius')
 # The same clock with another identity must change the operation.
 alternate=multiplication(clock,1,zero)
 require(alternate!=mult and all(alternate[a][1]==a for a in range(32)),'phase/identity guard')
 rows.append((g,minimal,chart[clock[one]],'all field axioms and polynomial products verified'))

require(len(addition_tables)==6,'six distinct source additions with the same multiplicative clock')

# A same-order clock with the linearity guard deliberately destroyed.
perm=list(range(32));perm[1],perm[2]=2,1
badclock=[perm[clock[perm[x]]] for x in range(32)]
badmult=multiplication(badclock,one,zero)
witness=next((a,b,c) for a in range(32) for b in range(32) for c in range(32)
 if badmult[a][first_add[b][c]]!=first_add[badmult[a][b]][badmult[a][c]])
a,b,c=witness
require(any(badclock[first_add[x][y]]!=first_add[badclock[x]][badclock[y]] for x in range(32) for y in range(32)),'nonlinear hostile clock')
# Same fixed XOR addition and identity1, but two distinct field multiplications.
# Degree-general polynomial multiplication for the separate F8 hostile.
def gf_multiply(a,b,f,m):
 z=0
 while b:
  if b&1:z^=a
  a<<=1;b>>=1
  if a&(1<<m):a^=f
 return z
for f in (11,13):
 require(f&1 and f.bit_count()%2==1,'irreducible F2 cubic: neither zero nor one is a root')
 for a in range(8):
  for b in range(8):
   for c in range(8):require(gf_multiply(a,b^c,f,3)==gf_multiply(a,b,f,3)^gf_multiply(a,c,f,3),'F8 distributivity control')
require(gf_multiply(2,4,11,3)==3 and gf_multiply(2,4,13,3)==5,'same addition different F8 fields')
print('Six four-port + Singer-clock + identity reconstructions:')
for row in rows:print(row)
print('Nonlinear order31 clock distributivity hostile:',witness,
 badmult[a][first_add[b][c]],first_add[badmult[a][b]][badmult[a][c]])
print('F8 without named clock: XOR addition and identity1 agree, but 2*4 equals3 modulo x^3+x+1 and5 modulo x^3+x^2+1.')
print('PASS: named zero, additive carrier, linear Singer clock, and multiplicative identity suffice; order and cycle size alone do not.')
