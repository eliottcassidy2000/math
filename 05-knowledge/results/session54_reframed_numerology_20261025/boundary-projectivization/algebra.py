#!/usr/bin/env python3
import json,collections as co,itertools as it,pathlib
D=pathlib.Path(__file__).resolve().parent
def ck(p,m):
 if not p:raise RuntimeError(m)
def mul(a,b):
 out=0
 while b:
  if b&1:out^=a
  b>>=1;a<<=1
  if a&64:a^=67 # X^6+X+1
 return out
def power(a,n):
 out=1
 while n:
  if n&1:out=mul(out,a)
  a=mul(a,a);n>>=1
 return out
powers=[power(2,e) for e in range(63)];ck(len(set(powers))==63 and power(2,63)==1,'primitive field element')
log={x:e for e,x in enumerate(powers)};F4={x for x in range(64) if power(x,4)==x};ck(len(F4)==4,'F4')
def inv(x):return power(x,62) if x else 0
packet=json.loads((D/'result.json').read_text());fibers={}
for row in packet['map']:fibers.setdefault(row['boundary'],[]).append(powers[row['exponent']])
for fiber in fibers.values():
 a=fiber[0];ck(set(fiber)=={mul(a,c) for c in F4 if c},'F4 fiber')
derivatives=[]
for a in range(1,64):
 dist=co.Counter(inv(x^a)^inv(x) for x in range(64));hist=co.Counter(dist.get(b,0) for b in range(64))
 ck(dict(hist)=={0:33,2:30,4:1},'derivative spectrum')
 exceptional=[x for x in range(64) if inv(x^a)^inv(x)==inv(a)]
 ck(set(exceptional)=={mul(a,c) for c in F4},'quartet')
 derivatives.append({'a':a,'exceptional_b':inv(a),'quartet':exceptional,'histogram':dict(hist)})
lines=set()
for a,b in it.combinations(range(1,64),2):
 if log[a]%21==log[b]%21:continue
 span={mul(a,c)^mul(b,d) for c in F4 for d in F4};ck(len(span)==16,'plane span')
 points=frozenset(log[x]%21 for x in span if x);ck(len(points)==5,'projective line size');lines.add(points)
ck(len(lines)==21,'PG line count')
preserved=[L for L in lines if frozenset(-e%21 for e in L) in lines]
ck(not preserved,'inversion unexpectedly preserves line')
L=min(sorted(L) for L in lines)
print(json.dumps({'status':'PASS','field_polynomial_mask':67,'primitive_element':2,'F4':sorted(F4),'projective_lines':[sorted(L) for L in sorted(lines,key=lambda L:sorted(L))],'inversion_line_images_that_are_lines':len(preserved),'line_noncollineation_witness':{'line':L,'inverse_image':sorted(-e%21 for e in L)},'derivatives':derivatives},sort_keys=True,indent=2))
