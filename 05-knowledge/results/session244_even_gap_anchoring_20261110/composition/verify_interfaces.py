"""Exact pure-gap higher interface verifier, independent subset/path enumeration.
No graph-class census or arbitrary-target coverage is asserted.
"""
import importlib.util,json,sys
from itertools import permutations,product
from pathlib import Path
s=importlib.util.spec_from_file_location('v',Path(__file__).with_name('verify.py'))
v=importlib.util.module_from_spec(s);s.loader.exec_module(v);p=v.p

def interface(a,cut):
 h,e=p.tables(a);_,start=p.tables(p.reverse(a));f=(1<<len(a))-1
 es=[sum(z[r] for r in range(len(a)) if cut>>r&1) for z in e]
 ss=[sum(z[r] for r in range(len(a)) if cut>>r&1) for z in start]
 return (h[f],sum(h[t]*h[f^t] for t in range(f+1)),
         sum(h[t]*ss[f^t] for t in range(f+1)),
         sum(es[t]*h[f^t] for t in range(f+1)),
         sum(es[t]*ss[f^t] for t in range(f+1)),ss[f],es[f])

def delete(a,r):
 ids=[i for i in range(len(a)) if i!=r]
 return [sum((a[i]>>j&1)<<k for k,j in enumerate(ids)) for i in ids],sum((a[r]>>j&1)<<k for k,j in enumerate(ids))

def response3(a,r):
 h,_=p.tables(a);_,start=p.tables(p.reverse(a));f=(1<<len(a))-1;out=0
 for mask in range(1,f+1):
  rem=f^mask;t=rem
  while True:
   out+=start[mask][r]*h[t]*h[rem^t]
   if t==0:break
   t=(t-1)&rem
 return out

def covers(a):
 n=len(a)
 for order in permutations(range(n)):
  for k in range(n+1):
   ps=(order[:k],order[k:])
   if all(all(a[x]>>y&1 for x,y in zip(path,path[1:])) for path in ps):yield ps

def hasarc(ps,x,y):return any((x,y) in tuple(zip(path,path[1:])) for path in ps)
def ends(ps,x):return bool(ps[0]) and ps[0][-1]==x

def mul(x,y):
 h,q,A,G,T,S,E=x;j,r,B,F,U,V,W=y
 return (h*j,q*r,A*r+h*B,F*q+j*G,A*F+T*j+h*U+E*V,S*j,h*W)
def close(z):return z[0]+2*z[2]-z[4]

def main():
 if '--fail' in sys.argv:p.need(False,'intentional failure')
 # K is vertex t=0 followed by a=1 -> b=2 -> c=3 -> a.
 K=[14,4,8,2]
 # Insert the nonempty r-start source segment as z=4; a->z is its only entry.
 Kz=[x|((1<<4) if i==1 else 0) for i,x in enumerate(K)]+[13]
 c=list(covers(Kz));added=[ps for ps in c if hasarc(ps,1,4)]
 p.need(len(added)==42,'two-path added skeleton coefficient')
 p.need(sum(ends(ps,1) for ps in added)==0 and sum(ends(ps,2) for ps in added)==12,'endpoint difference skeleton coefficient')
 small=[ps for ps in covers([10,4,1,6]) if hasarc(ps,0,3)]
 p.need(len(small)==14 and sum(ends(ps,1) for ps in small)==4,'bounded-error skeleton coefficients')
 removed=[]
 for j in range(2):
  q=[ps for ps in covers(K) if ps[j] and ps[j][0]==1]
  counts=[sum(ends(ps,x) for ps in q) for x in (1,2)]
  p.need(len(q)==3 and counts==([1,1] if j==0 else [0,0]),'removed-edge skeleton cancellation')
  removed.append(dict(path=j,total=len(q),endpoint_a_b=counts))
 literal=[];modules=[]
 sources=[(p.labelled(n,code),r) for n in range(1,4) for code in range(1<<(n*(n-1)//2)) for r in range(n)]
 # Additional nonisomorphic-looking labels exercise nonsymmetric root responses.
 sources += [(p.labelled(4,code),r) for code in (0,1,3,7,13,19,27,42,63) for r in range(4)]
 for a,r in sources:
  D,E,U,d=p.state(a,r);Q=p.pair(a)[1];V=p.state(p.reverse(a),r)[2];W=response3(a,r)
  states=[]
  for w in ('OJJJJOZ','OJOZJOZ'):
   graph,root=v.literal(a,r,w);deletion,cut=delete(graph,root);z=interface(deletion,cut);states.append(z)
   p.need(z[0]==3*(D+2*U) and z[1]==4*(24*Q+42*W-6*V),'module h/q formula')
   p.need(z[2]==18*D+44*U and z[6]==0,'module A/E formula')
   p.need(p.state(graph,root)[0]==z[0]+z[2]-z[4]-z[6],'one-root closure identity')
  p.need(tuple(y-x for x,y in zip(*states))==(0,0,0,-12*W,-2*U,0,0),'full seven-coordinate identity')
  raw=[interface(*delete(*v.literal(a,r,w))) for w in ('JJJJ','JOZJ')]
  p.need(tuple(y-x for x,y in zip(*raw))==(0,0,E,-4*W,0,0,-U),'raw singleton-cut interface')
  escape=[interface(*delete(*v.literal(a,r,w))) for w in ('JJJJOZ','JOZJOZ')]
  p.need(tuple(y-x for x,y in zip(*escape))==(0,0,4*E,-4*W,E-U,0,0),'bounded-error higher interface')
  p.need(escape[0][0]==3*D-E+3*U and escape[0][1]==4*(12*Q+14*W-6*V),'bounded-error common h/q')
  if len(a)<=2:modules.append((states,U,W))
  literal.append(dict(rows=a,root=r,U=U,W=W,V=V,deletion_h=states[0][0],deletion_q=states[0][1]))
 # Finite word checks include every word over all three letters, not only O/Z.
 iff_cases=0
 for n in range(9):
  for chars in product('OZJ',repeat=n):
   w=''.join(chars);delta=v.act((1,0,0,0),w)
   p.need((delta==(1,0,0,0))==('J' not in w),'pure channel iff no J')
   if 'J' in w:p.need(delta[3]>=1,'positive deletion after first J')
   iff_cases+=1
 # Every prescribed half-gap n>=2 has the E=1 bounded-response-error pair.
 escape_modules=[];escape_rows=[];short_rows=[]
 for n in range(2,18):
  N=2*n-3
  if N==1:a,r=[0],0
  else:
   letters=[];u=N
   while u>2:
    letters.append('O' if u%2==0 else 'Z');u=u//2 if u%2==0 else (u+1)//2
   a,r=v.literal([2,0],0,''.join(reversed(letters)))
  z=p.state(a,r);p.need(z[1]==1 and z[2]==N,'E-one response selector')
  raw=[v.act(z,w) for w in ('JJJJ','JOZJ')]
  p.need(tuple(y-x for x,y in zip(*raw))==(N+1,N,1,0),'raw response error one with endpoint retained')
  p.need(raw[1][0]<=12*N-1,'raw anchor bound')
  for w,out in zip(('JJJJ','JOZJ'),raw):
   b,t=v.literal(a,r,w);p.need(p.state(b,t)==out and p.strong(b),'literal strong raw pair')
  short_rows.append(dict(half_gap=(N+1)//2,lower=raw[0][0],upper=raw[1][0],endpoint_difference=N,response_difference=1))
  outs=[v.act(z,w) for w in ('JJJJOZ','JOZJOZ')]
  for w,out in zip(('JJJJOZ','JOZJOZ'),outs):
   b,t=v.literal(a,r,w);p.need(p.state(b,t)==out and p.strong(b),'literal strong bounded-error pair')
  p.need(tuple(y-x for x,y in zip(*outs))==(2*n,0,4,0),'bounded-error prescribed gap')
  p.need(outs[0][3]%3==2 and outs[1][0]<=76*n-117,'coprime deletion and linear anchor')
  for m in range(7):
   seen=set()
   for chars in product('OZ',repeat=m):
    w=''.join(chars);delta=tuple(y-x for x,y in zip(*(v.act(o,w) for o in outs)))
    digit=sum(2**i for i,c in enumerate(w) if c=='Z')
    p.need(delta==(2*n+4*digit,0,4*2**m,0),'exact bounded-error tail digits');seen.add(delta[0])
   p.need(seen==set(range(2*n,2*n+4*2**m,4)),'complete tail gap interval at step four')
  if n<=4:
   # Different tail lengths, with fixed words in each paired module.
   tail='OZ'[:n-2]
   ims=[interface(*delete(*v.literal(a,r,w+tail))) for w in ('JJJJOZ','JOZJOZ')]
   escape_modules.append((ims,len(tail)))
  escape_rows.append(dict(half_gap=n,response=N,lower=outs[0][0],upper=outs[1][0],deletion=outs[0][3]))
 def valuation2(x):return (x & -x).bit_length()-1
 address=[]
 for k in range(1,8):
  ms=[escape_modules[i%len(escape_modules)] for i in range(k)];coeff=[]
  for i,(zs,t) in enumerate(ms):
   c=4*2**t
   for j,(other,_) in enumerate(ms):
    if j<i:c*=other[0][0]
    if j>i:c*=other[0][1]
   coeff.append(c)
  vals=[valuation2(c) for c in coeff]
  p.need(all(a>b for a,b in zip(vals,vals[1:])),'strict dyadic address valuations')
  image=[]
  for bits in product((0,1),repeat=k):
   z=ms[0][0][bits[0]]
   for i in range(1,k):z=mul(z,ms[i][0][bits[i]])
   image.append(z[2])
  p.need(len(set(image))==2**k,'no complementary response cancellation')
  p.need(image==[image[0]+sum(b*c for b,c in zip(bits,coeff)) for bits in product((0,1),repeat=k)],'exact address weights')
  address.append(dict(modules=k,valuations=vals,distinct_responses=len(image)))
 # Independent multiplication recomputes every binary choice in fixed module banks.
 cubes=[]
 for k in range(1,7):
  ms=[modules[(2*i+1)%len(modules)] for i in range(k)]
  outs=[]
  for bits in product((0,1),repeat=k):
   z=ms[0][0][bits[0]]
   for i in range(1,k):z=mul(z,ms[i][0][bits[i]])
   outs.append((bits,close(z)))
  base=outs[0][1];steps=[]
  for i,(zs,U,W) in enumerate(ms):
   if i:
    left=ms[0][0][0]
    for j in range(1,i):left=mul(left,ms[j][0][0])
    hp,Ap=left[0],left[2]
   else:hp,Ap=1,0
   hs=1
   for j in range(i+1,k):hs*=ms[j][0][0][0]
   step=2*hs*(hp*U+6*Ap*W);steps.append(step)
   p.need(step%(2*3**(k-1))==0,'3-adic step divisor')
  for bits,out in outs:p.need(out==base+sum(b*t for b,t in zip(bits,steps)),'affine physical cube')
  cubes.append(dict(modules=k,counts=len(outs),base=base,steps=steps,guaranteed_divisor=2*3**(k-1)))
 print(json.dumps(dict(status='PASS',skeleton_added_total=42,skeleton_added_endpoint_b=12,skeleton_removed=removed,
   literal_interfaces=literal,literal_pairs=len(literal),word_iff_cases=iff_cases,cubes=cubes,bounded_error_gaps=escape_rows,dyadic_addresses=address,short_strong_pairs=short_rows,
   scope='Universal proved interface and continuation laws, finite exact corroboration; arbitrary-target coverage remains open.'),sort_keys=True,indent=2))
if __name__=='__main__':main()
