import importlib.util,json
from pathlib import Path
from itertools import product
from math import gcd
s=importlib.util.spec_from_file_location('verify',Path(__file__).with_name('verify.py'));v=importlib.util.module_from_spec(s);s.loader.exec_module(v);p=v.p

def response_word(n):
 w=[]
 while n>2:
  if n%2:w.append('Z');n=(n+1)//2
  else:w.append('O');n//=2
 return ''.join(reversed(w))

def efficient(N):
 a,r=([0],0) if N==1 else v.literal([2,0],0,response_word(N))
 g,eps,u,d=p.state(a,r);p.need(u==N and d==1,'response source')
 k=len(a);a,r=v.literal(a,r,'OZ');r=k
 z=p.state(a,r)
 p.need(z==(g+2*N,N,2*g+2*N-eps,g+N-eps),'efficient endpoint full state')
 return a,r,z,g,eps

def main():
 cheap=[]
 for N in range(1,33):
  a,r,z,g,eps=efficient(N);counts=[]
  for w in ('JJJJZ','JOZJZ'):
   b,s=v.literal(a,r,w);h=p.tables(b)[0][-1]
   p.need(h==v.act(z,w)[0] and p.strong(b),'efficient compiled physical count')
   counts.append(h)
  p.need(counts==[31*g+37*N-11*eps,31*g+39*N-11*eps] and counts[1]<=70*N,'efficient bound')
  p.need(len(b)<=(8 if N==1 else (N-1).bit_length()+8),'efficient size')
  cheap.append(dict(N=N,g=g,epsilon=eps,counts=counts,order=len(b)))
 mixed=[];lifts=0
 for n in range(7):
  for x in product('OZ',repeat=n):
   tail=''.join(x)
   if tail and tail[-1]!='Z':continue
   a,r=v.literal([2,0],0,'Z'+tail);r=2;z=p.state(a,r)
   p.need(z[1]==1 and p.strong(a),'mixed tail endpoint one')
   c,s=v.literal(a,r,'J');near=p.state(c,s)
   p.need(near==(z[0]+z[2]-1,z[0]-1,z[0]+z[2],z[0]) and p.strong(c),'near response full state')
   p.need(gcd(near[0],near[2])==1,'primitive response')
   mixed.append(dict(tail=tail,source=z,near=near))
   if n<=3:
    for b in (p.tt(1),p.tt(2),p.endpoint(2)):
     bh,bq=p.pair(b);ab=p.join(c,b);out=p.graft(ab,s,[0]);full=p.state(out,len(out)-1)
     expecth=(bh+bq)*near[0]+bq
     p.need(full==(expecth,bh*near[0],expecth,bh*near[0]) and p.strong(out),'affine lift full state')
     if bh==1:p.need(gcd(full[0],full[1])==1,'primitive saturated lift')
     lifts+=1
 special=[]
 for k in range(1,8):
  a,r=v.literal([2,0],0,'Z'+'O'*k+'Z');a,r=v.literal(a,2,'J');z=p.state(a,r)
  h=4*3**k+9*2**k+1
  p.need(z[0]==h and z[2]==h+1,'special family domain k positive')
  special.append(dict(k=k,H=h,R=h+1))
 print(json.dumps(dict(status='PASS',efficient_compiled=cheap,mixed_endpoint_one_sources=mixed,affine_lift_cases=lifts,special_near_response=special,scope='Exact families and conditional target descent, not arbitrary scalar or sidecar coverage.'),sort_keys=True,indent=2))
if __name__=='__main__':main()
