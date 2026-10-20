from math import gcd
from hashlib import sha256
import json
checks=0;edges=0;cases=0;digest=sha256()
def req(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def B(w):
 b=0;M=1
 for x in w:
  if x:b=3*b+M
  M*=2
 return b
def rot(w,r):return w[r:]+w[:r]
for q in range(4,41):
 for p in range(1,q):
  if gcd(q,p)>1:continue
  c=[(i+1)*p//q-i*p//q for i in range(q)];D=2**q-3**p;pinv=pow(p,-1,q)
  for i in range(q):
   if c[i]==c[(i+1)%q]:continue
   w=c[:];w[i],w[(i+1)%q]=w[(i+1)%q],w[i]
   eps=1 if c[i]==1 else -1;star=(i+1)%q;E={};Bs={}
   g=gcd(abs(D),B(w));good=0;upper=abs(D)
   for r in range(q):
    cr=rot(c,r);wr=rot(w,r);Bs[r]=B(wr)
    if r==star:continue
    j=(i-r)%q;s=sum(cr[j+2:]);E[r]=eps*2**j*3**s
    req(j<=q-2,'linear swap location')
    req(B(wr)==B(cr)+E[r],'single swap formula')
   for t in range(q-1):
    r=t*pinv%q;rr=(t+1)*pinv%q
    k=((q-1-t)*pinv%q)-1;h=p-1-((k+1)*p//q);M=2**k*3**h
    req(B(rot(c,r))-B(rot(c,rr))==M,'mechanical ladder formula')
    if r==star or rr==star:continue
    value=M+E[r]-E[rr]
    req(value==Bs[r]-Bs[rr] and value!=0,'trinomial difference')
    req(value%g==0,'content divides trinomial')
    upper=gcd(upper,abs(value));good+=1;edges+=1
   req(good>=q-3,'exceptional cut edge count')
   req(upper%g==0,'upper bound in divisibility')
   digest.update(f'{q}:{p}:{i}:{g}:{upper}:{good}\n'.encode());cases+=1
# Separate cheap q<=80 hostile probe, counting swap occurrences, not necklaces.
probe_cases=nonunit=0;first=None;maximum=None;integer_positive_nontrivial=[]
for q in range(3,81):
 for p in range(1,q):
  if gcd(q,p)>1:continue
  c=[(i+1)*p//q-i*p//q for i in range(q)];D=2**q-3**p
  for i in range(q):
   j=(i+1)%q
   if c[i]==c[j]:continue
   w=c[:];w[i],w[j]=w[j],w[i];bw=B(w);g=gcd(abs(D),bw);probe_cases+=1
   if g>1:
    nonunit+=1;row={'q':q,'p':p,'swap_start':i,'D':D,'B':bw,'g':g,'word':''.join(map(str,w))}
    if first is None:first=row
    if maximum is None or g>maximum['g']:maximum=row
   if D>1 and g==D:integer_positive_nontrivial.append((q,p,i))
req(first=={'q':7,'p':2,'swap_start':3,'D':119,'B':112,'g':7,'word':'0000101'},'minimal hostile')
print(json.dumps({'status':'PASS','exact_trinomial_gate_universe':'all swap occurrences of coprime lower mechanical words, 4<=q<=40','gate_cases':cases,'gate_edges':edges,'explicit_checks':checks,'gate_digest_sha256':digest.hexdigest(),'hostile_probe_universe':'all swap occurrences of coprime lower mechanical words, 3<=q<=80; repetitions of output necklaces retained','probe_cases':probe_cases,'nonunit_occurrences':nonunit,'first_nonunit':first,'maximum_content':maximum,'positive_integer_candidates_D_gt_1':integer_positive_nontrivial,'boundary':'This finite absence is not a proof of the all-q proper-divisor target.'},indent=2,sort_keys=True))
