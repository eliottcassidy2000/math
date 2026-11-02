"""Independent local-universality checks; imports own graph kernel only."""
from kernel import *
from itertools import product
import json,hashlib,pathlib

actual=0;strong_count=0
for size in range(8):
 for word in product((0,1),repeat=size):
  a=[sum(1<<j for j in range(i+1,size)) for i in range(size)]
  a=extend(tuple(a),sum(b<<i for i,b in enumerate(word)))
  c=t=ell=0
  for bit in word:
   if bit:c,t,ell=c,2*t+ell,1
   else:c,t,ell=c+t+ell,2*t+ell,0
  require(h_recursive(a)==2*c+1,('actual word',word))
  if word and word[0] and not word[-1]:
   require(len(scc(a))==1,('strong word',word));strong_count+=1
  actual+=1

full_moduli=[];state_count=0
for m in range(1,81):
 seed=(0,1%m);seen={seed};queue=[seed]
 for d,t in queue:
  for p in ((d,2*t%m),((d-t)%m,(4*t+1)%m)):
   if p not in seen:seen.add(p);queue.append(p)
 require({(d+t)%m for d,t in seen}==set(range(m)),('positive block orbit',m))
 full_moduli.append(m);state_count+=len(seen)

commutator=0
for n in range(3,64,2):
 i2=pow(2,-1,n);i4=pow(4,-1,n);i8=pow(8,-1,n)
 def C(p):d,t=p;return((d-t)%n,(t+i4)%n)
 def D(p):d,t=p;return((d-i2*t)%n,(t+i2)%n)
 def Ci(p):d,t=p;return((d+t-i4)%n,(t-i4)%n)
 def Di(p):d,t=p;return((d+i2*t-i4)%n,(t-i2)%n)
 for p in product(range(n),repeat=2):
  q=D(C(Di(Ci(p))))
  require(q==((p[0]+3*i8)%n,p[1]),('commutator',n,p))
  for fn,inv in ((C,Ci),(D,Di)):
   q=p
   for _ in range(n-1):q=fn(q)
   require(q==inv(p),('positive inverse',n,p))
  commutator+=1

for a in range(10):
 m=2**a;start=max(1,(a+1)//2)
 require({((3*r+1-4**r)//9)%m for r in range(start,start+m)}==set(range(m)),('dyadic start',a))

j={'status':'PASS','actual_binary_words_through_length7':actual,'strong_words_directly_checked':strong_count,'positive_block_orbit_all_moduli':[1,80],'positive_block_states_checked':state_count,'odd_commutator_and_positive_inverse_states':commutator,'odd_moduli_for_commutator':[3,63],'dyadic_initialization_exponents':[0,9]}
j['script_sha256']=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
print(json.dumps(j,indent=2,sort_keys=True))
