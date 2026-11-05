#!/usr/bin/env python3
from verify_extension import *

def induced(q,vs):return tuple(sum(((q[u]>>v)&1)<<j for j,v in enumerate(vs)) for u in vs)
def substitute(q,r,t):
 vs=[i for i in range(len(q)) if i!=r];N=len(vs);a=[0]*(N+len(t))
 for i,u in enumerate(vs):
  for j,v in enumerate(vs):
   if q[u]>>v&1:a[i]|=1<<j
  for j in range(len(t)):
   if q[u]>>r&1:a[i]|=1<<(N+j)
   else:a[N+j]|=1<<i
 for i,row in enumerate(t):a[N+i]|=row<<N
 validate(a);return tuple(a)
def triangles(q,r):return [(o,i) for o in range(len(q)) for i in range(len(q)) if q[r]>>o&1 and q[o]>>i&1 and q[i]>>r&1]
def cutoff(M,n,h):
 if M<h:return -1
 s=0;p=h
 while p*(n+2*(s+1))<=M:s+=1;p*=n+2*s
 return s

def run():
 contexts=virtual=physical=0
 inputs=[(0,),transitive(2),(2,4,1)]
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   q=tournament(n,code)
   for r in range(n):
    contexts+=1;ts=triangles(q,r);outside=induced(q,[v for v in range(n) if v!=r]);B=bin_polynomial(outside);W=[0]
    for o,i in ts:W=add(W,bin_polynomial(induced(q,[v for v in range(n) if v not in (r,o,i)])))
    d=value(B,1);w=value(W,1);require(d>=1 and w>=len(ts),'weighted context parameters')
    op=skeleton_operator(q,r)
    for par,k in product((0,1),range(5)):
     f=basis(k,par);low=add(mul(B,f),mul(mon(1),mul(W,delta(f))))
     cs=coordinates(add(apply_skeleton(op,f),scale(low,-1)),(par+n-1)%2)
     require(all(v>=0 for v in cs),'weighted positive residual');virtual+=1
    for t in inputs:
     h=h_subsets(t)[-1];b=value(bin_polynomial(t),2);got=h_subsets(substitute(q,r,t))[-1]
     require(got>=d*h+w*b,'physical weighted inequality')
     if ts:require(got>=(len(t)+2)*h and got>=h+2**len(t),'two active inequalities')
     physical+=1
 cutgraphs=cutpairs=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   t=tournament(n,code);words=[p for p in permutations(range(n)) if all(t[u]>>v&1 for u,v in zip(p,p[1:]))];images=set()
   for word in words:
    for cut in range(n+1):
     image=(word[:cut],word[cut:]);require(image not in images,'fixed-labelled-cut injectivity');images.add(image);cutpairs+=1
   b=value(bin_polynomial(t),2);require(len(images)==(n+1)*len(words) and b>=len(images) and b>=2**n,'bin lower bounds');cutgraphs+=1
 targets=0
 for n in range(1,11):
  for h in (1,3,5,13):
   for M in range(1,501):
    s=cutoff(M,n,h)
    if s<0:require(M<h,'impossible seed cutoff')
    else:
     p=h
     for j in range(1,s+1):p*=n+2*j
     require(p<=M and p*(n+2*(s+1))>M,'exact maximal active depth')
    targets+=1
 # Several genuinely nonuniform choices of contexts, but each whole-module at its own step.
 tri=(2,4,1);passive=join(tri,(0,));oddactive=join(tri,(0,));bank=[(tri,0),(passive,3),(oddactive,1),((2,0),0)]
 adaptive=0
 for seed in inputs:
  for codes in product(range(4),repeat=3):
   t=seed;initial=h_subsets(t)[-1];lower=initial;s=0
   for code in codes:
    q,r=bank[code];oldn=len(t);oldh=h_subsets(t)[-1];oldb=value(bin_polynomial(t),2)
    t=substitute(q,r,t);h=h_subsets(t)[-1]
    require(h>=oldh,'monotone context')
    if triangles(q,r):
     s+=1;lower*=len(seed)+2*s
     require(h>=lower and 2**oldn<=h-oldh,'adaptive count/order bound')
     d=h_subsets(induced(q,[v for v in range(len(q)) if v!=r]))[-1]
     w=sum(h_subsets(induced(q,[v for v in range(len(q)) if v not in (r,o,i)]))[-1] for o,i in triangles(q,r))
     require(w<=(h-d*oldh)//oldb,'exact context budget')
   adaptive+=1
 for n in range(1,11):require(h_subsets(wrapper(transitive(n)))[-1]==1+2**n,'padding must remain before active wrapper')
 p=[0,1];counts=[1]
 for _ in range(3):p=D(p);counts.append(value(p,1))
 require(counts==[1,3,15,123],'initial sharp depth controls')
 print(json.dumps({'status':'PASS','contexts':contexts,'weighted_virtual_columns':virtual,'physical_inequalities':physical,'cut_input_graphs':cutgraphs,'literal_labelled_cuts':cutpairs,'target_cutoffs':targets,'adaptive_words':adaptive,'initial_triangle_counts':counts},sort_keys=True,indent=2))
if __name__=='__main__':run()
