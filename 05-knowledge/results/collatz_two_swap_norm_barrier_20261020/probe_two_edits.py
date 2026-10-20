#!/usr/bin/env python3
"""Finite exact two-edit census; floating norm sizes are exploratory only."""
import argparse,cmath,json,math
from collections import Counter

ap=argparse.ArgumentParser();ap.add_argument('--max-p',type=int,default=60);ap.add_argument('--norm-max-p',type=int,default=24);args=ap.parse_args()
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise ArithmeticError(msg)
def carry(mask,q):
 B=0
 for i in range(q):
  bit=(mask>>i)&1;B=(3 if bit else 1)*B+(bit<<i)
 return B
def rotations(mask,q):
 full=(1<<q)-1
 return {((mask>>j)|(mask<<(q-j)))&full for j in range(q)}
def edits(mask,q):
 for i in range(q):
  j=(i+1)%q
  if ((mask>>i)^(mask>>j))&1:yield i,mask^(1<<i)^(1<<j)
def fmt(mask,q):return ''.join(str((mask>>i)&1) for i in range(q))
def move_labels(labels,i,q,disp):
 j=(i+1)%q
 labels=dict(labels);disp=dict(disp)
 if i in labels: label=labels.pop(i);labels[j]=label;epsilon=1
 else:label=labels.pop(j);labels[i]=label;epsilon=-1
 disp[label]=disp.get(label,0)+epsilon
 return labels,{a:d for a,d in disp.items() if d},label,epsilon

def canonical_coefficients(p,d):
 coeff={}
 def add(a,c):
  quotient,a=divmod(a,p);coeff[a]=coeff.get(a,0)+c*2**quotient
 if not d:add(0,1)
 elif len(d)==2 and all(delta==1 for a,delta in d):
  add(0,1)
  for a,delta in d:add(a+1,1);add(a,-1)
 elif len(d)==2 and all(delta==-1 for a,delta in d):
  a,b=sorted(a for a,delta in d)
  for exponent,c in [(p-a,1),(1,-1),(0,1),(b-a+1,-1),(b-a,1)]:add(exponent,c)
 elif len(d)==2:
  a=next(a for a,delta in d if delta==1);b=next(a for a,delta in d if delta==-1)
  for exponent,c in [(p-b,1),(p-b+a+1,1),(p-b+a,-1),(1,-1),(0,1)]:add(exponent,c)
 else:raise ArithmeticError(('unexpected displacement',p,d))
 return sorted((a,c) for a,c in coeff.items() if c)
def norm_estimate(p,d):
 # Root product for the canonical INTEGER polynomial. This float is only
 # a ranking heuristic; exact 2,3-free norms are computed separately.
 coeff=canonical_coefficients(p,d)
 total=0.0
 for j in range(p):
  z=cmath.rect(2**(1/p),2*math.pi*j/p)
  Q=sum(c*z**a for a,c in coeff)
  if abs(Q)<1e-14:return None
  total+=math.log(abs(Q))
 return math.exp(total/p)

rows=[];hits=[];classes=Counter();path_count=0;distinct_count=0;norm_patterns={};same_two=[];max_content=None;max_strict_primitive=None;max_strict_repeated=None
for p in range(1,args.max_p+1):
 for sign in ['negative','positive']:
  q=(3**p).bit_length()-(sign=='negative')
  if q<=p:continue
  source=sum((((i+1)*p)//q-(i*p)//q)<<i for i in range(q))
  clock=2**q-3**p;need((clock<0)==(sign=='negative'),'clock sign')
  k,r=divmod(q,p)
  if math.gcd(p,q)==1:
   labels={((j*k+(-(-j*r//p))-1)%q):j for j in range(p)}
  else:labels={i:j for j,i in enumerate(i for i in range(q) if(source>>i)&1)}
  need(len(labels)==p and set(labels)=={i for i in range(q) if(source>>i)&1},'labels')
  orbit0=rotations(source,q);orbit1=set()
  first=list(edits(source,q))
  for i,v in first:orbit1.update(rotations(v,q))
  seen={};shapeclasses=Counter()
  for i,v in first:
   labels1,d1,l1,e1=move_labels(labels,i,q,{})
   for j,w in edits(v,q):
    labels2,d2,l2,e2=move_labels(labels1,j,q,d1)
    need(sum(d2.values())==e1+e2,'displacement sum')
    cl='undo' if not d2 else 'same_one_two' if len(d2)==1 else ('two_plus' if e1==e2==1 else 'two_minus' if e1==e2==-1 else 'mixed')
    shapeclasses[cl]+=1;classes[cl]+=1;path_count+=1
    if cl=='same_one_two':same_two.append([p,q,i,j])
    seen.setdefault(w,(i,j,cl,d2))
  maxg=0;strict2=0;shapehits=[]
  for w,(i,j,cl,d2) in seen.items():
   B=carry(w,q);g=math.gcd(abs(clock),B);maxg=max(maxg,g)
   distance=0 if w in orbit0 else 1 if w in orbit1 else 2
   strict2+=distance==2
   if g==abs(clock):
    hit=dict(p=p,q=q,sign=sign,shape_gcd=math.gcd(p,q),word=fmt(w,q),first_swap=i,second_swap=j,path_class=cl,necklace_distance=distance,fixed_point=B//clock)
    hits.append(hit);shapehits.append(hit)
   if max_content is None or math.log(g)/p>max_content['log_content_per_p']:
    max_content=dict(p=p,q=q,g=g,log_content_per_p=math.log(g)/p,root_content=g**(1/p),word=fmt(w,q),path=[i,j],necklace_distance=distance)
   if distance==2:
    candidate=dict(p=p,q=q,g=g,log_content_per_p=math.log(g)/p,root_content=g**(1/p),word=fmt(w,q),path=[i,j],necklace_distance=distance)
    if math.gcd(p,q)==1:
     if max_strict_primitive is None or candidate['log_content_per_p']>max_strict_primitive['log_content_per_p']:max_strict_primitive=candidate
    elif max_strict_repeated is None or candidate['log_content_per_p']>max_strict_repeated['log_content_per_p']:max_strict_repeated=candidate
   if math.gcd(p,q)==1 and 2<=p<=args.norm_max_p:
    ds=tuple(sorted(((-lab*r)%p,delta) for lab,delta in d2.items()))
    key=(p,ds)
    if key not in norm_patterns:
     est=norm_estimate(p,ds)
     norm_patterns[key]=dict(p=p,q=q,displacements=ds,norm_root_estimate=est,path=[i,j],word=fmt(w,q),content=g,necklace_distance=distance)
  distinct_count+=len(seen)
  rows.append(dict(p=p,q=q,sign=sign,shape_gcd=math.gcd(p,q),paths=sum(shapeclasses.values()),distinct_outputs=len(seen),exact_necklace_distance_two_outputs=strict2,max_content=maxg,classes=dict(shapeclasses),integer_hits=len(shapehits)))
# Known radius-three negative cycle must not enter either zero-to-two orbit.
p,q=7,11;source=sum((((i+1)*p)//q-(i*p)//q)<<i for i in range(q));target=int('11110111000'[::-1],2)
b={source}
for unused in range(2):b|={v for u in list(b) for i,v in edits(u,q)}
need(not (rotations(target,q)&b),'minus17 two-edit hostile')
need(all(p==1 for p,q,i,j in same_two),'same-one displacement two should be impossible in critical 1/2 gap words with p>=2')
ranked=sorted(norm_patterns.values(),key=lambda row:row['norm_root_estimate'] or -1,reverse=True)
print(json.dumps(dict(status='FINITE_EXACT_CENSUS_WITH_EXPLORATORY_FLOAT_NORM_RANKING',max_p=args.max_p,norm_max_p=args.norm_max_p,checks=checks,shape_count=len(rows),two_actual_edit_paths=path_count,distinct_outputs_by_shape=distinct_count,path_classes=dict(classes),same_one_displacement_two=same_two,integer_hits=hits,maximum_content_root=max_content,maximum_strict_primitive_content=max_strict_primitive,maximum_strict_nonprimitive_content=max_strict_repeated,norm_patterns=len(norm_patterns),top_norm_root_estimates=ranked[:30],minus17_not_within_two=True,rows=rows),indent=2))
