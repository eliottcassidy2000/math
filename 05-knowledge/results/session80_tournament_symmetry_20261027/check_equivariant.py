#!/usr/bin/env python3
"""Independent finite checks of recovered tournament/C3 equivariant indices.

No saved receipts, supplied prime census, or original verifier scripts are read.
Run beside recover_converse.py and recover_c3.py; standard library only.
"""
import itertools,json
import recover_converse as r
import recover_c3 as c

def fixed_masks(part):
 n=sum(part);sigma=[];start=0
 for d in part:sigma.extend(list(range(start+1,start+d))+[start]);start+=d
 pairs=list(itertools.combinations(range(n),2));ix={p:i for i,p in enumerate(pairs)}
 seen=set();orbits=[]
 for i in range(len(pairs)):
  if i in seen:continue
  here=i;parity=0;orb={}
  while here not in orb:
   orb[here]=parity;seen.add(here);a,b=pairs[here];a,b=sigma[a],sigma[b]
   parity^=1^(a>b);here=ix[tuple(sorted((a,b)))]
  if parity!=orb[here]:return
  orbits.append(orb)
 for choice in range(1<<len(orbits)):
  rows=[0]*n
  for k,orb in enumerate(orbits):
   for i,parity in orb.items():
    a,b=pairs[i]
    if parity^((choice>>k)&1):a,b=b,a
    rows[a]|=1<<b
  yield rows

def prime_subsets(rows):
 n=len(rows)
 if n<3:return False
 full=(1<<n)-1
 for m in range(1,full):
  if m.bit_count()<2:continue
  outside=full^m
  while outside:
   bit=outside&-outside;outside-=bit;v=bit.bit_length()-1
   if rows[v]&m not in (0,m):break
  else:return False
 return True

def prime_closure(rows):
 n=len(rows)
 if n<3:return False
 full=(1<<n)-1
 for i in range(n):
  for j in range(i+1,n):
   m=(1<<i)|(1<<j)
   while True:
    old=m
    for k in range(n):
     if not (m>>k)&1 and (rows[k]&m) not in (0,m):m|=1<<k
    if m==old:break
   if m!=full:return False
 return True

def cycles(p):
 seen=set();out=[]
 for i in range(len(p)):
  if i in seen:continue
  j=i;k=0
  while j not in seen:seen.add(j);k+=1;j=p[j]
  out.append(k)
 return tuple(sorted(out))

def strong(rows):
 n=len(rows)
 if not n:return False
 full=(1<<n)-1
 for reverse in (False,True):
  reach=1
  while True:
   old=reach
   for i in range(n):
    if reach>>i&1:
     if reverse:reach|=sum(1<<j for j in range(n) if rows[j]>>i&1)
     else:reach|=rows[i]
   if old==reach:break
  if reach!=full:return False
 return True

def check_anti():
 r.N=18
 *_,pt=r.recover();tests=[]
 for n in range(8):
  for part in r.parts(n):
   count=pc=0
   for rows in fixed_masks(part):count+=1;pc+=prime_subsets(rows)
   expected=r.fixed(part,True);ep=pt.get(part,0)*r.z(part)
   if count!=expected or pc!=ep:raise RuntimeError((part,count,expected,pc,str(ep)))
   tests.append({'partition':part,'fixed':count,'prime_fixed':pc})
 return {'independent_xor_and_all_subsets_through':7,'partition_tests':len(tests),'tests':tests}

def check_full_cycle():
 r.N=18
 *_,pt=r.recover()
 f={m:2**((m-1)//2) for m in range(1,10,2)};h={}
 rows_out=[]
 for m in range(1,10,2):
  h[m]=2**m-2*f[m]-sum(h[d]*f[m//d] for d in range(1,m) if m%d==0)
  total=prime=0
  for rows in fixed_masks((2*m,)):
   total+=1;p=prime_closure(rows);prime+=p
   if 2*m<=10 and p!=prime_subsets(rows):raise RuntimeError('closure vs all subsets')
  if prime!=h[m] or prime!=2*m*pt.get((2*m,),0) or total!=2**m:raise RuntimeError((m,total,prime,h[m]))
  rows_out.append({'n':2*m,'anti_fixed':total,'prime_anti_fixed':prime,'recurrence':h[m]})
 for n in range(4,19,4):
  if list(fixed_masks((n,))):raise RuntimeError('divisible4 no-go')
 return {'full_cycle_checks':rows_out,'module_check':'minimal closure of each vertex pair; cross-checked against exhaustive proper subsets through10','divisible_by_four_full_cycles_through16':'zero, checked independently by XOR'}

def check_c3():
 r.N=6
 expected_c,expected_d,_=c.recover_c3();actual_c={};actual_d={};reports=[]
 for n in range(7):
  pairs=list(itertools.combinations(range(n),2));triples=list(itertools.combinations(range(n),3));ti={t:k for k,t in enumerate(triples)}
  image={}
  for mask in range(1<<len(pairs)):
   rows=[0]*n
   for k,(a,b) in enumerate(pairs):
    if mask>>k&1:a,b=b,a
    rows[a]|=1<<b
   h=0
   for k,(a,b,d) in enumerate(triples):
    if ((rows[a]>>b&1)==(rows[b]>>d&1)==(rows[d]>>a&1)):h|=1<<k
   s=strong(rows)
   if h in image and image[h]!=s:raise RuntimeError('strong realization disagreement')
   image[h]=s
  perms=[]
  for p in itertools.permutations(range(n)):
   trans=[ti[tuple(sorted(p[j] for j in t))] for t in triples]
   perms.append((trans,cycles(p)))
  unseen=set(image);classes=dclasses=0
  while unseen:
   h=min(unseen);orbit=set();auts={};edges=[k for k in range(len(triples)) if h>>k&1]
   for trans,ct in perms:
    mapped=sum(1<<trans[k] for k in edges);orbit.add(mapped)
    if mapped==h:auts[ct]=auts.get(ct,0)+1
   if not orbit<=set(image):raise RuntimeError('image not relabel invariant')
   unseen-=orbit;classes+=1;dclasses+=image[h];order=sum(auts.values())
   for ct,mul in auts.items():
    actual_c[ct]=actual_c.get(ct,0)+r.F(mul,order)
    if image[h]:actual_d[ct]=actual_d.get(ct,0)+r.F(mul,order)
  reports.append({'n':n,'labeled_tournaments':1<<len(pairs),'distinct_labeled_C3_hypergraphs':len(image),'C3_isotypes':classes,'strong_realizable_C3_isotypes':dclasses})
 if actual_c!=expected_c or actual_d!=expected_d:raise RuntimeError('actual hypergraph index mismatch')
 return {'checks':'PASS actual C3 image sets, relabeling orbits and hypergraph automorphism cycle indices; all coefficients C and D through6','census':reports,'actual_C_cycle_index':r.serial(actual_c),'actual_D_cycle_index':r.serial(actual_d)}

if __name__=='__main__':
 print(json.dumps({'anti_fixed':check_anti(),'anti_full_cycles':check_full_cycle(),'c3_image_sets':check_c3()},indent=2))
