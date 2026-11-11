import itertools,json
from pathlib import Path

def rows(n,code):
 r=[0]*n
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if code>>k&1:r[i]|=1<<j
  else:r[j]|=1<<i
 return r

def triangles(r):
 return [sum(1<<i for i in t) for t in itertools.combinations(range(len(r)),3) if all((r[i]&sum(1<<j for j in t)).bit_count()==1 for i in t)]

def inflate(r,sizes):
 types=[i for i,s in enumerate(sizes) for _ in range(s)]
 return [sum(1<<b for b,j in enumerate(types) if (i==j and a<b) or (i!=j and r[i]>>j&1)) for a,i in enumerate(types)]

def hp(r):
 n=len(r);d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for s in range(1,1<<n):
  for v,c in enumerate(d[s]):
   if c:
    t=r[v]&~s
    while t:
     b=t&-t;t-=b;d[s|b][b.bit_length()-1]+=c
 return sum(d[-1])

def cycle_masks(r):
 result=[]
 for v in range(len(r)):
  def rec(w,mask,length):
   if length>=3 and r[w]>>v&1:result.append(mask)
   t=r[w]&~mask&~((1<<(v+1))-1)
   while t:
    b=t&-t;t-=b;rec(b.bit_length()-1,mask|b,length+1)
  rec(v,1<<v,1)
 return result

def formula(r,sizes):
 from math import prod
 return 1+sum(prod((1<<sizes[i])-1 for i in range(len(r)) if mask>>i&1)-(-1)**mask.bit_count() for mask in cycle_masks(r))

def oddcycles(r):
 count=0
 for v in range(len(r)):
  def rec(w,mask,length):
   nonlocal count
   if length>=3 and length%2 and r[w]>>v&1:count+=1
   t=r[w]&~mask&~((1<<(v+1))-1)
   while t:
    b=t&-t;t-=b;rec(b.bit_length()-1,mask|b,length+1)
  rec(v,1<<v,1)
 return count

def main():
 # Complete labelled five-template universe, quotiented only by exact relabelling.
 seen=set();cores=[];pairs=list(itertools.combinations(range(5),2))
 for code in range(1024):
  if code in seen:continue
  r=rows(5,code);cores.append((code,r))
  for p in itertools.permutations(range(5)):
   seen.add(sum((r[p[i]]>>p[j]&1)<<k for k,(i,j) in enumerate(pairs)))
 witnesses={};cases=0;nonapex={};core_records=[]
 for code,r in cores:
  ts=triangles(r);common=31
  for t in ts:common&=t
  if ts and not common:
   safe=[s for s in range(32) if all(a&b&~s for a in ts for b in ts)]
   maximal=[s for s in safe if not any(s!=t and s&t==s for t in safe)]
   from collections import Counter
   core_records.append({"code":code,"rows":r,"triangles":ts,"maximal_expandable_masks":maximal,"cycle_support_multiplicities":dict(Counter(cycle_masks(r)))})
 for code,r in cores:
  ts=triangles(r)
  for sizes in itertools.product(range(1,3),repeat=5):
   singleton=sum(1<<i for i,s in enumerate(sizes) if s==1)
   predicted=all(a&b&singleton for a in ts for b in ts)
   actual_ts=triangles(inflate(r,sizes))
   actual=all(a&b for a in actual_ts for b in actual_ts)
   assert predicted==actual
 hostile=inflate([2,4,1],[2,2,2])
 assert hp(hostile)!=formula([2,4,1],[2,2,2])
 for code,r in cores:
  ts=triangles(r)
  for sizes in itertools.product(range(1,5),repeat=5):
   if sum(sizes)>11:continue
   singleton=sum(1<<i for i,s in enumerate(sizes) if s==1)
   if any(not(a&b&singleton) for a in ts for b in ts):continue
   cases+=1;t=inflate(r,sizes);h=hp(t)
   assert h==formula(r,sizes)
   w={'template_code':code,'sizes':sizes,'rows':t,'h':h}
   witnesses.setdefault(h,w)
   common=(1<<5)-1
   for a in ts:common&=a
   if ts and not common:nonapex.setdefault(h,w)
 for h,w in nonapex.items():
  assert hp(w['rows'])==1+2*oddcycles(w['rows'])
 result={'scope':'All 12 unlabelled five-vertex templates; sizes 1..4, total order <=11; safe by singleton intersection criterion.', 'templates':len(cores),'nonapex_core_classification':core_records,'criterion_controls':384,'hostile_C3_TT2':{'rows':hostile,'actual_H':hp(hostile),'invalid_formula_value':formula([2,4,1],[2,2,2])},'cases':cases,'h_values':sorted(witnesses),'missing_odd_through_101':[h for h in range(1,102,2) if h not in witnesses], 'nonapex_values':sorted(nonapex),'nonapex_witnesses':nonapex}
 Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
 print({k:v for k,v in result.items() if k!='nonapex_witnesses'})
if __name__=='__main__':main()
