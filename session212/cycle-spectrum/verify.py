#!/usr/bin/env python3
"""Independent exact counters and hostile controls. All checks survive -O."""
from itertools import combinations,permutations
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent

def require(ok,reason):
 if not ok:raise RuntimeError(reason)

def valid(a):
 n=len(a);return all(not(a[i]>>i&1) for i in range(n)) and all(((a[i]>>j&1)+(a[j]>>i&1))==1 for i in range(n) for j in range(i)) and all(0<=r<(1<<n) for r in a)

def endpoint_counts(a):
 n=len(a);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for u,c in enumerate(dp[mask]):
   if not c:continue
   rem=a[u]&~mask
   while rem:
    bit=rem&-rem;rem-=bit;dp[mask|bit][bit.bit_length()-1]+=c
 return dp[-1]

def count_cycles(a):
 n=len(a)
 if n<3:return 0
 dp=[[0]*n for _ in range(1<<n)];dp[1][0]=1
 for mask in range(1,1<<n,2):
  for u,c in enumerate(dp[mask]):
   if not c:continue
   rem=a[u]&~mask
   while rem:
    bit=rem&-rem;rem-=bit;dp[mask|bit][bit.bit_length()-1]+=c
 return sum(c for u,c in enumerate(dp[-1]) if a[u]&1)

def literal(a,paths=True):
 n=len(a);e=[0]*n;c=0
 if paths:
  for p in permutations(range(n)):
   if all(a[u]>>v&1 for u,v in zip(p,p[1:])):e[p[-1]]+=1
 if n>=3:
  for q in permutations(range(1,n)):
   p=(0,)+q
   if all(a[p[i]]>>p[(i+1)%n]&1 for i in range(n)):c+=1
 return e,c

def strong(a):
 n=len(a);r=[v|1<<i for i,v in enumerate(a)]
 for k in range(n):
  for i in range(n):
   if r[i]>>k&1:r[i]|=r[k]
 return all(v==(1<<n)-1 for v in r)

def transitive_after(a,deleted):
 v=[u for u in range(len(a)) if u not in deleted]
 return all(not((a[x]>>y&1) and (a[y]>>z&1) and (a[z]>>x&1)) and not((a[x]>>z&1) and (a[z]>>y&1) and (a[y]>>x&1)) for x,y,z in combinations(v,3))

def near(m):
 require(m>=3,'near order');a=[sum(1<<j for j in range(i+1,m)) for i in range(m)];a[0]^=1<<(m-1);a[-1]^=1;return a

def gate(a,B):
 m=len(a);z=m;v=m+1
 return [r|((1<<v) if i in B else 0) for i,r in enumerate(a)]+[(1<<m)-1,(1<<z)|sum(1<<i for i in range(m) if i not in B)]

def realize(k,n=None):
 require(isinstance(k,int) and k>=0,'target k')
 if n is None:n=max(5,4+(k-2).bit_length()) if k>=2 else 5
 require(n>=5 and k<=2**(n-4)+1,'target order')
 m=n-2;weights=[1]+[1<<(i-1) for i in range(1,m-1)]+[1];left=k;B=[]
 for i in sorted(range(m),key=lambda i:weights[i],reverse=True):
  if weights[i]<=left:left-=weights[i];B.append(i)
 require(left==0,'subset sum compiler');return gate(near(m),B)

def all_labelled(n):
 pairs=list(combinations(range(n),2))
 for code in range(1<<len(pairs)):
  a=[0]*n
  for k,(i,j) in enumerate(pairs):
   if code>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
  yield a

def main():
 receipt={};checks=0;hostile={}
 for n in range(1,6):
  for a in all_labelled(n):
   e,c=literal(a);require(endpoint_counts(a)==e,'endpoint literal');require(count_cycles(a)==c,'cycle literal');require((c>0)==(n>=3 and strong(a)),'Camion control');checks+=1
 receipt['literal_all_labelled_orders_1_to_5']=checks
 for m in range(3,13):
  a=near(m);require(endpoint_counts(a)==[1]+[1<<(i-1) for i in range(1,m-1)]+[1],'binary endpoint formula');require(count_cycles(a)==1,'unique cycle')
 receipt['endpoint_orders']=list(range(3,13));gates=0
 for m in range(3,7):
  a=near(m);e=endpoint_counts(a)
  for mask in range(1<<m):
   B=[i for i in range(m) if mask>>i&1];t=gate(a,B);k=sum(e[i] for i in B)
   require(valid(t),'tournament gate');require(count_cycles(t)==k,'gate law');require(strong(t)==(k>0),'strong gate');require(transitive_after(t,{0,m+1}),'feedback two upper')
   if k>=2:require(not any(transitive_after(t,{v}) for v in range(m+2)),'feedback two lower')
   gates+=1
 receipt['all_subset_gates_m3_to_6']=gates
 targets=sorted(set(range(66))|{127,128,129,255,256,257,511,512,513,962,970,971})
 for k in targets:
  a=realize(k);require(count_cycles(a)==k,'integer compiler');require(valid(a),'compiled tournament');require(strong(a)==(k>0),'compiled strongness')
 receipt['compiled_targets']=targets
 fixed={3:[1],4:[1],5:list(range(1,4)),6:list(range(1,6)),7:list(range(1,16))+[17,24],8:list(range(1,55))}
 witness_checks=0
 for n,expected in fixed.items():
  d=json.loads((ROOT/f'fixedcycle{n}.json').read_text());require(d['spectrum']==expected,'fixed-cycle spectrum')
  require(d['labelled_or_extension_presentations']==1<<(n*(n-1)//2-n),'fixed-cycle universe')
  for k,a in d['witnesses'].items():
   require(valid(a) and count_cycles(a)==int(k),'fixed-cycle witness');witness_checks+=1
 d=json.loads((ROOT/'spectrum9.json').read_text());require(d['labelled_or_extension_presentations']==6880*256,'extension universe')
 expected_gaps=[215,219,223,226,229]+list(range(232,243))+list(range(244,249));require(d['gaps']==expected_gaps,'order9 gaps')
 for k,a in d['witnesses'].items():
  require(valid(a) and count_cycles(a)==int(k),'order9 witness');witness_checks+=1
 receipt['independent_census_witness_checks']=witness_checks
 inherited=json.loads((ROOT/'inherited_summary.json').read_text())
 for n in (8,9):
  mine=json.loads((ROOT/('fixedcycle8.json' if n==8 else 'spectrum9.json')).read_text())['spectrum']
  if n==8:mine=[0]+mine
  require(inherited[str(n)]['spectrum']==mine,'independent inherited agreement')
 ten=inherited['10'];require(ten['spectrum']==list(range(962))+[971],'inherited n10 spectrum');require((ten['classes'],ten['strong_classes'])==(9733056,9355949),'inherited n10 universe')
 maxwit=json.loads((ROOT/'maximum10_witness.json').read_text());a=maxwit['rows'];require(valid(a),'n10 witness validity');require(sum(endpoint_counts(a))==15745,'n10 path maximum witness');require(count_cycles(a)==971,'n10 cycle maximum DP');require(literal(a,False)[1]==971,'n10 literal maximum cycle witness')
 receipt['inherited_order10']={'scope':'re-extraction, not independent full enumeration','classes':9733056,'strong_classes':9355949,'spectrum':'[0,961] union {971}','maximum_witness_independently_counted':True}
 receipt['hostile_order7_missing_values']=[16,18,19,20,21,22,23]
 receipt['minimum_feedback_number_for_k_ge_2']=2
 print(json.dumps(receipt,sort_keys=True,indent=2))
if __name__=='__main__':main()
