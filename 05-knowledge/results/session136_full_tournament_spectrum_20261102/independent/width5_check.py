"""Independent width<=5 capped fixed point; no producer imports.
Quotient-word kernel is grown forwards by length rather than end-recursion.
"""
from kernel import *
from itertools import product,combinations
from collections import defaultdict
import json,hashlib,pathlib
B=127;N=7

def prime(a):
 n=len(a)
 for size in range(2,n):
  for subset in combinations(range(n),size):
   if all(len({(a[v]>>u)&1 for u in subset})==1 for v in range(n)if v not in subset):return False
 return True

def canon(a):return min(induced(a,p)for p in permutations(range(len(a))))
classes=defaultdict(list)
for code in range(1024):
 a=tournament(5,code)
 if prime(a):classes[canon(a)].append(code)
counts=sorted(len(v)for v in classes.values());reps=[tournament(5,min(v))for v in classes.values()]
require(counts==[24,120,120] and sorted(min(v)for v in classes.values())==[10,41,76],'prime5 classification')
quotients=[(2,4,1)]+reps

def words_forward(q):
 n=len(q);level={}
 for i in range(n):
  c=tuple(int(j==i)for j in range(n));level[(c,i)]=1
 table={}
 for size in range(1,11):
  nxt=defaultdict(int)
  for (c,last),count in level.items():
   table[c]=table.get(c,0)+count
   if size<10:
    for v in range(n):
     if q[last]>>v&1 and c[v]<6:
      d=list(c);d[v]+=1;nxt[(tuple(d),v)]+=count
  level=nxt
 return table
Ws=[words_forward(q)for q in quotients]

def hp_from_profiles(table,ps):
 total=0
 for c in product(*(range(1,len(p))for p in ps)):
  coeff=table.get(c,0)
  if coeff:
   for p,k in zip(ps,c):coeff*=factorial(k)*p[k]
   total+=coeff
 return total

def state(a):return ('s',len(a),path_cover_literal(a))
def H(s):return s[2][1]if s[0]=='s'else s[1]

def budget_tuples(states,k,budget):
 if k==0:yield ();return
 for s in states:
  if s[1]<=budget-(k-1):
   for tail in budget_tuples(states,k-1,budget-s[1]):yield(s,)+tail

states={('s',1,(0,1))};graphs={next(iter(states)):(0,)};rounds=0;checks=0
while True:
 old=list(states);small=[s for s in old if s[0]=='s'];new={}
 for x,y in product(old,repeat=2):
  if x[0]==y[0]=='s'and x[1]+y[1]<N:
   a=ordinal(graphs[x],graphs[y]);new.setdefault(state(a),a)
  else:new.setdefault(('L',min(B+1,H(x)*H(y))),None)
 for q,W in zip(quotients,Ws):
  for ss in budget_tuples(small,len(q),10):
   n=sum(s[1]for s in ss)
   h=hp_from_profiles(W,[s[2]for s in ss])
   if n<N:
    a=substitute(q,[graphs[s]for s in ss]);s=state(a)
    require(s[2][1]==h,('kernel direct',q,ss));new.setdefault(s,a);checks+=1
   else:new.setdefault(('L',min(B+1,h)),None)
 # All omitted active outputs are saturated: a large module already yields
 # >=129; all-small total>=11 is strong and H>=f11=225.
 if any(s[0]=='L'for s in old):new.setdefault(('L',128),None)
 new={s:a for s,a in new.items()if s not in states}
 if not new:break
 states.update(new);graphs.update({s:a for s,a in new.items()if a is not None});rounds+=1
values=sorted({H(s)for s in states if H(s)<=B});missing=[h for h in range(1,B+1,2)if h not in values]
require(len(states)==86 and len(values)==49,'width5 closure size')
require(missing==[7,21,31,35,43,53,67,71,77,83,89,101,103,113,127],'width5 missing support')
a=tournament(6,146);require(a==(36,17,2,7,13,30) and h_literal(a)==31,'width6 H31')
print(json.dumps({'status':'PASS','prime5_representative_codes':sorted(min(v)for v in classes.values()),'prime5_class_labelled_sizes':counts,'states':len(states),'attained':values,'missing':missing,'small_by_order':{n:sum(s[0]=='s'and s[1]==n for s in states)for n in range(1,N)},'rounds_with_new_states':rounds,'literal_small_kernel_checks':checks,'H31_order6_witness':a,'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()},sort_keys=True,indent=2))
