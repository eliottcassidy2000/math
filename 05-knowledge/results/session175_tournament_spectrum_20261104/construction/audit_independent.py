import itertools,json,functools,importlib.util
from pathlib import Path
p=Path(__file__).parent
spec=importlib.util.spec_from_file_location('author',p/'check.py');author=importlib.util.module_from_spec(spec);spec.loader.exec_module(author)
def require(x,msg):
 if not x:raise RuntimeError(msg)
def count(a):
 @functools.lru_cache(None)
 def f(mask,last):
  if not mask:return 1
  return sum(f(mask^(1<<v),v) for v in range(len(a)) if mask>>v&1 and a[last]>>v&1)
 return sum(f(((1<<len(a))-1)^(1<<v),v) for v in range(len(a)))
def perms(a):
 return [w for w in itertools.permutations(range(len(a))) if all(a[w[i]]>>w[i+1]&1 for i in range(len(w)-1))]
words=0;bank={1};minbound=999999
for n in range(10):
 for w in itertools.product((0,1),repeat=n):
  a=[0]*(n+1)
  for i in range(n):
   for j in range(i+1,n):a[i]|=1<<j
   if w[i]:a[n]|=1<<i
   else:a[i]|=1<<n
  h=count(a);word=''.join(map(str,w))
  require(h==author.stream(word)==author.literal(word),'word '+word)
  if h<=199:bank.add(h)
  words+=1
own_missing=[h for h in range(1,200,2) if h not in bank]
solve,bound=author.census(199)
require(bound==9,'length bound')
require(bank==set(solve),'complete inverse')
require(all(h in own_missing for h in [7,13,21,29,35,41,43,49,55,57]),'gaps')
orphan_checks=0
for n in range(1,5):
 for code in range(1<<(n*(n-1)//2)):
  a=author.labelled(n,code);old=perms(a);oldset=set(old)
  for b in itertools.product((0,1),repeat=n):
   aug=a[:]+[0]
   for i,x in enumerate(b):
    if x:aug[n]|=1<<i
    else:aug[i]|=1<<n
   new=perms(aug)
   orphan=sum(tuple(v for v in w if v!=n) not in oldset for w in new)
   d=sum(sum(b[w[i]]==1 and b[w[i+1]]==0 for i in range(n-1)) for w in old)
   require(len(new)==len(old)+d+orphan,'orphan accounting')
   orphan_checks+=1
B=[8,1,3,6];fam=[]
for r in range(1,9):
 A=[sum(1<<j for j in range(i+1,r)) for i in range(r)]
 h=count(author.wrapper(A,B));require(h==27*2**r-17,'family')
 fam.append(h)
for k in range(100):
 h=27*2**(3*k+2)-17
 require(h%7==0,'seven')
 require((h//7)%7==(6+3*k)%7,'valuation')
result={'status':'PASS','independent_recursive_word_cases':words,'complete_missing_through_199':own_missing,'exhaustive_orphan_cases':orphan_checks,'asymmetric_family_r1_to8':fam,'valuation_cases':100}
print(json.dumps(result,indent=2))
