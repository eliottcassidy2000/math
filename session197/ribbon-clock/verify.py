import itertools,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'ribbon'))
from transfer import graph,held_karp

def clock(n):
 a=[[False]*n for _ in range(n)];a[0][1]=True
 for b in range(2,n):a[1][b]=True;a[b][0]=True
 for i in range(2,n):
  for j in range(i+1,n):a[i][j]=True
 return a

def span(a,order):return max([0]+[j-i for j in range(len(order))for i in range(j)if a[order[j]][order[i]]])
def balanced(n):
 w=(n+1)//3
 g1=max(0,n-1-2*w);g2=min(w,n-1-g1);g3=n-1-g1-g2
 assert 0<=g1<=w and 1<=g2<=w and 0<=g3<=w
 p=g1;q=g1+g2;it=iter(range(2,n))
 return [1 if i==p else 0 if i==q else next(it)for i in range(n)]

res={}
for n in range(3,9):
 a=clock(n);mn=min(span(a,p)for p in itertools.permutations(range(n)))
 assert mn==(n+1)//3
 res[str(n)]={'permutations':__import__('math').factorial(n),'minimum_span':mn}
for n in range(3,13):assert held_karp(clock(n))==2**(n-2)+1
for n in range(3,201):assert span(clock(n),balanced(n))==(n+1)//3
witness=list(map(int,'022333332'));assert held_karp(graph(witness,2))==513
result={'exhaustive_orderings':res,'held_karp_through_order':12,'balanced_construction_through_order':200,
 'equal_count_contrast':{'clock_order':11,'clock_minimum_span':4,'count':513,'span2_word':'022333332','span2_order':9}}
(HERE/'VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
