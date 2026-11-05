from functools import lru_cache
from pathlib import Path
import json
A=[2,8,3,5];B=[12,5,8,2]
def build(k):
 n=k+8;a=[[i<j for j in range(n)]for i in range(n)]
 for i in range(4):
  for j in range(4):a[i][j]=bool(A[i]>>j&1);a[k+4+i][k+4+j]=bool(B[i]>>j&1)
 spine=[2]+list(range(4,k+4))+[k+6]
 for x,y in zip(spine,spine[1:]):a[x][y]=False;a[y][x]=True
 return a

def all_cycles(a):
 n=len(a);cs=[0]*(1<<n)
 for root in range(n):
  @lru_cache(None)
  def ways(S,last):
   if S==1<<root:return int(last==root)
   if last==root:return 0
   R=S^(1<<last)
   return sum(ways(R,u)for u in range(root,n)if R>>u&1 and a[u][last])
  for S in range(1<<root,1<<n,1<<(root+1)):
   if S.bit_count()>=3:cs[S]=sum(ways(S,u)for u in range(root+1,n)if S>>u&1 and a[u][root])
 return cs

def hp(a):
 @lru_cache(None)
 def paths(S,last):
  if S==1<<last:return 1
  R=S^(1<<last)
  return sum(paths(R,u)for u in range(len(a))if R>>u&1 and a[u][last])
 return sum(paths((1<<len(a))-1,u)for u in range(len(a)))

def span2_order(k):return [3,0,1,2]+list(range(4,k+4))+[k+6,k+5,k+4,k+7]

def check_reorder(k):
 a=build(k);p=span2_order(k)
 assert sorted(p)==list(range(k+8))
 return max([0]+[j-i for j in range(len(p))for i in range(j)if a[p[j]][p[i]]])

rows=[];tot=0
for k in range(8):
 a=build(k);cs=all_cycles(a);h=hp(a)
 assert cs[-1]==4 and max(cs[:-1])==2
 assert check_reorder(k)==2
 rows.append({'k':k,'n':len(a),'H':h,'C':cs[-1],'proper_subsets':len(cs)-1,'maximum_proper_C':max(cs[:-1])})
 tot+=len(cs)-1
assert [x['H']for x in rows[:6]]==[105,187,349,641,1177,2167]
assert all(rows[k]['H']==sum(rows[j]['H']for j in range(k-3,k))for k in range(3,8))
assert all(check_reorder(k)==2 for k in range(101))
result={'independent_family_checks':rows,'all_proper_subset_checks':tot,'uniform_span2_order_checked_through_k':100}
Path(__file__).with_name('INDEPENDENT_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
