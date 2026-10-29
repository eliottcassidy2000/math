"""Independent exhaustive matrix audit; stdlib, no producer imports or asserts."""
from itertools import combinations, product
import json

def need(ok,msg):
 if not ok: raise RuntimeError(msg)
def matrix(n,bits):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=1 if bits>>k&1 else -1;a[j][i]=-a[i][j]
 return a
def switch(a,d):return [[a[i][j]*d[i]*d[j] for j in range(len(a))] for i in range(len(a))]
def module(a,S):return all(len({a[x][y] for x in S})==1 for y in range(len(a)) if y not in S)
def rankone(a,S):
 C=[i for i in range(len(a)) if i not in S]
 return all(a[x][y]*a[u][v]==a[x][v]*a[u][y] for x,u in combinations(S,2) for y,v in combinations(C,2))
def normal(a):
 d=[1]+[a[0][j] for j in range(1,len(a))]
 return switch(a,d)
counts=[]; hostile=None
for n in range(2,6):
 nc=ns=0
 for bits in range(1<<(n*(n-1)//2)):
  a=matrix(n,bits);b=normal(a)
  need(all(b[0][j]==1 for j in range(1,n)), 'anchor source normalization')
  for tail in product((-1,1),repeat=n-1):
   d=(1,)+tail; c=switch(a,d)
   need(normal(c)==b,'normalized representative changes inside switch orbit');ns+=1
  for mask in range(1,(1<<n)-1):
   S=[i for i in range(n) if mask>>i&1];C=[i for i in range(n) if not mask>>i&1]
   r=rankone(a,S)
   if 0 not in S:
    need(module(b,S)==r,'root-outside rankone cut versus normalized deletion module')
   # Construct switches only on S using first outside reference. No rank routine reuse.
   d=[a[x][C[0]] if x in S else 1 for x in range(n)]
   need(module(switch(a,d),S)==r,'rankone versus constructive normalization');nc+=1
   if hostile is None and len(S)>1 and len(C)>1 and r and not module(a,S):
    hostile={'n':n,'bits':bits,'S':S,'matrix':a,'switch':d}
 counts.append({'n':n,'tournaments':1<<(n*(n-1)//2),'proper_cut_checks':nc,'switch_normal_checks':ns})
# Distinguish operation masks independently of orientation; exhaustive cut-mask comparison.
operations=[]
for n in range(1,9):
 edges=list(combinations(range(n),2));path={(i,i+1) for i in range(n-1)}
 mask=frozenset(e for e in edges if e not in path)
 cuts={frozenset((i,j) for i,j in edges if ((s>>i)^(s>>j))&1) for s in range(1<<n)}
 is_switch=mask in cuts;is_converse_switch=frozenset(set(edges)-mask) in cuts
 need(is_switch==(n<=2),'native move/switching boundary')
 need(is_converse_switch==(n<=3),'native move/converse-switching boundary')
 operations.append({'n':n,'native_is_switch':is_switch,'native_is_converse_switch':is_converse_switch})
# Distinguish ordinary primality from nontrivial signed-cut primality.
c5=[[0 if i==j else (1 if (j-i)%5 in (1,2) else -1) for j in range(5)] for i in range(5)]
p7=[[0 if i==j else (1 if (j-i)%7 in (1,2,4) else -1) for j in range(7)] for i in range(7)]
source_p7=[[0]+[1]*7]+[[-1]+row for row in p7]
def countcuts(a):
 n=len(a);return sum(rankone(a,[i for i in range(n) if mask>>i&1]) for mask in range(1,1<<n) if 2<=mask.bit_count()<=n-2)
need(countcuts(c5)==10,'cyclic5 ten signed-cut sides')
need(not any(module(c5,[i for i in range(5) if mask>>i&1]) for mask in range(1,32) if 2<=mask.bit_count()<=4),'cyclic5 ordinary prime')
need(countcuts(source_p7)==0 and module(source_p7,list(range(1,8))),'source Paley7 signed prime / ordinary decomposable')
print(json.dumps({'status' :'FINITE-EXACT','cut_tests':counts,'rankone_nonmodule_hostile':hostile,'operation_masks':operations},sort_keys=True,indent=2))
