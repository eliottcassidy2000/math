"""Independent literal strong-wrapper controls and exact 5/17 arithmetic."""
from itertools import permutations
from pathlib import Path
import json

def require(ok,why):
 if not ok: raise RuntimeError(why)

def book(k):
 n=k+2;a=[[False]*n for _ in range(n)];a[0][1]=True
 for j in range(2,n):
  a[1][j]=True;a[j][0]=True
  for t in range(j+1,n):a[j][t]=True
 return a

def p2_brute(a):
 n=len(a);total=0
 for v in permutations(range(n)):
  path=lambda l,u:all(a[x][y] for x,y in zip(v[l:u],v[l+1:u]))
  total+=sum(path(0,c) and path(c,n) for c in range(n+1))
 return total

def ordinal(blocks):
 n=sum(map(len,blocks));a=[[False]*n for _ in range(n)];offset=0
 for b in blocks:
  for i in range(len(b)):
   for j in range(len(b)):a[offset+i][offset+j]=b[i][j]
   for j in range(offset+len(b),n):a[offset+i][j]=True
  offset+=len(b)
 return a

def strong_wrapper(M):
 # Physical names X=0,Y=1 remain fixed under all phase changes.
 n=len(M)+2;a=[[False]*n for _ in range(n)];a[0][1]=True
 for i in range(len(M)):
  a[1][i+2]=True;a[i+2][0]=True
  for j in range(len(M)):a[i+2][j+2]=M[i][j]
 return a

def count(a,start=None,deleted=None):
 n=len(a);incoming=[sum(1<<i for i in range(n) if a[i][j]) for j in range(n)]
 dp=[[0]*n for _ in range(1<<n)]
 for j in range(n):dp[1<<j][j]=int(start is None or j==start)
 for mask in range(1,1<<n):
  rest=mask
  while rest:
   bit=rest&-rest;rest-=bit;j=bit.bit_length()-1;before=mask^bit
   if not before:continue
   opts=before&incoming[j];value=0
   while opts:
    kb=opts&-opts;opts-=kb;value+=dp[before][kb.bit_length()-1]
   dp[mask][j]=value
 full=sum(dp[-1]);card=None if deleted is None else sum(dp[((1<<n)-1)^(1<<deleted)])
 return full,card

def strongly_connected(a):
 n=len(a)
 for start in range(n):
  seen={start};front=[start]
  while front:
   u=front.pop()
   for v in range(n):
    if a[u][v] and v not in seen:seen.add(v);front.append(v)
  if len(seen)!=n:return False
 return True

A=book(2);B=book(1);alpha=p2_brute(A);beta=p2_brute(B)
require((alpha,beta)==(32,12),'small P2 values')
phase=(2,1);controls=[]
for j in range(6):
 r,s=phase
 names=[('A',i) for i in range(r)]+[('B',i) for i in range(s)]+[('A',i) for i in range(r,2)]+[('B',i) for i in range(s,2)]
 M=ordinal([A if kind=='A' else B for kind,i in names]);T=strong_wrapper(M)
 require(strongly_connected(T),'literal strongness')
 require([i for i in range(len(T)) if T[i][1]]==[0],'same indegree-one root and port')
 full,card=count(T,deleted=1);start,_=count(T,start=1)
 require(full==225+(32*12)**2==147681,'full strong target')
 require(card==start==225,'root deletion/start identity, defect zero')
 controls.append({'coordinates':phase,'blocks':[kind+str(i+1) for kind,i in names],'vertices':len(T),'full_H':full,'root':1,'sole_in_neighbor':0,'root_start_H':start,'root_deletion_H':card,'defect':card-start})
 phase=(2-s,r+s-1)
require(phase==(2,1),'exact controller return')
B17=book(4);beta17=p2_brute(B17)
require(beta17==224,'H17 core P2')
N=85**2+(32*224)**2
require(N==51387449==29*1771981,'5/17 strong target and proper factor')
hostile=[]
for module in (ordinal([book(1),book(1)]),ordinal([book(3),[[False]]])):
 h,_=count(module);cover=p2_brute(module);outer,_=count(strong_wrapper(module))
 require(h==9 and outer==h+cover,'same scalar, different contextual response')
 hostile.append({'order':len(module),'H':h,'P2':cover,'wrapper_H':outer})
require([x['P2'] for x in hostile]==[144,168],'hostile P2 values')
require([x['wrapper_H'] for x in hostile]==[153,177],'hostile strong responses')
record={'equal_order_equal_H_hostile':hostile,'status':'PASS','literal_full_and_root_subset_DP_controls':controls,'H5_H17_wrapper':{'core_counts':[5,17],'P2_values':[32,224],'vertices':22,'full_H_by_proved_cover_identity':N,'factorization':[29,1771981],'root_start_and_deletion':7225,'defect':0},'Fermat_unit_branch':[{'m':m,'H':1+2**(4*m)} for m in (1,2,4,8)]}
Path(__file__).with_name('strong_wrapper_receipts.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'status':'PASS','literal_controls':6,'literal_target':147681,'H5_H17_target':N,'P2_H5_H17':[32,beta17]}))
