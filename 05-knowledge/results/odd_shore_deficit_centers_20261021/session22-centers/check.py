#!/usr/bin/env python3
"""Fresh scoped exact controls, stdlib only. Run python -B check.py."""
import itertools, json, hashlib, pathlib, datetime, time
ROOT=pathlib.Path(__file__).resolve().parent
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x: raise RuntimeError(msg)
def strong(A,S=None):
 S=set(range(len(A))) if S is None else set(S)
 if not S:return False
 def reach(rev):
  got={next(iter(S))}; todo=list(got)
  while todo:
   u=todo.pop()
   for v in S-got:
    if A[v][u] if rev else A[u][v]:got.add(v);todo.append(v)
  return got==S
 return reach(0) and reach(1)
def circ(m,mask=0):
 D={i if mask>>(i-1)&1 else m-i for i in range(1,(m+1)//2)}
 return [[int(i!=j and (j-i)%m in D) for j in range(m)] for i in range(m)]
def cycle(A,r,S=None):
 S=list(range(len(A))) if S is None else list(S)
 def dfs(p,left):
  if len(p)==r:return p if A[p[-1]][p[0]] else None
  for v in left:
   if A[p[-1]][v]:
    z=dfs(p+[v],[w for w in left if w!=v])
    if z:return z
 for v in S:
  z=dfs([v],[w for w in S if w!=v])
  if z:return z
 return None
def frames(m,ma,mb,cross):
 A=[[0]*(2*m) for _ in range(2*m)]
 C,D=circ(m,ma),circ(m,mb)
 for j in range(m):
  for k in range(m):
   A[2*j][2*k]=C[j][k];A[2*j+1][2*k+1]=D[j][k]
   bit=cross[(k-j)%m];A[2*j][2*k+1]=bit;A[2*k+1][2*j]=1-bit
 return A
def certificate(A,C,v):
 n=len(A); S=set(C); comp=set(range(n))-S
 direction='exit' if sum(A[x][v] for x in C)>=len(C)//2+1 else 'entry'
 B=A if direction=='exit' else [list(r) for r in zip(*A)]
 C=C if direction=='exit' else list(reversed(C))
 D=cycle(B,len(comp),comp);need(D is not None,'complement Hamilton cycle')
 i=D.index(v);W=D[i:]+D[:i]
 words=[]
 for i,x in enumerate(C):
  if B[x][v]:words.append(C[i+1:]+C[:i+1]+W)
 need(len(words)>=5,'five exits')
 for P in words:need(len(set(P))==n and all(B[x][y] for x,y in zip(P,P[1:])),'word')
 return {'source_rows':[''.join(map(str,r)) for r in A],'cycle':C,'terminal':v,'direction':direction,'words':words}
def audit_receipt(R):
 if type(R) is not dict or set(R)!={'source_rows','cycle','terminal','direction','words'}:return False
 rows=R['source_rows']
 if type(rows) is not list or len(rows)<16:return False
 n=len(rows);labels=set(range(n))
 if any(type(row) is not str or len(row)!=n or set(row)-{'0','1'} for row in rows):return False
 if type(R['direction']) is not str or R['direction'] not in ('exit','entry'):return False
 A=[[int(x) for x in row] for row in rows]
 if any(A[i][i] for i in range(n)) or any(A[i][j]+A[j][i]!=1 for i in range(n) for j in range(i+1,n)):return False
 C=R['cycle'];v=R['terminal'];words=R['words']
 if type(C) is not list or len(C)!=8 or any(type(x) is not int or x not in labels for x in C):return False
 S=set(C)
 if len(S)!=8 or type(v) is not int or v not in labels-S:return False
 if type(words) is not list or not 5<=len(words)<=8:return False
 if any(type(W) is not list or len(W)!=n or any(type(x) is not int or x not in labels for x in W) or set(W)!=labels for W in words):return False
 B=A if R['direction']=='exit' else [list(r) for r in zip(*A)]
 if not all(B[C[i]][C[(i+1)%8]] for i in range(8)) or not strong(A,labels-S):return False
 expected=[];suffix=None
 rotations={tuple(C[i:]+C[:i]) for i in range(8)}
 for W in words:
  if W[8]!=v or not all(B[x][y] for x,y in zip(W,W[1:])):return False
  if tuple(W[:8]) not in rotations:return False
  if suffix is None:suffix=W[8:]
  elif W[8:]!=suffix:return False
  expected.append(W[:8])
 return len({tuple(p) for p in expected})==len(expected)
def counts(A):
 n=len(A);perms=list(itertools.permutations(range(n)))
 aut=sum(all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(i+1,n)) for p in perms)
 H=B=0
 for p in perms:
  if all(A[x][y] for x,y in zip(p,p[1:])):
   H+=1
   B+=all(A[p[i]][p[j]]==A[p[n-1-j]][p[n-1-i]] for i in range(n) for j in range(i+1,n))
 return [H,aut,H//aut,B//aut]
start=datetime.datetime.now(datetime.timezone.utc).isoformat();t=time.monotonic();examples=[];framecounts={}
# Exhaust all two-shore m3/m5 marked cycle frames; complement deletion of
# any nonempty proper subset of first shore, not merely cycle subsets.
for m in (3,5):
 count=0
 for ma in range(1<<((m-1)//2)):
  for mb in range(1<<((m-1)//2)):
   for free in itertools.product((0,1),repeat=m-2):
    cross=[1]+list(free)+[0] # A_j->B_j->A_(j+1)
    A=frames(m,ma,mb,cross);need(strong(A),'frame strong')
    for mask in range(1,(1<<m)-1):
     S={2*j for j in range(m) if mask>>j&1}
     need(strong(A,set(range(2*m))-S),'proper shore deletion')
    count+=1
 framecounts[str(m)]=count
# m9: all16 first-shore circulants, one second-shore circulant,
# all128 compatible cross patterns. A complete stated subfamily, not allframes.
count=0
for ma in range(16):
 C0=cycle(circ(9,ma),8);need(C0 is not None,'eight cycle exists')
 C=[2*x for x in C0]
 for free in itertools.product((0,1),repeat=7):
  cross=[1]+list(free)+[0];A=frames(9,ma,15,cross)
  totals=[sum(A[x][v] for x in C) for v in range(1,18,2)]
  need(sum(totals)==8*sum(cross),'exact cross total')
  need(any(x!=4 for x in totals),'odd shore cannot balance')
  need(strong(A,set(range(18))-set(C)),'strong remainder')
  v=next(v for v,x in zip(range(1,18,2),totals) if x!=4)
  # Literal words on16 representative cases, plus a hostile.
  if count%128==0:
   R=certificate(A,C,v);need(audit_receipt(R),'receipt');examples.append(R)
  count+=1
framecounts['9_scoped']=count
bad=json.loads(json.dumps(examples[0]));bad['words'][0][0]=bad['words'][0][1]
need(not audit_receipt(bad),'duplicate word vertex rejected')
bad=json.loads(json.dumps(examples[0]));bad['cycle']*=2
need(not audit_receipt(bad),'doubled cycle rejected')
# Schema and Python-index alias hostiles: never allow values to acquire
# vertex meaning merely because indexing/equality accepts them.
for replacement in (-1, True, 17.0, 18):
 bad=json.loads(json.dumps(examples[0]))
 bad['words']=[[replacement if x==17 else x for x in W] for W in bad['words']]
 need(not audit_receipt(bad),'invalid literal Hamilton vertex rejected')
for replacement in (-1, True, 0.0, 18):
 bad=json.loads(json.dumps(examples[0]));bad['cycle'][0]=replacement
 need(not audit_receipt(bad),'invalid literal cycle vertex rejected')
for replacement in (-1, True, 1.0, 18, examples[0]['cycle'][0]):
 bad=json.loads(json.dumps(examples[0]));bad['terminal']=replacement
 need(not audit_receipt(bad),'invalid terminal rejected')
for direction in ('garbage',0,True,None):
 bad=json.loads(json.dumps(examples[0]));bad['direction']=direction
 need(not audit_receipt(bad),'invalid direction rejected')
for bad in (None,[],{}, {'source_rows': []}):need(not audit_receipt(bad),'invalid schema rejected')
bad=json.loads(json.dumps(examples[0]));bad['source_rows'][0]='2'+bad['source_rows'][0][1:]
need(not audit_receipt(bad),'nonbinary source rejected')
bad=json.loads(json.dumps(examples[0]));bad['words'][0][0]=[0]
need(not audit_receipt(bad),'nested label rejected without exception')
# Explicit q1 regular49 witness; no all-circulant census.
A=circ(49,(1<<24)-1);C=list(range(15))+[25]
need(all(A[C[i]][C[(i+1)%16]] for i in range(16)),'q1 sixteen cycle')
need(strong(A,set(range(49))-set(C)),'q1 remainder')
# Fresh small exact known positive/hostile controls, all Hamilton words and auts.
small={'C3':counts(circ(3,1)),'R5':counts(circ(5,3))}
A=[[int(i!=j and (j-i)%7 in {1,2,4}) for j in range(7)] for i in range(7)]
small['Paley7']=counts(A)
need(small['C3']==[3,3,1,1],'C3')
need(small['R5']==[15,5,3,3],'R5 pureblue hostile')
need(small['Paley7']==[189,21,9,3],'Paley7 deficit6 positive')
receipt={'status':'FINITE-EXACT scoped controls, not complete center classification','start_utc':start,'end_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_seconds':time.monotonic()-t,'checks':checks,'frame_counts':framecounts,'small_controls':small,'source_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'prefix_receipts':examples}
(ROOT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k!='prefix_receipts'},indent=2))
