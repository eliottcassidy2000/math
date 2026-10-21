#!/usr/bin/env python3
"""Literal insertion words for all marked q2 frames of shore sizes5 and7."""
import itertools,json,pathlib,hashlib,time,datetime
ROOT=pathlib.Path(__file__).resolve().parent;checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def circ(m,mask):
 D={i if mask>>(i-1)&1 else m-i for i in range(1,(m+1)//2)}
 return [[int(i!=j and (j-i)%m in D) for j in range(m)] for i in range(m)]
def frame(m,ma,mb,bits):
 A=[[0]*(2*m) for _ in range(2*m)];C=circ(m,ma);D=circ(m,mb)
 for j in range(m):
  for k in range(m):
   A[2*j][2*k]=C[j][k];A[2*j+1][2*k+1]=D[j][k]
   z=bits[(k-j)%m];A[2*j][2*k+1]=z;A[2*k+1][2*j]=1-z
 return A
def cycle(A,S):
 start=S[0]
 def rec(p,left):
  if not left:return p if A[p[-1]][start] else None
  for x in left:
   if A[p[-1]][x]:
    z=rec(p+[x],[y for y in left if y!=x])
    if z:return z
  return None
 return rec([start],S[1:])
def rotate(C,start):
 i=C.index(start);return C[i:]+C[:i]
def construct(A):
 n=len(A);m=n//2;L=list(range(0,n,2));R=list(range(1,n,2))
 if sum(A[L[0]][v] for v in R)<m/2:L,R=R,L
 v=R[0];C=cycle(A,L);D=cycle(A,R[1:]);need(C is not None and D is not None,'shore and punctured cycles')
 j=next(j for j in range(m) if A[C[j]][v] and A[v][C[(j+1)%m]])
 def append(P):
  w=next(w for w in R[1:] if A[P[-1]][w]);return P+rotate(D,w)
 words=[]
 for i in range(1,m-1):
  cut=C[(j-i+1)%m:]+C[:(j-i+1)%m]
  words.append(append(cut[:i]+[v]+cut[i:]))
 u=next(u for u in L if A[v][u]);words.append(append([v]+rotate(C,u)))
 patterns=[''.join('A' if x in L else 'B' for x in W) for W in words]
 for W in words:need(set(W)==set(range(n)) and len(W)==n and all(A[x][y] for x,y in zip(W,W[1:])),'literal Hamilton word')
 expected=['A'*i+'B'+'A'*(m-i)+'B'*(m-1) for i in range(1,m-1)]+['B'+'A'*m+'B'*(m-1)]
 need(patterns==expected,'claimed patterns')
 mates=[''.join('B' if x=='A' else 'A' for x in reversed(P)) for P in patterns]
 need(len(set(patterns+mates))==2*m-2,'all original and mate patterns distinct')
 bad='A'*(m-1)+'BA'+'B'*(m-1)
 need(bad==''.join('B' if x=='A' else 'A' for x in reversed(bad)),'omitted index hostile is fixed')
 return {'source_rows':[''.join(map(str,r)) for r in A],'high_cross_shore':L,'other_shore':R,'terminal':v,'words':words,'patterns':patterns,'mate_patterns':mates}
t=time.monotonic();start=datetime.datetime.now(datetime.timezone.utc).isoformat();examples=[];counts={};digest=hashlib.sha256()
for m in (5,7):
 count=0
 for ma in range(1<<((m-1)//2)):
  for mb in range(1<<((m-1)//2)):
   for free in itertools.product((0,1),repeat=m-2):
    bits=(1,)+free+(0,);A=frame(m,ma,mb,bits);R=construct(A);R['frame']=[m,ma,mb,list(bits)]
    digest.update((json.dumps(R,sort_keys=True,separators=(',',':'))+'\n').encode())
    if mb==0 and sum(free)==0:examples.append(R)
    count+=1
 counts[str(m)]=count
receipt={'status':'FINITE-EXACT literal construction controls; deficit consequence requires self-converseness','start_utc':start,'end_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_seconds':time.monotonic()-t,'checks':checks,'frame_counts':counts,'all_constructions_sha256':digest.hexdigest(),'source_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'sample_receipts':examples}
(ROOT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='sample_receipts'},indent=2))
