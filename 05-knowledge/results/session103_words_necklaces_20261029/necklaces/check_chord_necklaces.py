#!/usr/bin/env python3
import itertools,json,hashlib
from collections import Counter
CHECKS=0
def need(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(msg)
def key(A):return tuple(A[i][j] for i in range(len(A)) for j in range(i+1,len(A)))
def rotate(A,r):
 n=len(A);return [[A[(i+r)%n][(j+r)%n] for j in range(n)] for i in range(n)]
def canon(A):return min(key(rotate(A,r)) for r in range(len(A)))
def cycle_count_cap2(A):
 n=len(A);out=[sum(1<<j for j in range(n) if A[i][j]) for i in range(n)];full=(1<<n)-1
 def dfs(v,used):
  if used==full:return A[v][0]
  todo=out[v]&~used;ans=0
  while todo:
   bit=todo&-todo;todo-=bit;ans+=dfs(bit.bit_length()-1,used|bit)
   if ans>=2:return 2
  return ans
 return dfs(0,1)
def families(n):
 A=[[0]*n for _ in range(n)];fixed=set()
 for i in range(n):j=(i+1)%n;A[i][j]=1;A[j][i]=0;fixed.add(tuple(sorted((i,j))))
 free=[p for p in itertools.combinations(range(n),2) if p not in fixed]
 reps={};rooted=0
 for mask in range(1<<len(free)):
  for b,(i,j) in enumerate(free):A[i][j]=(mask>>b)&1;A[j][i]=1-A[i][j]
  if cycle_count_cap2(A)==1:
   rooted+=1;k=canon(A)
   if k not in reps:
    r=next(r for r in range(n) if key(rotate(A,r))==k);reps[k]=rotate(A,r)
 return reps,rooted,1<<len(free)
def observer(A):
 n=len(A);ports=list(itertools.permutations(range(n),2));ix={p:i for i,p in enumerate(ports)};M=[[0]*len(ports) for _ in ports]
 for p in itertools.permutations(range(n)):
  bad=[i for i in range(n-1) if not A[p[i]][p[i+1]]]
  if len(bad)>1:continue
  for k in range(1,n):
   if bad and bad[0]!=k-1:continue
   M[ix[p[0],p[k]]][ix[p[k-1],p[-1]]]+=1
 return M
def determinant(M):
 a=[r[:] for r in M];n=len(a);prev=1;sgn=1
 for k in range(n-1):
  p=next((i for i in range(k,n) if a[i][k]),None)
  if p is None:return 0
  if p!=k:a[p],a[k]=a[k],a[p];sgn=-sgn
  d=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    q=a[i][j]*d-a[i][k]*a[k][j]
    if q%prev:raise RuntimeError('Bareiss nonintegral')
    a[i][j]=q//prev
   a[i][k]=0
  prev=d
 return sgn*a[-1][-1]
def dual_restored(A):
 n=len(A);return [[A[(-j)%n][(-i)%n] for j in range(n)] for i in range(n)]
def extend(A,bits):
 n=len(A);B=[r+[0] for r in A]+[[0]*(n+1)]
 B[n][0]=1;B[0][n]=0;B[n-1][n]=1;B[n][n-1]=0
 for j in range(1,n-1):B[n][j]=(bits>>(j-1))&1;B[j][n]=1-B[n][j]
 return B

def main():
 allreps={};dets={};rows=[]
 for n,expected,units in [(3,1,1),(4,1,1),(5,3,3),(6,8,7),(7,21,16)]:
  reps,rooted,total=families(n);allreps[n]=reps;need(len(reps)==expected,'isotype census')
  hist=Counter();selfdual=0;orbit_sizes=Counter()
  for k,A in reps.items():
   M=observer(A);d=determinant(M);dets[n,k]=d;hist[d]+=1
   B=dual_restored(A);need(canon(B) in reps,'dual closure');selfdual+=canon(B)==k
   need(observer([[A[j][i] for j in range(n)] for i in range(n)])==[list(r) for r in zip(*M)],'observer transpose')
   orbit_sizes[len({key(rotate(A,r)) for r in range(n)})]+=1
   word=[[A[i][(i+d)%n] for d in range(1,n)] for i in range(n)]
   need(all(word[i][d-1]+word[(i+d)%n][n-d-1]==1 for i in range(n) for d in range(1,n)),'reciprocity')
  need(sum(v for d,v in hist.items() if abs(d)==1)==units,'unit census')
  rows.append({'n':n,'forced_cycle_candidates':total,'rooted_unique_cycle_tournaments':rooted,'isotypes':len(reps),'det_histogram':dict(sorted(hist.items())),'self_converse_isotypes':selfdual,'rotation_orbit_sizes':dict(sorted(orbit_sizes.items()))})
 candidates={}
 for k,A0 in allreps[5].items():
  for r in range(5):
   A=rotate(A0,r)
   for bits in range(8):
    B=extend(A,bits)
    if cycle_count_cap2(B)!=1:continue
    d=dets[6,canon(B)]
    candidates.setdefault(bits,[]).append((k,r,A,B,d))
 hostile=None
 for bits,cs in sorted(candidates.items()):
  for x in cs:
   for y in cs:
    if x[0]!=y[0] and x[-1]==0 and abs(y[-1])==1:
     hostile={'orientation_bits':bits,'rule':'append vertex5 after4 before0; bit(j-1)=1 means5->j forj1..3','singular_parent':x[2],'unimodular_parent':y[2],'parent_canonical_keys':[list(x[0]),list(y[0])],'parent_determinants':[dets[5,x[0]],dets[5,y[0]]],'singular_child':x[3],'unimodular_child':y[3],'child_determinants':[x[-1],y[-1]]};break
   if hostile:break
  if hostile:break
 need(hostile is not None,'same insertion hostile exists')
 need(sorted(map(sum,hostile['singular_parent']))==sorted(map(sum,hostile['unimodular_parent']))==[1,1,2,3,3],'same parent scores')
 need(sorted(map(sum,hostile['singular_child']))==sorted(map(sum,hostile['unimodular_child']))==[1,1,2,3,4,4],'same child scores')
 A=hostile['singular_parent'];n=len(A);word=[[A[i][(i+d)%n] for d in range(1,n)] for i in range(n)];plain=[word[(-i)%n] for i in range(n)]
 need(any(plain[i][d-1]+plain[(i+d)%n][n-d-1]!=1 for i in range(n) for d in range(1,n)),'plain reversal is not a tournament word')

 A=hostile['singular_child'];repair=[r+[0] for r in A]+[[0]*7]
 for j in range(6):repair[6][j]=int(j in (0,1,5));repair[j][6]=1-repair[6][j]
 need(cycle_count_cap2(repair)==1,'repair unique HC')
 repair_det=determinant(observer(repair));need(repair_det==-1,'repair unimodular')
 hostile['repair_child']=repair;hostile['repair_determinant']=repair_det
 hostile['repair_rule']='adjoin vertex6 with6->0,1,5 and2,3,4->6; existing6-vertex child retained'
 R5=[[int(i!=j and (j-i)%5 in (1,2)) for j in range(5)] for i in range(5)]
 other=[[R5[(2*i)%5][(2*j)%5] for j in range(5)] for i in range(5)]
 need(cycle_count_cap2(R5)==2 and canon(R5)!=canon(other),'nonunique cycle marking hostile')
 mark_hostile={'tournament':'cyclic R5 differences1,2','first_cycle':[0,1,2,3,4],'second_cycle':[0,2,4,1,3],'constant_chord_rows':['1100','1010'],'necklaces_distinct':True,'unmarked_tournament_same':True}

 print(json.dumps({'status':'FINITE-EXACT','checks':CHECKS,'census':rows,'hostile':hostile,'nonunique_marking_hostile':mark_hostile,'scope':'Faithful chord-row necklace model; no identification of missing b/e insertion semantics.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
