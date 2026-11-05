#!/usr/bin/env python3
"""Exact all-subset cycle checks for the infinite serial-seam obstruction family."""
from pathlib import Path
from itertools import combinations
import json
from verify import require,count_paths,count_cycles
from verify_rooted import converse,profile,glue,induced

def all_subset_cycles(rows):
 n=len(rows);N=1<<n;out=[0]*N
 for root in range(n):
  dp=[[0]*n for _ in range(N)];dp[1<<root][root]=1
  for mask in range(1<<root,N,1<<(root+1)):
   # Exactly the masks whose least set bit is root.
   for u,value in enumerate(dp[mask]):
    if not value:continue
    if mask.bit_count()>=3 and rows[u]>>root&1:out[mask]+=value
    todo=rows[u]&~mask&~((1<<root)-1)
    while todo:
     bit=todo&-todo;todo-=bit;dp[mask|bit][bit.bit_length()-1]+=value
 return out

def main():
 A=[2,8,3,5];a=2;B=converse(A);b=2
 require(profile(A,a)[:6]==(5,1,9,2,1,3),'left endpoint profile')
 require(profile(B,b)[:6]==(5,2,6,1,1,3),'right endpoint profile')
 # Literal small-root property used by the structural proof.
 for rows,root,wanted in [(A,a,2),(B,b,2)]:
  for mask in range(1<<4):
   if not mask>>root&1 or mask==15:continue
   vs=[i for i in range(4) if mask>>i&1];P=profile(induced(rows,vs),vs.index(root))
   require(P[3 if rows==A else 1]<=1,'proper endpoint response')
 left=A;mark=a;reports=[];hseq=[];subset_checks=0
 for k in range(6):
  T,r=glue(left,mark,B,b);n=len(T);cycles=all_subset_cycles(T)
  h=count_paths(T);hseq.append(h)
  require(cycles[-1]==count_cycles(T)==4,'full cycle count')
  require(all(c<=2 for c in cycles[:-1]),'proper induced obstruction failure')
  subset_checks+=len(cycles)-1
  # All non-singleton block roots form the backward spine.
  spine=[2]+list(range(4,4+k))+[4+k+2]
  for x,y in zip(spine,spine[1:]):require(T[y]>>x&1,'spine direction')
  triangles=[[0,1,3],[4+k,5+k,7+k]]+[list(range(4+j,7+j)) for j in range(0,k-2,3)]
  require(len(set(sum(triangles,[])))==3*len(triangles),'triangle packing disjointness')
  require(all(count_cycles(induced(T,t))==1 for t in triangles),'triangle packing cycles')
  order=[3,0,1,2]+list(range(4,4+k))+[4+k+j for j in [2,1,0,3]]
  ordered=induced(T,order)
  require(max(j-i for j in range(n) for i in range(j) if ordered[j]>>i&1)==2,'span-two order')
  reports.append({'span_two_order':order,'minimum_backward_span':2,'k':k,'order':n,'rows':T,'root':r,'H':h,'C':4,'maximum_proper_induced_C':max(cycles[:-1]),'proper_C_histogram':{str(c):cycles[:-1].count(c) for c in range(3)},'feedback_vertex_lower_bound':len(triangles)})
  left,mark=glue(left,mark,[0],0)
 require(hseq==[105,187,349,641,1177,2167],'path initial values')
 require(all(hseq[k]==sum(hseq[k-3:k]) for k in range(3,len(hseq))),'tribonacci recurrence checks')
 S=[[1,-1,1],[1,-1,0],[1,0,1]]
 def mm(X,Y):return [[sum(X[i][k]*Y[k][j] for k in range(3)) for j in range(3)]for i in range(3)]
 S2=mm(S,S);S3=mm(S2,S)
 require(S3==[[S2[i][j]+S[i][j]+int(i==j)for j in range(3)]for i in range(3)],'matrix recurrence certificate')
 payload={'left_rows':A,'left_root':a,'right_rows':B,'right_root':b,'all_proper_subset_checks':subset_checks,'family':reports,'path_recurrence':'h[k+3]=h[k+2]+h[k+1]+h[k]'}
 text=json.dumps(payload,sort_keys=True,indent=2)+'\n';(Path(__file__).resolve().parent/'serial_obstructions_certificate.json').write_text(text);print(text,end='')
if __name__=='__main__':main()
