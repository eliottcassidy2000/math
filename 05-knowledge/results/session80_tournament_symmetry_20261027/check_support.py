#!/usr/bin/env python3
"""Standalone exact support/depth checks with direct modules and automorphisms."""
def tournament(n,mask):
 a=[[0]*n for _ in range(n)]; b=0
 for i in range(n):
  for j in range(i+1,n):
   a[i][j]=(mask>>b)&1; a[j][i]=1-a[i][j]; b+=1
 return a

def auts(a):
 n=len(a); ds=list(map(sum,a)); order=sorted(range(n),key=lambda v:ds.count(ds[v])); p=[-1]*n; used=set(); out=[]
 def rec(k):
  if k==n: out.append(tuple(p)); return
  v=order[k]
  for w in range(n):
   if w in used or ds[v]!=ds[w]: continue
   if any(a[v][u]!=a[w][p[u]] for u in order[:k]): continue
   p[v]=w; used.add(w); rec(k+1); used.remove(w); p[v]=-1
 rec(0); return out

def induced(a,B): return [[a[i][j] for j in B] for i in B]

def canonical(a):
 n=len(a); full=(1<<n)-1; mods=[]
 for s in range(1,full+1):
  B=[i for i in range(n) if s>>i&1]
  if all(len({a[v][u] for u in B})==1 for v in range(n) if not(s>>v&1)): mods.append(s)
 strong=[s for s in mods if not any((s&t) and (s&~t) and (t&~s) for t in mods)]
 maximal=[s for s in strong if s!=full and not any(s!=t and t!=full and s&t==s for t in strong)]
 return [[i for i in range(n) if s>>i&1] for s in maximal]

def inflate(q,children):
 tags=[(i,v) for i,c in enumerate(children) for v in range(len(c))]
 return [[children[i][v][w] if i==j else q[i][j] for j,w in tags] for i,v in tags]

from collections import Counter
import json

def depth(a):
 if len(a)<2: return 0
 B=canonical(a); q=[[a[b[0]][c[0]] for c in B] for b in B]
 # A tournament is transitive exactly when its scores are 0,...,k-1.
 prime=sorted(map(sum,q))!=list(range(len(q)))
 return int(prime)+max(depth(induced(a,b)) for b in B)

def verify(a):
 n=len(a); h=depth(a); G=auts(a); support=sum(any(p[i]!=i for p in G) for i in range(n)); r=n-2*h-1
 if r<0 or support>3*(r+1) or (h>=1 and support>(3*r)//2+3): raise RuntimeError((n,h,support,r))
 return n,h,r,support,len(G)

def main():
 histogram=Counter(); total=0
 for n in range(1,7):
  for mask in range(1<<(n*(n-1)//2)):
   row=verify(tournament(n,mask)); histogram[row]+=1; total+=1
 C=[[0,1,0],[0,0,1],[1,0,0]]; a=[[0]]; spines=[]
 for h in range(1,9):
  a=inflate(C,[a,[[0]],[[0]]]); row=verify(a)
  if row!=(2*h+1,h,0,3,3): raise RuntimeError(('spine',row))
  spines.append(row)
 sharp=[]
 for k in range(4):
  a=C
  for _ in range(k): a=inflate(C,[a,C,[[0]]])
  for extra in range(3):
   row=verify(a); h=k+1+extra
   if row!=(2*h+1+2*k,h,2*k,3*k+3,3**(k+1)): raise RuntimeError(('even sharp',row,k,extra))
   sharp.append(row); a=inflate(C,[a,[[0]],[[0]]])
 print(json.dumps({'sharp_even_excess':sharp,'all_labelled_n1_to6':total,'columns':['n','h','r','support','aut_order'],'histogram':[{ 'parameters':k,'count':v} for k,v in sorted(histogram.items())],'tight_spines':spines},indent=2))
if __name__=='__main__': main()
