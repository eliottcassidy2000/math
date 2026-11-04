#!/usr/bin/env python3
# Independent implementation: no imports from check.py.
import itertools,json,math

def need(b,m):
 if not b:raise RuntimeError(m)

def words(g,vertices=None):
 if vertices is None:vertices=range(len(g))
 return sum(all((g[w[i]]>>w[i+1])&1 for i in range(len(w)-1)) for w in itertools.permutations(vertices))

def covers(g):
 n=len(g);q=[0]*n
 for w in itertools.permutations(range(n)):
  for cuts in range(1<<(n-1)):
   if all(cuts>>i&1 or g[w[i]]>>w[i+1]&1 for i in range(n-1)):q[cuts.bit_count()]+=1
 return q

def attach(g):
 n=len(g);return [x|1<<n for x in g]+[1<<(n+1),(1<<n)-1]

def skeleton(r):
 m=2*r;out=[]
 for k in range(1,r+2):
  total=0
  for w in itertools.permutations(range(m)):
   # Distinct gaps of exterior word contain the indistinguishable Xs.
   for gaps in itertools.combinations(range(m+1),k):
    seq=[]
    for j in range(m+1):
     if j in gaps:seq.append(-1)
     if j<m:seq.append(w[j])
    def arc(a,b):
     if a==-1:return b%2==0
     if b==-1:return a%2==1
     return a<b
    total+=all(arc(a,b) for a,b in zip(seq,seq[1:]))
  out.append(total)
 return out

def triangle_sum(n,j):
 blocks=[[3*i,3*i+1,3*i+2] for i in range(j)]+[[v] for v in range(3*j,n)]
 g=[0]*n
 for ii,B in enumerate(blocks):
  for C in blocks[ii+1:]:
   for a in B:
    for b in C:g[a]|=1<<b
  if len(B)==3:
   for a,b in zip(B,B[1:]+B[:1]):g[a]|=1<<b
 return g

def color_count(g,t):
 n=len(g);h={s:words(g,[v for v in range(n) if s>>v&1]) for s in range(1<<n)};total=0
 for coloring in itertools.product(range(t),repeat=n):
  masks=[0]*t
  for v,c in enumerate(coloring):masks[c]|=1<<v
  total+=math.prod(h[s] for s in masks)
 return total

def main():
 A=[16,9,3,5,14];B=[24,17,3,6,12]
 need(covers(A)==[15,78,198,240,120],'profile A')
 need(covers(B)==[15,90,210,240,120],'profile B')
 need(words(attach(A))==123 and words(attach(B))==135,'probe collision')
 need(words([2,4,1])==3 and words([0,5,1])==1,'nonuniform witness')
 K=[skeleton(r) for r in range(4)]
 need(K==[[1],[3,1],[9,10,2],[31,83,50,6]],'literal skeleton')
 checks=0
 for n in range(1,7):
  for j in range(n//3+1):
   g=triangle_sum(n,j)
   for t in range(4):
    need(color_count(g,t)==t**(n-3*j)*(t**3+2*t)**j,'actual triangle sum');checks+=1
 print(json.dumps({'status':'PASS','independent_literal_profile_examples':2,'independent_seven_vertex_H_examples':[123,135],'literal_skeleton_rows':K,'actual_color_assignment_checks':checks},indent=2))
if __name__=='__main__':main()
