import sys,random
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from verify import *
from math import gcd
from functools import reduce

def cycle(n):return [[0<(j-i)%n<=(n-1)//2 for j in range(n)] for i in range(n)]
def data(a,mask):
 n=len(a);full=(1<<n)-1;paths=[]
 for s in range(full+1):
  ids=[i for i in range(n) if s>>i&1]
  paths.append([p for p in permutations(ids) if all(a[u][v] for u,v in zip(p,p[1:]))])
 H=len(paths[-1]);Q=sum(len(paths[s])*len(paths[full^s]) for s in range(full+1));M=[[0]*6 for _ in range(6)]
 for s in range(full+1):
  for P in paths[s]:
   for R in paths[full^s]:
    for c,(last,begun) in enumerate(states):
     if R and not begun and not(mask>>R[0]&1):continue
     new=(2 if mask>>P[-1]&1 else 1) if P else last
     row=states.index((new,int(begun or bool(R))));M[row][c]+=1
 return H,Q,M

def bound(gs):
 if len(gs)<=2:return 1
 return reduce(gcd,(reduce(lambda a,b:a*b,(g for k,g in enumerate(gs) if k not in (i,j)),1) for i in range(len(gs)) for j in range(i+1,len(gs))))
def compose(items):
 n=1+sum(len(a) for a,mask in items);out=[[False]*n for _ in range(n)];off=1
 for a,mask in items:
  for i in range(len(a)):
   u=off+i;out[0][u]=bool(mask>>i&1);out[u][0]=not out[0][u]
   for j in range(len(a)):out[u][off+j]=a[i][j]
   for q in range(off+len(a),n):out[u][q]=True
  off+=len(a)
 return out

def run(items,cache):
 x=v;gs=[]
 for index in items:
  a,mask,H,Q,M=cache[index];x=mv(M,x);gs.append(gcd(H,Q))
 D=sum(a*b for a,b in zip(ell,x));ck(D%bound(gs)==0,'general gcd divisor');return D

def main():
 blocks=[[[False]],cycle(3),graph([0,1],0),cycle(5),cycle(7)];cache=[];signatures=[]
 for a in blocks:
  n=len(a)
  masks=range(1<<n) if n<=5 else (0,1,21,127)
  for mask in masks:
   H,Q,M=data(a,mask);cache.append((a,mask,H,Q,M))
  signatures.append([n,H,Q,gcd(H,Q)])
 rng=random.Random(228);checks=literal=nontrivial=0
 for z in range(6000):
  inds=[rng.randrange(len(cache)) for _ in range(1+z%8)]
  D=run(inds,cache);checks+=1
  gs=[gcd(cache[i][2],cache[i][3]) for i in inds];nontrivial+=bound(gs)>1
  if sum(len(cache[i][0]) for i in inds)<=10 and literal<100:
   ck(profile(compose([(cache[i][0],cache[i][1]) for i in inds]))[0]==D,'literal general profile');literal+=1
 # Exhaustive all length<=5 words over singleton bit0/1 and triangle outdegree0..3.
 menu=[]
 for a in blocks[:2]:
  for k in range(len(a)+1):
   mask=(1<<k)-1
   menu.append(next(i for i,x in enumerate(cache) if x[0]==a and x[1]==mask))
 extra=[]
 for specs in [((3,1),(3,2),(3,4)),((3,1),(5,3),(5,5)),((5,1),(5,3),(5,5))]:
  inds=[next(i for i,x in enumerate(cache) if len(x[0])==n and x[1]==mask) for n,mask in specs]
  D=run(inds,cache);K=bound([gcd(cache[i][2],cache[i][3]) for i in inds])
  ck(profile(compose([(cache[i][0],cache[i][1]) for i in inds]))[0]==D,'nontrivial literal gcd control')
  extra.append({'blocks':specs,'D':D,'K':K})
 tri_checks=0
 for m in range(1,6):
  for inds in product(menu,repeat=m):
   D=run(inds,cache);t=sum(len(cache[i][0])==3 for i in inds)
   ck(D%3**max(t-2,0)==0,'triangle bound');tri_checks+=1
 c0=next(i for i,x in enumerate(cache) if len(x[0])==3 and x[1]==0)
 c1=next(i for i,x in enumerate(cache) if len(x[0])==3 and x[1]==1)
 for t in range(2,41):ck(run([c0]*(t-2)+[c1,c1],cache)==3**(t-2)*41,'sharp triangle construction')
 print(json.dumps({'status':'PASS','block_signatures_n_H_Q_g':signatures,'nontrivial_literal_controls':extra,'general_rooted_block_words':checks,'nontrivial_gcd_bound_cases':nontrivial,'literal_general_counts':literal,'exhaustive_triangle_singleton_words_lengths_1_to_5':tri_checks,'sharp_triangle_t_range':[2,40]},sort_keys=True,indent=2))
if __name__=='__main__':main()
