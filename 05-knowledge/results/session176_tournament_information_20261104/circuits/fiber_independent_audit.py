#!/usr/bin/env python3
import itertools,json

def need(p,s):
 if not p:raise RuntimeError(s)
def cycle(reverse=False):
 a=[[0,1,0],[0,0,1],[1,0,0]]
 return list(map(list,zip(*a))) if reverse else a

def concat(a,b):
 n=len(a);m=len(b)
 return [r+[1]*m for r in a]+[[0]*n+r for r in b]
def literal_k(a):
 n=len(a);out=[[0]*n for _ in range(n)]
 for p in itertools.permutations(range(n)):
  bad=[i for i in range(n-1) if not a[p[i]][p[i+1]]]
  for split in range(n-1):
   if all(i==split for i in bad):out[p[split]][p[split+1]]+=1
 return out

def paths(a):
 n=len(a)
 return sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in itertools.permutations(range(n)))
def extend(a,b):
 n=len(a)
 return [r+[1-b[i]] for i,r in enumerate(a)]+[list(b)+[0]]
def main():
 checked=0
 for k in (1,2):
  signatures=[]
  for bs in itertools.product((0,1),repeat=k):
   a=[]
   for b in bs:a=concat(a,cycle(b))
   K=literal_k(a);n=3*k;ws=[];degrees=[0]*n
   for u in range(n):
    for v in range(u+1,n):
     i,j=u//3,v//3
     if i==j:need(K[u][v]+K[v][u]==2*3**(k-1),'internal')
     else:
      t=j-i-1
      need(K[u][v]==(3**(k-2) if t==0 else 0),'forward')
      need(K[v][u]==3**k*4**t,'reverse')
     w=(K[u][v]+K[v][u])//2;ws.append(w)
     if a[u][v]:degrees[u]+=w
     else:degrees[v]+=w
   signatures.append((paths(a),ws,degrees));checked+=1
  need(all(s==signatures[0] for s in signatures),'fiber discrepancy')
 # Two added vertices can separate the seed opposite triangles even though
 # every single-vertex insertion response agrees. No producer code is used.
 a=cycle();b=cycle(True);witness=None;tests=0
 for x in itertools.product((0,1),repeat=3):
  aa=extend(a,x);bb=extend(b,x)
  need(paths(aa)==paths(bb),'single insertion should collide')
  for y in itertools.product((0,1),repeat=4):
   ha=paths(extend(aa,y));hb=paths(extend(bb,y));tests+=1
   if ha!=hb and witness is None:witness={'first_apex_outbits':x,'second_apex_outbits_to_base_then_first_apex':y,'counts':[ha,hb]}
 need(witness is not None,'expected continuation separation')
 print(json.dumps({'status':'PASS','literal_cover_fiber_tournaments':checked,'two_step_extension_pairs':tests,'two_step_separation':witness,'scope':'fiber bases are labelled and all isomorphic as unlabelled tournaments'},sort_keys=True,indent=2))
if __name__=='__main__':main()
