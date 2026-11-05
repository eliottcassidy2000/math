#!/usr/bin/env python3
"""Independent audit of the newly proposed anti-automorphism and fixed words."""
import json
from check_graph import graph,check

TAU=[7,13,12,11,10,9,8,0,6,5,4,3,2,1]

def main():
 rows=graph();n=len(rows);tau=TAU
 check(sorted(tau)==list(range(n)) and all(tau[tau[v]]==v and tau[v]!=v for v in range(n)),'fixed-point-free involution')
 check(all(((rows[u]>>v)&1)==((rows[tau[v]]>>tau[u])&1)
           for u in range(n) for v in range(n) if u!=v),'anti-automorphism')
 pairs=[[i,tau[i]] for i in range(7)];pair_of={v:i for i,P in enumerate(pairs) for v in P}
 final=[v for v in range(n) if rows[v]>>tau[v]&1]
 dp=[[0]*n for _ in range(128)];parent={}
 for v in range(n):dp[1<<pair_of[v]][v]=1
 for mask in range(1,128):
  for v in range(n):
   if not dp[mask][v]:continue
   for w in range(n):
    bit=1<<pair_of[w]
    if mask&bit or not rows[v]>>w&1:continue
    dp[mask|bit][w]+=dp[mask][v];parent.setdefault((mask|bit,w),(mask,v))
 F=sum(dp[127][v] for v in final)
 check(F==6111,'independent paired-subset DP fixed count')
 # A second independent count: inclusion-exclusion over omitted pair colors.
 IE=0
 for mask in range(1,128):
  verts=[v for v in range(n) if mask>>pair_of[v]&1]
  counts=[int(v in verts) for v in range(n)]
  edges=[(v,w) for v in verts for w in verts if rows[v]>>w&1]
  for _ in range(6):
   nxt=[0]*n
   for v,w in edges:nxt[w]+=counts[v]
   counts=nxt
  IE+=(-1 if (7-mask.bit_count())&1 else 1)*sum(counts[v] for v in final)
 check(IE==F,'pair-color walk inclusion-exclusion')
 last=next(v for v in final if dp[127][v]);mask=127;half=[]
 while True:
  half.append(last)
  if mask.bit_count()==1:break
  mask,last=parent[mask,last]
 half=half[::-1];word=half+[tau[v] for v in half[::-1]]
 check(sorted(word)==list(range(n)),'fixed word permutation')
 check(all(rows[v]>>w&1 for v,w in zip(word,word[1:])),'fixed word Hamilton arcs')
 check(word==[tau[v] for v in word[::-1]],'fixed word symmetry')
 rotations=[[7*(i//7)+(i%7+k)%7 for i in range(n)] for k in range(7)]
 rho=rotations[1];rho_inv=rotations[-1]
 check([tau[rho[tau[v]]] for v in range(n)]==rho_inv,'dihedral conjugation')
 # Every anti-map is a fixed automorphism times tau; the prior graph audit gave exactly C7.
 anti=[[rot[tau[v]] for v in range(n)] for rot in rotations]
 check(len({tuple(x) for x in anti})==7,'seven distinct anti-automorphisms')
 p=2163467;burnside=(p+F)//2
 check((p+F)%2==0 and burnside==1084789,'dihedral path-orbit count consequence')
 print(json.dumps({'status':'PASS','tau':tau,'pairs':pairs,'central_allowed_endpoints':final,
                   'endpoint_fixed_counts':{str(v):dp[127][v] for v in final},
                   'fixed_Hamilton_words':F,'independent_pair_color_walk_count':IE,
                   'explicit_fixed_word':word,'dihedral_path_orbit_count':burnside,
                   'scope':'Separate audit of proposed self-converse/fixed-word extension; original frozen packet unchanged.'},
                  indent=2,sort_keys=True))

if __name__=='__main__':main()
