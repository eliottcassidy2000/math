from itertools import combinations,permutations
from collections import Counter
from functools import lru_cache
from pathlib import Path
import json

def req(v,msg):
 if not v:raise RuntimeError(msg)
V=range(7);edges=list(combinations(V,2));idx={e:i for i,e in enumerate(edges)}
base=sum(1<<i for i,(a,b) in enumerate(edges) if (b-a)%7 in (1,2,4))
def arc(m,a,b):return bool(m>>idx[min(a,b),max(a,b)]&1)==(a<b)
def toggle(C):return sum(1<<idx[a,b] for a,b in combinations(sorted(C),2))
@lru_cache(None)
def hp(m):
 # Independent subset dynamic programming, rather than permutation count.
 d={(1<<a,a):1 for a in V}
 for mask in range(1,128):
  for a in V:
   n=d.get((mask,a),0)
   for b in V:
    if n and not mask>>b&1 and arc(m,a,b):d[mask|1<<b,b]=d.get((mask|1<<b,b),0)+n
 return sum(d.get((127,a),0) for a in V)
Ls=[sorted({tuple(sorted((t+s*a)%7 for a in (0,1,3))) for t in V}) for s in (1,-1)]
cycles=[C for C in combinations(V,3) if all(sum(arc(base,a,b) for b in C if b!=a)==1 for a in C)]
req(len(cycles)==14 and set(cycles)==set(Ls[0])|set(Ls[1]) and not set(Ls[0])&set(Ls[1]),'two plane partition')
# All one-line flips are regular, H171.
req(hp(base)==189,'base H')
for C in cycles:
 m=base^toggle(C);req(hp(m)==171,'triangle H171')
 req(all(sum(arc(m,a,b) for b in V if b!=a)==3 for a in V),'triangle regular')
# Entire fixed-plane cube: locate every H189 mask, independently.
cube_summary=[]
for L in Ls:
 dist=Counter();pal=[]
 for x in range(128):
  m=base
  for j,C in enumerate(L):
   if x>>j&1:m^=toggle(C)
  dist[hp(m)]+=1
  if hp(m)==189:pal.append(x)
 pencils=[sum(1<<j for j,C in enumerate(L) if a in C) for a in V]
 req(set(pal)=={0,127,*pencils,*(127^x for x in pencils)},'Paley mask structural classification')
 cube_summary.append({'H_distribution':dict(dist),'H189_masks':pal})
profiles=Counter();moves=Counter();reps=[];cube_checks=0;nonzero=0
for P in permutations(V):
 if not all(arc(base,a,b) for a,b in zip(P,P[1:])):continue
 pe={tuple(sorted(e)) for e in zip(P,P[1:])}
 U=[[j for j,C in enumerate(L) if not pe & set(combinations(C,2))] for L in Ls]
 blue=all(arc(base,P[a],P[b])==arc(base,P[6-b],P[6-a]) for a,b in combinations(V,2))
 ds=[(b-a)%7 for a,b in zip(P,P[1:])]
 req(blue==(ds==ds[::-1]),'palindrome')
 u,v=map(len,U);req(blue==(u==v),'balance');profiles[u,v,blue]+=1
 adm=[C for C in cycles if not pe & set(combinations(C,2))]
 req(len(adm)==u+v,'all admissible');moves[blue,len(adm)]+=1
 # Consecutive step ratio 2 marks one family, ratio 4 the other.
 ratio=Counter(ds[j+1]*pow(ds[j],-1,7)%7 for j in range(5))
 req((u,v)==(1+ratio[2],1+ratio[4]),'ratio count')
 if P[:2]==(0,1):reps.append({'path':P,'increments':ds,'ratios_2_4':[ratio[2],ratio[4]],'u+':u,'u-':v,'blue':blue})
 for L,J in zip(Ls,U):
  for bits in range(1<<len(J)):
   m=base
   for k,j in enumerate(J):
    if bits>>k&1:m^=toggle(L[j])
   req(all(arc(m,a,b) for a,b in zip(P,P[1:])),'path survives')
   req(all(sum(arc(m,a,b) for b in V if b!=a)==3 for a in V),'cube regular')
   req(hp(m)==(189 if bits==0 else 171),'cube H');cube_checks+=1;nonzero+=bool(bits)
req(len(reps)==9 and sum(profiles.values())==189,'orbit coverage')
autos=[p for p in permutations(V) if all(arc(base,a,b)==arc(base,p[a],p[b]) for a,b in edges)]
req(len(autos)==21,'automorphism count')
for r in reps:req(len({tuple(p[a] for a in r['path']) for p in autos})==21,'free orbit')
outs=[98,84,41,19,37,74,28]
hostile=sum(1<<i for i,(a,b) in enumerate(edges) if outs[a]>>b&1)
req(all(sum(arc(hostile,a,b) for b in V if b!=a)==3 for a in V) and hp(hostile)==171,'hostile regular H171')
hostile_rows=[]
for P,want in [((0,1,2,3,4,5,6),True),((0,1,4,2,5,6,3),False)]:
 req(all(arc(hostile,a,b) for a,b in zip(P,P[1:])),'hostile HP')
 blue=all(arc(hostile,P[a],P[b])==arc(hostile,P[6-b],P[6-a]) for a,b in combinations(V,2))
 pe={tuple(sorted(e)) for e in zip(P,P[1:])}
 cs=[C for C in combinations(V,3) if all(sum(arc(hostile,a,b) for b in C if b!=a)==1 for a in C) and not pe & set(combinations(C,2))]
 req(blue==want and len(cs)==5,'hostile blue black same moves')
 hostile_rows.append({'path':P,'blue':blue,'admissible_triangle_count':len(cs)})
out={'hostile_regular_H171':hostile_rows,'status':'PASS','profiles':[{'u+':u,'u-':v,'blue':b,'count':n} for (u,v,b),n in sorted(profiles.items())],'admissible_moves':[{'blue':b,'moves':k,'paths':n} for (b,k),n in sorted(moves.items())],'full_plane_cubes':cube_summary,'path_cube_states_tested':cube_checks,'nonzero_path_cube_states_tested':nonzero,'nine_orbit_certificate':reps}
Path(__file__).with_name('fano_paths_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
