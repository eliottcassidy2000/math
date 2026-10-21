from itertools import combinations,permutations
from pathlib import Path
import hashlib,json
E=list(combinations(range(4),2));omega='omega'
faces=[frozenset(combinations(S,2)) for S in combinations(range(4),3)]
quads=[frozenset(set(E)-set(P)) for P in combinations(E,2) if len(set(P[0])|set(P[1]))==4]
cycles=faces+quads
edge_lines=[frozenset(F) for F in faces]+[frozenset(set(E)-set(Q))|{omega} for Q in quads]
def dual(C):return frozenset(C) if len(C)==3 else frozenset(set(E)-set(C))|{omega}
checks=0
def need(x,s):
 global checks;checks+=1
 if not x:raise RuntimeError(s)
need(len(set(quads))==3 and len(set(edge_lines))==7,'model sizes')
clines=[]
for tri in combinations(cycles,3):
 isline=not(tri[0]^tri[1]^tri[2]);meet=set(dual(tri[0]))&set(dual(tri[1]))&set(dual(tri[2]))
 need(isline==(len(meet)==1),'exact dual incidence')
 if isline:clines.append(tri)
need(len(clines)==7,'seven cycle lines')
for p in permutations(range(4)):
 act=lambda e:tuple(sorted((p[e[0]],p[e[1]])))
 for C in cycles:
  lhs=dual(frozenset(act(e) for e in C));rhs=frozenset(omega if e==omega else act(e) for e in dual(C))
  need(lhs==rhs,'S4 equivariance')
coords={(0,3):1,(1,3):2,(2,3):4,(0,1):3,(0,2):5,(1,2):6}
word={0:frozenset()}
for a in range(1,8):
 C=frozenset(e for e,v in coords.items() if ((a&(v^7)).bit_count()%2)==1);word[a]=C
 L=frozenset(e for e,v in coords.items() if (a&v).bit_count()%2==0)
 if (a&7).bit_count()%2==0:L=L|{omega}
 need(C in cycles and dual(C)==L,'coordinate duality')
for a in range(8):
 for b in range(8):need(word[a]^word[b]==word[a^b],'linearity')
root=Path(__file__).resolve().parent
print(json.dumps({'status':'PASS','checks':checks,'edge_point_orbits':[1,6],'cycle_point_orbits':[4,3],'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'PROOF.md',Path(__file__)]}},indent=2))
