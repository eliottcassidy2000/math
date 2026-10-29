from itertools import combinations
from collections import Counter
import json

def need(x,msg):
 if not x:raise RuntimeError(msg)
lines=[frozenset(c) for c in combinations(range(1,8),3) if c[0]^c[1]^c[2]==0]
faces=[frozenset(i+1 for i in range(7) if m>>i&1) for m in range(128)]
faces=[S for S in faces if not any(L<=S for L in lines)]
f=Counter(map(len,faces));need(len(lines)==7,'lines');need(len(faces)==64,'faces')
need(dict(f)=={0:1,1:7,2:21,3:28,4:7},'f-vector')
caps=[S for S in faces if len(S)==4]
need(set(caps)=={frozenset(range(1,8))-L for L in lines},'caps are line complements')
simplex={frozenset()}|set(caps)
need(all(a^b in simplex for a in simplex for b in simplex),'simplex code')
need(any(a^b not in set(faces) for a in faces for b in faces),'face family is not linear code')
need(3**7-len(faces)==2123,'rank formula')
print(json.dumps({'Fano_lines':len(lines),'allowed_active_subsets':len(faces),'face_sizes':dict(sorted(f.items())),'rank_by_audited_hypergraph_theorem':2123,'four_faces_plus_zero_simplex_code_size':8}))
