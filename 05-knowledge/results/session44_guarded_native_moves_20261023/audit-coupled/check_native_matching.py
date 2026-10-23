"""Independent source enumeration and literal-pin matching audit."""
import contextlib,io,itertools,json,runpy
from pathlib import Path
with contextlib.redirect_stdout(io.StringIO()):H=runpy.run_path(str(Path(__file__).with_name('check_native_packet.py')))
checks=0

def need(c,msg):
 global checks
 checks+=1
 if not c:raise RuntimeError(msg)
A=H['matrix'](int(H['nodes'][2242]['canon_a']));rows=[];allpaths=0
for P in itertools.permutations(range(8)):
 if not all(A[u][v] for u,v in zip(P,P[1:])):continue
 allpaths+=1;word=sum(A[P[i]][P[j]]<<b for b,(i,j) in enumerate(H['chords']))
 need(H['merged'][H['owner'][word]]==2242,'literal source')
 target=H['merged'][H['owner'][H['full']^word]]
 if target==3442:rows.append((P,word))
need(allpaths==407 and len(rows)==3,'complete source and target fibre')
expected={(1,0,5,3,2,4,6,7),(1,5,4,3,2,0,6,7),(4,1,5,3,2,0,6,7)}
need({P for P,w in rows}==expected,'exact actual path inventory')
for P,w in rows:
 need((H['x'][2242]-H['x'][3442])%8==5,'original electrical phase')
 K=H['kmatrix'](A,P)
 need(any(H['isomorphism'](K,H['matrix'](int(H['nodes'][3442][key]))) is not None for key in ('canon_a','canon_b')),'direct target isomorphism')
 need(P[3:5]==(3,2) and P[6:8]==(6,7),'fixed shared actual arcs')
pins=[((0,3,4,6,7),(1,3,2,6,7)),((3,4,5,6,7),(3,2,0,6,7)),((2,3,4,6,7),(5,3,2,6,7))]
eligible=[[P for P,w in rows if tuple(P[j] for j in clocks)==values] for clocks,values in pins]
need([len(s) for s in eligible]==[2,2,2],'literal pin eligibilities')
states=[a for a in itertools.product(*eligible) if len(set(a))==3]
need(len(states)==2,'exact injection fibre')
need(sum(a!=b for a,b in zip(*states))==3,'sharp three-slot support')
need(set(states[0])==set(states[1]),'global occupied resource set unchanged')
# Independently test the finite abelian obstruction for all local cosets in Z/3 and Z/6.
from math import gcd
systems=0
for m in (3,6):
 subgroups=[tuple(x for x in range(m) if x%d==0) for d in range(1,m+1) if m%d==0]
 for K1,K2 in itertools.product(subgroups,repeat=2):
  Dset={(a-b)%m for a,b in itertools.product(K1,K2)}
  ker=sum((a-b)%m==0 for a,b in itertools.product(K1,K2))
  for r,s,delta in itertools.product(range(m),repeat=3):
   actual=sum((r+a-s-b)%m==delta for a,b in itertools.product(K1,K2))
   need(actual==(ker if (delta-r+s)%m in Dset else 0),'coset obstruction and retained joint lifts')
   systems+=1
print(json.dumps({'status':'PASS','checks':checks,'source_paths':allpaths,'target_paths':[P for P,w in rows],'eligible':eligible,'states':states,'minimum_slot_budget':3,'finite_abelian_systems':systems},indent=2))
