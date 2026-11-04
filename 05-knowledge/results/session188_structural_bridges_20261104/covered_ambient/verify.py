"""Exact ordered-cover reduction checks. Imports frozen companion polynomial primitives."""
import importlib.util,itertools,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('covered',ROOT.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(sp);sp.loader.exec_module(v)
def product_num(G):
 def left(I):return v.add(v.R,*[G[u][w] for u in range(3) for w in I])
 def right(I):return v.add(v.R,*[G[u][w] for u in I for w in range(3)])
 return v.prod(left([0]),right([1,2]),left([1]),right([0,2]))
D2=v.add(product_num(v.G),v.scale(product_num([list(x) for x in zip(*v.G)]),-1))
v.need(v.add(*[v.permute(D2,pi) for pi in itertools.permutations(range(3))])=={},'separated-slot symmetrization')
v.need(v.add(D2,v.permute(D2,(1,0,2)))=={},'separated-slot exchange antisymmetry')
E=v.add(v.X,v.Y,v.Z,v.prod(v.X,v.Y),v.prod(v.Y,v.Z),v.prod(v.X,v.Z),v.scale(v.P,3))
qs=sorted({tuple(v.profile(v.tournament(n,c))) for n in range(1,6) for c in range(1<<(n*(n-1)//2))})
num={};p=v.ONE
for h in range(1,18):num[h]=v.mul(v.N,p);p=v.mul(p,E)
def contrast(q,h):
 z=h*v.functional(num[h],h+1,list(q));v.need(z%3==0,'exact coordinate integrality');return z//3
records=[]
for q in qs:
 ds=[0]+[contrast(q,h) for h in range(1,3*(len(q)-1)+3)]
 v.need(all(x<=0 for x in ds) and ds[1]<0,'finite actual-profile domination')
 records.append({'q':q,'T_minus_U_profile':ds})
literal=[]
for name,a in [('TT1',[0]),('TT2',[2,0]),('TT3',[6,4,0]),('C3',[2,4,1])]:
 m=len(a);mask=(1<<m)-1;q=v.profile(a);T=v.extend(v.blowup(a),mask,mask<<m);U=v.extend(v.blowup(a,True),mask,mask<<m);p=v.profile(T);r=v.profile(U)
 v.need([x-y for x,y in zip(p,r)]==[0]+[contrast(q,h) for h in range(1,len(T)+1)],'literal all-coordinate ordered-cover control')
 literal.append({'child':name,'q':q,'T_plus_rows':T,'U_plus_rows':U,'q_T_plus':p,'q_U_plus':r})
artificial=[0,9,0,1,1,1,0];opposite=contrast(artificial,2);v.need(opposite==146,'opposite-sign artificial obstruction')
print(json.dumps({'status':'PASS','scope':'Exact all-order ordered-cover formula; arbitrary-ambient domination is conditional, not proved. Artificial sequence is not a tournament profile.','distinct_actual_profiles_through_five':len(qs),'actual_profile_records':records,'literal_cover_controls':literal,'artificial_obstruction':{'q':artificial,'h':2,'T_minus_U':opposite,'h1_T_minus_U':contrast(artificial,1)},'symbolic_separated_slot_numerator':[[*e,c] for e,c in sorted(D2.items())],'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'verify.py',ROOT.parent/'covered_response'/'verify.py']}},sort_keys=True,indent=2))
