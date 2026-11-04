"""Independent physical reconstruction and predecessor-DP ambient audit."""
from pathlib import Path
import itertools,json,hashlib,importlib.util
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'even_cyclic_ambient'
s=importlib.util.spec_from_file_location('covered',ROOT.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
data=json.loads((SRC/'normal.json').read_text());v.need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer modes match')
def module(a,M):return all(a[u]&M in (0,M) for u in range(len(a)) if not M>>u&1)
def flip(a,fs):
 r=list(a)
 for M in fs:
  for u,w in itertools.combinations(range(len(a)),2):
   if M>>u&1 and M>>w&1:r[u]^=1<<w;r[w]^=1<<u
 return r
def blow(Q,blocks):
 sizes=list(map(len,blocks));off=[sum(sizes[:i]) for i in range(len(sizes))];out=[0]*sum(sizes)
 for i,A in enumerate(blocks):
  for u,row in enumerate(A):
   out[off[i]+u]=row<<off[i]
   for j,n in enumerate(sizes):
    if Q[i]>>j&1:out[off[i]+u]|=((1<<n)-1)<<off[j]
 return out
def put(Q,r,A):return blow(Q,[A if i==r else [0] for i in range(len(Q))])
counts={'families':0,'scalar_equal':0,'scalar_strict':0};checks=[]
for record in data['cases']:
 C,Q,r=record['C'],record['Q'],record['root'];M=blow([2,4,1],[[2,0],[2,0],C]);full=(1<<len(M))-1;mods=[m for m in range(1,full) if m.bit_count()>1 and module(M,m)];v.need(mods==record['local_original_modules'],'independent original module enumeration');T=put(Q,r,M);mm=full<<r
 E=sum(1<<u for u in range(len(T)) if not mm>>u&1 and T[u]>>r&1);query=[E|(3<<r),E|(12<<r)];v.need(query==record['query'],'physical root-incidence lift');base=[v.hp(T),v.hp(v.extend(T,*query))];v.need(base==record['baseline'],'literal source counts');n_equal=n_strict=0
 for bits,expected in enumerate(record['all_values_by_family_bits']):
  selected=[full]+[m for i,m in enumerate(mods) if bits>>i&1];local=flip(M,selected);U=put(Q,r,local);actual=[v.hp(U),v.hp(v.extend(U,*query))];v.need(actual==expected,'independent target literal counts');v.need(module(v.extend(T,*query),mm|(3<<len(T))) and module(v.extend(U,*query),mm|(3<<len(U))),'augmented port core stays a module')
  if actual[0]==base[0]:n_equal+=1;v.need((base[1]-actual[1])%8==4,'scalar-equal residue')
  else:n_strict+=1;v.need(actual[0]>base[0],'scalar-strict case')
 v.need((n_equal,n_strict)==(record['scalar_equal'],record['scalar_strict']),'case counts');counts['families']+=n_equal+n_strict;counts['scalar_equal']+=n_equal;counts['scalar_strict']+=n_strict;checks.append([len(T),r,n_equal,n_strict])
r=data['multiple_highest_prime_control'];T,V,U=r['T'],r['V'],r['U'];mods=[m for m in range(1,1<<len(T)) if module(T,m)];strong=[m for m in mods if all(m&n in (0,m,n) for n in mods)];changed=[]
for M in strong:
 if M.bit_count()<2:continue
 children=[n for n in strong if n!=M and n&M==n and not any(n!=k and k!=M and n&k==n and k&M==k for k in strong)];reps=[(n&-n).bit_length()-1 for n in children]
 if any(((T[a]>>b)^(U[a]>>b))&1 for a in reps for b in reps if a!=b):changed.append(M)
tops=[m for m in changed if not any(m!=n and m&n==m for n in changed)];v.need(set(tops)==set(r['highest_nodes']),'independent highest changed nodes');values=[v.hp(v.extend(a,*r['query'])) for a in (T,V,U)];v.need(values==r['values'] and values[1]==values[2] and (values[0]-values[2])%8==4,'nuisance highest node exact invariance');v.need(v.hp(T)==v.hp(U)==r['scalar'],'highest control scalar equality')
print(json.dumps({'status':'PASS','contexts':len(checks),**counts,'context_controls':checks,'highest_nodes':tops,'highest_query':r['query'],'highest_values':values,'source_sha256':{p:hashlib.sha256((SRC/p).read_bytes()).hexdigest() for p in ('PROOF.md','verify.py','normal.json','optimized.json')}},sort_keys=True,indent=2))
