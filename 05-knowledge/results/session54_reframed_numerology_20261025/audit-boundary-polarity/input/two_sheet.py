#!/usr/bin/env python3
import json,pathlib,collections as co
D=pathlib.Path(__file__).resolve().parent
def ck(p,m):
 if not p:raise RuntimeError(m)
def mul(a,b):
 out=0
 while b:
  if b&1:out^=a
  b>>=1;a<<=1
  if a&64:a^=67
 return out
def power(a,n):
 out=1
 for _ in range(n):out=mul(out,a)
 return out
def trace(x):return x^power(x,4)^power(x,16)
powers=[power(2,e) for e in range(63)]
def inc(p,l):return trace(powers[(p-l)%63])==0
base={(p,l) for p in range(21) for l in range(21) if inc(p,l)}
ck(len(base)==105,'flag count')
ck(all(sum((p,l) in base for l in range(21))==5 for p in range(21)),'point degree')
ck(all(sum((p,l) in base for p in range(21))==5 for l in range(21)),'line degree')
ck(all(inc(p,l)==inc(-l%21,-p%21) for p in range(21) for l in range(21)),'polarity')
for k in range(3):
 for b in range(7):
  ck(all(inc(p,l)==inc((16**k*p+9*b)%21,(16**k*l+9*b)%21) for p in range(21) for l in range(21)),'G collineations')
rows=json.loads((D/'result.json').read_text())['map'];exps={r['blocker']:r['exponent'] for r in rows};boundary={r['blocker']:r['boundary'] for r in rows}
lift={(a,b) for a in exps for b in exps if inc(exps[a],exps[b])}
ck(len(lift)==945,'lift flags')
ck(all(sum((a,b) in lift for b in exps)==15 for a in exps),'lift degree')
for a in exps:
 for b in exps:
  ck(({c for c in exps if (a,c) in lift}=={c for c in exps if (b,c) in lift})==(boundary[a]==boundary[b]),'twins equal fibers')
hist={};witness={}
for a in exps:
 for b in exps:
  t=(a&b).bit_count();value=(a,b) in lift
  hist.setdefault(t,co.Counter())[value]+=1;witness.setdefault((t,value),(a,b))
t=next(t for t in hist if len(hist[t])==2)
# A trace polarity can have absolute points; record rather than infer none from bare inversion.
absolute=[p for p in range(21) if inc(p,-p%21)]
print(json.dumps({'status':'PASS','base_vertices':42,'base_flags':105,'base_degree':5,'lifted_vertices':126,'lifted_flags':945,'lifted_degree':15,'base_incidence':sorted(base),'absolute_points':absolute,'native_intersection_hostile':{'shared_arc_count':t,'incident_pair':witness[t,True],'nonincident_pair':witness[t,False]},'native_intersection_histogram':{k:{str(a):b for a,b in v.items()} for k,v in hist.items()}},sort_keys=True,indent=2))
