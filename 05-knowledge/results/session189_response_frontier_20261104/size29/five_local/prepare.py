#!/usr/bin/env python3
import itertools,json,pathlib,collections,hashlib
P=pathlib.Path(__file__).resolve().parent
N=5; PAIRS=list(itertools.combinations(range(N),2)); PERMS=list(itertools.permutations(range(N)))
def rows(c):
 r=[0]*N
 for b,(i,j) in enumerate(PAIRS):
  if c>>b&1:r[i]|=1<<j
  else:r[j]|=1<<i
 return r
def code(r):return sum(((r[i]>>j)&1)<<b for b,(i,j) in enumerate(PAIRS))
def canonical(r):return min(sum(((r[p[i]]>>p[j])&1)<<b for b,(i,j) in enumerate(PAIRS)) for p in PERMS)
def profile(r):
 q=[0]*N
 for p in PERMS:
  bad=sum((not(r[p[i]]>>p[i+1]&1))<<i for i in range(N-1))
  for cuts in range(16):
   if cuts&bad==bad:q[cuts.bit_count()]+=1
 return q
def modules(r):return [m for m in range(1,32) if m.bit_count()>1 and all((r[v]&m) in (0,m) for v in range(N) if not(m>>v&1))]
def toggle(c,ms):
 for m in ms:c^=sum(1<<b for b,(i,j) in enumerate(PAIRS) if m>>i&1 and m>>j&1)
 return c
def aug(r,a,c):return [z+((a>>i&1)<<5)+((c>>i&1)<<6) for i,z in enumerate(r)]+[(31^a)|64,31^c]
def hamilton(r):
 n=len(r);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for i,v in enumerate(d[s]):
   todo=r[i]&~s
   while v and todo:
    b=todo&-todo;todo-=b;d[s|b][b.bit_length()-1]+=v
 return sum(d[-1])
def control(r):
 v=[hamilton(aug(r,m&31,m>>5)) for m in range(1024)];c=v[:]
 for j in range(10):
  for m in range(1024):
   if m>>j&1:c[m]-=c[m^(1<<j)]
 return dict(rows=r,values=v,coefficients=c)
def main():
 rs=[rows(c) for c in range(1024)];cs=[canonical(r) for r in rs];qs=[profile(r) for r in rs];classes=[];bad=[];allcount=noncount=0
 for rep in sorted(set(cs)):
  r=rs[rep];ms=modules(r);targets=[]
  for f in range(1<<len(ms)):
   fam=[m for j,m in enumerate(ms) if f>>j&1];t=toggle(rep,fam)
   targets.append(dict(family=fam,code=t,canonical_code=cs[t],profile=qs[t]))
  q=qs[rep];multiplicity=cs.count(rep);allcount+=multiplicity*len(targets)
  if q[0]>1:noncount+=multiplicity*len(targets)
  changed=[t for t in targets if t['profile']!=q]
  if changed:bad.append(dict(source=rep,targets=changed))
  classes.append(dict(code=rep,rows=r,sorted_degrees=sorted(x.bit_count() for x in r),labeled_count=multiplicity,nontransitive=q[0]>1,profile=q,normalized_mod8=[x*pow(q[0],-1,8)%8 for x in q],normalized_mod16=[x*pow(q[0],-1,16)%16 for x in q],modules=ms,targets=targets,all_family_profiles_equal=not changed,converse_class=cs[rep^1023]))
 summary=dict(classes=len(classes),nontransitive_classes=sum(c['nontransitive'] for c in classes),labeled_inputs=1024,nontransitive_labeled=sum(c['labeled_count'] for c in classes if c['nontransitive']),family_outputs=allcount,nontransitive_family_outputs=noncount,profile_changing_classes=[b['source'] for b in bad])
 (P/'five_child_records.json').write_text(json.dumps(classes,indent=2)+'\n');(P/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(json.dumps(summary));print('\n'.join(str({k:c[k] for k in ('code','rows','profile','normalized_mod8','normalized_mod16','modules','all_family_profiles_equal','converse_class')}) for c in classes),flush=True)
 controls=[]
 for c in classes:
  controls.append(dict(type='five_code_'+str(c['code']),variable_order=['a'+str(i) for i in range(5)]+['c'+str(i) for i in range(5)],source=control(c['rows']),converse=control(rs[c['code']^1023])))
 (P/'named_controls.json').write_text(json.dumps(controls,indent=2)+'\n')
 stats=collections.Counter()
 for ctl in controls:
  for side in ('source','converse'):
   for m,v in enumerate(ctl[side]['coefficients']):
    if v:stats[m.bit_count()]+=1
 print('nonzero coefficients by degree',dict(stats))
 (P/'control_summary.json').write_text(json.dumps(dict(vectors=24,coefficients=24576,nonzero_by_degree=dict(stats)),indent=2)+'\n')
if __name__=='__main__':main()
