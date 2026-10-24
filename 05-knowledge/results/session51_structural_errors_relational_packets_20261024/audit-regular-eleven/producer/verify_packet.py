#!/usr/bin/env python3
"""Portable explicit checks; no assertions, network, or original scratch paths."""
import collections,csv,hashlib,json,math,pathlib
ROOT=pathlib.Path(__file__).resolve().parent
CHECKS=collections.Counter()
def require(x,m):
 CHECKS[m]+=1
 if not x:raise RuntimeError(m)
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read_atlas(path,n):
 words=iter((ROOT/path).read_text().split());size=int(next(words));out=[]
 for _ in range(size):
  r=tuple(int(next(words)) for _ in range(n));g=tuple(tuple(int(next(words)) for _ in range(n)) for _ in range(int(next(words))));out.append((r,g))
 require(next(words,None) is None,'atlas consumes exactly input');return out
def main():
 man=read('MANIFEST.json')
 for f in man['files']:
  p=ROOT/f['path'];require(p.stat().st_size==f['bytes'] and sha(f['path'])==f['sha256'],'manifest file hash')
 pins=read('SOURCE_PINS.json')
 for f in pins['sources']:require(sha(f['path'])==f['sha256'],'source hash')
 require(sha('sources/parent-nine-manifest.json')==pins['parent_nine_producer_manifest'],'frozen source manifest')
 require(sha('sources/parent-nine-audit-manifest.json')==pins['parent_nine_audit_manifest'],'frozen source audit')
 require(sha('sources/nonregular-eleven-manifest.json')==pins['nonregular_eleven_manifest'],'frozen nonregular manifest')
 nonreg_manifest=read('sources/nonregular-eleven-manifest.json')
 for original,archived in [('THEOREM.md','sources/nonregular-eleven-theorem.md'),('normal.json','sources/nonregular-eleven-receipt.json')]:
  require(sha(archived)==next(x['sha256'] for x in nonreg_manifest['files'] if x['path']==original),'nonregular payload matches archived manifest')
 source_manifest=read('sources/parent-nine-manifest.json');entry=next(x for x in source_manifest['files'] if x['path']=='parents-nine.txt')
 require(sha('sources/parents-nine.txt')==entry['sha256'],'nine atlas matches authenticated manifest')
 small=read('sources/parents-through-seven.json')
 for m in (5,7):
  actual=read_atlas(f'sources/parents-{m}.txt',m);source=next(x['parents'] for x in small['levels'] if x['order']==m)
  require(actual==[(tuple(p['rows']),tuple(map(tuple,p['actions']))) for p in source],'control source text equals certified JSON')
 dp=read('degree-mass.json');require(dp==read('degree-mass-optimized.json'),'normal optimized degree DP equality')
 require(read('hostile-normal.json')==read('hostile-optimized.json'),'normal optimized hostile equality')
 require(all(x['rejected'] for x in read('hostile-normal.json')['controls']),'all hostile mutations rejected')
 for label in ('seven','nine','eleven'):
  for suffix in ('.json','-parents.txt','-states.tsv','-components.tsv'):
   require((ROOT/f'{label}-O2{suffix}').read_bytes()==(ROOT/f'{label}-O3{suffix}').read_bytes(),'all O2 O3 artifacts byte identical')
  r=read(f'{label}-O2.json');n=r['parent_order'];k=n+1;atlas=read_atlas(f'{label}-O2-parents.txt',n)
  require(len(atlas)==r['almost_regular_parent_classes'],'parent count matches receipt')
  require(len({p for p,g in atlas})==len(atlas),'distinct literal parent normal forms')
  for rows,g in atlas:
   require(collections.Counter(x.bit_count() for x in rows)=={n//2-1:n//2,n//2:n//2},'almost regular parent scores')
   require(len(g)==len(set(g)) and tuple(range(n)) in g,'distinct actions and identity')
   for a in g:
    require(sorted(a)==list(range(n)),'parent action permutation')
    require(all(sum(((rows[i]>>j)&1)<<a[j] for j in range(n))==rows[a[i]] for i in range(n)),'parent action all arcs')
  states=list(csv.DictReader((ROOT/f'{label}-O2-states.tsv').open(),delimiter='\t'))
  comps=list(csv.DictReader((ROOT/f'{label}-O2-components.tsv').open(),delimiter='\t'))
  require(len(states)==len(atlas),'one state row per parent')
  require(len(comps)==r['regular_classes'],'component count matches receipt')
  D=[];conv=[];comp_of=[]
  for i,s in enumerate(states):
   require(int(s['state'])==i,'consecutive state ids');D.append({int(j):int(v) for j,v in (z.split(':') for z in s['reroot_destinations'].split(','))});conv.append(int(s['converse_state']));comp_of.append(int(s['component']))
   require(sum(D[-1].values())==k and i in D[-1],'reroot row sum and self root')
   require(int(s['parent_aut_order'])==len(atlas[i][1]),'state parent automorphism count')
  eh=collections.Counter();ah=collections.Counter();seh=collections.Counter();oh=collections.Counter();coverage=set();mass=0
  for cid,c in enumerate(comps):
   require(int(c['component'])==cid,'consecutive component ids')
   members=list(map(int,c['states'].split(',')));orbits=sorted(map(int,c['root_orbit_sizes'].split(',')));group=int(c['aut_order']);error=int(c['error']);mate=int(c['converse_component'])
   require(not (set(members)&coverage),'components disjoint');coverage.update(members)
   require(sum(orbits)==k and len(orbits)==len(members),'component orbit sizes')
   require(error==k*(len(members)-1),'component defect formula')
   require(sum(D[i][i] for i in members)==k,'component diagonal trace')
   for i in members:
    require(comp_of[i]==cid and set(D[i])==set(members),'exact component support')
    require(D[i]==D[members[0]],'all component rows identical')
    require(sorted(D[i].values())==orbits,'component orbit multiplicities')
    require(D[i][i]*len(atlas[i][1])==group,'full group order from parent stabilizer')
    require(comp_of[conv[i]]==mate,'converse maps components')
   require(math.factorial(k)%group==0,'group order factorial divisor');mass+=math.factorial(k)//group;eh[error]+=1;ah[group]+=1;oh[tuple(orbits)]+=1
   if mate==cid:seh[error]+=1
  require(coverage==set(range(len(atlas))),'component coverage')
  for i in range(len(atlas)):
   require(conv[conv[i]]==i,'parent converse involution')
   require({conv[j]:v for j,v in D[i].items()}==D[conv[i]],'converse matrix intertwining')
  require(sum(D[i][i] for i in range(len(atlas)))==r['regular_trace']==k*len(comps),'global trace identity')
  require(sum(e*c for e,c in eh.items())==r['regular_error']==k*(len(atlas)-len(comps)),'global defect identity')
  require(sorted(map(list,eh.items()))==r['regular_error_histogram'],'error histogram')
  require(sorted(map(list,ah.items()))==r['regular_automorphism_histogram'],'full group histogram')
  require(sorted(map(list,seh.items()))==r['self_converse_regular_error_histogram'],'self converse error histogram')
  require([[list(a),b] for a,b in sorted(oh.items())]==r['root_orbit_size_histogram'],'orbit multiset histogram')
  require(sum(i==j for i,j in enumerate(conv))==r['self_converse_parents']==r['self_converse_regular_classes']==sum(seh.values()),'self converse pointing equality')
  independent=next(x['labelled_regular_count'] for x in dp['results'] if x['regular_order']==k)
  require(mass==sum(math.factorial(n)//len(g) for p,g in atlas)==independent==r['labelled_regular_mass']==r['labelled_almost_regular_mass'],'three independent labeled masses')
  require(r['target_canonicalizations']==0,'no target canonicalization')
 old=read('sources/regular-nine-reference.json');new=read('nine-O2.json')
 for a,b in [('pointed_regular_states','almost_regular_parent_classes'),('regular_trace','regular_trace'),('regular_classes','regular_classes'),('regular_error','regular_error')]:require(old[a]==new[b],'inherited regular nine control')
 require(old['root_orbit_size_multisets']==new['root_orbit_size_histogram'],'inherited nine full orbit histogram')
 regular=read('eleven-O2.json');nonreg=read('sources/nonregular-eleven-receipt.json');combined=read('combined-correction.json')
 for inputname,outname in [('nonregular_strong_positive_histogram','strong_positive_error_histogram'),('known_total_positive_histogram','total_positive_error_histogram')]:
  h=collections.Counter(dict(nonreg[inputname]));h.update({e:c for e,c in regular['regular_error_histogram'] if e})
  require(sorted(map(list,h.items()))==combined[outname],'combined histogram structural addition')
  require(sum(e*c for e,c in h.items())==combined['strong_error' if inputname.startswith('nonregular') else 'total_error'],'combined exact first moment')
 require(combined['strong_error']==regular['regular_error']+nonreg['nonregular_strong_error'],'combined strong correction')
 require(combined['total_error']==regular['regular_error']+nonreg['known_total_error'],'combined total correction')
 result={'status':'PASS','manifest_sha256':sha('MANIFEST.json'),'checks':dict(CHECKS),'scope':'frozen source pins; parent matrices and actions; all state/component arithmetic and converse; independent degree mass; two-mode byte equality'}
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
