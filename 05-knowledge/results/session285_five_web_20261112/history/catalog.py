#!/usr/bin/env python3
"""Pinned lexical inventory; hits are leads, not mathematical/session counts."""
import argparse, collections, gzip, hashlib, json, pathlib, re, subprocess
P=pathlib.Path(__file__).resolve(); ROOT=pathlib.Path(subprocess.check_output(['git','rev-parse','--show-toplevel'],cwd=P.parent,text=True).strip())
p=argparse.ArgumentParser(); p.add_argument('--ref',default='7510e38f9c777a61a1b3210334077e0ebe0df2ad'); p.add_argument('--frozen-scope',type=pathlib.Path,help='Replay only recorded path/blob scope'); a=p.parse_args()
out=P.parent
seed={
 'five_word':r'\bfive\b|\bfivefold\b|\bfive[- ](?:core|port|state|parameter|dimensional|fold)',
 'five_set':r'(?:set|group|collection|family|families|list|types?|classes|elements|parameters?)\s+(?:of\s+)?(?:exactly\s+)?(?:five|5)\b|\b(?:five|5)[- ](?:element|parameter|dimensional|member|fold|core|port|state)',
 'pentagonal':r'pentagon|pentagonal|pentagram|pentad|quintic|quintuple',
 'platonic':r'platonic|icosahed|dodecahed|golden triangle',
 'fermat':r'Fermat (?:prime|number)|P[eé]pin',
 'ogg_moonshine':r'\bOgg\b|moonshine|supersingular',
 'stress':r'isothermal|isometric net|stress traject|Airy stress|Nem[eé]nyi',
 'prior_search':r'(?:search|histor|recover|collect|catalog|atlas|inventory).{0,70}(?:five|Platonic|Fermat)|(?:five|Platonic|Fermat).{0,70}(?:search|histor|recover|collect|catalog|atlas|inventory)'}
# Expansion vocabulary is explicitly chosen after reading inherited seed sources.
topics={
 'involution':r'involution|self.converse|anti.automorphism',
 'curvature':r'curvature|Gauss.Bonnet|angle.defect|hyperbolic|isothermal',
 'finite_state':r'automat|finite.state|transfer matrix|Hankel|realization rank',
 'golden':r'Fibonacci|Zeckendorf|golden|Lucas|Cassini',
 'modular_forms':r'modular form|moonshine|supersingular|Atkin.Lehner|Gleason',
 'symplectic':r'symplectic|Lagrangian|Richelot|Pauli|Hamming',
 'arithmetic_dynamics':r'Collatz|FRACTRAN|PRIMEGAME|register machine|cellular automat',
 'interface':r'interface|endpoint|receiver|rooted|port tensor',
 'coherence':r'associahed|pentagon|coherence|holonomy|flat connection|cocycle',
 'cancellation':r'cancellation|moment|signed|kernel|positive cone',
 'geometry':r'conformal|elastic|stress|harmonic|biharmonic|Airy|potential',
 'number_theory':r'aliquot|idoneal|divisor|totient|Carmichael|Pythagorean',
 'graph_spectrum':r'square.sum|metagraph|Hamilton|tournament|Cayley',
 'prime_core':r'prime core|indecomposable|modular decomposition|substitution|critical tournament'}
regex={k:re.compile(v,re.I) for k,v in seed.items()}; treg={k:re.compile(v,re.I) for k,v in topics.items()}
raw=subprocess.check_output(['git','ls-tree','-rz',a.ref],cwd=ROOT)
entries=[]
for e in raw.split(b'\0'):
 if not e: continue
 meta,path=e.split(b'\t',1); mode,kind,sha=meta.decode().split(); path=path.decode()
 if kind=='blob' and re.match(r'0[0-7]-',path) and pathlib.Path(path).suffix.lower() in {'.md','.txt','.tex','.rst'}: entries.append((path,sha))
# One cat-file process, repeated blobs read once. Disable missing-object fetches:
# omitted blobs are listed, not silently called searched.
available={line.split()[0]:int(line.split()[2]) for line in subprocess.check_output(['git','cat-file','--batch-all-objects','--batch-check=%(objectname) %(objecttype) %(objectsize)'],cwd=ROOT,text=True).splitlines() if line.split()[1]=='blob'}
frozen=None if a.frozen_scope is None else {r['path']:r['blob'] for r in json.loads(gzip.decompress(a.frozen_scope.read_bytes()) if a.frozen_scope.suffix=='.gz' else a.frozen_scope.read_text())}
proc=subprocess.Popen(['git','cat-file','--batch'],cwd=ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE,env={**__import__('os').environ,'GIT_NO_LAZY_FETCH':'1'})
cache={}; unavailable=[]; rows=[]; counts=collections.Counter(); expanded=[]
for ix,(path,sha) in enumerate(entries):
 if ix%1000==0: (out/'progress.tmp').write_text(str(ix)+' / '+str(len(entries)))
 if frozen is not None and frozen.get(path)!=sha:
  unavailable.append({'path':path,'sha':sha,'reason':'outside_frozen_scope'}); continue
 if sha not in available or available[sha]>2000000:
  unavailable.append({'path':path,'sha':sha,'reason':'missing_local_blob' if sha not in available else 'over_2MB'}); continue
 if sha not in cache:
  proc.stdin.write((sha+'\n').encode()); proc.stdin.flush(); header=proc.stdout.readline().decode().strip().split()
  if len(header)==2 and header[1]=='missing': cache[sha]=None
  else:
   size=int(header[2]); content=proc.stdout.read(size); end=proc.stdout.read(1);
   if end!=b'\n': raise RuntimeError('invalid batch framing')
   cache[sha]=content.decode('utf8','replace')
 text=cache[sha]
 if text is None: unavailable.append({'path':path,'sha':sha,'reason':'missing_local_blob' if sha not in available else 'over_2MB'}); continue
 tags=[k for k,r in treg.items() if r.search(text)]
 hits={}; excerpts={}
 for k,r in regex.items():
  lines=[(i,line) for i,line in enumerate(text.splitlines(),1) if r.search(line)]
  if lines: hits[k]=len(lines); excerpts[k]=[[i,line[:350]] for i,line in lines[:3]]; counts[k]+=1
 if hits:
  rows.append({'path':path,'blob':sha,'queries':hits,'topics':tags,'excerpts':excerpts})
 elif len(tags)>=2:
  expanded.append({'path':path,'blob':sha,'topics':tags})
proc.stdin.close(); proc.wait()
co=collections.Counter()
for row in rows:
 for i,x in enumerate(row['topics']):
  for y in row['topics'][i+1:]: co[(x,y)]+=1
summary={'ref':a.ref,'scope':'tracked prose under 00–07 at pinned ref; no agents mirrors, untracked work, binary attachments or deleted versions; skip blobs over2MB and unavailable local objects','candidate_paths':len(entries),'distinct_candidate_blobs':len(set(s for _,s in entries)),'read_paths':len(entries)-len(unavailable),'unavailable_paths':len(unavailable),'seed_hit_paths':len(rows),'seed_hit_distinct_blobs':len(set(r['blob'] for r in rows)),'query_path_counts':dict(counts),'expansion_only_paths_with_at_least_two_topics':len(expanded),'expansion_vocabulary_provenance':'manually selected from deeply read inherited five-core, aliquot, Platonic, PRIMEGAME and coherence sources; co-occurrence edges are not mathematical bridges','seed_regex':seed,'topic_regex':topics}
for name,obj in [('inventory.json',rows),('expansion.json',expanded),('unavailable.json',unavailable),('summary.json',summary),('cooccurrence.json',[{'a':a,'b':b,'paths':v,'status':'LEXICAL_ONLY'} for (a,b),v in co.most_common()])]: (out/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(summary,indent=2)); print('output_sha256',hashlib.sha256((out/'inventory.json').read_bytes()).hexdigest())
