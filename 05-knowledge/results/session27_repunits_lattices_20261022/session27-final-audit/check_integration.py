#!/usr/bin/env python3
"""Source-bound checks, read only except stdout. Optional ROOT SHARED arguments."""
import ast,hashlib,itertools,json,pathlib,re,subprocess,sys
ROOT=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else pathlib.Path('/workspace/scratch/07d5532f2fce/session27-root')
SHARED=pathlib.Path(sys.argv[2]) if len(sys.argv)>2 else pathlib.Path('/workspace/shared')
checks={}
def require(test,label):
 if not test: raise RuntimeError(label)
 checks[label]=True
sha=lambda b:hashlib.sha256(b).hexdigest()
for n,h in {'SYNTHESIS.md':'df71f5b5e0792695a1d732a6bd92d93ee940a173b23ef99f2e0ab8999ae3302a','EXPLORATION.md':'587f98e51f63e875ed85db73867e34df13bb050d9bb9030af4c30f5f223bf271','LORENTZIAN_SCORE_BRIDGE.md':'662cd827e092b4ecc513c351c841ed4810167ad9847b590ead25f3f5dfd7b95c'}.items():
 require(sha((ROOT/n).read_bytes())==h,'frozen '+n)
c=SHARED/'session27-corrections'; manifest=json.loads((c/'manifest.json').read_text())
require(len(manifest['files'])==6,'six correction targets')
for f in manifest['files']:
 a=(c/'originals'/f['path']).read_bytes();b=(c/'proposed'/f['path']).read_bytes()
 require(sha(a)==f['original_sha256'] and sha(b)==f['proposed_sha256'],'snapshot hashes '+f['path'])
 require(hashlib.sha1(b'blob '+str(len(a)).encode()+b'\0'+a).hexdigest()==f['git_blob_oid'],'original Git blob '+f['path'])
 if f['path'].endswith('.py'):
  require(ast.dump(ast.parse(a.decode('utf8',errors='surrogateescape')))==ast.dump(ast.parse(b.decode('utf8',errors='surrogateescape'))),'unchanged AST '+f['path'])
 if '/07-reflections/' in '/'+f['path'] or f['path'].endswith('.out') or f['path'].endswith('repunit_gap_synthesis.py'):
  require(b.endswith(a),'unchanged suffix '+f['path'])
 if f['path'].endswith('ramanujan_near_miss_kpc1.py'):
  require(b.splitlines()[0]==a.splitlines()[0] and b.endswith(a.split(b'\n',1)[1]),'shebang and old body preserved')
p=subprocess.run(['git','apply','--check',str(c/'corrections.patch')],cwd=c/'originals',capture_output=True,text=True)
require(p.returncode==0,'patch dry run: '+p.stderr.strip())
r=SHARED/'session27-repunits';idx=json.loads((r/'remote_index.json').read_text());cat=json.loads((r/'catalogue.json').read_text());summary=json.loads((r/'remote_summary.json').read_text())
allrows=[];counts={}
for q,num,absent in [('repunit',100,48),('Mersenne',100,58),('1729',100,36),('Wagstaff',38,8),('Fibonacci fan',100,32)]:
 rows=json.loads((r/('remote_'+q.replace(' ','_')+'.json')).read_text())['results'];indexed=[x for x in idx if x['query']==q.replace(' ','_')]
 require(len(rows)==len(indexed)==num,'raw count '+q)
 require(sum(not x['present_locally'] for x in indexed)==absent,'saved absence count '+q)
 require([{k:v for k,v in x.items() if k not in ('query','present_locally')} for x in indexed]==rows,'raw-index identity '+q)
 require(all('/54bf96642934c6cdf01f02e364f3bba00892ae22/' in x['url'] for x in rows),'pinned URLs '+q)
 allrows+=rows;counts[q]=num
paths={x['path'] for x in idx};absent={x['path'] for x in idx if not x['present_locally']};new=paths-{x['path'] for x in cat['entries']}
require(len(idx)==438 and len(paths)==415 and len(absent)==174 and len(new)==193,'remote totals')
require(summary['unique_paths']==len(paths) and summary['unique_absent_local']==len(absent) and set(summary['new_paths'])==new,'summary agrees with raw index')
for f in json.loads((r/'remote_fetch_manifest.json').read_text()):
 b=(r/'remote_fetched'/pathlib.Path(f['path']).name).read_bytes()
 require(hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==f['sha'],'remote fetched exact blob '+f['path'])
require(cat['tracked_paths']==113808 and cat['text_files_scanned']==72894 and cat['missing_worktree_files']==40488 and len(cat['entries'])==7290,'local catalogue summary')
require(72894+40488+424+2==113808,'local inventory partition')
links=[]
for n in ['SYNTHESIS.md','EXPLORATION.md','TOPOLOGY_SCOPE.md']:
 for target in re.findall(r'\]\(([^)]+)\)',(ROOT/n).read_text()):
  if '://' in target: continue
  target=target.split('#')[0]
  p=(SHARED/target) if target.startswith('session27-') else ROOT/target
  require(p.is_file(),'packet link '+n+' -> '+target);links.append(target)
arcs=[{(0,1),(0,2),(1,2)},{(0,1),(1,2),(2,0)}]
h=[sum(all((p[i],p[i+1]) in a for i in range(2)) for p in itertools.permutations(range(3))) for a in arcs]
require(h==[1,3] and {frozenset(e) for e in arcs[0]}=={frozenset(e) for e in arcs[1]},'topology minimal hostile')
require(1+2*(13-13)+4*3+8==2*13-5 and 1+2*3+4*3+8==27 and 1**3+1**3-2**3==-6,'correction arithmetic')
print(json.dumps({'status':'PASS','checks':checks,'count':len(checks),'scope':'No source mutations, no new exhaustive census; absence flags refer to saved discovery snapshot.'},indent=2))
