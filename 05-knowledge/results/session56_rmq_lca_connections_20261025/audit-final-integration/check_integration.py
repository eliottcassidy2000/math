from pathlib import Path, PurePosixPath
import hashlib,json,re,posixpath,subprocess
root=Path('/tmp/session56-root/final'); mp=root/'PREAUDIT_MAP.json'
data=json.loads(mp.read_text()); rows=data['files']; base=data['destination']
repo='/tmp/session56-publication/clean-repo'
prior=set(subprocess.check_output(['git','-C',repo,'ls-tree','-r','--name-only','HEAD','--',base],text=True).splitlines())
new={r['destination']:r for r in rows}
if len(new)!=len(rows):raise RuntimeError('duplicate destination')
for r in rows:
 b=Path(r['source']).read_bytes()
 if len(b)!=r['bytes'] or hashlib.sha256(b).hexdigest()!=r['sha256']:raise RuntimeError(('hash mismatch',r['source']))
targets=prior|set(new); links=[]
for name in ['README.md','CONNECTIONS.md','SYNTHESIS.md']:
 s=(root/name).read_text()
 for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',s):
  if '://' in target or target.startswith('#'):continue
  d=posixpath.normpath(posixpath.join(base,target.split('#')[0]))
  if d not in targets:raise RuntimeError(('broken root link',name,target))
  links.append([name,target])
if len(re.findall(r'^## \d+\.',(root/'CONNECTIONS.md').read_text(),re.M))!=17:raise RuntimeError('not 17 cards')
for directory, key in [('session56-diagrams','path'),('session56-cayley-mirror-bridge','local')]:
 p=Path('/tmp')/directory
 pins=json.loads((p/'SOURCE_PINS.json').read_text())
 if isinstance(pins,dict):pins=pins['sources']
 for pin in pins:
  if hashlib.sha256((p/pin[key]).read_bytes()).hexdigest()!=pin['sha256']:raise RuntimeError(('source pin',pin))
p=Path('/tmp/session56-debt-range'); manifest=json.loads((p/'MANIFEST.json').read_text())['files']
for f,h in manifest.items():
 if hashlib.sha256((p/f).read_bytes()).hexdigest()!=h:raise RuntimeError(('debt manifest',f))
print(json.dumps({'status':'PASS','preaudit_map_sha256':hashlib.sha256(mp.read_bytes()).hexdigest(),'new_files':len(rows),'new_bytes':sum(r['bytes'] for r in rows),'prior_files':len(prior),'root_relative_links':len(links),'connection_cards':17,'bundled_source_pins':'PASS','root_sha256':{n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in ['README.md','CONNECTIONS.md','SYNTHESIS.md']}},indent=2))
