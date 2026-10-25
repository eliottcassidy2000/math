from pathlib import Path
import json,hashlib,subprocess,sys,re,posixpath,tempfile,shutil
A=Path(__file__).resolve().parent
def need(x,s):
 if not x:raise RuntimeError(s)
def sha(b):return hashlib.sha256(b).hexdigest()
def validate(r):
 fs=json.loads((r/'MANIFEST.json').read_text())['files'];need({str(p.relative_to(r)) for p in r.rglob('*') if p.is_file()}==set(fs)|{'MANIFEST.json'},'exact closure')
 for rel,h in fs.items():
  p=r/rel;need(not p.is_symlink() and p.resolve().is_relative_to(r.resolve()),'payload path');need(sha(p.read_bytes())==h,'digest '+rel)
validate(A)
mp=json.loads((A/'inputs/pre-audit-map.json').read_text());base=mp['destination'];by={x['destination']:x for x in mp['files']};need(len(by)==len(mp['files'])==514 and len(mp['packet_manifests'])==28,'map counts')
for summary in ['README.md','SYNTHESIS.md']:
 b=(A/'inputs'/summary).read_bytes();need(sha(b)==by[base+'/'+summary]['sha256'],'summary hash')
 for link in re.findall(r'\]\(([^)]+)\)',b.decode()):need(posixpath.normpath(base+'/'+link.split('#')[0]) in by,'summary link')
for rel,h in json.loads((A/'inputs/source-doc-pins.json').read_text()).items():need(sha((A/'inputs/source-docs'/rel).read_bytes())==h==by[base+'/'+rel]['sha256'],'source prose '+rel)
for p in (A/'inputs/figures').glob('*.svg'):need(sha(p.read_bytes())==by[base+'/diagrams/'+p.name]['sha256'],'figure bytes')
for row in json.loads((A/'inputs/endpoint-paths.json').read_text()):need(row['source_sha256']==by[base+'/larger-paley-certificates/certificate-'+str(row['q'])+'.json']['sha256'],'endpoint witness provenance')
for row in json.loads((A/'inputs/progress-routing-receipt.json').read_text())['targets']:
 for mode,key in [('before','source_sha256'),('after','result_sha256')]:need(sha((A/'inputs/progress'/mode/row['path']).read_bytes())==row[key],'annotation source hash')
for opt,fixture in [([], 'arithmetic.normal.json'),(['-O'],'arithmetic.optimized.json')]:need(subprocess.check_output([sys.executable,'-B',*opt,str(A/'check_arithmetic.py')],cwd='/tmp')==(A/fixture).read_bytes(),'fresh arithmetic')
if '--self-test' in sys.argv:
 with tempfile.TemporaryDirectory() as t:
  r=Path(t)/'audit';shutil.copytree(A,r);validate(r)
  for kind in ['missing','corrupt','extra','nested-extra']:
   p=r/('inputs/extra' if kind=='nested-extra' else 'extra' if kind=='extra' else 'inputs/README.md');old=p.read_bytes() if p.exists() else None
   if kind=='missing':p.unlink()
   else:p.write_bytes((old or b'')+b'!')
   try:validate(r)
   except RuntimeError:pass
   else:raise RuntimeError('accepted '+kind)
   if old is None:p.unlink()
   else:p.write_bytes(old)
  for opt in [[],['-O']]:subprocess.run([sys.executable,'-B',*opt,str(r/'replay.py')],cwd='/tmp',check=True,capture_output=True)
print('PASS standalone final synthesis audit: exact closure, frozen source/map/link checks, normal/-O arithmetic, optional relocation and mutation hostiles')
