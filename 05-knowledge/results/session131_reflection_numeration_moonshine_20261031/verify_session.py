#!/usr/bin/env python3
"""Portable Session131 release gate; exact hashes, isolated replays, hostiles."""
import hashlib,json,pathlib,shutil,subprocess,sys,tempfile
ROOT=pathlib.Path(__file__).resolve().parent
CASES=[('integral_normal_form/check_normal_form.py','integral_normal_form/normal.json'),('moonshine/check.py','moonshine/results.json'),('binary_carries/check_binary_carries.py','binary_carries/normal.json'),('binary_carries/verify_hostiles.py','binary_carries/hostile.normal.json'),('zeckendorf/check_bridge.py','zeckendorf/normal.json'),('square_sum/check.py','square_sum/verification.normal.json')]+[(f'independent/check_{script}.py',f'independent/{stem}.normal.json') for script,stem in [('module','module'),('moonshine_polynomials','moonshine'),('numeration_certificates','numeration'),('tournament_carries','carries'),('square_sum','square_sum')]]
def req(c,m):
 if not c:raise RuntimeError(m)
report={'status':'PASS','manifests':[],'replays':[],'negative_gates':[]}
for manifest in sorted(ROOT.glob('*/SHA256SUMS')):
 count=0
 for line in manifest.read_text().splitlines():
  sha,name=line.split('  ',1);p=manifest.parent/name
  req(hashlib.sha256(p.read_bytes()).hexdigest()==sha,'hash mismatch '+str(p));count+=1
 report['manifests'].append({'path':str(manifest.relative_to(ROOT)),'entries':count})
req(len(report['manifests'])==5,'missing SHA manifest')
sq=json.loads((ROOT/'square_sum/MANIFEST.json').read_text())
for name,sha in sq.items():req(hashlib.sha256((ROOT/'square_sum'/name).read_bytes()).hexdigest()==sha,'square manifest mismatch')
report['manifests'].append({'path':'square_sum/MANIFEST.json','entries':len(sq)})
for pin in json.loads((ROOT/'independent/INPUT_PINS.json').read_text())['inputs']:
 req(hashlib.sha256((ROOT/'independent'/pin['path']).read_bytes()).hexdigest()==pin['sha256'],'audited input changed '+pin['path'])
with tempfile.TemporaryDirectory(prefix='s131-release-') as tmp:
 local=pathlib.Path(tmp)/'packet';shutil.copytree(ROOT,local,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 for script,expected in CASES:
  for opt in (False,True):
   p=subprocess.run([sys.executable,'-I',*(['-O'] if opt else []),str(local/script)],cwd=tmp,capture_output=True)
   req(p.returncode==0,script+' failed '+p.stderr.decode());req(p.stdout==(local/expected).read_bytes(),script+' output changed')
   report['replays'].append({'script':script,'optimized':opt,'stdout_sha256':hashlib.sha256(p.stdout).hexdigest(),'status':'PASS'})
 for script in ['integral_normal_form/check_normal_form.py','zeckendorf/check_bridge.py']+[x for x,y in CASES if x.startswith('independent/')]:
  for opt in (False,True):
   p=subprocess.run([sys.executable,'-I',*(['-O'] if opt else []),str(local/script),'--test-failing-gate'],cwd=tmp,capture_output=True)
   req(p.returncode==1,'negative gate inactive '+script);report['negative_gates'].append({'script':script,'optimized':opt,'exit_code':p.returncode})
 p=subprocess.run([sys.executable,'-I',str(local/'square_sum/verify_packet.py')],cwd=tmp,capture_output=True)
 req(p.returncode==0,'square whitelist verifier failed '+p.stderr.decode());report['square_sum_release']=json.loads(p.stdout)
print(json.dumps(report,sort_keys=True,indent=2))
