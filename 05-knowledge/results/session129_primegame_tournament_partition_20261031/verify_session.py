#!/usr/bin/env python3
"""Portable Session129 release gate: hashes, isolated normal/-O, hostile gates."""
import hashlib,json,pathlib,shutil,subprocess,sys,tempfile
ROOT=pathlib.Path(__file__).resolve().parent
CASES=[('encoding/check_encoding.py','encoding/normal.json'),('bridge_audit/check_bridge.py','bridge_audit/results.json'),('tournaments/check_tournaments.py','tournaments/census.json'),('partition/check_partition.py','partition/results.json'),('figurate/check_figurate.py','figurate/normal.json'),('independent/check_independent.py','independent/normal.json'),('periodic_density/check_density.py','periodic_density/normal.json')]
def require(c,m):
 if not c:raise RuntimeError(m)
report={'status':'PASS','manifests':[],'replays':[],'negative_gates':[]}
for manifest in sorted(ROOT.glob('*/SHA256SUMS')):
 count=0
 for line in manifest.read_text().splitlines():
  sha,name=line.split('  ',1);p=manifest.parent/name
  require(hashlib.sha256(p.read_bytes()).hexdigest()==sha,'hash mismatch '+str(p));count+=1
 report['manifests'].append({'path':str(manifest.relative_to(ROOT)),'entries':count})
require(len(report['manifests'])==7,'missing manifest')
with tempfile.TemporaryDirectory(prefix='s129-release-') as tmp:
 local=pathlib.Path(tmp)/'packet';shutil.copytree(ROOT,local,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 for script,expected in CASES:
  for opt in (False,True):
   args=[sys.executable,'-I',*(['-O'] if opt else []),str(local/script)]
   p=subprocess.run(args,cwd=tmp,capture_output=True)
   require(p.returncode==0,script+' failed: '+p.stderr.decode())
   require(p.stdout==(local/expected).read_bytes(),script+' changed output')
   report['replays'].append({'script':script,'optimized':opt,'stdout_sha256':hashlib.sha256(p.stdout).hexdigest(),'status':'PASS'})
 for script in ('bridge_audit/check_bridge.py','independent/check_independent.py'):
  for opt in (False,True):
   p=subprocess.run([sys.executable,'-I',*(['-O'] if opt else []),str(local/script),'--test-failing-gate'],cwd=tmp,capture_output=True)
   require(p.returncode==1,'inactive failure gate '+script)
   report['negative_gates'].append({'script':script,'optimized':opt,'exit_code':p.returncode})
 p=subprocess.run([sys.executable,'-I',str(local/'figurate/verify_release.py')],cwd=tmp,capture_output=True)
 require(p.returncode==0,'figurate mutation gate failed '+p.stderr.decode())
 require(p.stdout==(local/'figurate/REPLAY.json').read_bytes(),'figurate replay mismatch')
 report['figurate_mutation_release']=json.loads(p.stdout)
print(json.dumps(report,sort_keys=True,indent=2))
