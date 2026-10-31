#!/usr/bin/env python3
"""Portable Session130 release: manifests, isolated normal/O, negative gates."""
import hashlib,json,pathlib,shutil,subprocess,sys,tempfile
ROOT=pathlib.Path(__file__).resolve().parent
CASES=[('gain_cohomology/check_lattice.py',[],'gain_cohomology/normal.json'),('primegame_lines/check_primegame_lines.py',[],'primegame_lines/normal.json'),('colored_incidence/check_colored_incidence.py',[],'colored_incidence/normal.json'),('inheritance/check_equivariant_ocf.py',[],'inheritance/normal.json'),('hamiltonian/check_hamiltonian_bridges.py',['--max-blue-n','8'],'hamiltonian/normal.json'),('hamiltonian/verify_hostile.py',[],'hamiltonian/hostile.normal.json'),('independent/check_independent.py',[],'independent/normal.json'),('independent/check_structural.py',[],'independent/structural.json')]
def req(c,m):
 if not c:raise RuntimeError(m)
report={'status':'PASS','manifests':[],'replays':[],'negative_gates':[]}
for manifest in sorted(ROOT.glob('*/SHA256SUMS')):
 count=0
 for line in manifest.read_text().splitlines():
  sha,name=line.split('  ',1);p=manifest.parent/name
  req(hashlib.sha256(p.read_bytes()).hexdigest()==sha,'hash mismatch '+str(p));count+=1
 report['manifests'].append({'path':str(manifest.relative_to(ROOT)),'entries':count})
req(len(report['manifests'])==6,'missing manifest')
with tempfile.TemporaryDirectory(prefix='s130-release-') as tmp:
 local=pathlib.Path(tmp)/'packet';shutil.copytree(ROOT,local,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 for script,options,expected in CASES:
  for opt in (False,True):
   p=subprocess.run([sys.executable,'-I',*(['-O'] if opt else []),str(local/script),*options],cwd=tmp,capture_output=True)
   req(p.returncode==0,script+' failed: '+p.stderr.decode())
   req(p.stdout==(local/expected).read_bytes(),script+' changed exact output')
   report['replays'].append({'script':script,'options':options,'optimized':opt,'stdout_sha256':hashlib.sha256(p.stdout).hexdigest(),'status':'PASS'})
 for script in ('gain_cohomology/check_lattice.py','independent/check_independent.py','independent/check_structural.py'):
  for opt in (False,True):
   p=subprocess.run([sys.executable,'-I',*(['-O'] if opt else []),str(local/script),'--test-failing-gate'],cwd=tmp,capture_output=True)
   req(p.returncode==1,'inactive negative gate '+script)
   report['negative_gates'].append({'script':script,'optimized':opt,'exit_code':p.returncode})
print(json.dumps(report,sort_keys=True,indent=2))
