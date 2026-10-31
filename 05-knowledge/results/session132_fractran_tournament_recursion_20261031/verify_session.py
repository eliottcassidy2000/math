#!/usr/bin/env python3
"""Session132 portable release; isolated stdlib with explicit bundled imports."""
import hashlib,json,pathlib,shutil,subprocess,sys,tempfile
ROOT=pathlib.Path(__file__).resolve().parent
CASES=[
 ('counting/check_counts.py',[],'counting/normal.json',[]),
 ('counting/verify_hostile.py',[],'counting/hostile.normal.json',[]),
 ('converse_folding/check_folding.py',[],'converse_folding/normal.json',[]),
 ('pure_blue_fractran/verify.py',[],'pure_blue_fractran/verification.out',['verification.json','first_transition_trace.json']),
 ('pure_blue_fractran/analyze_geometry.py',[],None,['geometry.json']),
 ('fractran/verify.py',[],'fractran/verification-normal.log',['VERIFIED.json','a000568.ir.json','a000568.fractions.00.txt']),
 ('prime_graphs/check.py',[],'prime_graphs/verification.normal.json',[]),
 ('independent/audit_fractions.py',['../fractran/a000568','--through','5'],'independent/a000568_audit.log',[]),
 ('independent/audit_fractions.py',['../fractran/literal_lines','--through','3'],'independent/literal_lines_audit.log',[]),
 ('independent/audit_compiler.py',[],'independent/compiler_audit.json',[]),
 ('independent/audit_small_programs.py',[],'independent/small_programs_audit.json',[]),
 ('independent/audit_geometry.py',[],'independent/geometry_audit.json',[]),
 ('independent/audit_counts.py',[],'independent/counts_audit.json',[]),
 ('independent/audit_prime_graphs.py',[],'independent/prime_graphs_audit.json',[])]
SHIM='import runpy,sys;from pathlib import Path;p=Path(sys.argv[1]);sys.path.insert(0,str(p.parent));sys.argv=sys.argv[1:];runpy.run_path(str(p),run_name="__main__")'
def req(c,m):
 if not c:raise RuntimeError(m)
def run(local,script,args,opt=False):
 return subprocess.run([sys.executable,'-I',*(['-O'] if opt else []),'-c',SHIM,str(local/script),*args],cwd=(local/script).parent,capture_output=True)
report={'status':'PASS','manifests':[],'replays':[],'negative_gates':[]}
for m in sorted(ROOT.glob('*/SHA256SUMS')):
 count=0
 for line in m.read_text().splitlines():
  h,f=line.split('  ',1);req(hashlib.sha256((m.parent/f).read_bytes()).hexdigest()==h,'hash mismatch '+str(m.parent/f));count+=1
 report['manifests'].append({'path':str(m.relative_to(ROOT)),'entries':count})
req(len(report['manifests'])==5,'missing SHA manifest')
pm=json.loads((ROOT/'prime_graphs/MANIFEST.json').read_text())['sha256']
for f,h in pm.items():req(hashlib.sha256((ROOT/'prime_graphs'/f).read_bytes()).hexdigest()==h,'graph manifest mismatch '+f)
report['manifests'].append({'path':'prime_graphs/MANIFEST.json','entries':len(pm)})
with tempfile.TemporaryDirectory(prefix='s132-release-') as tmp:
 local=pathlib.Path(tmp)/'packet';shutil.copytree(ROOT,local,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 for script,args,reference,generated in CASES:
  first=None
  for opt in (False,True):
   p=run(local,script,args,opt);req(p.returncode==0,script+' failed: '+p.stderr.decode())
   if reference:req(p.stdout==(ROOT/reference).read_bytes(),script+' changed stdout')
   elif first is not None:req(p.stdout==first,script+' differs under optimization')
   first=p.stdout
   for f in generated:
    relative=pathlib.Path(script).parent/f
    req((local/relative).read_bytes()==(ROOT/relative).read_bytes(),'regenerated artifact changed '+str(relative))
   report['replays'].append({'script':script,'args':args,'optimized':opt,'stdout_sha256':hashlib.sha256(p.stdout).hexdigest(),'status':'PASS'})
 for opt in (False,True):
  p=run(local,'converse_folding/check_folding.py',['--test-failing-gate'],opt);req(p.returncode==1,'fold negative gate inactive')
  report['negative_gates'].append({'control':'intentional fold gate','optimized':opt,'exit_code':p.returncode})
 for name,script,args,file in [
  ('damaged A000568 row','independent/audit_fractions.py',['../fractran/a000568','--through','0'],'fractran/a000568.fractions.00.txt'),
  ('damaged pure-blue row','pure_blue_fractran/verify.py',[],'pure_blue_fractran/pure_blue.fractran')]:
  f=local/file;original=f.read_bytes();lines=f.read_text().splitlines();lines[0]='1/1';f.write_text('\n'.join(lines)+'\n')
  for opt in (False,True):
   p=run(local,script,args,opt);req(p.returncode==1,name+' survived');report['negative_gates'].append({'control':name,'optimized':opt,'exit_code':p.returncode})
  f.write_bytes(original)
 p=run(local,'prime_graphs/verify_packet.py',[]);req(p.returncode==0,'graph whitelist verifier failed '+p.stderr.decode());report['graph_whitelist_stdout']=p.stdout.decode().strip()
print(json.dumps(report,indent=2,sort_keys=True))
