from pathlib import Path
import hashlib,json,shutil,subprocess,tempfile
P=Path(__file__).resolve().parent
def need(x,s):
 if not x:raise RuntimeError(s)
def run(args,cwd):return subprocess.run([str(x)for x in args],cwd=cwd,check=True,capture_output=True)
with tempfile.TemporaryDirectory(prefix='s197-feedback533-replay-')as tmp:
 T=Path(tmp);shutil.copytree(P/'source',T,dirs_exist_ok=True)
 L=T/'session197-lowcycle';B=T/'session197-bridges';A=B/'feedback533-audit'
 frozen=(L/'feedback533_certificate.json').read_bytes()
 for opts in ([],['-O']):
  r=run(['python3',*opts,L/'verify_feedback533.py'],T)
  need(r.stdout==frozen and (L/'feedback533_certificate.json').read_bytes()==frozen,'source exact replay')
 source=T/'source.csv';ind=B/'c3-independent-records.csv'
 for src,exe in [(L/'verify_c3_complete.cpp',T/'census'),(B/'check_c3_complete.cpp',T/'independent'),(A/'filter.cpp',T/'filter')]:
  run(['g++','-O2','-std=c++17',src,'-o',exe],T)
 run([T/'census',source],T)
 r=run([T/'independent',ind,source],T);need(r.stdout==(B/'c3-independent.txt').read_bytes(),'independent bank receipt')
 r=run([T/'filter',ind],T)
 need(r.stdout==(A/'independent-hits.csv').read_bytes() and r.stderr==(A/'independent-counts.txt').read_bytes(),'complete uncapped filter')
 for opts in ([],['-O']):
  r=run(['python3',*opts,A/'check.py'],T);need(r.stdout==(A/'controls.json').read_bytes(),'independent graph controls')
 print(json.dumps({'status':'PASS','source_replays':2,'independent_bank_classes':27945,'C1_C2_graphs':12077,'exact_H533_cores':6,'all_feedback_numbers':3},sort_keys=True))
