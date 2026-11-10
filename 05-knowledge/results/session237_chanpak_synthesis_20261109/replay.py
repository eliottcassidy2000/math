#!/usr/bin/env python3
"""Portable compact/full Session237 exact replay, with live failure controls."""
import argparse,hashlib,json,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
BASE=Path(__file__).resolve().parent
P=Path('05-knowledge/results/session237_chanpak_synthesis_20261109')
SCRIPTS=[
'04-computation/experiments/session237_arithmetic_states/verify.py',
'04-computation/experiments/session237_arithmetic_states/verify_endpoint.py',
'04-computation/experiments/session237_cycle_tilings/verify.py',
'04-computation/experiments/session237_symmetry_grafts/witness_audit.py',
'04-computation/experiments/session237_endpoint_slack/verify.py',
str(P/'arithmetic_states_audit/independent_verify.py'),
str(P/'cycle_tilings_audit/check.py'),
str(P/'symmetry_grafts_audit/check.py'),
str(P/'endpoint_slack_audit/check.py'),
str(P/'endpoint_slack_audit/endpoint_dfa.py')]
def require(test,msg):
 if not test:raise RuntimeError(msg)
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');ap.add_argument('--fail',action='store_true');args=ap.parse_args()
 require(not args.fail,'deliberate replay failure')
 scripts=SCRIPTS+(['04-computation/experiments/session237_symmetry_grafts/probe.py'] if args.full else [])
 rows=[]
 for script in scripts:
  print('replay '+script,flush=True);outs=[]
  for opt in [[],['-O']]:
   r=subprocess.run([sys.executable,*opt,script],cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=1200)
   require(r.returncode==0,script+' failed '+r.stderr.decode());outs.append(r.stdout)
   f=subprocess.run([sys.executable,*opt,script,'--fail'],cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=15)
   require(f.returncode!=0 and b'deliberate' in f.stderr,script+' failed sentinel')
  require(outs[0]==outs[1],script+' optimized mismatch')
  rows.append({'script':script,'source_sha256':sha((ROOT/script).read_bytes()),'output_sha256':sha(outs[0]),'normal_and_optimized':'PASS','both_failure_controls':'PASS'})
 a=(ROOT/'04-computation/experiments/session237_symmetry_grafts/probe.out').read_bytes()
 b=(BASE/'symmetry_grafts_audit/probe_optimized.out').read_bytes()
 require(a==b,'full symmetry retained output mismatch')
 require(b'202013' in b and b'4275' in b,'full symmetry counts missing')
 result={'status':'PASS','mode':'full' if args.full else 'compact','checks':rows,'full_symmetry_preserved_normal_optimized_output_sha256':sha(a),'full_symmetry_freshly_replayed':args.full,'scopes':'See REPRODUCE.md; no all-target coverage claimed.'}
 (BASE/'REPLAY.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
