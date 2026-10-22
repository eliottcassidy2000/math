"""Portable representative replay. Copies the packet before checks may write outputs."""
import argparse,hashlib,json,os,pathlib,shutil,subprocess,tempfile,time
p=argparse.ArgumentParser();p.add_argument('--packet-root',type=pathlib.Path,required=True);p.add_argument('--repository-root',type=pathlib.Path,required=True);p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args()
packet=a.packet_root.resolve();repo=a.repository_root.resolve();out=a.out.resolve()
commands=[
 ['python3','session30-root/check_centering.py'],
 ['python3','session30-root/check_isospectral_metagraph.py'],
 ['python3','session30-root/check_random_controls.py'],
 ['python3','session30-root/check_selfconverse_completion.py'],
 ['python3','session30-root/check_completion_blue_gates.py'],
 ['python3','session30-root/check_blue_flux_census.py'],
 ['python3','session30-deficit-six/check.py'],
 ['python3','session30-nonregular/prove_closed_kernels.py'],
 ['python3','session30-nonregular/compare_symbolic.py','session30-decorated-engine/symbolic-kernels.json'],
 ['python3','session30-decorated-engine/check_prediction_stability.py'],
 ['python3','session30-deleted-drt/derive_a1.py'],
 ['python3','session30-deleted-drt/derive_a2_leading.py'],
 ['python3','session30-flat-universality/verify_formal.py'],
 ['python3','session30-flat-audit/audit.py'],
 ['python3','session30-collatz/verify.py','--repo-root',str(repo),'--out','base-replay.json'],
 ['python3','session30-collatz/verify_extension.py','--repo-root',str(repo),'--out','extension-replay.json'],
]
result={'scope':'Representative portable replay, not a new exhaustive mathematical audit','checks':[]}
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
with tempfile.TemporaryDirectory(prefix='session30-smoke-') as tmp:
 dst=pathlib.Path(tmp)/'packet';shutil.copytree(packet,dst)
 for cmd in commands:
  start=time.monotonic()
  r=subprocess.run(cmd,cwd=dst,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)
  result['checks'].append({'command':cmd,'exit_code':r.returncode,'seconds':round(time.monotonic()-start,4),'output_sha256':hashlib.sha256(r.stdout).hexdigest(),'output_tail':r.stdout.decode(errors='replace')[-1200:]})
  out.write_text(json.dumps(result,indent=2)+'\n')
  print(r.returncode,' '.join(cmd),flush=True)
  if r.returncode:raise RuntimeError(r.stdout.decode(errors='replace'))
result['status']='PASS';out.write_text(json.dumps(result,indent=2)+'\n')
