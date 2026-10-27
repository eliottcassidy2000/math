import concurrent.futures as cf, subprocess, pathlib, time, json, sys
root=pathlib.Path(__file__).resolve().parent
which=sys.argv[1] if len(sys.argv)>1 else 'normal'
jobs=sorted((root/f'jobs-{which}').glob('*.job'))
def one(job):
 t=time.monotonic(); p=subprocess.run([str(root/'tc_window_scan_portable'),str(job),'check'],capture_output=True,text=True)
 return {'job':job.name,'returncode':p.returncode,'seconds':time.monotonic()-t,'stdout':p.stdout,'stderr':p.stderr}
with (root/f'replay-{which}.jsonl').open('w') as f, cf.ThreadPoolExecutor(max_workers=2) as ex:
 for r in ex.map(one,jobs):
  f.write(json.dumps(r)+'\n');f.flush(); print(r['job'],r['returncode'],round(r['seconds'],2),flush=True)
