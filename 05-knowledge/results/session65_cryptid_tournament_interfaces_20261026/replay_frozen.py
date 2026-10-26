#!/usr/bin/env python3
"""Replay a fixed packet in separate normal/optimized copies; no bare assertions."""
import argparse,ast,hashlib,json,os,shutil,subprocess,sys,tempfile,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('packet');p.add_argument('tasks');p.add_argument('receipt');p.add_argument('--scratch',required=True);a=p.parse_args()
packet=Path(a.packet).resolve();tasks=json.loads(Path(a.tasks).read_text());scratch=Path(a.scratch).resolve();scratch.mkdir(parents=True,exist_ok=True)
def hashes(root):
 symlinks=[p for p in root.rglob('*') if p.is_symlink()]
 if symlinks:raise RuntimeError('Symlinks: '+str(symlinks))
 return {f.relative_to(root).as_posix():hashlib.sha256(f.read_bytes()).hexdigest() for f in root.rglob('*') if f.is_file() and '__pycache__' not in f.parts}
source=hashes(packet);runs=[];captured={}
for mode,flag in [('normal',[]),('optimized',['-O'])]:
 with tempfile.TemporaryDirectory(prefix=mode+'-',dir=scratch) as d:
  frozen=Path(d)/'packet';shutil.copytree(packet,frozen,ignore=shutil.ignore_patterns('__pycache__'))
  if hashes(frozen)!=source:raise RuntimeError('Frozen source mismatch')
  for task in tasks:
   name=task['script'];f=frozen/name
   tree=ast.parse(f.read_text());bare=sum(isinstance(n,ast.Assert) for n in ast.walk(tree))
   if bare:raise RuntimeError(f'Primary script has {bare} optimization-removed assertions: {name}')
   command=[sys.executable,*task.get('python_flags',[]),*flag,str(f),*task.get('args',[])]
   env=os.environ.copy();env.pop('PYTHONPATH',None);env.pop('PYTHONHOME',None);env['PYTHONDONTWRITEBYTECODE']='1';env['PYTHONHASHSEED']='0';env.update(task.get('environment',{}))
   start=time.monotonic();r=subprocess.run(command,cwd=frozen/task.get('cwd','.'),env=env,capture_output=True,timeout=task.get('timeout',600));elapsed=time.monotonic()-start
   if r.returncode:raise RuntimeError(f'{name} {mode} exit {r.returncode}\n{r.stderr.decode(errors="replace")}')
   if mode=='normal':captured[name]=(r.stdout,r.stderr)
   elif captured[name]!=(r.stdout,r.stderr):
    Path(a.receipt+'.failure.json').write_text(json.dumps({'status':'FAIL','script':name,'prior_runs':runs,'normal_stdout':captured[name][0].decode(errors='replace'),'normal_stderr':captured[name][1].decode(errors='replace'),'optimized_stdout':r.stdout.decode(errors='replace'),'optimized_stderr':r.stderr.decode(errors='replace')},indent=2)+'\n')
    raise RuntimeError(f'Normal/optimized output mismatch: {name}')
   stored=task.get('stored_stdout')
   if stored and r.stdout!=(frozen/stored).read_bytes():raise RuntimeError(f'Stored stdout mismatch: {name} -> {stored}')
   after=hashes(frozen)
   if after!=source:raise RuntimeError(f'Replay mutates packet bytes: {name}: '+str({k for k in set(after)|set(source) if after.get(k)!=source.get(k)}))
   runs.append(dict(script=name,mode=mode,source_sha256=source[name],command=command,exit_code=0,bare_assertions=bare,elapsed_seconds=round(elapsed,6),stdout_sha256=hashlib.sha256(r.stdout).hexdigest(),stderr_sha256=hashlib.sha256(r.stderr).hexdigest(),stored_stdout=stored,environment_overrides=task.get('environment',{}),source_unchanged=True))
if hashes(packet)!=source:raise RuntimeError('Original packet changed during replay')
receipt=dict(status='PASS',python=sys.version,script_count=len(tasks),successful_runs=len(runs),normal_optimized_outputs_match=True,original_sources_unchanged=True,runs=runs)
Path(a.receipt).write_text(json.dumps(receipt,indent=2)+'\n');print(f'PASS: {len(tasks)} scripts, {len(runs)} normal/optimized executions; exact outputs and unchanged packet')
