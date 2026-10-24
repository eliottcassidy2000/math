#!/usr/bin/env python3
"""Compile and replay physical frozen production copies, O2 and O3/NDEBUG."""
from pathlib import Path
import subprocess,json,time,hashlib
B=Path(__file__).resolve().parent
rows=[]
source_hashes={name:hashlib.sha256((B/'producer'/name).read_bytes()).hexdigest() for name in ('build_nine.cpp','trace_ten.cpp')}
for label,flags in [('O2',['-O2']),('O3-ndebug-full',['-O3','-DNDEBUG','-DFULL_REFINEMENT'])]:
 for source,exe in [('build_nine.cpp','build'),('trace_ten.cpp','trace')]:
  command=['g++','-std=c++17',*flags,str(B/'producer'/source),'-o',str(B/f'{exe}-{label}')]
  subprocess.run(command,check=True,capture_output=True)
 t=time.monotonic(); cmd=[str(B/f'build-{label}'),str(B/'producer/parents-eight.txt'),str(B/f'parents-nine-{label}.txt'),str(B/f'build-{label}.json')]
 r=subprocess.run(cmd,check=True,capture_output=True,text=True); (B/f'build-{label}.out').write_text(r.stdout+r.stderr)
 if (B/f'parents-nine-{label}.txt').read_bytes()!=(B/'producer/parents-nine.txt').read_bytes():raise RuntimeError('atlas mismatch')
 if json.loads((B/f'build-{label}.json').read_text())!=json.loads((B/'producer/build-normal.json').read_text()):raise RuntimeError('build receipt mismatch')
 rows.append({'mode':label,'operation':'build','elapsed_seconds':time.monotonic()-t,'atlas_sha256':hashlib.sha256((B/f'parents-nine-{label}.txt').read_bytes()).hexdigest()})
 t=time.monotonic(); cmd=[str(B/f'trace-{label}'),str(B/'producer/parents-nine.txt'),str(B/f'trace-{label}.json'),str(B/f'trace-{label}.tsv')]
 r=subprocess.run(cmd,check=True,capture_output=True,text=True);(B/f'trace-{label}.out').write_text(r.stdout+r.stderr)
 if (B/f'trace-{label}.tsv').read_bytes()!=(B/'producer/trace-normal-per-parent.tsv').read_bytes():raise RuntimeError('per-parent trace mismatch')
 actual=json.loads((B/f'trace-{label}.json').read_text());expected=json.loads((B/'producer/trace-normal.json').read_text())
 for k,v in expected.items():
  if k in ('refinement_rounds','singleton_early_exits'):continue
  if actual[k]!=v:raise RuntimeError('trace receipt mismatch: '+k)
 rows.append({'mode':label,'operation':'trace','elapsed_seconds':time.monotonic()-t,'per_parent_sha256':hashlib.sha256((B/f'trace-{label}.tsv').read_bytes()).hexdigest()})
 print(label,'PASS',flush=True)
(B/'producer-replays.json').write_text(json.dumps({'status':'PASS','physical_copy':True,'source_hashes':source_hashes,'runs':rows},indent=2)+'\n')
