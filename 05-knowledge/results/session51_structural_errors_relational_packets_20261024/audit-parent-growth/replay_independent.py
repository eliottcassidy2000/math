#!/usr/bin/env python3
from pathlib import Path
import subprocess,json,time
B=Path(__file__).resolve().parent
results=[]
for label,flags in [('O2',['-O2']),('O3',['-O3','-DNDEBUG','-DFULL_REFINEMENT'])]:
 for src,exe in [('independent.cpp','independent'),('incoming_trace.cpp','incoming')]:
  subprocess.run(['g++','-std=c++17',*flags,str(B/src),'-o',str(B/f'{exe}-{label}')],check=True,capture_output=True)
 t=time.monotonic();r=subprocess.run([str(B/f'independent-{label}'),str(B/'producer/parents-eight.txt'),str(B/'producer/parents-nine.txt'),str(B/f'independent-{label}.json')],check=True,capture_output=True,text=True);(B/f'independent-{label}.out').write_text(r.stdout+r.stderr);results.append({'mode':label,'operation':'independent_atlas','elapsed_seconds':time.monotonic()-t})
 t=time.monotonic();r=subprocess.run([str(B/f'incoming-{label}'),str(B/'producer/parents-nine.txt'),str(B/'producer/trace-normal-per-parent.tsv'),str(B/f'incoming-{label}.tsv'),str(B/f'incoming-{label}.json')],check=True,capture_output=True,text=True);(B/f'incoming-{label}.out').write_text(r.stdout+r.stderr);results.append({'mode':label,'operation':'incoming_trace','elapsed_seconds':time.monotonic()-t})
 print(label,'PASS',flush=True)
if (B/'independent-O2.json').read_bytes()!=(B/'independent-O3.json').read_bytes():raise RuntimeError('independent mode mismatch')
if (B/'incoming-O2.tsv').read_bytes()!=(B/'incoming-O3.tsv').read_bytes():raise RuntimeError('incoming per-parent mode mismatch')
(B/'independent-replays.json').write_text(json.dumps({'status':'PASS','runs':results},indent=2)+'\n')
