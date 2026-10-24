#!/usr/bin/env python3
"""Replay both independent exhaustive algorithms, native inheritance and checks."""
import argparse,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);args=p.parse_args();repo=args.repo.resolve()
def run(cmd,output=None):
 print('Running: '+' '.join(map(str,cmd)),flush=True)
 if output:
  with output.open('w')as out:subprocess.run(list(map(str,cmd)),stdout=out,check=True)
 else:subprocess.run(list(map(str,cmd)),check=True)
build=HERE/'replay-build';build.mkdir(exist_ok=True)
for name in ('census','census_support','independent_bfs','independent_depth'):
 run(['g++','-std=c++17','-O3','-o',build/name,HERE/(name+'.cpp')])
run([build/'census',HERE/'input.txt',0,14,22],HERE/'budget22.jsonl')
run([build/'census_support',HERE/'input.txt',0,14,22,8],HERE/'budget22-support8.jsonl')
run([build/'independent_bfs',HERE/'input.txt',0,14,22],HERE/'independent-budget22-support8.jsonl')
run([build/'independent_depth',HERE/'input.txt',0,14,22],HERE/'independent-depth-budget22-support8.jsonl')
for name in ('construct_route.py','move_bank.py','deletion_probe.py'):
 run([sys.executable,'-B',HERE/name,'--repo',repo],HERE/(name+'.replay.log'))
for mode in ('normal','optimized'):
 cmd=[sys.executable]+(['-O']if mode=='optimized'else[])+['-B']
 run(cmd+[repo/'05-knowledge/results/native_packet_decomposition_20261023/native/check.py'],HERE/('inherited_native.'+mode+'.log'))
 run(cmd+[HERE/'check.py','--repo',repo,'--output',HERE/(mode+'.json')])
 run(cmd+[HERE/'audit_incoming.py','--repo',repo,'--output',HERE/('incoming.'+mode+'.json')])
 run(cmd+[HERE/'audit_deletions.py','--output',HERE/('deletions.'+mode+'.json')])
for stem in ('','incoming.','deletions.'):
 if (HERE/(stem+'normal.json')).read_bytes()!=(HERE/(stem+'optimized.json')).read_bytes():raise RuntimeError('normal/optimized checker receipts differ: '+stem)
print('PASS complete28 joint8slot22clock exhaustive replay and independent audit',flush=True)
