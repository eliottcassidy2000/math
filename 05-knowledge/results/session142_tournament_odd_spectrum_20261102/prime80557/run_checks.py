#!/usr/bin/env python3
import argparse,json,subprocess,sys,hashlib
from pathlib import Path
P=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--work-dir',type=Path,required=True);args=p.parse_args();W=args.work_dir.resolve();W.mkdir(parents=True,exist_ok=True)
def run(cmd,path):
 r=subprocess.run([str(x) for x in cmd],capture_output=True,check=True);path.write_bytes(r.stdout);return r.stdout
for script,receipt,extra in [('verify_source_pins.py','source_pins',[]),('freeze_and_verify.py','verify',[]),('sources/square_sum_spectrum_overnight_20261031d/source_inverse.py','inverse',[]),('sources/square_sum_spectrum_overnight_20261031d/bank_query.py','query',[80557,100000007,'--cap',56])]:
 ext='jsonl' if receipt=='query' else 'json'
 a=run([sys.executable,'-B',P/script,*extra],W/f'{receipt}.normal.{ext}');b=run([sys.executable,'-B','-O',P/script,*extra],W/f'{receipt}.optimized.{ext}')
 if a!=b or a!=(P/f'{receipt}.normal.{ext}').read_bytes():raise ValueError(receipt)
for opt in ('O0','O2'):
 exe=W/('subset_dp_'+opt);run(['g++','-'+opt,'-std=c++17','-Wall','-Wextra','-pedantic',P/'subset_dp.cpp','-o',exe],W/f'compile.{opt}.txt')
 r=run([exe,P/'adjacency_80557.txt'],W/f'subset_dp.{opt}.json')
 if r!=(P/'subset_dp.O2.json').read_bytes():raise ValueError(opt)
 for mode,flags in [('normal',[]),('optimized',['-O'])]:
  r=run([sys.executable,'-B',*flags,P/'hostile_controls.py',exe],W/f'hostile.{opt}.{mode}.json')
  if r!=(P/'hostile.normal.json').read_bytes():raise ValueError('hostile')
for t in (80557,100000007):
 for mode,flags in [('normal',[]),('optimized',['-O'])]:
  dest=W/f'reconstructed_{t}.{mode}.txt'
  run([sys.executable,'-B',*flags,P/'construct.py',t,'--adjacency',dest,'--recipe',W/f'reconstructed_{t}.{mode}.json'],W/f'construct_{t}.{mode}.json')
  if dest.read_bytes()!=(P/f'adjacency_{t}.txt').read_bytes():raise ValueError('reconstruction')
for line in (P/'SHA256SUMS').read_text().splitlines():
 digest,name=line.split('  ',1)
 if hashlib.sha256((P/name).read_bytes()).hexdigest()!=digest:raise ValueError(name)
print(json.dumps({'status':'passed','portable_replay':True,'full_subset_DP_orders':[19],'direct_graph_partition_orders':[19,29],'source_inverse_queries':80498},sort_keys=True))
