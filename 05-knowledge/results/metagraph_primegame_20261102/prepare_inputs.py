#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def main():
 pin=json.loads((HERE/'SOURCES.json').read_text())['public_commit']
 target=Path(sys.argv[1]);target.mkdir(parents=True,exist_ok=True)
 refs=json.loads((HERE/'INPUTS.json').read_text())
 for row in refs:
  n=row['n'];raw=subprocess.check_output(['git','show',pin+':'+row['source']],cwd=HERE)
  codes=raw.decode().splitlines() if n<9 else [l.split()[-1] for l in raw.decode().splitlines()]
  b=('\n'.join(codes)+'\n').encode();need(hashlib.sha256(b).hexdigest()==row['input_sha256'],'source code bank hash')
  (target/f'quotients_n{n}.txt').write_bytes(b)
  source=f'05-knowledge/results/spectrum_analogies_nearjoin_20261102/data/census2_{n}.txt'
  old=subprocess.check_output(['git','show',pin+':'+source],cwd=HERE).decode().splitlines()
  d={l.split()[-1]:l.split()[0] for l in old}
  need(set(d)==set(codes) and len(d)==row['count'],'independent inherited census universe')
  (target/f'expected_n{n}.txt').write_text('\n'.join(d[c] for c in codes)+'\n')
 print('PASS restored pinned class banks and inherited independent H counts')
if __name__=='__main__':main()
