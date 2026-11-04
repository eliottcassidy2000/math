"""Verify frozen merge evidence; optionally rerun exact normal/-O producers."""
import argparse
import gzip
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def need(ok,message):
 if not ok:raise RuntimeError(message)


def read(path):
 if os.name=='nt':path=Path('\\\\?\\'+str(path.resolve()))
 return path.read_bytes()


def digest(path,lf=False):
 b=read(path)
 return sha256(b.replace(b'\r\n',b'\n') if lf else b).hexdigest()


def main():
 p=argparse.ArgumentParser();p.add_argument('--replay',action='store_true');args=p.parse_args()
 m=json.loads(read(HERE/'manifest.json'));files={x.name for x in HERE.iterdir() if x.is_file() and x.name!='manifest.json'}
 need(files==set(m['files']),'packet inventory')
 for name,h in m['files'].items():need(digest(HERE/name)==h,'packet '+name)
 for name,h in m['dependencies_lf_sha256'].items():need(digest(ROOT/name,True)==h,'dependency '+name)
 receipt=json.loads(read(HERE/'receipt.json'))
 need(read(HERE/'receipt.json')==read(HERE/'receipt.optimized.json'),'normal/-O receipt')
 table=read(HERE/'full_response.json.gz')
 need(table==read(HERE/'full_response.optimized.json.gz'),'normal/-O full tables')
 need(sha256(table).hexdigest()==receipt['tables_sha256'],'archive hash')
 raw=gzip.decompress(table);need(sha256(raw).hexdigest()==receipt['uncompressed_sha256'],'lossless data hash')
 for name,h in receipt['sources'].items():need(digest(HERE.parent/name)==h,'exact source '+name)
 data=json.loads(raw)
 need(data['denominator']==1024 and len(data['masks'])==16 and len(data['queries'])==259,'named basis')
 need(len(data['coefficients'])==len(data['values'])==len(data['old_edge_walsh_numerators'])==1024,'labelled universe')
 need(all(len(v)==256 for v in data['coefficients']+data['old_edge_walsh_numerators']),'complete polynomial and transform')
 need(all(len(v)==259 for v in data['values']),'actual counts')
 need(receipt['actual_queries']==265216 and receipt['quartic_normal_form_checks']==102400,'verification universe')
 need(receipt['odd_old_support_quartic_zeros']==51200 and receipt['literal_full_polynomials']==12,'independent controls')
 if args.replay:
  with tempfile.TemporaryDirectory(prefix='even-merge-replay-') as folder:
   for opt in (False,True):
    output=Path(folder)/'receipt.json';tables=Path(folder)/'tables.gz'
    cmd=[sys.executable]+(['-O'] if opt else [])+[str(HERE/'exact.py'),'--output',str(output),'--tables',str(tables)]
    subprocess.run(cmd,check=True)
    need(read(output)==read(HERE/'receipt.json'),'receipt replay')
    need(read(tables)==table,'full named table replay')
 print(json.dumps(dict(status='PASS',files=len(files),dependencies=len(m['dependencies_lf_sha256']),
                       full_named_responses=1024,actual_queries=265216,replayed=args.replay),indent=2))


if __name__=='__main__':main()
