"""Complete inherited class bank; matrix identities, not a complete response census."""
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BANK=ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz'


def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'census.json');args=p.parse_args()
 pins={f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['census.py','matrix.cpp']}
 data=json.loads(gzip.decompress(BANK.read_bytes()))
 lines=[f"{n} {rec['code']}" for n in range(2,10) for rec in data['members'][str(n)]]
 if len(lines)!=198947:raise RuntimeError('complete bank count')
 with tempfile.TemporaryDirectory(prefix='even-obstruction-') as temp:
  exe=Path(temp)/'matrix.exe'
  subprocess.run(['g++','-O3','-std=c++17','-Wall','-Wextra','-Werror',str(HERE/'matrix.cpp'),'-o',str(exe)],check=True)
  out=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,capture_output=True,check=True)
 records=[]
 for line in out.stdout.splitlines():
  n,count,b0,d0,h0,joint0,all0=map(int,line.split())
  if count!=len(data['members'][str(n)]):raise RuntimeError('bank order count')
  records.append(dict(n=n,classes=count,endpoint_balanced=b0,cover_skew_zero=d0,
                      shared_cubics_zero=h0 if n%2==0 else None,
                      cover_and_shared_zero=joint0 if n%2==0 else None,
                      all_three_zero=all0 if n%2==0 else None))
 if pins!={f:sha256((HERE/f).read_bytes()).hexdigest() for f in pins}:raise RuntimeError('source changed')
 result=dict(status='PASS finite complete bank; all-order claim rests on E1-E4 proof',
             classes=len(lines),counts=records,source_sha256=pins,bank_sha256=sha256(BANK.read_bytes()).hexdigest())
 args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
 print(json.dumps(result,indent=2))


if __name__=='__main__':main()
