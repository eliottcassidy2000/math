"""Frozen proof packet verification and nonmutating optional exact replay."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def need(ok,message):
 if not ok:raise RuntimeError(message)


def digest(p,lf=False):
 b=p.read_bytes()
 return sha256(b.replace(b'\r\n',b'\n') if lf else b).hexdigest()


def main():
 p=argparse.ArgumentParser();p.add_argument('--replay',action='store_true');args=p.parse_args()
 m=json.loads((HERE/'manifest.json').read_text())
 files={f.name for f in HERE.iterdir() if f.is_file() and f.name!='manifest.json'}
 need(files==set(m['files']),'file inventory')
 for name,value in m['files'].items():need(digest(HERE/name)==value,'packet hash '+name)
 for name,value in m['dependencies_lf_sha256'].items():need(digest(ROOT/name,True)==value,'dependency '+name)
 for stem in ['census','certificate']:
  receipt=json.loads((HERE/(stem+'.json')).read_text())
  need((HERE/(stem+'.json')).read_bytes()==(HERE/(stem+'.optimized.json')).read_bytes(),'normal/-O '+stem)
  for name,value in receipt['source_sha256'].items():need(digest(HERE/name)==value,'producer '+name)
 census=json.loads((HERE/'census.json').read_text());cert=json.loads((HERE/'certificate.json').read_text())
 need(digest(ROOT/m['bank_path'])==census['bank_sha256']==m['bank_sha256'],'bank')
 need(census['classes']==198947,'bank universe')
 need(all(r['cover_and_shared_zero']==0 for r in census['counts'] if r['n'] in (4,6,8)),'finite even obstruction')
 need(cert['small_labelled_integer_controls']==1099 and cert['literal_hamilton_counts']==708,'literal controls')
 need(cert['actual_compiled_counts']==1268 and cert['full_cubic_actual_queries']==312,'actual query controls')
 tables=HERE/'full_named_cubic.json'
 need(tables.read_bytes()==(HERE/'full_named_cubic.optimized.json').read_bytes(),'complete cubic normal/-O')
 need(digest(tables)==cert['tables_sha256'],'cubic receipt')
 need(all(len(x['coefficients'])==1753 for x in json.loads(tables.read_text())['cases']),'full named cubic coordinates')
 if args.replay:
  with tempfile.TemporaryDirectory(prefix='parity-replay-') as folder:
   for opt in [False,True]:
    for stem in ['census','certificate']:
     output=Path(folder)/(stem+'.json')
     cmd=[sys.executable]+(['-O'] if opt else [])+[str(HERE/(stem+'.py')),'--output',str(output)]
     if stem=='certificate':tmp=Path(folder)/'tables.json';cmd+=['--tables',str(tmp)]
     subprocess.run(cmd,check=True)
     need(output.read_bytes()==(HERE/(stem+'.json')).read_bytes(),'replay '+stem)
     if stem=='certificate':need(tmp.read_bytes()==tables.read_bytes(),'complete cubic replay')
 print(json.dumps(dict(status='PASS',scope='E1-E4: every single nontrivial module converse is visible; general injectivity not asserted',
                       files=len(files),dependencies=len(m['dependencies_lf_sha256']),
                       classes=198947,actual_compiled_counts=1268,full_named_cubics=2,replayed=args.replay),indent=2))


if __name__=='__main__':main()
