"""Frozen packet integrity and optional nonmutating normal/-O exact replay."""
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


def digest(path,lf=False):
 b=path.read_bytes()
 return sha256(b.replace(b'\r\n',b'\n') if lf else b).hexdigest()


def main():
 p=argparse.ArgumentParser();p.add_argument('--replay',action='store_true');args=p.parse_args()
 m=json.loads((HERE/'manifest.json').read_text())
 files={f.name for f in HERE.iterdir() if f.is_file() and f.name!='manifest.json'}
 need(files==set(m['files']),'packet inventory')
 for name,value in m['files'].items():need(digest(HERE/name)==value,'packet '+name)
 for name,value in m['dependencies_lf_sha256'].items():need(digest(ROOT/name,True)==value,'dependency '+name)
 for stem in ('certificate','full_response','audit'):
  need((HERE/(stem+'.json')).read_bytes()==(HERE/(stem+'.optimized.json')).read_bytes(),'normal/-O '+stem)
 cert=json.loads((HERE/'certificate.json').read_text());audit=json.loads((HERE/'audit.json').read_text())
 for name,value in cert['source_sha256'].items():need(digest(HERE/name)==value,'producer '+name)
 for name,value in audit['source_sha256'].items():
  path=HERE/name if name=='audit.py' else HERE.parent/name
  need(digest(path)==value,'audit source '+name)
 need(cert['tables_sha256']==audit['tables_sha256']==digest(HERE/'full_response.json'),'complete tables')
 need(cert['actual_compiled_queries']==37210 and cert['full_response_states']==40,'actual query universe')
 need(cert['complete_labelled_quotients']==75 and cert['small_reversal_states']==1098,'complete finite universe')
 need(cert['small_actual_queries']==4306 and cert['literal_complete_polynomials']==8,'decoder and word checks')
 need(audit['full_quartic_coordinates']==16432 and audit['whole_converse_pairs']==20,'quartic audit universe')
 need(audit['mixed_quartic']['coefficient']==-2 and audit['nested_separator']['counts']==[251,249],'hostile controls')
 tables=json.loads((HERE/'full_response.json').read_text())['states']
 for s in tables:
  n=len(s['rows']);size=1+n+n*(n-1)//2
  need(len(s['masks'])==size,'basis')
  need(len(s['coefficients'])==len(s['values'])==len(s['contrast_from_top0_epsilon0'])==size**2,'full named response inventory')
 if args.replay:
  with tempfile.TemporaryDirectory(prefix='module-cube-replay-') as folder:
   for opt in (False,True):
    prefix=[sys.executable]+(['-O'] if opt else [])
    output=Path(folder)/'certificate.json';table=Path(folder)/'full_response.json'
    subprocess.run(prefix+[str(HERE/'certify.py'),'--output',str(output),'--tables',str(table)],check=True)
    need(output.read_bytes()==(HERE/'certificate.json').read_bytes(),'certificate replay')
    need(table.read_bytes()==(HERE/'full_response.json').read_bytes(),'complete response replay')
    output=Path(folder)/'audit.json'
    subprocess.run(prefix+[str(HERE/'audit.py'),'--output',str(output)],check=True)
    need(output.read_bytes()==(HERE/'audit.json').read_bytes(),'independent audit replay')
 print(json.dumps(dict(status='PASS',scope='S1-S3 disjoint modules; nested generality OPEN',
                       files=len(files),dependencies=len(m['dependencies_lf_sha256']),
                       actual_compiled_queries=37210,full_named_responses=40,
                       independent_quartic_coordinates=16432,replayed=args.replay),indent=2))


if __name__=='__main__':main()
