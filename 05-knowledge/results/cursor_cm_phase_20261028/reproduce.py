"""Replay exact receipts in a temporary mirror, without changing the saved results."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
REL=HERE.relative_to(ROOT)
DEPS=[Path('05-knowledge/results/joint_port_germs_20261028/port_core.py'),
      Path('05-knowledge/results/relation_tower_frobenius_ports_20261028/check_root_switch.py')]

def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--direct-eight',action='store_true')
 args=parser.parse_args()
 tasks=[('endpoint_kernel.py',[]),('cm_series.py',[]),('unit_clocks.py',[]),
        ('run_probe.py',[]),('cursor_symmetry.py',[])]
 if args.direct_eight:tasks.append(('cursor_probe.py',['--max-order','8']))
 inputs=list(HERE.glob('*.py'))+[ROOT/p for p in DEPS]
 report={'status':'FINITE-EXACT replay','python':sys.version.split()[0],
         'source_sha256_LF':{p.relative_to(ROOT).as_posix():digest(p) for p in inputs},
         'runs':[]}
 with tempfile.TemporaryDirectory(prefix='math_cursor_') as tmp:
  base=Path(tmp).resolve();target=base/REL;target.mkdir(parents=True)
  for p in HERE.glob('*.py'):shutil.copyfile(p,target/p.name)
  for p in DEPS:
   dest=base/p;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/p,dest)
  for name,extra in tasks:
   expected=json.loads((HERE/Path(name).with_suffix('.json')).read_text(encoding='utf-8'))
   for optimized in (False,True):
    command=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(target/name)]+extra
    result=subprocess.run(command,cwd=base,text=True,capture_output=True)
    if result.returncode:raise RuntimeError(result.stdout+'\n'+result.stderr)
    produced=json.loads((target/Path(name).with_suffix('.json')).read_text(encoding='utf-8'))
    if produced!=expected:raise RuntimeError('receipt mismatch: '+name)
    row={'script':name,'optimized':optimized,'receipt_sha256_LF':digest(HERE/Path(name).with_suffix('.json'))}
    report['runs'].append(row);print(json.dumps(row),flush=True)
 (HERE/'audit_receipt.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
 print('PASS: exact receipts agree in normal and optimized execution.',flush=True)
if __name__=='__main__':main()
