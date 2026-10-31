"""Run all exact proof checks normally and optimized against frozen receipts."""
from pathlib import Path
import subprocess,sys,json,hashlib,tempfile,argparse
ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path);args=parser.parse_args()
expected=json.loads((ROOT/'receipt_hashes.json').read_text())
results=[]
with tempfile.TemporaryDirectory(prefix='bandwidth-check-') as temp:
 outdir=args.output_dir or Path(temp)
 outdir.mkdir(parents=True,exist_ok=True)
 for script,stem,suffix in [('transfer.py','transfer','out'),('verify.py','verify','json'),('independent_check.py','independent','json')]:
  outputs=[]
  for mode,flags in [('normal',[]),('optimized',['-O'])]:
   name=f'{stem}.{mode}.{suffix}'
   result=subprocess.run([sys.executable,*flags,'-B',str(ROOT/script)],cwd=temp,capture_output=True,check=True)
   if result.stderr:raise RuntimeError(result.stderr.decode())
   digest=hashlib.sha256(result.stdout).hexdigest()
   if digest!=expected[f'{stem}.normal.{suffix}']:raise RuntimeError(('receipt mismatch',name,digest))
   (outdir/name).write_bytes(result.stdout);outputs.append(result.stdout)
   results.append(dict(receipt=name,sha256=digest))
  if outputs[0]!=outputs[1]:raise RuntimeError('optimization changed output')
 print(json.dumps(dict(status='PASS',executions=6,standard_library_only=True,normal_optimized_equal=True,receipts=results),sort_keys=True))
