"""Replay producer and independent checks in an isolated temporary directory."""
from pathlib import Path
import json,shutil,subprocess,tempfile
from reconstruct import reconstruct
D=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='three-core-replay-') as tmp:
 t=Path(tmp)
 for p in D.iterdir():
  if p.is_file():shutil.copy2(p,t/p.name)
 reconstruct(t/'RECIPES.bin')
 result=subprocess.check_output(['python3',str(t/'verify.py')],text=True,cwd=t)
 if json.loads(result)!=json.loads((D/'VERIFY.json').read_text()):raise RuntimeError('Producer verification changed')
 subprocess.run(['g++','-O2','-std=c++17',str(D/'audit/audit.cpp'),'-o',str(t/'audit')],check=True,cwd=t)
 result=subprocess.check_output([str(t/'audit'),str(t/'CORE_STATES.txt'),str(t/'RECIPES.bin')],text=True,cwd=t)
 if [json.loads(s) for s in result.splitlines()]!=[json.loads(s) for s in (D/'audit/RESULT.jsonl').read_text().splitlines()]:raise RuntimeError('Independent verification changed')
 result=subprocess.check_output(['python3',str(D/'audit/check_exceptions.py'),str(t/'EXCEPTIONS.json')],text=True,cwd=t)
 if json.loads(result)!=json.loads((D/'audit/EXCEPTIONS_AUDIT.json').read_text()):raise RuntimeError('Exception verification changed')
 print(json.dumps({'status':'PASS','certificate_reconstructed_and_hash_verified':True,'producer_receipt_matches':True,'independent_receipt_matches':True,'exception_receipt_matches':True},indent=2))
