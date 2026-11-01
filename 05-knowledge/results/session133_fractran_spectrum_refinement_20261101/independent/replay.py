"""Portable independent-audit driver; normal and -O receipts must match."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
TESTS=[('compiler','audit_compiler.py'),('path','audit_path_counter.py'),
       ('linear','audit_linear_stream.py'),('direct','audit_direct_programs.py'),
       ('basis','audit_forbidden_basis.py'),('seven','audit_seven_counters.py'),
       ('fibres','audit_pointed_fibres.py')]

def main():
    pins=HERE/'INPUT_PINS.json'
    if pins.exists():
        for relative,expected in json.loads(pins.read_text())['inputs'].items():
            source=HERE.parent/relative
            if hashlib.sha256(source.read_bytes()).hexdigest()!=expected:
                raise RuntimeError(f'input hash changed: {relative}')
    entries=[]
    for stem,script in TESTS:
        receipts=[]
        for mode,flags in [('normal',[]),('optimized',['-O'])]:
            result=subprocess.run([sys.executable,*flags,str(HERE/script)],
                                  cwd=HERE,capture_output=True,check=True)
            json.loads(result.stdout)
            path=HERE/f'{stem}.{mode}.json';path.write_bytes(result.stdout)
            receipts.append(result.stdout)
        if receipts[0]!=receipts[1]:raise RuntimeError(f'{stem}: normal/-O differ')
        entries.append({'test':stem,'script':script,'status':'PASS',
                        'script_sha256':hashlib.sha256((HERE/script).read_bytes()).hexdigest(),
                        'receipt_sha256':hashlib.sha256(receipts[0]).hexdigest()})
        print(stem+': PASS',flush=True)
    (HERE/'REPLAY.json').write_text(json.dumps({'status':'PASS','tests':entries},indent=2)+'\n')

if __name__=='__main__':main()
