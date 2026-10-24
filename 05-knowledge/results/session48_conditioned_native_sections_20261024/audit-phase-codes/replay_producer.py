#!/usr/bin/env python3
"""Replay the byte-frozen producer; do not import its code."""
from pathlib import Path
from hashlib import sha256
import argparse,json,subprocess,sys
HERE=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    results={}
    for label,flags in [('normal',['-I','-B']),('optimized',['-I','-O','-B'])]:
        done=subprocess.run([sys.executable,*flags,str(HERE/'inputs/producer/check.py')],capture_output=True,text=True,check=True)
        results[label]={'stdout':done.stdout,'stderr':done.stderr,'exit_code':done.returncode}
    if results['normal']!=results['optimized']:raise RuntimeError('Producer modes differ')
    done=subprocess.run([sys.executable,'-I','-B',str(HERE/'inputs/producer/verify_manifest.py')],capture_output=True,text=True,check=True)
    results['manifest']={'stdout':done.stdout,'stderr':done.stderr,'exit_code':done.returncode}
    results['producer_manifest_sha256']=sha256((HERE/'inputs/producer/MANIFEST.json').read_bytes()).hexdigest()
    results['producer_receipt_sha256']=sha256((HERE/'inputs/producer/receipt.json').read_bytes()).hexdigest()
    if args.write:(HERE/'PRODUCER_REPLAY.json').write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')
    elif json.loads((HERE/'PRODUCER_REPLAY.json').read_text())!=results:raise RuntimeError('Producer replay differs from frozen result')
    print('PASS frozen producer replay: normal, optimized, and closed manifest')

if __name__=='__main__':main()
