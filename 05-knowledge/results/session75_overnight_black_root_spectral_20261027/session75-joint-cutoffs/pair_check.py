#!/usr/bin/env python3
"""All-rectangle source-pair birth check at a declared completed scale."""
from pathlib import Path
import argparse,hashlib,json
import literal_verify as verifier
ROOT=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('N',type=int);ap.add_argument('--suffix',default='');args=ap.parse_args();N=args.N
    rpath=ROOT/'results'/f'native_sources_N{N}.json';rhash=hashlib.sha256(rpath.read_bytes()).hexdigest()
    receipt=json.loads((ROOT/'results'/f'receipt_N{N}.json').read_text());literal=json.loads((ROOT/'results'/f'literal_receipt_N{N}.json').read_text())
    verifier.require(rhash==receipt['source_ledger_sha256']==literal['source_ledger_sha256'],'source ledger bound to fresh full reconstruction')
    result=verifier.pair_birth_verify(json.loads(rpath.read_text()),N,receipt)
    result.update(N=N,status='PASS',source_ledger_sha256=rhash)
    (ROOT/'results'/f'pair_birth_receipt_N{N}{args.suffix}.json').write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
