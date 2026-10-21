"""Pack producer output: keep path digits, retain the full-target commitment.

The no-search consumer reconstructs every target, in a separate implementation,
and must reproduce the SHA256 of the original path-plus-target byte stream.
"""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import zlib

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def sha(b):
    return hashlib.sha256(b).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('raw_bank',type=Path)
    args=parser.parse_args()
    raw=args.raw_bank.read_bytes()
    if len(raw)!=18124808 or raw[:8]!=b'D12SCR01':
        raise ValueError('complete producer bank required')
    words=raw[:8]+b''.join(raw[i:i+6] for i in range(8,len(raw),15))
    paths=[
        '05-knowledge/results/deficit_six_centers_20261021b/COMBINED_CUTOFF.md',
        '05-knowledge/results/deficit_six_centers_20261021b/combined_receipt.json',
        '05-knowledge/results/paired_black_leaf_reduction_20261020/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md',
    ]
    result={'schema':'D12SCR01','frame_count':262144,'balanced_count':151040,
            'record_bytes':15,'bank_bytes':len(raw),'bank_sha256':sha(raw),
            'word_bank_bytes':len(words),'word_bank_sha256':sha(words),
            'producer_sha256':sha((HERE/'census.cpp').read_bytes().replace(b'\r\n',b'\n')),
            'pins':{p:sha((ROOT/p).read_bytes().replace(b'\r\n',b'\n')) for p in paths},
            'bank_base85_zlib':base64.b85encode(zlib.compress(words,9)).decode('ascii')}
    (HERE/'inputs.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print('Packed',len(raw),'expanded bytes into',len(result['bank_base85_zlib']),'base85 characters.')


if __name__=='__main__':
    main()
