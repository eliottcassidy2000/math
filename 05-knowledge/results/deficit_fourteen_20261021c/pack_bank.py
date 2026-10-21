"""Save the independent complete path-and-target producer stream."""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import zlib
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def sha(x):return hashlib.sha256(x).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('raw_bank',type=Path);a=p.parse_args()
    raw=a.raw_bank.read_bytes()
    if len(raw)!=8+2048*8*19 or raw[:8]!=b'D14SCR01':raise ValueError('complete bank required')
    paths=['05-knowledge/results/deficit_twelve_20261021c/check.py',
           '05-knowledge/results/deficit_twelve_20261021c/receipt.json',
           '05-knowledge/results/deficit_six_centers_20261021b/combined_receipt.json',
           '05-knowledge/results/paired_black_leaf_reduction_20261020/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md']
    data={'schema':'D14SCR01','frames':2048,'bytes':len(raw),'sha256':sha(raw),
          'producer_sha256':sha((HERE/'census.cpp').read_bytes().replace(b'\r\n',b'\n')),
          'pins':{p:sha((ROOT/p).read_bytes().replace(b'\r\n',b'\n')) for p in paths},
          'bank_base85_zlib':base64.b85encode(zlib.compress(raw,9)).decode('ascii')}
    (HERE/'inputs.json').write_text(json.dumps(data,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print('saved',len(raw),'raw bytes')


if __name__=='__main__':main()
