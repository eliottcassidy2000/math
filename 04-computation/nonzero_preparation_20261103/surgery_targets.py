#!/usr/bin/env python3
"""Eight literal positive joint-gap2 targets recovered from the full core bank.
No absence or complete spectrum claim is made by this small witness checker.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import surgery_ready as ready
import surgery_native as native

WITNESSES = [{'H': 38609, 'gap': 2, 'order': 15, 'core': [2, 4, 201, 211, 231, 207, 131, 3], 'ports': [0, 2], 'word': [0, 1, 1, 0, 1], 'rows': [24828, 16633, 24626, 16420, 16392, 16400, 16444, 16508, 24831, 16895, 17407, 26623, 20479, 5882, 8192], 'root': 13, 'aux': 14}, {'H': 38611, 'gap': 2, 'order': 14, 'core': [6, 4, 440, 499, 483, 451, 391, 259, 3], 'ports': [0], 'word': [0, 0, 1], 'rows': [12792, 8697, 8259, 8196, 8204, 8220, 8248, 8316, 8444, 12799, 13311, 10239, 2558, 4096], 'root': 12, 'aux': 13}, {'H': 38615, 'gap': 2, 'order': 15, 'core': [2, 4, 169, 147, 167, 203, 159, 3], 'ports': [0], 'word': [0, 0, 0, 0, 1], 'rows': [24828, 16633, 16466, 16484, 16456, 16404, 16416, 16508, 24831, 25087, 25599, 26623, 20479, 4350, 8192], 'root': 13, 'aux': 14}, {'H': 38617, 'gap': 2, 'order': 14, 'core': [2, 4, 73, 467, 487, 463, 387, 263, 7], 'ports': [0, 1], 'word': [1, 1, 1], 'rows': [12796, 12793, 8626, 8228, 8200, 8208, 8252, 8312, 8440, 8703, 9215, 10239, 4092, 4096], 'root': 12, 'aux': 13}, {'H': 38627, 'gap': 2, 'order': 15, 'core': [2, 4, 89, 83, 99, 79, 3], 'ports': [0], 'word': [1, 0, 0, 0, 1, 1], 'rows': [24700, 16505, 16418, 16420, 16396, 16400, 16444, 16511, 24831, 25087, 25599, 18431, 20479, 6398, 8192], 'root': 13, 'aux': 14}, {'H': 38629, 'gap': 2, 'order': 14, 'core': [198, 452, 488, 499, 487, 451, 384, 256, 1], 'ports': [0], 'word': [0, 1, 1], 'rows': [12600, 8249, 8211, 8196, 8200, 8220, 8255, 8319, 8446, 12799, 9215, 10239, 3582, 4096], 'root': 12, 'aux': 13}, {'H': 38633, 'gap': 2, 'order': 14, 'core': [86, 244, 248, 243, 224, 193, 128, 1], 'ports': [0], 'word': [1, 1, 1, 1], 'rows': [12456, 8201, 8195, 8196, 8207, 8222, 8255, 8318, 8447, 8703, 9215, 10239, 4094, 4096], 'root': 12, 'aux': 13}, {'H': 38635, 'gap': 2, 'order': 14, 'core': [2, 4, 9, 19, 39, 463, 927, 799, 543, 63], 'ports': [0, 2], 'word': [0, 1], 'rows': [13308, 9209, 13298, 9188, 9160, 8720, 8224, 8288, 8416, 8640, 13311, 10239, 3066, 4096], 'root': 12, 'aux': 13}]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();records=[]
    for witness in WITNESSES:
        rows,r,a=ready.ready_target(tuple(witness['core']),tuple(witness['ports']),tuple(witness['word']))
        ready.require(rows==tuple(witness['rows']) and r==witness['root'] and a==witness['aux'],'literal target reproduction')
        ready.require(ready.strong(rows),'actual strong target')
        obs=native.joint_readout(rows,r,a)
        direct=ready.counts(rows,(r,))
        ready.require(direct['H']==witness['H'] and direct['F'][r]==obs['F'],'independent full target H/F')
        ready.require(obs['sigma']==2 and obs['tau']==0,'actual positive joint C2 guard')
        records.append(dict(witness,observers=obs))
    total=ready.GATES+native.GATES
    result={'scope':'eight actual positive strong sole-exit gap2 targets; no all-order absence claim','records':records,'gates':total}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();args.output.write_bytes(raw)
    print(json.dumps({'gates':total,'targets':len(records),'receipt_bytes':len(raw),'receipt_sha256':sha256(raw).hexdigest()},sort_keys=True))


if __name__=='__main__':main()
