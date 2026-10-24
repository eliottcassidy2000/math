#!/usr/bin/env python3
"""Replay unchanged frozen producer code with SHA-pinned filesystem IO only.
The adapter replaces the producer's git-show reader; it executes no Git command.
"""
from pathlib import Path
import argparse,hashlib,importlib,json,sys
HERE=Path(__file__).resolve().parent
P=HERE/'inputs/native';O=HERE/'inputs/original'
sys.dont_write_bytecode=True;sys.path.insert(0,str(P))
args=argparse.ArgumentParser();args.add_argument('which',choices=['native570','audit570']);a=args.parse_args()
import common
expected=json.loads((P/'native570.normal.json').read_text())['inputs']
def read(path):
 data=(O/path).read_bytes();h=hashlib.sha256(data).hexdigest()
 if h!=expected[path]:raise RuntimeError('Original input hash mismatch: '+path)
 common.INPUTS[path]=h;return data
common.read=read
bank='compatible_charged_helpers_20261023/supply/selected_supply.json'
def forbidden_git_adapter(argv,**kwargs):
 if argv!=['git','show','1d94c1d1d9:05-knowledge/results/'+bank]:raise RuntimeError('Unexpected subprocess; no subprocess/Git execution allowed')
 data=(O/bank).read_bytes()
 if hashlib.sha256(data).hexdigest()!=json.loads((P/'audit570.normal.json').read_text())['incoming_supply_sha256']:raise RuntimeError('Helper input hash mismatch')
 return data
common.subprocess.check_output=forbidden_git_adapter
r=importlib.import_module(a.which).main()
text=json.dumps(r,indent=2,sort_keys=True)+'\n'
if text!=(P/(a.which+'.normal.json')).read_text():raise RuntimeError('Replay differs from saved producer receipt')
print(text,end='')
