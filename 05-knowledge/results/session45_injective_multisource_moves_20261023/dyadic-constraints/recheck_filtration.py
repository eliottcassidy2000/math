#!/usr/bin/env python3
"""Rerun the pinned original order-eight filtration, without absent old controls."""
import argparse,csv,gzip,importlib.util,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);a=p.parse_args()
base=a.repo/'05-knowledge/results';path=base/'metagraph_order8_electrical_20261023/check.py'
spec=importlib.util.spec_from_file_location('original_electrical_checker',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
L=[];B=[]
for row in csv.DictReader(gzip.decompress((base/'metagraph_order8_20261023/merged_edges.tsv.gz').read_bytes()).decode().splitlines(),delimiter='\t'):
 i,j,b,k=(int(row[c]) for c in ['source','target','blue','black'])
 if k%2:raise RuntimeError('black weight must be paired')
 if b:L.append((i,j,b))
 if k:B.append((i,j,k//2))
result=module.calculate(3528,L,B)
expected=json.loads((base/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
if result!=expected:raise RuntimeError('original filtration rerun differs')
print(json.dumps({'status':'VERIFIED exact original order8 filtration rerun','checks':module.CHECKS,'matches_entire_original_order8_receipt':True,'scope':'Original implementation rerun; old order≤7 control input absent from this checkout and not rerun.','two_primary_factor_exponents':[1]*3358+[2]*4+[8]},sort_keys=True,indent=2))
