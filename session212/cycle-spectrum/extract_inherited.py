#!/usr/bin/env python3
"""Re-extract cycle support from the pinned THM-4889 (H,C) census.
Does not regenerate or independently certify the order-ten tournament universe.
"""
from collections import Counter
from pathlib import Path
import argparse,hashlib,json,subprocess
ROOT=Path(__file__).resolve().parent
REF='f9de668e6caaec77b829c13e3d87cf1ad8db61a3'
PREFIX='05-knowledge/results/path_spectrum_two_angles_metagraph_20261106/'

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--repo',required=True);args=parser.parse_args()
 manifest=subprocess.check_output(['git','-C',args.repo,'show',f'{REF}:{PREFIX}SHA256SUMS']).decode();expected=dict((line.split()[1],line.split()[0]) for line in manifest.splitlines() if line.strip())
 summary={};maximum=None
 for n in (8,9,10):
  rel=f'data/census_{n}.txt';blob=subprocess.check_output(['git','-C',args.repo,'show',f'{REF}:{PREFIX}{rel}']);digest=hashlib.sha256(blob).hexdigest()
  if expected[rel]!=digest:raise RuntimeError('inherited source hash mismatch')
  hist=Counter()
  for line in blob.decode().splitlines():
   f=line.split()
   if f and f[0]=='C':hist[int(f[2])]+=int(f[3])
   if n==10 and f[:2]==['W','15745']:
    bits=f[2];rows=[0]*10;k=0
    for i in range(10):
     for j in range(i+1,10):
      if bits[k]=='1':rows[i]|=1<<j
      else:rows[j]|=1<<i
      k+=1
    maximum={'rows':rows,'H':15745,'C':971,'source':'THM-4889 census_10.txt W 15745','upper_triangle_bits':bits}
  summary[str(n)]={'order':n,'classes':sum(hist.values()),'strong_classes':sum(v for k,v in hist.items() if k),'spectrum':sorted(hist),'gaps':[k for k in range(max(hist)+1) if k not in hist],'multiplicities':dict(sorted(hist.items())),'source_sha256':digest}
 (ROOT/'inherited_summary.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
 if maximum is None:raise RuntimeError('maximum witness absent')
 (ROOT/'maximum10_witness.json').write_text(json.dumps(maximum,sort_keys=True,indent=2)+'\n')
 print(json.dumps({str(n):{'source_sha256':summary[str(n)]['source_sha256'],'classes':summary[str(n)]['classes'],'strong_classes':summary[str(n)]['strong_classes'],'spectrum_size':len(summary[str(n)]['spectrum'])} for n in (8,9,10)},sort_keys=True,indent=2))
if __name__=='__main__':main()
