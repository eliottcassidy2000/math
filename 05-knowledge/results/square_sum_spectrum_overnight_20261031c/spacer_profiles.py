#!/usr/bin/env python3
"""Exact fixed-parent profiles for D + transitive k-spacer + D."""
import argparse,json
from pathlib import Path
from one_block_formula import block_profiles,formula
from audit_feedback import subset_dp,require
HERE=Path(__file__).resolve().parent

def compose(b,c):
 H=b['H']*c['H'];E=b['E']*c['E']
 return {'H':H,'E':E,'p':c['H']*b['p']+c['p']*b['E'],
         'q':c['E']*b['q']+c['q']*b['H'],
         's':b['s']*c['H'],'t':b['H']*c['t'],
         'F':c['H']*b['F']+c['p']*b['q']+(c['H']-c['s'])*(b['H']-b['t'])+(c['F']-c['H'])*b['H']}

def rows(k,mask):
 D=[8,1,3,6];n=8+k;out=[]
 blocks=[D]+[[0] for _ in range(k)]+[D]
 start=0
 for block in blocks:
  after=((1<<n)-1)^((1<<(start+len(block)))-1)
  out.extend((r<<start)|after for r in block);start+=len(block)
 for v in range(n):
  if not(mask>>v&1):out[v]|=1<<n
 return out+[mask]

def main():
 p=argparse.ArgumentParser();p.add_argument('k',type=int);a=p.parse_args();k=a.k
 require(0<=k<=5,'small declared spacer range')
 ds=block_profiles([8,1,3,6]);ss=block_profiles([0]);records=[];checks=0
 for mask in range(1<<(8+k)):
  b=ds[mask&15]
  for j in range(k):b=compose(b,ss[(mask>>(4+j))&1])
  b=compose(b,ds[mask>>(4+k)]);b['mask']=mask;records.append(b)
  if k<=1 or mask%max(1,(1<<(8+k))//64)==0:
   require(subset_dp(rows(k,mask))==b['F'],'literal adjacency profile');checks+=1
 require(all(b['H']==25 and b['E']==1024*2**k for b in records),'ordinal block base counts')
 (HERE/f'DT{k}D_profiles.tsv').write_text('\n'.join(' '.join(str(b[x]) for x in ('H','E','p','q','s','t','F','mask')) for b in records)+'\n')
 print(json.dumps({'status':'passed','spacer':k,'profiles':len(records),'literal_checks':checks,
                   'scope':'Profile composition is exact; literal controls are exhaustive only for k<=1.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
