#!/usr/bin/env python3
"""Bounded-memory reuse of the proved source-specific marked-block inverse."""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'square_sum_spectrum_overnight_20261031b'))
import diamond_inverse as inv

def main():
 p=argparse.ArgumentParser();p.add_argument('targets',nargs='*',type=int);p.add_argument('--cache',type=int,default=500000);a=p.parse_args()
 if a.cache<1:raise ValueError('positive cache limit required')
 original=inv.inverse.__wrapped__;records=[]
 for F in a.targets or (594793,774073,949777):
  inv.inverse=lru_cache(maxsize=a.cache)(original)
  inv.unary_F.cache_clear()
  record=inv.target(F);record['cache']=inv.inverse.cache_info()._asdict();records.append(record)
 print(json.dumps({'status':'passed','cache_limit':a.cache,'results':records,
                   'scope':'Exact source-specific inverse; cache eviction repeats work but never removes candidate branches.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
