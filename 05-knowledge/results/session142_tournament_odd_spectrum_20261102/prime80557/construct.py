#!/usr/bin/env python3
"""Recover an actual witness from the inherited x<=56 seven-parent bank.
A failed query is bank-relative, never a global nonexistence certificate.
"""
import argparse,json,sys
from pathlib import Path
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'sources/square_sum_spectrum_overnight_20261031d'))
from bank_query import query
from freeze_and_verify import adjacency,two_lane_count,validate,strong
p=argparse.ArgumentParser(description=__doc__);p.add_argument('target',type=int);p.add_argument('--cap',type=int,default=56);p.add_argument('--adjacency',type=Path,required=True);p.add_argument('--recipe',type=Path,required=True);args=p.parse_args()
if args.target<1 or args.target%2==0 or args.cap<0:p.error('positive odd target and nonnegative cap required')
r=query(args.target,args.cap)
if r['x'] is None:
 print(json.dumps({'status':'bank_relative_miss','query':r,'global_membership':'OPEN or needs a different certificate'},sort_keys=True));sys.exit(2)
rows,blocks=adjacency(r['word']);validate(rows);actual=two_lane_count(rows,blocks)
if actual!=args.target:raise ValueError('literal adjacency count mismatch')
n=len(rows);args.adjacency.write_text(f'{n} {actual}\n'+'\n'.join(''.join(str(row>>j&1) for j in range(n)) for row in rows)+'\n')
r.update({'order':n,'rows':rows,'blocks':blocks,'strong':strong(rows),'verified_H':actual});args.recipe.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
print(json.dumps({'status':'verified_existing_grammar_witness','order':n,'H':actual,'strong':r['strong']},sort_keys=True))
