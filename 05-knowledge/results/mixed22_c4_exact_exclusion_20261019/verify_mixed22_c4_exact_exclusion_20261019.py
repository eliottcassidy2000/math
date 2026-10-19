#!/usr/bin/env python3
"""Standalone standard-library replay of the complete mixed-(2,2), c=4 packet.

The placeholder hashes are filled only after every construction shard finishes.
Without --part, independently replay the whole packet. The four --part runs
may instead be run in parallel and consolidated with --assemble.
"""
from pathlib import Path
import argparse,gzip,hashlib,heapq,importlib.util,json,time
ROOT=Path(__file__).resolve().parent
RESIDUAL_SHAS={'residual_certificates_part0_20261019.jsonl.gz': '838f5d80504adad2ceabab9764947c0be8b27cd7881c829851ee9105265f96fa', 'residual_certificates_part1_20261019.jsonl.gz': 'f1b63ea57b8f859a444f20b7639bd7a32beb41c6ad11b4c7efdaa643eaba1a83', 'residual_certificates_part2_20261019.jsonl.gz': '5f988a6bd9387c6bfefcb6b61538a4fbd3ff2d99fd68880ace07a63a8d4c325b', 'residual_certificates_part3_20261019.jsonl.gz': '28bf42336c857a034963c52dd7da81512bc54fafd63daba27b43020193b1a5c5'}
INVENTORY_SHA='63697acce82106e50e52b56834821dd64799302fe4cb7966b4c83e334a7b71c9'
ANALYTIC_VERIFIER_SHA='2e3e1b3640449af276bf7fc1ba5a14369195b4edf77dcbe148018fa8d8fb9fc1'
REPLAY_HELPER_SHA='2dafdbfdaeaacc2ca28a65ffd61d915b9077cb127f68eb15e4e7be4c8e287226'
EXPECTED={'residual_cores': 215267, 'generator_orders': 2902976, 'directed_screen_orders': 113423189, 'distinct_certificates': 184897, 'one_row_certificates': 72131, 'two_row_certificates': 112616, 'larger_certificates': 150, 'maximum_certificate_support': 12, 'maximum_integer_weight': 164426573285821101279024000}

def check(ok,message):
 if not ok:raise RuntimeError(message)

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''):h.update(block)
 return h.hexdigest()

def load_module(name,filename):
 s=importlib.util.spec_from_file_location(name,ROOT/filename)
 m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def fresh_stats():
 return {'residual_cores':0,'generator_orders':0,'directed_screen_orders':0,
         'distinct_certificates':0,'one_row_certificates':0,'two_row_certificates':0,
         'larger_certificates':0,'maximum_certificate_support':0,'maximum_integer_weight':0}

def accumulate(s,r):
 s['residual_cores']+=1;s['generator_orders']+=r['residual_orders'];s['directed_screen_orders']+=r['analytic_directed_orders']
 for c in r['certificates']:
  terms=c['terms'];k=len(terms);s['distinct_certificates']+=1
  s['one_row_certificates' if k==1 else 'two_row_certificates' if k==2 else 'larger_certificates']+=1
  s['maximum_certificate_support']=max(s['maximum_certificate_support'],k)
  s['maximum_integer_weight']=max([s['maximum_integer_weight']]+[t[2] for t in terms])

def write(path,value):
 temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');temp.replace(path)

def main():
 ap=argparse.ArgumentParser();group=ap.add_mutually_exclusive_group();group.add_argument('--part',type=int,choices=range(4));group.add_argument('--assemble',action='store_true');args=ap.parse_args()
 start=time.time();inventory=ROOT/'residual_inventory_20261019.json.gz'
 for name,expected in RESIDUAL_SHAS.items():check(sha(ROOT/name)==expected,'residual stream hash')
 check(sha(inventory)==INVENTORY_SHA,'inventory hash')
 check(sha(ROOT/'verify_c4_analytic_and_probes_20261019.py')==ANALYTIC_VERIFIER_SHA,'analytic verifier hash')
 check(sha(ROOT/'exact_cone_replay_20261019.py')==REPLAY_HELPER_SHA,'replay helper hash')
 V=load_module('analytic_control','verify_c4_analytic_and_probes_20261019.py');R=V.R
 inv=json.loads(gzip.decompress(inventory.read_bytes()));codes=inv['codes'];check(len(codes)==215267 and inv['classification_sha256']==V.A_SHA,'inventory scope')
 self_sha=sha(Path(__file__));stats=fresh_stats();part_stats=[fresh_stats() for _ in range(4)]
 if args.part in (None,0) and not args.assemble:V.main()
 seen=0
 def records(k):
  with gzip.open(ROOT/f'residual_certificates_part{k}_20261019.jsonl.gz','rt') as f:
   previous=-1
   for line in f:
    r=json.loads(line);i=r['residual_index'];check(i%4==k and i>previous,'stream partition/order');previous=i;yield r
 for r in heapq.merge(*(records(k) for k in range(4)),key=lambda r:r['residual_index']):
   i=r['residual_index'];check(i==seen,'residual indices not exactly range(215267)');seen+=1
   pred=V.decode(codes[i][0]);check(tuple(r['pred'])==pred,'residual core correspondence')
   check(not r.get('short_span_orders',0),'unexpected short-span dependency')
   accumulate(stats,r);accumulate(part_stats[i%4],r)
   if not args.assemble and (args.part is None or i%4==args.part):
    R.replay(pred,r)
    if (part_stats[i%4]['residual_cores'] if args.part is not None else seen)%1000==0:
     print(json.dumps({'part':args.part,'last_residual_index':i,'seconds':time.time()-start}),flush=True)
 check(seen==215267 and stats==EXPECTED,'complete residual coverage/counts')
 base={'residual_streams_sha256':RESIDUAL_SHAS,'inventory_sha256':INVENTORY_SHA,'classification_sha256':V.A_SHA,'analytic_verifier_sha256':ANALYTIC_VERIFIER_SHA,'replay_helper_sha256':REPLAY_HELPER_SHA,'verifier_sha256':self_sha}
 if args.part is not None:
  out=dict(base,status='INDEPENDENT_EXACT_REPLAY_PART_PASSED',part=args.part,statistics=part_stats[args.part],seconds=time.time()-start)
  write(ROOT/f'independent_replay_part{args.part}_20261019.json',out);print(json.dumps(out),flush=True);return
 part_seconds=[]
 if args.assemble:
  for k in range(4):
   rec=json.loads((ROOT/f'independent_replay_part{k}_20261019.json').read_text())
   check(rec['status']=='INDEPENDENT_EXACT_REPLAY_PART_PASSED' and rec['part']==k,'part status')
   check(all(rec[key]==value for key,value in base.items()) and rec['statistics']==part_stats[k],'part identity/statistics');part_seconds.append(rec['seconds'])
 census=json.loads((ROOT/'independent_census_replay_receipt_20261019.json').read_text())
 check(census['status']=='COMPLETE_ANALYTIC_CENSUS_AND_SPARSE_PROBES_INDEPENDENTLY_REPLAYED' and census['classification_sha256']==V.A_SHA and census['inventory_sha256']==INVENTORY_SHA and census['verifier_sha256']==ANALYTIC_VERIFIER_SHA,'analytic replay identity')
 out=dict(base,status='ALL_1165395_DEGREE_VALID_CORES_INDEPENDENTLY_EXACT_REPLAYED',scope='Seven-thin mixed-(2,2), c=4 maximum-span root morphology only.',natural_nonroot_posets=2800472,degree_invalid=1635077,degree_valid_cores=1165395,analytic_balanced_cores=733492,analytic_cycle_cores=216636,statistics=stats,short_span_orders=0,parallel_part_seconds=part_seconds,replay_mode='four exact parts plus complete set-equality assembly' if args.assemble else 'single complete exact replay',seconds=time.time()-start)
 write(ROOT/'independent_full_replay_receipt_20261019.json',out);print(json.dumps(out),flush=True)
if __name__=='__main__':main()
