#!/usr/bin/env python3
"""Portable catalog/base/constructor replay. Independent audits have separate CLIs."""
import argparse,hashlib,itertools,json
from pathlib import Path
from verify_integrity import verify as verify_integrity
from core import Constructor,DEFAULT_BASE,DEFAULT_SMALL,DEFAULT_Q49,HERE,threshold
p=argparse.ArgumentParser();p.add_argument('--base-dir',type=Path,default=DEFAULT_BASE);p.add_argument('--small-base-file',type=Path,default=DEFAULT_SMALL);p.add_argument('--q49-dir',type=Path,default=DEFAULT_Q49);a=p.parse_args()
def need(x,msg):
 if not x:raise RuntimeError(msg)
verify_integrity()
meta=json.loads((HERE/'CERTIFICATES.json').read_text());need(len(meta)==200,'exact Q25 catalog count')
for f,h in meta.items():need(hashlib.sha256((HERE/f).read_bytes()).hexdigest()==h,'certificate hash:'+f)
for symbol in['A','B']:
 for directory,stem in[('paired','B'),('cycles','C')]:need({f.name for f in(HERE/symbol/directory).glob('*.json')}=={f'{stem}_p{p}_d{d}.json' for p in[0,1] for d in range(12,37)},'exact state coverage')
def make():return Constructor(a.base_dir,a.small_base_file,a.q49_dir)
c=make()
for n in range(42,2082):c.bases.base(n)
base_report={'pairs':2040,'checks':c.bases.checks,'chunks':len(c.bases.used_base_chunks)}
controls=[]
for word in map(''.join,itertools.product('ABC',repeat=2)):
 for extra in[0,1]:controls.append(make().cycle(threshold(word)+extra,word))
for word in list(map(''.join,itertools.permutations('ABC')))+['AAA','BBB','CCC']:controls.append(make().cycle(threshold(word),word))
rejected=0
for word,n in[('',1062),('D',1062),('B',1061),('AB',threshold('AB')-1)]:
 try:make().cycle(n,word)
 except(ValueError,RuntimeError):rejected+=1
need(rejected==4,'hostile constructor domains')
# The A/B-only path must not access any Q49 dependency.
no49=Constructor(a.base_dir,a.small_base_file,HERE/'NONEXISTENT_Q49');no49.cycle(1062,'B');need(not no49.used_q49,'Q25 independence')
print(json.dumps({'status':'passed','Q25_certificates':200,'base_report':base_report,'controls':controls,'hostile_domains_rejected':rejected,'Q25_requires_Q49':False},indent=2))
