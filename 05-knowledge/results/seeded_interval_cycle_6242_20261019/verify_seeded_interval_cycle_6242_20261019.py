#!/usr/bin/env python3
"""Reconstruct the exact optimal-endpoint seeded interval cycle. Python stdlib only."""
import argparse,gzip,hashlib,json
from math import isqrt
from pathlib import Path
ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--write-cycle',action='store_true');args=parser.parse_args()
checks=0

def check(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)

def edge(x,y):return(min(x,y),max(x,y))
def square(n):return n>=0 and isqrt(n)**2==n
def edges(path):return{edge(a,b)for a,b in zip(path,path[1:]+path[:1])}
def bad_edges(path):return{e for e in edges(path)if not square(sum(e))}
def canonical(path):
 n=len(path);i=path.index(482);step=1 if path[(i+1)%n]==3362 else-1
 out=[path[(i+step*k)%n]for k in range(n)]
 check(out[:3]==[482,3362,6242],'seed orientation')
 return out

base=json.loads((ROOT/'single_cut_input_127.json').read_text());path=base['path']
check(base['status']=='witness' and len(path)==127,'wrong input size')
check(sorted(path)==list(range(1,128)),'input interval coverage')
check(path[:2]==[69,127] and path[-1]==10,'input prefix and endpoint')
for a,b in zip(path,path[1:]):check(square(a+b),'input nonsquare edge')
i7=path.index(7);cut=path.index(2)
check(cut==i7+1 and (cut+1)%2==1,'input protected cut state')

f=json.loads((ROOT/'final_split_p1_d24_cutodd_reverse.json').read_text())
check(f['status']=='witness' and f['parity']==1 and f['residue']==24,'wrong split certificate')
check(f['cutparity']==1 and f['order0']==-1,'wrong cut convention')
seen_segments=set();seen_singles=set();lifted=[]
for atom in f['word']:
 if atom[0]=='single':
  _,v=atom;check(1<=v<=24 and v not in seen_singles,'singleton coverage')
  seen_singles.add(v);lifted.append(v)
 else:
  _,c,i,seg,orientation=atom
  check(-24<=c<=24 and i==0 and seg in(0,1)and orientation in(-1,1),'segment atom')
  check((c,seg)not in seen_segments,'repeated segment');seen_segments.add((c,seg))
  block=[49*v+(-1)**(j+1)*c for j,v in enumerate(path)]
  part=block[:cut]if seg==0 else block[cut:]
  lifted.extend(part if orientation==1 else part[::-1])
check(seen_segments=={(c,s)for c in range(-24,25)for s in(0,1)},'missing segment')
check(seen_singles==set(range(1,25)),'missing singleton')
check(len(lifted)==6247 and sorted(lifted)==list(range(1,6248)),'lifted interval coverage')
check(lifted[0]==1 and lifted[-1]==8,'lifted closure endpoints')
for a,b in zip(lifted,lifted[1:]+lifted[:1]):check(square(a+b),'lifted nonsquare edge')
lifted=canonical(lifted)
lift_raw=('\n'.join(map(str,lifted))+'\n').encode()
check(hashlib.sha256(lift_raw).hexdigest()=='72afa24b30cdd4d774e1f048f079b0f4808a19ba8bbffc241807e40976ee8078','lifted hash')

trace=json.loads((ROOT/'seeded_cycle_6242_repair_trace.json').read_text())
check(trace['source_cycle_length']==6247 and trace['deleted_vertices']==list(range(6243,6248)),'deletion list')
cycle=[v for v in lifted if v<=6242]
protected={edge(482,3362),edge(3362,6242)}
check(len(bad_edges(cycle))==5,'initial broken edge count')
for move in trace['moves']:
 l=move['left_cut_index'];r=move['right_cut_index'];n=len(cycle)
 check(0<=l<r<n and r-l>1 and not(l==0 and r==n-1),'invalid two-edge reversal indices')
 removed={edge(cycle[l],cycle[(l+1)%n]),edge(cycle[r],cycle[(r+1)%n])}
 check(removed=={tuple(e)for e in move['removed_edges']},'removed-edge trace mismatch')
 check(not(removed&protected),'protected seed edge cut')
 before=len(bad_edges(cycle));check(before==move['bad_before'],'bad-before mismatch')
 cycle[l+1:r+1]=reversed(cycle[l+1:r+1])
 added={edge(cycle[l],cycle[(l+1)%n]),edge(cycle[r],cycle[(r+1)%n])}
 check(added=={tuple(e)for e in move['added_edges']},'added-edge trace mismatch')
 after=len(bad_edges(cycle));check(after==move['bad_after'] and after<=before,'bad-after mismatch')
check(len(trace['moves'])==24,'repair trace length')
check(sorted(cycle)==list(range(1,6243)),'final interval coverage')
for a,b in zip(cycle,cycle[1:]+cycle[:1]):check(square(a+b),'final nonsquare edge')
check(protected<=edges(cycle),'seed edges missing')
cycle=canonical(cycle)
raw=('\n'.join(map(str,cycle))+'\n').encode();digest=hashlib.sha256(raw).hexdigest()
check(digest=='a892ae08c5304c35fae319dfcd9429faa5ea37b2f21e89d3a715e17912ceca6d','canonical final hash')
if args.write_cycle:(ROOT/'seeded_interval_cycle_6242.txt.gz').write_bytes(gzip.compress(raw,mtime=0))
print(json.dumps({'status':'passed','checks':checks,'input_length':127,'cut_positions_1_based':[i7+1,cut+1],
 'lifted_cycle_length':6247,'deleted_vertices':list(range(6243,6248)),'repair_moves':len(trace['moves']),
 'cycle_length':6242,'vertex_set':[1,6242],'canonical_prefix':cycle[:20],
 'canonical_plaintext_sha256':digest,'canonical_plaintext_bytes':len(raw),
 'scope':'one exact interval cycle with prescribed consecutive seed; interval endpoint is minimal because the seed contains6242'},sort_keys=True,indent=2))
