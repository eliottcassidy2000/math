#!/usr/bin/env python3
"""Literal graph recovery and independent two-lane enumeration; no inverse imports."""
import argparse,json,itertools,hashlib
from pathlib import Path
from math import isqrt
BASE=Path(__file__).resolve().parent

def require(ok,msg):
 if not ok: raise ValueError(msg)

def adjacency(word):
 rows=[];cuts=[];blocks=[]
 for tag,mask in word:
  local=[0] if tag=='S' else [8,1,3,6] if tag=='D' else None
  require(local is not None,'unknown block')
  k=len(local);require(0<=mask<1<<k,'mask out of range');old=len(rows)
  rows=[r|(((1<<k)-1)<<old) for r in rows]+[r<<old for r in local]
  cuts.extend((mask>>i)&1 for i in range(k));blocks.append(list(range(old,old+k)))
 n=len(rows)
 rows=[r|((1-cuts[i])<<n) for i,r in enumerate(rows)]+[sum(b<<i for i,b in enumerate(cuts))]
 return rows,blocks

def validate(rows):
 n=len(rows)
 require(all(type(r)==int and 0<=r<1<<n and not r>>i&1 for i,r in enumerate(rows)),'invalid rows')
 require(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i)),'not tournament')

def path(rows,seq):return all(rows[a]>>b&1 for a,b in zip(seq,seq[1:]))

def literal_count(rows):
 return sum(path(rows,p) for p in itertools.permutations(range(len(rows))))

def two_lane_count(rows,blocks):
 """Place each non-apex vertex before/after the apex, retaining both paths."""
 validate(rows);a=len(rows)-1
 require(sorted(sum(blocks,[]))==list(range(a)),'partition')
 require(all(rows[v]>>w&1 for i,B in enumerate(blocks) for C in blocks[i+1:] for v in B for w in C),'not ordinal blocks')
 require(all(len(B)<=4 for B in blocks),'local factorial cap')
 # B-state0 empty,1 last vertex -> apex,2 not. A-state0 empty,1 apex -> first,2 not.
 counts={(0,0):1}
 for block in blocks:
  options={}
  for mask in range(1<<len(block)):
   left=[v for j,v in enumerate(block) if mask>>j&1];right=[v for j,v in enumerate(block) if not mask>>j&1]
   for P in itertools.permutations(left):
    if not path(rows,P):continue
    for Q in itertools.permutations(right):
     if not path(rows,Q):continue
     b=0 if not P else 1 if rows[P[-1]]>>a&1 else 2
     c=0 if not Q else 1 if rows[a]>>Q[0]&1 else 2
     options[b,c]=options.get((b,c),0)+1
  new={}
  for (b,c),v in counts.items():
   for (d,e),w in options.items():
    key=(d or b,c or e);new[key]=new.get(key,0)+v*w
  counts=new
 return sum(v for (b,c),v in counts.items() if b!=2 and c!=2)

def singleton_formula(rows):
 m=len(rows)-1;a=m
 require(all(rows[i]>>j&1 for i in range(m) for j in range(i+1,m)),'nontransitive parent')
 if m==0:return 1
 b=[(rows[a]>>i)&1 for i in range(m)]
 return b[0]+1-b[-1]+sum((1-b[i])*b[i+1] for i in range(m-1))+sum((1-b[i])*b[j]*(1<<(i-j-1)) for i in range(m) for j in range(i))

def strong(rows):
 full=(1<<len(rows))-1
 for start in range(len(rows)):
  seen=1<<start
  while True:
   nxt=seen
   for i,row in enumerate(rows):
    if seen>>i&1:nxt|=row
   if nxt==seen:break
   seen=nxt
  if seen!=full:return False
 return True

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def controls():
 tests=0
 alphabet=[('S',0),('S',1)]+[('D',i) for i in range(16)]
 words=[[]]+[[x] for x in alphabet]+[[x,y] for x in alphabet for y in alphabet if x[0]=='S' or y[0]=='S']
 for word in words:
  rows,blocks=adjacency(word);require(two_lane_count(rows,blocks)==literal_count(rows),'literal control mismatch');tests+=1
 for length in range(9):
  for bits in itertools.product(range(2),repeat=length):
   rows,blocks=adjacency([['S',b] for b in bits]);require(two_lane_count(rows,blocks)==singleton_formula(rows),'formula mismatch');tests+=1
 return tests

def main():
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();records=[]
 for line in (BASE/'query.normal.jsonl').read_text().splitlines():
  q=json.loads(line);rows,blocks=adjacency(q['word']);n=len(rows);h=two_lane_count(rows,blocks)
  require(h==q['target'],'count mismatch');require(strong(rows),'not strong');require(prime(h),'not prime')
  if q['parent']=='singleton': require(singleton_formula(rows)==h,'closed formula mismatch')
  cert={'target':h,'order':n,'rows':rows,'blocks':blocks,'word':q['word'],'strong':True,'prime':True,'scope':'recovered existing grammar witness; no novelty/minimum-order/global-completeness claim'}
  dest=BASE/f'witness_{h}.json'
  if args.write:
   dest.write_text(json.dumps(cert,indent=2)+'\n');(BASE/f'adjacency_{h}.txt').write_text(f'{n} {h}\n'+'\n'.join(''.join(str(row>>j&1) for j in range(n)) for row in rows)+'\n')
  else:require(json.loads(dest.read_text())==cert,'frozen certificate changed')
  records.append({'target':h,'order':n,'strong':True,'prime':True,'count_two_lane':h,'count_singleton_formula':h if q['parent']=='singleton' else None,'certificate_sha256':hashlib.sha256(dest.read_bytes()).hexdigest()})
 print(json.dumps({'status':'passed','controls':controls(),'records':records},sort_keys=True,indent=2))
if __name__=='__main__':main()
