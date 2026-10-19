#!/usr/bin/env python3
"""Construct and independently expand any requested certified N>=128453. Stdlib only."""
import argparse,gzip,hashlib,json,struct
from functools import lru_cache
from math import isqrt
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;BASE=ROOT/'square_sum_B_base_interval_20261019';FINAL=ROOT/'seeded_final_split_table_20261019'
checks=0;used=[]
def check(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def check_pair(paths,family):
 for path in paths:
  n=len(path);check(sorted(path)==list(range(1,n+1)),'pair interval coverage')
  check(path[:2]==([1,3] if family=='B' else [69,127]),'pair prefix')
  terminal=(10 if n%2==0 else 8) if family=='Aplus' else (8 if n%2==0 else 10)
  check(path[-1]==terminal,'pair endpoint')
  for x,y in zip(path,path[1:]):check(isqrt(x+y)**2==x+y,'pair square edge')
  if family!='B':
   p2=path.index(2);p7=path.index(7)
   check(abs(p2-p7)==1 and p2%2==0 and p7%2==1,'protected2–7 position')
 check(len(paths[1])==len(paths[0])+1,'pair consecutive lengths')
 pos={v:(j+1)%2 for j,v in enumerate(paths[0])}
 check(all(pos[v]==(j+1)%2 for j,v in enumerate(paths[1]) if v in pos),'pair parity compatibility')
def base_pair(n):
 manifest=json.loads((BASE/'B_base_chunks.json').read_text());chunk=next(c for c in manifest if c['first']<=n<=c['last'])
 rawzip=(BASE/chunk['file']).read_bytes();check(hashlib.sha256(rawzip).hexdigest()==chunk['sha256'],'base chunk hash')
 raw=gzip.decompress(rawzip);lo,hi=struct.unpack('<II',raw[8:16]);check(raw[:8]==b'SSBAS1\0\0' and lo<=n<=hi,'base header')
 offset=16+(n-lo)*(n+lo-2);paths=[]
 for length in (n,n+1):
  roots=raw[offset:offset+length-1];offset+=length-1;path=[1]
  for r in roots:path.append(r*r-path[-1])
  paths.append(path)
 check_pair(paths,'B');return paths
def transition(paths,family,n):
 m=len(paths[0]);d=n-49*m
 witness=json.loads((HERE/'certified-pair-transitions'/f'{family}_p{m%2}_d{d}.json').read_text());used.append([family,m%2,d,witness.get('split_prefix_blocks',[])])
 outpair=[]
 for word in witness['words']:
  out=[]
  for atom in word:
   if atom[0]=='single':out.append(atom[1])
   elif atom[0]=='head':out.append(49-atom[1])
   else:
    kind,c,i,o=atom;block=[49*v+(-1)**(j+1)*c for j,v in enumerate(paths[i])]
    if kind=='tail':block=block[1:]
    out.extend(block if o==1 else block[::-1])
  outpair.append(out)
 check(len(outpair[0])==n,'transition length');check_pair(outpair,family);return outpair
@lru_cache(None)
def build_B(n):
 check(n>=53,'B range')
 if n<=2620:return base_pair(n)
 m=(n-24)//49;return transition(build_B(m),'B',n)
def cycle(N):
 check(N>=128453,'certified final range');n=(N-24)//49;d=N-49*n;j=int(d>=57)
 family='Aplus' if (n+j)%2==0 else 'Aminus';m=(n-24)//49
 A=transition(build_B(m),family,n)
 orders=tuple(1 if path.index(2)<path.index(7) else -1 for path in A)
 table=json.loads((FINAL/'final_split_glue_certificates_20261019.json').read_text())
 f=next(x for x in table if (x['parity'],x['residue'],x['order0'],x['order1'])==(n%2,d,*orders))
 cuts=[min(path.index(2),path.index(7))+1 for path in A];out=[]
 for atom in f['word']:
  if atom[0]=='single':out.append(atom[1])
  else:
   _,c,i,seg,o=atom;block=[49*v+(-1)**(j+1)*c for j,v in enumerate(A[i])]
   part=block[:cuts[i]] if seg==0 else block[cuts[i]:];out.extend(part if o==1 else part[::-1])
 check(len(out)==N,'cycle length');check(sorted(out)==list(range(1,N+1)),'cycle interval coverage')
 for x,y in zip(out,out[1:]+out[:1]):check(isqrt(x+y)**2==x+y,'cycle square edge')
 idx=out.index(3362);check({out[(idx-1)%N],out[(idx+1)%N]}=={482,6242},'seed adjacency')
 idx=out.index(482);step=1 if out[(idx+1)%N]==3362 else -1;canonical=[out[(idx+step*k)%N] for k in range(N)]
 check(canonical[:3]==[482,3362,6242],'oriented seed prefix')
 raw=('\n'.join(map(str,canonical))+'\n').encode();name=f'seeded_cycle_{N}.txt.gz';(HERE/name).write_bytes(gzip.compress(raw,mtime=0))
 return {'N':N,'A_index':n,'A_family':family,'final_residue':d,'cut_orders':list(orders),'cycle_file':name,'plaintext_sha256':hashlib.sha256(raw).hexdigest(),'canonical_prefix':canonical[:12]}
p=argparse.ArgumentParser();p.add_argument('N',nargs='*',type=int,default=[128453,132063,132128]);args=p.parse_args()
# Additional control directly exercises the repaired B-to-B transition.
control=build_B(2694);check(len(control[0])==2694,'B48 repaired expansion')
reports=[cycle(N) for N in args.N]
print(json.dumps({'status':'passed','checks':checks,'B_repaired_control_lengths':list(map(len,control)),'used_pair_transitions':used,'cycles':reports},indent=2,sort_keys=True))
