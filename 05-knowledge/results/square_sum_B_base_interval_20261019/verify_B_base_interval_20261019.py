#!/usr/bin/env python3
"""Independent stdlib decoder and exhaustive check of the full B base interval."""
import gzip,hashlib,json,struct
from pathlib import Path
ROOT=Path(__file__).resolve().parent
checks=0

def check(ok,message):
 global checks
 checks+=1
 if not ok:raise RuntimeError(message)

manifest=json.loads((ROOT/'B_base_chunks.json').read_text())
expected=53;pairs=0;total_vertices=0;total_roots=0
for chunk in manifest:
 check(chunk['first']==expected,'nonconsecutive chunk range')
 compressed=(ROOT/chunk['file']).read_bytes()
 check(len(compressed)==chunk['bytes'],'compressed length')
 check(hashlib.sha256(compressed).hexdigest()==chunk['sha256'],'compressed digest')
 raw=gzip.decompress(compressed)
 check(len(raw)==chunk['decoded_bytes'],'decoded length')
 check(hashlib.sha256(raw).hexdigest()==chunk['decoded_sha256'],'decoded digest')
 check(raw[:8]==b'SSBAS1\0\0','format version')
 lo,hi=struct.unpack('<II',raw[8:16])
 check(lo==chunk['first'] and hi==chunk['last'] and lo<=hi<=2620,'header range')
 pos=16
 for n in range(lo,hi+1):
  common_parity=bytearray(n+1)
  for member,length in enumerate((n,n+1)):
   roots=raw[pos:pos+length-1];pos+=length-1
   check(len(roots)==length-1,'truncated path')
   seen=bytearray(length+1);seen[1]=1;v=1
   if member==0:common_parity[1]=1
   for position,r in enumerate(roots,2):
    check(2<=r<=72,'root byte range')
    v=r*r-v
    check(1<=v<=length and not seen[v],'vertex range or repeated vertex')
    seen[v]=1
    if position==2:check(v==3,'common initial edge 1-3')
    if member==0:common_parity[v]=position%2
    elif v<=n:check(common_parity[v]==position%2,'common positional parity')
   check(sum(seen)==length,'complete interval coverage')
   check(v==(8 if length%2==0 else 10),'required endpoint')
   total_vertices+=length;total_roots+=length-1
  pairs+=1
 check(pos==len(raw),'trailing bytes')
 expected=hi+1
check(expected==2621 and pairs==2568,'complete base interval')
print(json.dumps({'status':'passed','first':53,'last':2620,'pairs':pairs,'chunks':len(manifest),
 'checked_vertices':total_vertices,'checked_roots':total_roots,'checks':checks,
 'state':'Both paths start 1,3; even length ends 8; odd length ends 10; shared vertices have equal positional parity.'},sort_keys=True,indent=2))
