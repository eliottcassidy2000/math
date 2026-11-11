#!/usr/bin/env python3
"""Direct-address witness compiler for the certified finite interval."""
import argparse,base64,hashlib,json,struct,zlib
from pathlib import Path
from spine import graph,value
P=Path(__file__).parent

def load():
 raw=zlib.decompress(base64.b85decode(''.join((P/'atlas.b85').read_text().split())))
 metadata=json.loads((P/'atlas.json').read_text())
 if hashlib.sha256(raw).hexdigest()!=metadata['raw_sha256']:raise RuntimeError('certificate checksum')
 return raw,metadata

def compile_target(N,raw=None):
 if N<1 or N%2==0 or N>1000000:raise ValueError('certified target range: positive odd N <= 1000000')
 if raw is None:raw,_=load()
 code=struct.unpack_from('<I',raw,4*(N//2))[0]
 if code==0xffffffff:return {'target':N,'attained_in_bounded_language':False}
 m=code>>28;bits=code&((1<<28)-1);word=[(bits>>(2*(m-1-i)))&3 for i in range(m)]
 if value(word)!=N:raise RuntimeError('certificate transfer mismatch')
 return {'target':N,'vertices':m+2,'word':word,'rows':graph(word),
 'minimum_scope':'minimum order among all words with fixed two-vertex core; not all tournaments',
 'strongness':'not required or asserted'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('target',type=int,nargs='?');p.add_argument('--decode',metavar='PATH');a=p.parse_args()
 if a.decode:
  raw,meta=load();Path(a.decode).write_bytes(raw);print(json.dumps(meta,indent=2))
 else:print(json.dumps(compile_target(a.target),indent=2))
