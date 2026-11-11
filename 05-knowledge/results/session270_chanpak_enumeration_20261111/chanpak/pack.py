#!/usr/bin/env python3
"""Pack census output into one uint32 word/length record per positive odd."""
import json,struct,zlib,base64,sys,hashlib
from pathlib import Path
P=Path(__file__).parent
x=json.loads(Path(sys.argv[1]).read_text());codes=[0xffffffff]*((x['limit']+1)//2)
for h,n,w in x['witnesses']:codes[h//2]=((n-2)<<28)|w
raw=struct.pack('<'+'I'*len(codes),*codes);packed=base64.b85encode(zlib.compress(raw,9)).decode()
(P/'atlas.b85').write_text('\n'.join(packed[i:i+100] for i in range(0,len(packed),100))+'\n')
Path('/tmp/s270-atlas.bin').write_bytes(raw)
summary={'scope':'marked two-apex G(2,0) words, all lengths 0..14; not tournament class counts',
 'target_limit':x['limit'],'max_spine':x['max_spine'],'max_vertices':x['max_spine']+2,
 'attained':len(x['witnesses']),'missing':[2*i+1 for i,c in enumerate(codes) if c==0xffffffff],
 'vertex_histogram':{str(n):sum(c!=0xffffffff and (c>>28)+2==n for c in codes) for n in range(2,17)},
 'certificate_format':'little-endian uint32 indexed by (N-1)/2, high4 bits spine length, low28 bits base4 word; FFFFFFFF absent',
 'raw_sha256':hashlib.sha256(raw).hexdigest(),'raw_bytes':len(raw),'base85_characters':len(packed),
 'visited_prefixes_with_repetition':x['visited_prefixes_with_repetition']}
(P/'atlas.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
