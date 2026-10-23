"""Freeze native producer outputs with deterministic gzip and source hashes."""
from pathlib import Path
import gzip,hashlib,json
P=Path(__file__).resolve().parent
names=['owner_map.u16','ports.bin','ordinary_edges.tsv','merged_edges.tsv']
out={}
for name in names:
 data=(P/name).read_bytes();path=P/(name+'.gz')
 with path.open('wb') as raw:
  with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=9) as f:f.write(data)
 out[name]={'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'gzip_bytes':path.stat().st_size,'gzip_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
(P/'artifact_hashes.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
