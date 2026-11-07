"""Reconstruct the losslessly chunked certificate outside the published package."""
from pathlib import Path
import base64,gzip,hashlib,json,sys
D=Path(__file__).resolve().parent
def reconstruct(output):
 m=json.loads((D/'certificate/MANIFEST.json').read_text());parts=[]
 for p in m['parts']:
  data=(D/p['path']).read_bytes()
  if len(data)!=p['bytes'] or hashlib.sha256(data).hexdigest()!=p['sha256']:raise ValueError('Invalid part '+p['path'])
  parts.append(data)
 compressed=base64.b64decode(b''.join(parts),validate=True)
 if hashlib.sha256(compressed).hexdigest()!=m['compressed_sha256']:raise ValueError('Compressed hash mismatch')
 raw=gzip.decompress(compressed)
 if len(raw)!=m['original_bytes'] or hashlib.sha256(raw).hexdigest()!=m['original_sha256']:raise ValueError('Original hash mismatch')
 output=Path(output)
 if output.exists():raise FileExistsError(output)
 output.write_bytes(raw)
 return m['original_sha256']
if __name__=='__main__':
 if len(sys.argv)!=2:raise SystemExit('Usage: python3 reconstruct.py /external/path/RECIPES.bin')
 print(reconstruct(sys.argv[1]))
