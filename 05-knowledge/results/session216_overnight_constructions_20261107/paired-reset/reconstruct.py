from pathlib import Path
import json,hashlib,base64,argparse
p=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true',help='verify all parts and reconstructed payload without writing');args=parser.parse_args()
m=json.loads((p/'RECONSTRUCTION.json').read_text())
for rec in m['files']:
 target=p/rec['path'];data=bytearray()
 if target.resolve().is_relative_to(p) is False:raise RuntimeError('unsafe path')
 for part in rec['parts']:
  src=p/part['file'];raw=src.read_bytes()
  if hashlib.sha256(raw).hexdigest()!=part['sha256']:raise RuntimeError('part hash '+part['file'])
  data.extend(base64.b64decode(raw.strip(),validate=True))
 if len(data)!=rec['bytes'] or hashlib.sha256(data).hexdigest()!=rec['sha256']:raise RuntimeError('payload hash '+rec['path'])
 if args.check:
  if target.exists() and target.read_bytes()!=data:raise RuntimeError('existing payload differs')
 else:
  target.parent.mkdir(parents=True,exist_ok=True)
  if target.exists() and target.read_bytes()!=data:raise RuntimeError('refusing to replace changed file')
  target.write_bytes(data)
for line in (p/'RECONSTRUCTED_SHA256SUMS').read_text().splitlines():
 digest,name=line.split('  ',1);target=p/name
 if not target.exists() and args.check and any(x['path']==name for x in m['files']):continue
 if hashlib.sha256(target.read_bytes()).hexdigest()!=digest:raise RuntimeError('reconstructed manifest '+name)
print('PASS: '+str(len(m['files']))+' split files verified'+(' without writing.' if args.check else ' and reconstructed.'))
