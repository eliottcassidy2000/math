#!/usr/bin/env python3
"""Read-only SHA256 check of the frozen packets used by the synthesis."""
from pathlib import Path
import argparse,hashlib,json
parser=argparse.ArgumentParser()
parser.add_argument("--root-dir",type=Path,default=Path(__file__).resolve().parent.parent,help="Workspace root or publication bundle root")
ROOT=parser.parse_args().root_dir.resolve()
ALIASES={"session47-helper-resilience/final-v3":"helper-resilience-v3","session47-audit-helper-resilience/third-audit":"audit-helper-resilience-v3","session47-audit-suspension-obstruction":"audit-five-to-seven","session47-audit-fresh570":"fresh570-audit"}
def packet_path(name):
 p=ROOT/name
 if p.exists():return p
 return ROOT/ALIASES.get(name,name.removeprefix("session47-"))
packets=['session47-helper-resilience/final-v3','session47-audit-helper-resilience/third-audit','session47-clock-connectivity','session47-owner-cycle-primes','session47-carrier-transfer','session47-carrier-transfer/five-to-seven','session47-audit-carrier-transfer','session47-audit-suspension-obstruction','session47-suspension-obstruction','session47-audit-wiring-theory','session47-global-observer','session47-audit-global-observer','session47-audit-fresh570','session47-bottleneck-theory']
reports=[]
for packet in packets:
 base=packet_path(packet);manifest=base/'MANIFEST.json';b=manifest.read_bytes();m=json.loads(b)
 entries=m.get('files_sha256',m.get('files',m.get('artifacts')))
 if isinstance(entries,dict):
  entries=[dict(path=k,sha256=v) if isinstance(v,str) else dict(path=k,**v) for k,v in entries.items()]
 if not isinstance(entries,list):raise RuntimeError(f'Unknown manifest format: {packet}')
 checked=[]
 for e in entries:
  name=e['path'];p=(base/name).resolve()
  if not p.is_relative_to(base.resolve()):raise RuntimeError(f'Manifest path escapes packet: {p}')
  raw=p.read_bytes();h=hashlib.sha256(raw).hexdigest()
  if h!=e['sha256']:raise RuntimeError(f'Manifest hash mismatch: {p}')
  if 'bytes' in e and len(raw)!=e['bytes']:raise RuntimeError(f'Byte count mismatch: {p}')
  checked.append(name)
 report={'packet':packet,'manifest_sha256':hashlib.sha256(b).hexdigest(),'checked_files':len(checked),'paths':checked}
 if packet=='session47-audit-fresh570':
  pins=json.loads((base/'SOURCE_PINS.json').read_text())['files']
  for pin in pins:
   name=pin['snapshot']
   if name not in checked:raise RuntimeError('Source-pinned snapshot omitted from delivery: '+name)
   raw=(base/name).read_bytes()
   if len(raw)!=pin['bytes'] or hashlib.sha256(raw).hexdigest()!=pin['sha256']:raise RuntimeError('Source-pin snapshot mismatch: '+name)
  report['source_pinned_snapshots_checked']=len(pins)
 reports.append(report)
print(json.dumps({'status':'PASS','packet_count':len(reports),'files_checked':sum(r['checked_files'] for r in reports),'packets':reports},indent=2,sort_keys=True))
