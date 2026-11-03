#!/usr/bin/env python3
"""Verify all pinned bytes and cross-check the independently produced receipts."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent

def demand(ok,message):
 if not ok:raise RuntimeError(message)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((root/'MANIFEST.json').read_text())
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts}
demand(actual==set(manifest['files']),'artifact whitelist mismatch')
for path,expected in manifest['files'].items():demand(sha(root/path)==expected,'checksum mismatch: '+path)
for record in json.loads((root/'SOURCE_PINS.json').read_text()):
 data=(root/record['local_copy']).read_bytes()
 demand(hashlib.sha256(data).hexdigest()==record['sha256'],'source SHA mismatch')
 demand(hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==record['git_blob'],'source Git blob mismatch')
for stem in ['verification','independent']:
 demand((root/(stem+'.normal.json')).read_bytes()==(root/(stem+'.optimized.json')).read_bytes(),'optimization changed checks')
a=json.loads((root/'verification.normal.json').read_text());b=json.loads((root/'independent.normal.json').read_text())
demand(a['status']==b['status']=='PASS','receipt status')
demand(a['census'][-1]['no_red_or_blue_triangle']==b['K5_avoiding']==29197,'census disagreement')
demand(a['paley11']['H']==b['paley11']['H']==95095,'Hamilton count disagreement')
demand(a['paley11']['endpoint_row']==b['paley11']['endpoint_row_0'],'endpoint disagreement')
demand(a['paley11']['specified_arc_incidence']==b['paley11']['arc_0_1']==17290,'arc count disagreement')
print(json.dumps({'status':'PASS','files':len(actual),'producer_checks':sum(a['checks'].values()),'independent_checks':sum(b['checks'].values()),'normal_optimized_equal':True},sort_keys=True))
