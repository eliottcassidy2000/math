#!/usr/bin/env python3
"""Check exact sibling-link targets and cross-packet audit receipt agreement."""
from pathlib import Path
import json,re
ROOT=Path(__file__).resolve().parent
text=(ROOT/'NOTE.md').read_text()
links=re.findall(r'\]\((\.\./[^)]+)\)',text)
for p in links:
 if not (ROOT/p.split('#')[0]).is_file():raise RuntimeError('Missing sibling link: '+p)
pairs=[('information','certificate.json','certificate.optimized.json'),('spectrum','certificate.json','certificate.optimized.json'),('gates','certificate.json','certificate.optimized.json'),('gates','independent_bit_audit.json','independent_bit_audit.optimized.json'),('circuits','certificate.json','certificate.optimized.json'),('circuits','fiber_independent_audit.json','fiber_independent_audit.optimized.json'),('resource','certificate.json','certificate.optimized.json')]
checked=[]
for folder,a,b in pairs:
 aa=json.loads((ROOT.parent/folder/a).read_text());bb=json.loads((ROOT.parent/folder/b).read_text())
 if aa!=bb:raise RuntimeError('Mode mismatch: '+folder+'/'+a)
 checked.append(folder+'/'+a)
print(json.dumps({'status':'PASS','sibling_links':len(links),'distinct_sibling_links':len(set(links)),'matched_receipt_pairs':checked,'scope':'link and selected receipt consistency, not an independent proof audit'},sort_keys=True,indent=2))
