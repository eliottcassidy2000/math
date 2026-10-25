"""Portable source-prefix, historical-snapshot and SVG integrity controls."""
from pathlib import Path
import hashlib,json,xml.etree.ElementTree as ET
root=Path(__file__).resolve().parent
def require(test,message):
 if not test:raise RuntimeError(message)
def sha(b):return hashlib.sha256(b).hexdigest()
rows=json.loads((root/'THIRD_SOURCE_VERIFICATION.json').read_text())['files']
prefixes=[]
for row in rows:
 b=(root/row['packaged']).read_bytes()
 require(sha(b)==row['packaged_sha256'],'packaged source hash '+row['packaged'])
 if row.get('packaging_changes'):
  n=row['source_bytes'];require(sha(b[:n])==row['source_sha256'],'audited prefix '+row['packaged'])
  require(b[n:].startswith(b'\n## Packaged navigation (presentation-only)\n'),'declared suffix '+row['packaged'])
  prefixes.append(row['packaged'])
 else:require(sha(b)==row['source_sha256'],'exact source copy '+row['packaged'])
require(sorted(prefixes)==['PROPOSED_UNIFYING_THEORY.md','SYNTHESIS.md'],'only approved navigation suffixes')
pins=json.loads((root/'recovery/SOURCE_PINS.json').read_text())
for pin in pins:
 relative=Path(pin['local']).relative_to('session63-recovery')
 b=(root/'recovery'/relative).read_bytes()
 require(sha(b)==pin['sha256'],'historical SHA256 '+str(relative))
 require(hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==pin['blob'],'historical Git blob '+str(relative))
svg=root/'figures/integral-observer-map.svg';tree=ET.fromstring(svg.read_bytes())
require(tree.tag=='{http://www.w3.org/2000/svg}svg','SVG namespace')
require(not any(x.tag.rsplit('}',1)[-1]=='script' for x in tree.iter()),'unexpected SVG script')
require(not (root/'figures/integral-observer-map.png').exists(),'PNG is a separate visual-QA artifact, not published')
print(json.dumps({'status':'PASS','copied_artifacts':len(rows),'audited_body_prefixes':prefixes,'historical_snapshots':len(pins),'svg_sha256':sha(svg.read_bytes()),'png_published':False},indent=2,sort_keys=True))
