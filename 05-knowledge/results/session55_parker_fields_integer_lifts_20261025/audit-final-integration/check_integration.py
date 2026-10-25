#!/usr/bin/env python3
"""Check root links against staged map plus earlier published local snapshot."""
from pathlib import Path
import hashlib,json,re
P=Path(__file__).resolve().parent;W=P.parent
M=W/'session55-root/final/PREAUDIT_MAP.json';m=json.loads(M.read_text())
base=W/'session50-publication/clean-repo';previous={}
for checkpoint in ('first','second'):
 for row in json.loads((W/f'session55-publication/{checkpoint}-plan.json').read_text()):previous[row['path']]=row
root=m['destination'];mapped={r['destination']:r for r in m['files']}
def need(t,s):
 if not t:raise RuntimeError(s)
for dst,r in mapped.items():
 p=Path(r['source']);need(p.stat().st_size==r['bytes'],'size '+dst);need(hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],'digest '+dst)
links=[]
for filename in ('README.md','SYNTHESIS.md'):
 source=Path(mapped[root+'/'+filename]['source'])
 for url in re.findall(r'\]\(([^)]+)\)',source.read_text()):
  if '://' in url:continue
  target=root+'/'+url.split('#')[0]
  need(target in mapped or target in previous,'missing link '+target)
  if target not in mapped:
   old=previous[target];need(Path(old['source']).is_file(),'missing frozen source '+target);need(hashlib.sha256(Path(old['source']).read_bytes()).hexdigest()==old['sha256'],'frozen digest '+target)
  links.append({'document':filename,'url':url,'target':target})
# Explicitly validate the new proof packets' dependency pins, not just manifests.
pins=[]
for r in json.loads((W/'session55-dyadic-reentry/SOURCE_PINS.json').read_text())['dependencies']:
 pins.append((r['path'],r['sha256']))
pins.extend(json.loads((W/'session55-ordinary-count-gf/source-pins.json').read_text()).items())
for path,h in pins:need(hashlib.sha256((W/path).read_bytes()).hexdigest()==h,'dependency '+path)
print(json.dumps({'status':'PASS','map_sha256':hashlib.sha256(M.read_bytes()).hexdigest(),'mapped_files':len(mapped),'root_links':links,'new_dependency_pins_checked':len(pins)},indent=2))
