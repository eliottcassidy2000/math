#!/usr/bin/env python3
"""Check a published Git tree against this packet's manifest and local files."""
import argparse,hashlib,json,posixpath,re,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PFX='05-knowledge/results/session103_words_necklaces_20261029/'
p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--commit',required=True);a=p.parse_args()
def git(*args):return subprocess.check_output(['git','-C',a.repo,*args])
rows=[];links=0
manifest=json.loads((ROOT/'MANIFEST.json').read_text())
for path,sha in manifest.items():
 local=(ROOT/path).read_bytes();remote=git('show',a.commit+':'+PFX+path)
 if local!=remote or hashlib.sha256(remote).hexdigest()!=sha:raise RuntimeError('byte/hash '+path)
 rows.append(path)
 if path.endswith('.md'):
  in_code=False;kept=[]
  for line in remote.decode().splitlines():
   if line.lstrip().startswith('```'):in_code=not in_code;continue
   if not in_code and not line.startswith('    '):kept.append(line)
  for target in re.findall(r'\[[^\]]+\]\(([^)]+)\)','\n'.join(kept)):
   if '://' in target or target.startswith('#'):continue
   target=target.split('#',1)[0]
   resolved=posixpath.normpath(posixpath.join(PFX,posixpath.dirname(path),target))
   git('cat-file','-e',a.commit+':'+resolved);links+=1
remote=git('show',a.commit+':'+PFX+'MANIFEST.json')
if remote!=(ROOT/'MANIFEST.json').read_bytes():raise RuntimeError('manifest mismatch')
print(json.dumps({'status':'VERIFIED','commit':a.commit,'files':len(rows)+1,'manifest_hashes':len(rows),'relative_links':links},indent=2,sort_keys=True))
