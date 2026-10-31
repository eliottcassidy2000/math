"""Verify the explicit publication whitelist and all non-self-referential hashes."""
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
def verify():
 files=(HERE/'PUBLIC_FILES').read_text().splitlines();meta=json.loads((HERE/'PUBLIC_FILES.json').read_text());manifest=json.loads((HERE/'MANIFEST.json').read_text())
 if len(files)!=len(set(files)) or files!=meta['files']:raise RuntimeError('publication whitelist')
 expected=set(files)-{'PUBLIC_FILES','PUBLIC_FILES.json'}
 if set(meta['sha256'])!=expected:raise RuntimeError('publication hash coverage')
 for name,h in meta['sha256'].items():
  if Path(name).is_absolute() or '..' in Path(name).parts:raise RuntimeError('unsafe relative path')
  if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=h:raise RuntimeError('publication hash:'+name)
 for name,h in manifest['files'].items():
  if name not in expected or hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=h:raise RuntimeError('proof manifest:'+name)
 return {'status':'passed','files':len(files),'hashed_files':len(expected)}
if __name__=='__main__':print(json.dumps(verify(),indent=2))
