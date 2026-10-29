"""Hash the packet, inherited executable dependencies, and primary proof inputs."""
from pathlib import Path
import sys,json,hashlib
import ports,paley,crt,counting,arithmetic,audit
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if '--write' in sys.argv:
 files={p for p in HERE.iterdir() if p.is_file() and p.name!='seal.json'}
 for module in list(sys.modules.values()):
  name=getattr(module,'__file__',None)
  if name:
   p=Path(name).resolve()
   if p.is_relative_to(ROOT) and p.suffix=='.py':files.add(p)
 for name in ('cursor_cm_phase_20261028/endpoint_kernel.json','cursor_cm_phase_20261028/ENDPOINTS.md',
              'lucas_prime_reflection_towers_20261028/FULL_HOLONOMY.md',
              'lucas_prime_reflection_towers_20261028/PORT_MOMENTS.md'):
  files.add(HERE.parent/name)
 files.add(ROOT/'01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md')
 data={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}
 (HERE/'seal.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
 print('Sealed',len(data),'files')
else:
 data=json.loads((HERE/'seal.json').read_text())
 bad=[name for name,h in data.items() if not (ROOT/name).is_file() or sha(ROOT/name)!=h]
 if bad:raise RuntimeError('Changed pinned files: '+', '.join(bad))
 print('Verified',len(data),'file hashes')
