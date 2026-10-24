#!/usr/bin/env python3
"""Portable normal/optimized replay. No network, git, producer or sibling imports."""
from pathlib import Path
from hashlib import sha256
import json,subprocess,sys
P=Path(__file__).resolve().parent
M=json.loads((P/'MANIFEST.json').read_text())
for item in M['files']:
 b=(P/item['path']).read_bytes()
 if len(b)!=item['bytes'] or sha256(b).hexdigest()!=item['sha256']:raise RuntimeError('manifest mismatch '+item['path'])
results=[]
for script,stem in [('verify_native.py','native'),('replay_filtration.py','filtration'),('audit.py','audit')]:
 normal=subprocess.check_output([sys.executable,'-B',str(P/script)],cwd='/tmp')
 optimized=subprocess.check_output([sys.executable,'-O','-B',str(P/script)],cwd='/tmp')
 if normal!=optimized:raise RuntimeError('normal/optimized divergence '+script)
 if normal!=(P/(stem+'.normal.json')).read_bytes() or optimized!=(P/(stem+'.optimized.json')).read_bytes():raise RuntimeError('frozen concise receipt mismatch '+script)
 results.append(json.loads(normal))
print(json.dumps({'status':'PASS','manifest_files':len(M['files']),'network':False,'git':False,'sibling_workspace_needed':False,'results':results},sort_keys=True))
