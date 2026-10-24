#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
import json
P=Path(__file__).resolve().parent
manifest=json.loads((P/'MANIFEST.json').read_text())
actual={str(f.relative_to(P)):sha256(f.read_bytes()).hexdigest() for f in P.rglob('*') if f.is_file() and f!=P/'MANIFEST.json'}
if actual!=manifest['sha256']:raise RuntimeError('manifest differs, including file-set closure')
print('PASS exact standalone manifest; root MANIFEST.json alone excluded; nested manifests included')
