"""Replay both paths in both modes against deterministic UTF-8/LF receipts."""
import hashlib
from pathlib import Path
import subprocess
import sys

p = Path(__file__).resolve().parent
manifest = p/'SHA256SUMS'
if manifest.exists():
    for line in manifest.read_text(encoding='utf-8').splitlines():
        expected, name = line.split('  ',1)
        if hashlib.sha256((p/name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('manifest mismatch: '+name)
    print('SHA256SUMS: all delivered packet files match')
for script, normal, optimized in [('check.py','normal.json','optimized.json'),
                                  ('audit.py','audit-normal.json','audit-optimized.json')]:
    outputs = []
    for flags, name in [([],normal), (['-O'],optimized)]:
        run = subprocess.run([sys.executable,'-B',*flags,str(p/script)], check=True,
                             capture_output=True, text=True, encoding='utf-8')
        output = run.stdout.replace('\r\n','\n')
        expected = (p/name).read_text(encoding='utf-8')
        if output != expected:
            raise RuntimeError('frozen receipt differs: '+name)
        outputs.append(output)
    if outputs[0] != outputs[1]:
        raise RuntimeError('normal/optimized mismatch: '+script)
    print(script+': normal = optimized; sha256 '+hashlib.sha256(outputs[0].encode()).hexdigest())
