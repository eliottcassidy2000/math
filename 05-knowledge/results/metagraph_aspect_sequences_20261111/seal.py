"""Fingerprint the packet using committed LF text and unchanged gzip bytes."""
from pathlib import Path
import hashlib
here=Path(__file__).resolve().parent
lines=[]
for p in sorted(here.iterdir()):
    if not p.is_file() or p.name=='SHA256SUMS':continue
    data=p.read_bytes()
    if p.suffix!='.gz':data=data.replace(b'\r\n',b'\n')
    lines.append(hashlib.sha256(data).hexdigest()+'  '+p.name)
(here/'SHA256SUMS').write_text('\n'.join(lines)+'\n',encoding='ascii')
print('Sealed',len(lines),'payload files; text hashes normalize CRLF to LF.')
