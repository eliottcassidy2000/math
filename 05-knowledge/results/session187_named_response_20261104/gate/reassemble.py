import pathlib, json, hashlib, base64, gzip
p=pathlib.Path(__file__).resolve().parent
m=json.loads((p/'REASSEMBLY.json').read_text())
chunks=[]
for part in m['parts']:
 b=(p/part['name']).read_bytes()
 if hashlib.sha256(b).hexdigest()!=part['sha256']: raise ValueError('chunk hash mismatch')
 chunks.append(b)
b=gzip.decompress(base64.b64decode(b''.join(chunks),validate=True))
if len(b)!=m['bytes'] or hashlib.sha256(b).hexdigest()!=m['sha256']: raise ValueError('reassembled hash mismatch')
(p/'full_response.json').write_bytes(b)
print('Verified full_response.json:',len(b),'bytes; SHA256',m['sha256'])
