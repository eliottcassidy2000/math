"""Lossless receipt transport. No third-party packages, network, or inference."""
from pathlib import Path
import argparse,base64,gzip,hashlib,json
ROOT=Path(__file__).resolve().parent

def receipt_bytes(name):
    meta=json.loads((ROOT/'receipt_archives.json').read_text())[name]
    raw=gzip.decompress(base64.b64decode(''.join((ROOT/meta['archive']).read_text().split()),validate=True))
    if len(raw)!=meta['bytes'] or hashlib.sha256(raw).hexdigest()!=meta['sha256']:
        raise RuntimeError('receipt archive commitment mismatch: '+name)
    return raw

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--destination',type=Path,default=ROOT);a=ap.parse_args()
    for name in json.loads((ROOT/'receipt_archives.json').read_text()):
        p=a.destination/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(receipt_bytes(name));print(p)
if __name__=='__main__':main()
