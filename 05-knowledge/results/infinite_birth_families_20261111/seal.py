"""Create or check canonical LF payload hashes (no executable artifacts)."""
from pathlib import Path
import argparse,hashlib,json

HERE=Path(__file__).resolve().parent

def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    paths=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='SHA256SUMS')
    actual={p.name:digest(p) for p in paths}
    manifest=HERE/'SHA256SUMS'
    if args.check:
        saved={line.split('  ',1)[1]:line.split('  ',1)[0] for line in manifest.read_text().splitlines()}
        if actual!=saved:
            changed=[k for k in sorted(set(actual)|set(saved)) if actual.get(k)!=saved.get(k)]
            raise RuntimeError('payload hash mismatch: '+', '.join(changed))
        print(json.dumps(dict(status='PASS',payload_files=len(paths),normalization='CRLF to LF')))
    else:
        manifest.write_text(''.join(v+'  '+k+'\n' for k,v in actual.items()),encoding='utf-8')
        print('sealed',len(paths),'payload files')

if __name__=='__main__':
    main()
