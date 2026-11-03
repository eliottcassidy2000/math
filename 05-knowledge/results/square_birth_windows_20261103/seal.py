"""Portable text seal for this frozen packet; excludes only its own manifest."""
from pathlib import Path
import hashlib,json,sys

HERE=Path(__file__).resolve().parent
def main():
    files={}
    for p in sorted(HERE.iterdir()):
        if not p.is_file() or p.name=='MANIFEST.json':continue
        data=p.read_bytes().replace(b'\r\n',b'\n')
        files[p.name]={'sha256_lf':hashlib.sha256(data).hexdigest(),'bytes_lf':len(data)}
    result={'normalization':'CRLF to LF for this text-only packet','files':files}
    dest=HERE/'MANIFEST.json'
    if sys.argv[1:]==['--write']:
        dest.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    elif sys.argv[1:]==['--check']:
        if json.loads(dest.read_text(encoding='utf-8'))!=result:raise RuntimeError('packet seal mismatch')
    else:raise RuntimeError('use --write or --check')
    print(json.dumps({'sealed_files':len(files),'mode':sys.argv[1]},sort_keys=True))
if __name__=='__main__':main()
