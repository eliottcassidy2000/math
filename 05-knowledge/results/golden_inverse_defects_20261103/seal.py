"""Check the exact packet file set, local links and LF-normalized hashes."""
import argparse
import hashlib
import json
import re
from pathlib import Path

HERE=Path(__file__).resolve().parent
FILES=tuple(sorted(('README.md','THEORY.md','TRANSPORT.md','RECOVERY.md','AUDIT.md',
                    'inverse.py','verify.py','seal.py','receipt.json')))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    present={p.name for p in HERE.iterdir() if p.is_file()}
    if present-{'MANIFEST.json'}!=set(FILES):
        raise RuntimeError(('packet file set',sorted((present-{'MANIFEST.json'})^set(FILES))))
    links=0
    for name in FILES:
        if not name.endswith('.md'):continue
        for target in re.findall(r'\[[^\]]+\]\(([^)\s]+)\)',(HERE/name).read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            if not (HERE/target.split('#',1)[0]).exists():
                raise RuntimeError(('missing local link',name,target))
            links+=1
    hashes={name:hashlib.sha256((HERE/name).read_text(encoding='utf-8').encode()).hexdigest() for name in FILES}
    data={'hash_basis':'UTF-8 text with universal-newline LF normalization','sha256':hashes}
    path=HERE/'MANIFEST.json'
    if args.write:path.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    elif json.loads(path.read_text(encoding='utf-8'))!=data:raise RuntimeError('packet content hash mismatch')
    print(json.dumps({'status':'PASS','sealed_files':len(FILES),'local_links':links}))


if __name__=='__main__':main()
