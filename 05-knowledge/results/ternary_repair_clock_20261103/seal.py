"""Exact packet membership, local links, and LF-normalized content hashes."""
import argparse,hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
FILES=tuple(sorted(('README.md','THEORY.md','CARRIES.md','AUDIT.md','RECOVERY.md',
    'research.py','native.py','seal.py','atlas7.cpp','audit_atlas.cpp','saturate.cpp',
    'check_closure.cpp','root_probe.cpp','receipt.json')))
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    actual={p.name for p in HERE.iterdir() if p.is_file()}-{'MANIFEST.json'}
    if actual!=set(FILES):raise RuntimeError(('file set',sorted(actual^set(FILES))))
    links=0
    for name in FILES:
        if not name.endswith('.md'):continue
        for target in re.findall(r'\[[^\]]+\]\(([^)\s]+)\)',(HERE/name).read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            if not(HERE/target.split('#',1)[0]).exists():raise RuntimeError(('missing link',name,target))
            links+=1
    data={'hash_basis':'UTF-8 text with universal-newline LF normalization',
          'sha256':{name:hashlib.sha256((HERE/name).read_text(encoding='utf-8').encode()).hexdigest() for name in FILES}}
    path=HERE/'MANIFEST.json'
    if args.write:path.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    elif json.loads(path.read_text())!=data:raise RuntimeError('content hash mismatch')
    print(json.dumps({'status':'PASS','sealed_files':len(FILES),'local_links':links}))
if __name__=='__main__':main()
