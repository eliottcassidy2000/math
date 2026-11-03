"""Check exact packet files, local links and normalized content hashes."""
import argparse,hashlib,json,os,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
FILES=tuple(sorted(('README.md','PROJECTIVE.md','TRANSPORT.md','LIVE_FIBRES.md',
 'RECOVERY.md','AUDIT.md','verify.py','audit.cpp','monoid.cpp','receipt.json','seal.py')))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    actual={p.name for p in HERE.iterdir() if p.is_file()}-{'MANIFEST.json'}
    if actual!=set(FILES):raise RuntimeError(('packet membership',actual^set(FILES)))
    links=0
    for name in FILES:
        if not name.endswith('.md'):continue
        for target in re.findall(r'\[[^\]]+\]\(([^)\s]+)\)',(HERE/name).read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            absolute=os.path.normpath(str(HERE/target.split('#')[0]))
            checked=Path('\\\\?\\'+absolute) if os.name=='nt' else Path(absolute)
            if not checked.exists():raise RuntimeError(('local link',name,target))
            links+=1
    data={'hash_basis':'UTF-8 with universal-newline LF normalization',
          'sha256':{n:hashlib.sha256((HERE/n).read_text(encoding='utf-8').encode()).hexdigest() for n in FILES}}
    path=HERE/'MANIFEST.json'
    if args.write:path.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    elif json.loads(path.read_text())!=data:raise RuntimeError('packet hash mismatch')
    print(json.dumps({'status':'PASS','sealed_files':len(FILES),'local_links':links}))
if __name__=='__main__':main()
