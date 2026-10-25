from pathlib import Path
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent

def verify(root):
    manifest=json.loads((root/'manifest.json').read_text())
    files={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
    if files!=set(manifest)|{'manifest.json'}:raise RuntimeError('exact file membership')
    for name,digest in manifest.items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('digest '+name)

def main():
    verify(ROOT)
    if '--verify-only' not in sys.argv:
        for opt,fixture in [([], 'check.out'),(['-O'],'check-O.out')]:
            got=subprocess.check_output([sys.executable,*opt,str(ROOT/'check.py')],cwd=ROOT)
            if got!=(ROOT/fixture).read_bytes():raise RuntimeError('fixture '+fixture)
    print('PASS closed file set, SHA256 manifest'+(' and normal/-O exact replay' if '--verify-only' not in sys.argv else ''))
if __name__=='__main__':main()
