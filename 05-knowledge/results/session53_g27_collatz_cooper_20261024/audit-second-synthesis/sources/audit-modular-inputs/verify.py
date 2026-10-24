"""Closed-manifest validation and exact independent replay."""
import hashlib,json,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parent

def require(ok,msg):
    if not ok: raise ValueError(msg)
def unique_pairs(pairs):
    d={}
    for k,v in pairs:
        require(k not in d,'duplicate JSON key')
        d[k]=v
    return d

def validate(root):
    m=json.loads((root/'MANIFEST.json').read_text(),object_pairs_hook=unique_pairs)
    require(set(m)=={'schema','files'},'manifest keys')
    require(m['schema']=='session53-modular-inputs-v1','schema')
    declared=set()
    for rec in m['files']:
        require(set(rec)=={'path','sha256','bytes'},'record keys')
        p=rec['path']; pp=pathlib.PurePosixPath(p)
        require(isinstance(p,str) and p and not pp.is_absolute() and '..' not in pp.parts and '.' not in pp.parts and '\\' not in p and str(pp)==p,'unsafe path')
        require(p!='MANIFEST.json' and p not in declared,'duplicate/self path')
        declared.add(p)
        f=root/p
        require(f.is_file() and not f.is_symlink(),'missing/unsafe payload')
        require(all(not part.is_symlink() for part in [f,*f.parents] if part!=root.parent),'symlink path')
        b=f.read_bytes()
        require(len(b)==rec['bytes'] and hashlib.sha256(b).hexdigest()==rec['sha256'],'changed payload: '+p)
    actual={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() or f.is_symlink()}
    require(actual==declared|{'MANIFEST.json'},'extra/missing payload')
    return m

def main():
    validate(ROOT)
    if '--manifest-only' not in sys.argv:
        for mode,name in [([], 'normal.json'),(['-O'],'optimized.json')]:
            out=subprocess.check_output([sys.executable,*mode,str(ROOT/'certificate.py')],cwd=ROOT)
            require(out==(ROOT/name).read_bytes(),'replay differs: '+name)
        require((ROOT/'normal.json').read_bytes()==(ROOT/'optimized.json').read_bytes(),'optimization changes result')
    print('PASS: manifest closure'+('' if '--manifest-only' in sys.argv else ' and independent normal/-O exact replay'))
if __name__=='__main__': main()
