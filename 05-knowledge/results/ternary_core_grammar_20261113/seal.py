"""Pin the local packet and directly used repository inputs; verify by default."""
from pathlib import Path
import sys,re,json,hashlib,importlib

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
MANIFEST=HERE/'manifest.json'

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def local_links():
    result=set()
    for doc in HERE.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text(encoding='utf8')):
            if '://' in target or target.startswith('#'):continue
            p=(doc.parent/target.split('#')[0]).resolve()
            if not p.exists():raise RuntimeError('broken link '+str(doc)+': '+target)
            if p.is_file() and p.parent!=HERE:result.add(p)
    return result

def inputs():
    for name in ('triangle_probe','local_obstruction','peeling','atomic_ears','neutral_entry','audit'):
        importlib.import_module(name)
    import triangle_probe
    triangle_probe.bank(7)
    paths=local_links()
    for module in list(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if not name:continue
        p=Path(name).resolve()
        if p.is_relative_to(ROOT) and p.parent!=HERE and p.suffix=='.py':paths.add(p)
    paths.update([HERE.parent/'multicore_metagraph_20261113/rooted_growth.py',
                  HERE.parent/'multicore_metagraph_20261113/rooted_growth.json.gz',
                  HERE.parent/'strong_H_spectrum_m9_values_kps_S134.out'])
    return paths

def main():
    local_links()
    if '--write' in sys.argv:
        own=[p for p in HERE.iterdir() if p.is_file() and p!=MANIFEST]
        data=dict(status='FINITE-EXACT receipts; scoped proofs in notes; Chan-Pak OPEN',
                  packet={p.relative_to(ROOT).as_posix():digest(p) for p in sorted(own)},
                  inputs={p.relative_to(ROOT).as_posix():digest(p) for p in sorted(inputs())})
        MANIFEST.write_text(json.dumps(data,indent=2)+'\n',encoding='utf8')
    data=json.loads(MANIFEST.read_text(encoding='utf8'))
    for group in ('packet','inputs'):
        for name,h in data[group].items():
            if digest(ROOT/name)!=h:raise RuntimeError('changed '+name)
    print('PASS',len(data['packet']),'packet files;',len(data['inputs']),'pinned inputs; local links valid')

if __name__=='__main__':main()
