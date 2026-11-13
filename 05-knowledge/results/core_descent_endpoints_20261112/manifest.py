"""Pin the exact research packet and inherited code/data; check by default."""
from pathlib import Path
import sys,json,hashlib,importlib.metadata,platform,subprocess
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
MANIFEST=HERE/'manifest.json'

def digest(path):
    data=path.read_bytes()
    if path.suffix in ['.py','.md','.json','.svg']:data=data.replace(b'\r\n',b'\n')
    return hashlib.sha256(data).hexdigest()

def write():
    from reproduce import JOBS
    paths={p.resolve() for p in HERE.iterdir() if p.is_file() and p!=MANIFEST and (p.suffix in ['.py','.md','.json','.gz','.svg','.png'])}
    # Each producer runs in its own interpreter, as in reproduce.py. Several
    # inherited folders have files called verify/common; never merge their
    # import namespaces just to discover dependencies.
    probe="""import sys,json,importlib.util
from pathlib import Path
p=Path(sys.argv[1]);sys.path.insert(0,str(p.parent))
s=importlib.util.spec_from_file_location('__pin_source__',p)
m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m)
print(json.dumps([str(Path(x.__file__).resolve()) for x in list(sys.modules.values()) if getattr(x,'__file__',None)]))
"""
    for script in dict.fromkeys(script for script,_,_ in JOBS):
        raw=subprocess.check_output([sys.executable,'-X','utf8','-c',probe,str(HERE/script)],text=True)
        for f in json.loads(raw):
            p=Path(f)
            if p.suffix=='.py' and p.is_relative_to(ROOT):paths.add(p)
    inherited=['hamilton_clocks_9','even_valuations','repaired_rays',
               'birth_extensions_11','birth_extensions_12','birth_extensions_13',
               'birth_extensions_13_synchrony','birth_extensions_14_rank5',
               'birth_sources_synchrony_12','mixed_grammar']
    paths.update(HERE.parent/'birth_grammar_repairs_20261112'/(n+'.json.gz') for n in inherited)
    paths.add(HERE.parent/'adaptive_state_induction_20261103'/'two_port_envelope.json')
    # Loaded during connection_checks.run via importlib without registration.
    paths.add(HERE.parent/'joint_incidence_clocked_frames_20261102'/'path_kernel.py')
    entries=[{'path':p.relative_to(ROOT).as_posix(),'sha256':digest(p)} for p in sorted(paths)]
    packages={n:importlib.metadata.version(n) for n in ['numpy','numba','sympy','networkx','mpmath']}
    out={'scope':'SHA256 pins with CRLF normalized to LF for .py/.md/.json/.svg; binary files exact. Packet and runtime-imported repository engines/data; logs and manifest excluded. Plot packages optional.','python':platform.python_version(),'packages':packages,'files':entries}
    MANIFEST.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('WROTE',len(entries),'file pins')

def check():
    data=json.loads(MANIFEST.read_text(encoding='utf8'));bad=[]
    for e in data['files']:
        p=ROOT/e['path']
        if not p.is_file() or digest(p)!=e['sha256']:bad.append(e['path'])
    if bad:raise RuntimeError('Manifest differences: '+repr(bad))
    print('PASS',len(data['files']),'file pins')

if __name__=='__main__':
    if sys.argv[1:]==['--write']:write()
    elif sys.argv[1:] in [[],['--check']]:check()
    else:raise SystemExit('Use --write to create pins or --check to validate them.')
