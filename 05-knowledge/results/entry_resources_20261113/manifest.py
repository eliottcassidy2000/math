"""Pin this packet and its actual imported repository dependencies."""
from pathlib import Path
import sys,json,hashlib,importlib.metadata,platform,subprocess
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];MANIFEST=HERE/'manifest.json'

def digest(path):
    data=path.read_bytes()
    if path.suffix in ['.py','.md','.json','.svg']:data=data.replace(b'\r\n',b'\n')
    return hashlib.sha256(data).hexdigest()

def write():
    from reproduce import JOBS
    paths={p.resolve() for p in HERE.iterdir() if p.is_file() and p!=MANIFEST and p.suffix in ['.py','.md','.json','.gz','.svg','.png']}
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
    inherited=[
      'endpoint_coverage_compiler_20261111/inflation_bank_7.json',
      'endpoint_coverage_compiler_20261111/wheel.json',
      'graceful_fixed_color_counterexample_20261019/witness.json',
      'graceful_connected_color_contacts_20261019/CONNECTED_COLOR_CONTACTS.md',
      'joint_incidence_clocked_frames_20261102/path_kernel.py',
      'joint_stability_forbidden_20261112/README.md',
      'session288_twin_lengths_20261112/audit/REVIEW.md',
      'session289_joint_births_20261112/absence/PROOF.md',
      'session291_cluster_ports_20261113/SYNTHESIS.md',
      'session292_triple_bridges_20261113/arithmetic/PROOF.md',
      'session293_cycle_returns_20261113/control/MODULE_IDEAL.md',
      'session293_cycle_returns_20261113/audit/REVIEW_IDEAL.md',
      'core_descent_endpoints_20261112/PROOFS.md']
    paths.update(HERE.parent/name for name in inherited)
    entries=[{'path':p.relative_to(ROOT).as_posix(),'sha256':digest(p)} for p in sorted(paths)]
    data={'scope':'CRLF normalized for text; binary hashes exact. Own packet, actual imported repository engines, consumed inherited data and principal inherited proof statements. Logs excluded.','python':platform.python_version(),'packages':{n:importlib.metadata.version(n) for n in ['numpy','numba']},'optional_figure_runtime':{'matplotlib':'3.11.2','scope':'Separate temporary installation; not needed for mathematical replay.'},'files':entries}
    MANIFEST.write_text(json.dumps(data,indent=2)+'\n',encoding='utf8',newline='\n');print('WROTE',len(entries),'pins')

def check():
    data=json.loads(MANIFEST.read_text(encoding='utf8'));bad=[e['path'] for e in data['files'] if not (ROOT/e['path']).is_file() or digest(ROOT/e['path'])!=e['sha256']]
    if bad:raise RuntimeError('Manifest differences: '+repr(bad))
    print('PASS',len(data['files']),'pins')

if __name__=='__main__':
    if sys.argv[1:]==['--write']:write()
    elif sys.argv[1:] in [[],['--check']]:check()
    else:raise SystemExit('Use --write or --check.')
