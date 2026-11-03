"""Replay exact certificates in normal and optimized Python; --full reruns discovery."""
from pathlib import Path
import subprocess,sys,tempfile,json,hashlib
HERE=Path(__file__).resolve().parent
def run(args):
    print('RUN',*map(str,args),flush=True);subprocess.run(list(map(str,args)),check=True,cwd=HERE)
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def main():
    temp=Path(tempfile.gettempdir());bins={}
    for stem in ('ears','audit_ready','search_ready'):
        exe=temp/('return6325-'+stem+'.exe');run(['g++','-std=c++20','-O3',HERE/(stem+'.cpp'),'-o',exe]);bins[stem]=exe
    if '--supplement-rank' in sys.argv:
        report=json.loads((HERE/'reproduction.json').read_text())
        if any(digest(HERE/name)!=value for name,value in report['output_sha256'].items()):raise ValueError('prior checked output changed')
        first=None
        for flags in ([],['-O']):
            run([sys.executable,'-B','-X','utf8',*flags,HERE/'rank_depth.py',bins['ears']])
            value=digest(HERE/'rank_depth.json')
            if first is None:first=value
            elif first!=value:raise ValueError('supplemental rank output drift')
        report['output_sha256']['rank_depth.json']=first
        report['input_sha256']['adaptive_four_returns_20261103/audit.json']=digest(HERE.parent/'adaptive_four_returns_20261103/audit.json')
        report['scripts_sha256']={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.cpp')}
        report['supplement']='Rank-depth theorem added after the full replay; its controls separately passed normal/O, earlier output hashes rechecked unchanged.'
        (HERE/'reproduction.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS supplemental rank normal/O',flush=True);return
    if '--full' in sys.argv:
        for name,args in (('search_ready.py',[bins['search_ready'],'12']),('one_ready_search.py',[bins['search_ready']]),('return6325.py',[]),('return6325_plus2.py',[])):
            run([sys.executable,'-B','-X','utf8',HERE/name,*args])
    files=['audit_ready.json','completion.json','halving_response.json','rank_depth.json','clone_probe.json','transport.json','audit.json'];first=None
    for optimized in (False,True):
        py=[sys.executable,'-B','-X','utf8']+(['-O'] if optimized else [])
        for name,args in (('audit_ready.py',[bins['audit_ready']]),('completion.py',[bins['ears']]),('halving_response.py',[bins['ears']]),('rank_depth.py',[bins['ears']]),('clone_probe.py',[]),('transport.py',[]),('audit.py',[])):
            run([*py,HERE/name,*args])
        hashes={name:digest(HERE/name) for name in files}
        if first is None:first=hashes
        elif hashes!=first:raise ValueError('optimized-mode result drift')
    inputs=['digit_cycle_ternary_20261103/return2107.json','two_entry_funnels_20261103/READY.md','ready_return_pentagon_20261103/ready_inverse.py','polarized_path_calculus_20261103/recognize.py','adaptive_four_returns_20261103/lift_search.cpp','adaptive_four_returns_20261103/audit.json','overnight_rank_return_20261029/atom_precision.py']
    payload=dict(status='PASS_NORMAL_AND_OPTIMIZED_IDENTICAL',hash_format='UTF-8 text bytes with CRLF normalized to LF',output_sha256=first,
        source_pin='31b13dd9a2c8',input_sha256={name:digest(HERE.parent/name) for name in inputs},
        scripts_sha256={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.cpp')},
        scope='Exact saved words and bounded universes; same-agent proof review, independent algorithms, no external proof audit')
    (HERE/'reproduction.json').write_text(json.dumps(payload,indent=2)+'\n');print('PASS normal/O identical',flush=True)
if __name__=='__main__':main()
