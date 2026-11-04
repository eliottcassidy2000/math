"""Check source pins, complete finite heads and independent literal receipts."""
from pathlib import Path
import subprocess,sys,json,tempfile,hashlib
HERE=Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def run(args,output=None):
    print('RUN',*map(str,args),flush=True);p=subprocess.run(list(map(str,args)),cwd=HERE,text=True,capture_output=True,check=True)
    if output:(HERE/output).write_text(p.stdout,encoding='utf8',newline='\n')
    else:print(p.stdout[-1600:],end='',flush=True)
def main():
    inputs=['return18979_palindrome_20261103/returns.json','return18979_palindrome_20261103/STABLE.md','return18979_palindrome_20261103/endpoint.cpp','return18979_palindrome_20261103/phases_bridge.py',
        'return18979_ports_20261103/PROOFS.md','return18979_ports_20261103/SQUARES.md','return6325_completion_20261103/halving_response.json','return6325_completion_20261103/transport.py',
        'two_entry_funnels_20261103/READY.md','ready_return_pentagon_20261103/OBSTRUCTION.md','polarized_path_calculus_20261103/recognize.py','polarized_path_calculus_20261103/FUNNELS.md',
        'session87_four_vertices_enumeration_20261028/fast_counts.py','tournament_enumeration_20261021a/check.py','apex_wedge_constructor_20261103/CARRIER.md']
    pins={p:digest(HERE.parent/p) for p in inputs};pinfile=HERE/'source_pins.json'
    if pinfile.exists():
        if json.loads(pinfile.read_text())['sha256']!=pins:raise ValueError('inherited source drift')
    else:pinfile.write_text(json.dumps(dict(inherited_commit='9af43540fd6e',hash_format='CRLF normalized to LF',sha256=pins),indent=2)+'\n')
    bins={}
    for stem in ('funnel_bank','endpoint_checked','audit_feedback'):
        exe=Path(tempfile.gettempdir())/('return56941-'+stem+'.exe');run(['g++','-std=c++20','-O3',HERE/(stem+'.cpp'),'-o',exe]);bins[stem]=exe
    if '--full' in sys.argv:
        for name,args in [('search_funnels.py',[]),('search_two_entry.py',[]),('targets.py',[]),('moving_returns.py',[bins['endpoint_checked']]),('births.py',[]),('feedback.py',[])]:run([sys.executable,'-B',HERE/name,*args])
        for H in (56933,56935,56939,56941):run([sys.executable,'-B',HERE/'two_births.py',H])
    run([sys.executable,'-B',HERE/'search_words.py'])
    outputs=['audit_bank.json','audit_moving.json','audit_clones.json','audit_feedback.json','count_controls.json','feedback_ports.json','feedback_transport.json'];first=None
    for flags in ([],['-O']):
        for name,args in [('audit_bank.py',[bins['funnel_bank'],bins['endpoint_checked']]),('audit_moving.py',[bins['endpoint_checked']]),('audit_clones.py',[bins['endpoint_checked']]),('check_counts.py',[]),('feedback_ports.py',[bins['endpoint_checked']])]:run([sys.executable,'-B','-X','utf8',*flags,HERE/name,*args])
        run([sys.executable,'-B','-X','utf8',*flags,HERE/'feedback_transport.py']);run([bins['audit_feedback']],output='audit_feedback.json');hashes={p:digest(HERE/p) for p in outputs}
        if first is None:first=hashes
        elif first!=hashes:raise ValueError('normal/optimized mathematical drift')
    result=dict(status='PASS_NORMAL_AND_OPTIMIZED_IDENTICAL',source_pin='9af43540fd6e',incoming_pin='83d4954e15',hash_format='CRLF normalized to LF',output_sha256=first,
        discovery_sha256={p:digest(HERE/p) for p in ('one_start_bank.json','two_entry_bank.json','targets.json','moving_returns.json','births.json','feedback.json','feedback_search.json')},
        scripts_sha256={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.cpp')},
        scope='Independent algorithms and same-agent proof review; not an external human audit, a spectrum proof, or a minimal closed carrier.')
    (HERE/'reproduction.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS complete finite heads and normal/O identical')
if __name__=='__main__':main()
