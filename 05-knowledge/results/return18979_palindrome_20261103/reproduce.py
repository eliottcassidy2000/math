"""Replay certificates and finite universes in normal and optimized Python."""
from pathlib import Path
import subprocess,sys,json,tempfile,hashlib
HERE=Path(__file__).resolve().parent
def run(args):
    print('RUN',*map(str,args),flush=True);subprocess.run(list(map(str,args)),cwd=HERE,check=True)
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def main():
    bins={}
    for stem in ('funnel_bank','endpoint'):
        exe=Path(tempfile.gettempdir())/('return18979-'+stem+'.exe');run(['g++','-std=c++20','-O3',HERE/(stem+'.cpp'),'-o',exe]);bins[stem]=exe
    if '--full' in sys.argv:
        for name in ('search_funnels.py','search_two_entry.py','targets.py','returns.py'):
            run([sys.executable,'-B','-X','utf8',HERE/name])
    outputs=['audit_bank.json','audit_returns.json','stable_words.json','arithmetic.json'];first=None
    for flags in ([],['-O']):
        for name,args in [('audit_bank.py',[bins['funnel_bank'],bins['endpoint']]),('audit_returns.py',[bins['endpoint']]),('stable_words.py',[bins['endpoint']]),('arithmetic.py',[])]:
            run([sys.executable,'-B','-X','utf8',*flags,HERE/name,*args])
        hashes={p:digest(HERE/p) for p in outputs}
        if first is None:first=hashes
        elif first!=hashes:raise ValueError('normal/optimized output drift')
    inputs=['return6325_completion_20261103/return6325.json','return6325_completion_20261103/return6325_plus2.json','return6325_completion_20261103/halving_response.json','return6325_completion_20261103/transport.json','return6325_completion_20261103/CLONE_LAW.md','two_entry_funnels_20261103/READY.md','polarized_path_calculus_20261103/recognize.py','polarized_path_calculus_20261103/FUNNELS.md','ready_return_pentagon_20261103/OBSTRUCTION.md','golden_inverse_defects_20261103/inverse.py','digit_memory_all_bases_20261103/TOURNAMENT.md']
    result=dict(status='PASS_NORMAL_AND_OPTIMIZED_IDENTICAL',source_pin='d53992645e',hash_format='CRLF normalized to LF',output_sha256=first,
        discovery_sha256={p:digest(HERE/p) for p in ('one_start_bank.json','two_entry_bank.json','targets.json','returns.json')},
        input_sha256={p:digest(HERE.parent/p) for p in inputs},scripts_sha256={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.cpp','.py')},
        scope='Same-agent proof review, independent matrix/DP algorithms; no external audit or all-size spectrum claim')
    incoming=['return18979_ports_20261103/PROOFS.md','session174_palindrome_symmetry_20261103/palindromes/REPORT.md','session174_palindrome_symmetry_20261103/counts7417/README.md','session174_palindrome_symmetry_20261103/synthesis/NOTE.md','all_prime_tournament_geometry_20261103/PROOFS.md']
    result['incoming_pin']='a421cb1ac1';result['incoming_sha256']={p:digest(HERE.parent/p) for p in incoming}
    (HERE/'reproduction.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS normal/O identical',flush=True)
if __name__=='__main__':main()
