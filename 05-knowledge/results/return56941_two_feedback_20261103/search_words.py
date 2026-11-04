"""Complete, monotone-pruned two-feedback search with an explicit receipt."""
from pathlib import Path
import subprocess,tempfile,json,re,hashlib
HERE=Path(__file__).resolve().parent
def main():
    exe=Path(tempfile.gettempdir())/'return56941-words.exe'
    subprocess.run(['g++','-std=c++20','-O3',str(HERE/'feedback_words.cpp'),'-o',str(exe)],check=True)
    p=subprocess.run([str(exe),'13','all'],capture_output=True,text=True,check=True)
    m=re.fullmatch(r'all_positions visits (\d+) pruned (\d+) prefixes (\d+) hits (\d+)\s*',p.stderr)
    if m is None:raise ValueError('missing complete-census receipt')
    counts=list(map(int,m.groups()));hits=[json.loads(x) for x in p.stdout.splitlines()]
    if counts!=[97648782,62397763,10833829,0] or hits:raise ValueError('full search differs')
    (HERE/'feedback_all_positions.jsonl').write_text(p.stdout,encoding='utf8',newline='\n')
    report=dict(status='FINITE_EXACT_COMPLETE_SEARCH_EXIT_ZERO',exterior_chain_bound=13,bulk_positions='every gap',bulk_letters=[1,3],other_letters='covered by global converse and apex exchange',apex_edges=[0,1],target_samples=[[4277,11735,26651],[4279,11737,26653]],visits=counts[0],pruned=counts[1],prefixes=counts[2],hits=hits,pruning='Hamilton insertion is injective; every sample can only increase',source_sha256=hashlib.sha256((HERE/'feedback_words.cpp').read_bytes().replace(b'\r\n',b'\n')).hexdigest())
    (HERE/'feedback_search.json').write_text(json.dumps(report,indent=2)+'\n');print(p.stderr.strip(),flush=True)
if __name__=='__main__':main()
