"""Rebuild the complete research packet; standard Python plus g++ on PATH."""
from pathlib import Path
import subprocess,sys,tempfile,time
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'session67_metagraph_enumeration_information_20261026/session67-root/output'
def py(name,*args):subprocess.run([sys.executable,'-B',str(HERE/name),*map(str,args)],check=True)
with tempfile.TemporaryDirectory(prefix='tournament-aspects-') as td:
    t=Path(td);exe={}
    for name in ('census','screen9','audit9'):
        exe[name]=t/(name+'.exe');subprocess.run(['g++','-O3','-std=c++17',str(HERE/(name+'.cpp')),'-o',str(exe[name])],check=True)
    raw=t/'records.jsonl'
    with raw.open('w',encoding='utf-8') as out:subprocess.run([str(exe['census']),str(SOURCE),'8'],stdout=out,check=True)
    py('analyze.py',raw);py('verify.py',7);py('hostiles.py');py('prepare9.py',t/'input9.txt')
    start=time.perf_counter()
    with (HERE/'screen9.json').open('w') as out:subprocess.run([str(exe['screen9']),str(t/'input9.txt')],stdout=out,check=True)
    (HERE/'screen9-seconds.txt').write_text(str(time.perf_counter()-start)+'\n')
    with (HERE/'audit9.json').open('w') as out:subprocess.run([str(exe['audit9']),str(t/'input9.txt'),str(HERE/'seed_witnesses9.tsv')],stdout=out,check=True)
    py('integrate.py')
    py('seal.py')
