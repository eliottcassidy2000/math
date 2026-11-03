"""Independent complete count/rank replay of the three sharp query banks."""
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
import json,shutil,subprocess
HERE=Path(__file__).resolve().parent

def main():
    data=json.loads((HERE/'forest_observations_five.json').read_text())
    last=data['actual_growth_ranks'][-1]
    banks={3:list(product(range(4),repeat=0)),4:[],5:[]}
    banks[3]=[w for k in range(3) for w in product(range(4),repeat=k)]
    banks[4]=[w for k in range(4) for w in product(range(5),repeat=k)]
    banks[5]=[w for k in range(4) for w in product(range(6),repeat=k)]+last['depth4_probe']['tested_words']
    chunks=[]
    for n,words in banks.items():
        chunks.append(f'{n} {len(words)}\n')
        chunks.extend(str(len(w))+' '+ ' '.join(map(str,w))+'\n' for w in words)
    compiler=shutil.which('g++')
    if not compiler:raise RuntimeError('A C++17 g++ compiler is required')
    with TemporaryDirectory(prefix='math-forest-rank-') as scratch:
        binary=Path(scratch)/'forest_rank_audit.exe'
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'forest_rank_audit.cpp'),'-o',str(binary)],check=True,capture_output=True)
        raw=subprocess.run([str(binary)],input=''.join(chunks),text=True,capture_output=True,check=True).stdout
    result=json.loads(raw)
    for row,target in zip(result['audits'],(7,28,126)):
        if row['rank_mod101']!=target or row['rank_mod103']!=target:raise RuntimeError('Independent rank mismatch')
    if len(result['audits'])!=3:raise RuntimeError('Missing audit bank')
    (HERE/'forest_rank_independent.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result))

if __name__=='__main__':main()
