"""Complete small-prefix quotient bank, with independent literal counts."""
from pathlib import Path
from itertools import permutations
import json, subprocess, tempfile

HERE = Path(__file__).resolve().parent

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def literal(rows):
    n = len(rows)
    paths = [0] * (1 << n)
    paths[0] = 1
    for mask in range(1, 1 << n):
        vertices = [i for i in range(n) if mask >> i & 1]
        paths[mask] = sum(all(rows[a] >> b & 1 for a,b in zip(word, word[1:]))
                          for word in permutations(vertices))
    full = (1 << n)-1
    return paths[-1], sum(paths[m]*paths[full ^ m] for m in range(full+1))

with tempfile.TemporaryDirectory(prefix='c4-prefix-') as directory:
    executable = Path(directory)/'probe.exe'
    subprocess.run(['g++','-O3','-std=c++20',str(HERE/'probe.cpp'),'-o',str(executable)],check=True)
    data = json.loads(subprocess.check_output([str(executable)],text=True))
need([row['classes'] for row in data['orders']] == [1,1,2,4,12,56,456], 'complete class controls')
summary = {h: set() for h in (5,29,35)}
checked = 0
for order in data['orders']:
    for item in order['selected']:
        need(literal(item['rows']) == (item['h'],item['p']), 'independent H/P2 count')
        odd = item['p']
        while odd % 2 == 0:
            odd //= 2
        item['odd_part_p'] = odd
        summary[item['h']].add(odd)
        checked += 1
data['summary_odd_parts'] = {h: sorted(values) for h,values in summary.items()}
data['independent_permutation_checks'] = checked
data['scope'] = 'All strong classes through order seven; larger-order prefix reduction requires the stated Busch bound and SCC factorization.'
(HERE/'receipt.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'status':'PASS','classes':[row['classes'] for row in data['orders']],
                  'independent_checks':checked,'odd_parts':data['summary_odd_parts']}))
