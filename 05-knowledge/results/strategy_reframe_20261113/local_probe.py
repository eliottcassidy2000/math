"""Full finite modular orbits, not truncated integer-word samples."""
from pathlib import Path
from collections import deque
import gzip,json

HERE=Path(__file__).resolve().parent
PRIOR=HERE.parent/'endpoint_context_recursion_20261113'

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def orbit(m,matrices):
    start=(5%m,9%m);seen={start};queue=deque([start])
    while queue:
        h,E=queue.popleft()
        for ((a,b),(c,d)) in matrices:
            nxt=((a*h+b*E)%m,(c*h+d*E)%m)
            if nxt not in seen:seen.add(nxt);queue.append(nxt)
    return seen

def run():
    d=json.loads(gzip.decompress((PRIOR/'renewal_catalog.json.gz').read_bytes()))
    mats=[row['target_matrix'] for row in d['records']]
    mods=list(range(2,101))+[128,210,256,315,420]
    rows=[]
    for m in mods:
        states=orbit(m,mats);hs={h for h,E in states}
        allowed=set(range(m)) if m%2 else set(range(1,m,2))
        missing=sorted(allowed-hs)
        rows.append(dict(modulus=m,states=len(states),target_residues=len(hs),missing_allowed=missing,allows_7=7%m in hs,allows_21=21%m in hs,allows_53084231=53084231%m in hs))
    # Independent exhaustive words through depth5 stay in each complete orbit.
    smallmods=(3,4,5,7,8,9,11)
    controls=0
    for m in smallmods:
        states=orbit(m,mats);level={(5,9)}
        for depth in range(6):
            need(all((h%m,E%m) in states for h,E in level),'integer words versus finite orbit')
            controls+=len(level)
            if depth<5:level={ (a*h+b*E,c*h+d*E) for h,E in level for ((a,b),(c,d)) in mats }
    result=dict(status='FINITE-EXACT',scope='Complete modular orbit under the nine inherited matrices, from (5,9), for each listed modulus. No all-moduli claim.',rows=rows,word_controls=controls)
    suffix='' if __debug__ else '.optimized'
    (HERE/('local'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(moduli=len(mods),not_full=[r for r in rows if r['missing_allowed']],word_controls=controls,examples=[r for r in rows if r['modulus'] in (8,45,97,420)])),flush=True)

if __name__=='__main__':run()
