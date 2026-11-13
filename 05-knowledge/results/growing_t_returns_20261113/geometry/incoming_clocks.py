"""Literal readiness barriers with immutable incoming Git-blob provenance."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def table(rows):
    n=len(rows); full=(1<<n)-1
    e=[[0]*n for _ in range(full+1)]; h=[1]+[0]*full
    for mask in range(1,full+1):
        for v in range(n):
            if mask>>v&1:
                rest=mask^(1<<v)
                e[mask][v]=sum(e[rest][u] for u in range(n) if rows[u]>>v&1) if rest else 1
        h[mask]=sum(e[mask])
    return e,h

def strong(rows):
    full=(1<<len(rows))-1
    for r in range(len(rows)):
        reached=1<<r
        while True:
            newer=reached
            for v,row in enumerate(rows):
                if reached>>v&1: newer|=row
            if newer==reached: break
            reached=newer
        if reached!=full: return False
    return True

pi=[462,508,120,496,481,449,384,260,4]
psi=[974,1020,568,1008,993,961,644,772,68,256]
from_order=[sum(1<<j for j in range(i+1,9)) for i in range(9)]
for i,j in [(0,4),(0,5),(2,7),(2,8)]:
    from_order[i]^=1<<j
    from_order[j]^=1<<i
need(from_order==pi,'incoming forward-order Pi constructor')
records=[]
for name,rows,expected in [('Pi',pi,(361,6,29)),('Psi',psi,(1077,18,85))]:
    need(strong(rows),'strong source')
    n=len(rows); full=(1<<n)-1
    e,h=table(rows)
    s,_=table([full^(1<<v)^rows[v] for v in range(n)])
    states=[]
    for r in range(n):
        f=sum(e[m][r]*h[full^m] for m in range(1,full+1))
        states.append([h[full],s[full][r],h[full^(1<<r)],f])
    observed=(h[full],min(t[1] for t in states),min(t[3]-t[0] for t in states))
    need(observed==expected,('literal profile',name,observed))
    need(observed[1]>=2 and observed[2]>=2,'all-additions criterion')
    records.append({'name':name,'rows':rows,'strong':True,'states':states,
                    'H':observed[0],'minimum_starts':observed[1],'minimum_readiness_defect':observed[2]})
path='05-knowledge/results/clone_clocks_negative_collatz_20261113/CLOCKS.md'
revision='94cb3ea'
commit=subprocess.check_output(['git','rev-parse',revision],cwd=ROOT).decode().strip()
blob_id=subprocess.check_output(['git','rev-parse',revision+':'+path],cwd=ROOT).decode().strip()
blob=subprocess.check_output(['git','show',revision+':'+path],cwd=ROOT)
out={'status':'FINITE-EXACT literal barriers; all-order consequence uses ROOTED_HEREDITY H5',
     'universe':'All roots of the two explicit incoming strong base graphs',
     'git_source':{'commit':commit,'path':path,'blob_sha1':blob_id,'sha256':hashlib.sha256(blob).hexdigest()},
     'records':records,'proof_sha256':hashlib.sha256((HERE/'ROOTED_HEREDITY.md').read_bytes()).hexdigest()}
suffix='.optimized' if sys.flags.optimize else ''
(HERE/('incoming_clocks'+suffix+'.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps([(r['name'],r['H'],r['minimum_starts'],r['minimum_readiness_defect']) for r in records]))
