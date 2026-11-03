"""Pinned 65-symbol seam bank; words are in PROCESSING order."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
ATLAS_PATH=HERE.parent/'golden_seam_recurrence_20261103/atom_atlas.json'
pins=json.loads((HERE/'source_pins.json').read_text())['sha256_lf']
atlas_key=ATLAS_PATH.relative_to(HERE.parents[2]).as_posix()
if hashlib.sha256(ATLAS_PATH.read_bytes().replace(b'\r\n',b'\n')).hexdigest()!=pins[atlas_key]:
    raise RuntimeError('Inherited atom atlas changed; audit the bank before repinning')
ATLAS=json.loads(ATLAS_PATH.read_text())
SIX=('s14','s268','s330','s248','s261','s85','s262','s453')
LABELS=tuple(k for k,v in ATLAS.items() if len(v['source_rows'])<=5)+SIX
MATS={k:tuple(map(tuple,ATLAS[k]['matrix'])) for k in LABELS}
COSTS={k:len(ATLAS[k]['source_rows']) for k in LABELS}

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

need(len(LABELS)==65 and all(l in LABELS for l in ('q0','q1','q2','q4','q5')),'fixed bank with all Q4 charts')

def step(v,m):return tuple(sum(v[i]*m[i][j] for i in range(3)) for j in range(3))

def vector(word):
    v=(1,0,0)
    for l in word:v=step(v,MATS[l])
    return v

def closure(cap):
    states=[(1,0,0)];ids={states[0]:0};parents=[None]
    for i,v in enumerate(states):
        for l,m in MATS.items():
            w=step(v,m)
            if w[0]>cap or w in ids:continue
            ids[w]=len(states);states.append(w);parents.append((i,l))
    return states,parents,ids

def recover(i,parents):
    w=[]
    while parents[i] is not None:i,l=parents[i];w.append(l)
    return tuple(reversed(w))

def actual_target(word):
    rows=[]
    for l in reversed(word):
        source=ATLAS[l]['source_rows'];n=len(rows);m=len(source)
        rows=[r|(((1<<m)-1)<<n) for r in rows]+[r<<n for r in source]
    n=len(rows)
    for i in range(n):
        for j in range(i+2,n):rows[i]^=1<<j;rows[j]^=1<<i
    return tuple(sum(1<<(n-1-j) for j in range(n) if rows[n-1-i]>>j&1) for i in range(n))
