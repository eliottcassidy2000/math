"""Frozen original native data; no remote or current-file dependencies."""
from pathlib import Path
import subprocess, json, csv, gzip, sys, hashlib
from array import array
from itertools import combinations
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
REV='cf87f98a581089b3f84d3fee16324195ab976e1e'
INPUTS={}
def read(path):
    raw=subprocess.check_output(['git','show',REV+':05-knowledge/results/'+path],cwd=ROOT)
    INPUTS[path]=hashlib.sha256(raw).hexdigest()
    return raw
def inventory():
    return sorted(json.loads(read('session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json'))['records'],key=lambda r:r['source'])
def data():
    nodes=list(csv.DictReader(read('metagraph_order8_20261023/merged_nodes.tsv').decode().splitlines(),delimiter='\t'))
    owner=array('H');owner.frombytes(gzip.decompress(read('metagraph_order8_20261023/owner_map.u16.gz')))
    if sys.byteorder!='little':owner.byteswap()
    merged={int(n[k]):int(n['merged_id']) for n in nodes for k in ('class_a','class_b')}
    d=json.loads(read('metagraph_order8_electrical_20261023/receipt.json'))['order8']
    return nodes,owner,merged,d['higher_layers'][-1]['lifts'][0],d['radical_lifts_mod4']
PAIRS=list(combinations(range(8),2))
CHORDS=[(i,j) for i,j in PAIRS if j>=i+2]
FULL=(1<<21)-1
def matrix(code):
    B=[[0]*8 for _ in range(8)]
    for k,(i,j) in enumerate(PAIRS):B[i][j]=(int(code)>>(27-k))&1;B[j][i]=1-B[i][j]
    return B
def paths(B):
    def dfs(p,left):
        if not left:yield p;return
        for v in left:
            if not p or B[p[-1]][v]:yield from dfs(p+(v,),left-{v})
    return dfs((),set(range(8)))
def word(B,p):return sum(B[p[i]][p[j]]<<k for k,(i,j) in enumerate(CHORDS))
def mirror(w):
    idx={p:i for i,p in enumerate(CHORDS)}
    return sum(((w>>i)&1)<<idx[7-v,7-u] for i,(u,v) in enumerate(CHORDS))
def orbit(w):return sorted({w,FULL^w,mirror(w),FULL^mirror(w)})
def label(X,Y,s):return sum((Y[j][s]%2)<<j for j in range(4)),X[s]%128
def need(v,msg):
    if not v:raise RuntimeError(msg)
