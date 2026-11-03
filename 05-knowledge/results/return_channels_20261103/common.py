"""Use the inherited, hash-guarded actual 65-symbol bank without reselection."""
from pathlib import Path
import sys,importlib.util
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'inverse65_zeckendorf_20261103'
sys.path.insert(0,str(OLD))
from bank import ATLAS,MATS,COSTS,need,step,vector,actual_target,closure,recover
from inverse import decode_state,physical_cone,source_certificate

def matrix(word):
    p=((1,0,0),(0,1,0),(0,0,1))
    for l in word:
        m=MATS[l];p=tuple(tuple(sum(p[i][k]*m[k][j] for k in range(3)) for j in range(3)) for i in range(3))
    return p

def old_graph_module():
    path=HERE.parent/'golden_hamilton_transfer_20261103/transfer.py'
    spec=importlib.util.spec_from_file_location('old_independent_graph',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def literal_frontier(rows):
    """Count directly from graph arcs using the inherited endpoint-forest DP."""
    from collections import Counter
    mod=old_graph_module();n=len(rows)
    if not n:return 1,0,1
    width=max([v-u for v in range(n) for u in range(v) if rows[v]>>u&1]+[1])
    states=Counter({():1});peak=1
    for v,row in enumerate(rows):
        mask=sum(1<<(v-u-1) for u in range(v) if row>>u&1)
        nxt=Counter()
        for state,count in states.items():
            for target,weight in mod.step(state,mask,width).items():nxt[target]+=count*weight
        states=nxt;peak=max(peak,len(states))
    return sum(c for s,c in states.items() if len(s)==1),width,peak
