"""Exploratory ready-return probes; not a completeness certificate."""
from pathlib import Path
import sys,json
from itertools import combinations
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'uniform_return_metagraph_20261111'/'reset'))
import classify_modules as m
from atlas_two_birth import construct
sys.path.insert(0,str(RESULTS/'ready_return_pentagon_20261103'))
import ready_inverse as inv

P=(34,244,249,243,225,192,129,1)
R=m.inherited.birth_beaten_only_by(P,0)

def embeddings(A,B,fixed):
    n=len(A); assignment=dict(fixed);used=set(fixed.values())
    def rec():
        if len(assignment)==n:
            yield tuple(assignment[v] for v in range(n));return
        todo=[v for v in range(n) if v not in assignment]
        domains={v:[w for w in range(len(B)) if w not in used and all(
            (A[v]>>u&1)==(B[w]>>z&1) for u,z in assignment.items())] for v in todo}
        v=min(todo,key=lambda x:len(domains[x]))
        for w in domains[v]:
            assignment[v]=w;used.add(w);yield from rec();used.remove(w);del assignment[v]
    yield from rec()

def main():
    S=construct((1,1,2),0,4,14,1)[-1]
    retained=[]
    for loss in range(len(S)):
        for drop in combinations(range(1,len(S)),loss):
            kept=[v for v in range(len(S)) if v not in drop]
            A=m.inherited.induced(S,kept)
            for mapping in embeddings(A,R,{0:5}):
                retained.append(dict(drop=drop,kept=kept,mapping=mapping))
        if retained:break
    print('RETENTION',len(S),len(R),retained,flush=True)
    hits=[]
    for n in range(3,12):
        for j in range(inv.class_counts(n)[n]):
            W=inv.unrank(n,j);h,s,q,F=inv.rooted_data(W,0)
            if h+F-1==665:
                T=m.inherited.birth_beaten_only_by(W,0)
                for v in range(1,len(W)):
                    co=m.two_birth.module_runs(m.inherited.induced(T,[v]+[w for w in range(len(T)) if w!=v]))
                    if co[2]==297 and not any(co[3:]):
                        hits.append(dict(n=n,rank=j,parent=W,module=v,runs=co))
                        print('UNIFORM665',hits[-1],flush=True)
            inv.tables.cache_clear()
        print('DONE',n,'hits',len(hits),flush=True)
    (HERE/'probe_ready.json').write_text(json.dumps(dict(retention=retained,next_hits=hits),indent=2)+'\n')

if __name__=='__main__':main()
