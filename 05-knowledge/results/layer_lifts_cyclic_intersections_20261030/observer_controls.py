#!/usr/bin/env python3
"""Full affine reachability and future-query separation with letters one/two."""
from collections import deque
from pathlib import Path
import argparse,json

checks=0
def need(ok,msg):
    global checks
    checks+=1
    if not ok:raise RuntimeError(msg)
def control(E):
    q=5**E
    def step(s,a):
        m,b=s;d=pow(2**a,-1,q);return 3*m*d%q,(3*b+1)*d%q
    def compose(s,t):
        m,b=s;n,c=t;return m*n%q,(m*c+b)%q
    paths={(1,0):()};queue=deque([(1,0)])
    while queue:
        s=queue.popleft()
        for a in (1,2):
            t=step(s,a)
            if t not in paths:paths[t]=paths[s]+(a,);queue.append(t)
    need(len(paths)==q*(q-q//5),'complete affine group via positive letters')
    for s in paths:
        m,b=s;mi=pow(m,-1,q);inverse=(mi,-b*mi%q)
        need(inverse in paths,'positive inverse word exists')
        actual=s
        for a in paths[inverse]:actual=step(actual,a)
        need(actual==(1,0),'future inverse word accepts its named source')
    G1=(3*pow(2,-1,q)%q,pow(2,-1,q))
    G2=(3*pow(4,-1,q)%q,pow(4,-1,q))
    inverse2=(pow(G2[0],-1,q),-G2[1]*pow(G2[0],-1,q)%q)
    need(compose(G1,inverse2)==(2,0),'dilation two from literal generators')
    A=compose((2,0),G1);B=compose(G1,(2,0));Bi=(pow(B[0],-1,q),-B[1]*pow(B[0],-1,q)%q)
    need(compose(A,Bi)==(1,pow(2,-1,q)),'unit translation commutator')
    pair_count=0
    if E==1:
        for s in paths:
            for t in paths:
                inv=(pow(s[0],-1,q),-s[1]*pow(s[0],-1,q)%q)
                actual=t
                for a in paths[inv]:actual=step(actual,a)
                need((actual==(1,0))==(s==t),'identity Hankel submatrix')
                pair_count+=1
    return dict(E=E,states=len(paths),max_positive_witness_length=max(map(len,paths.values())),
                exhaustive_distinguishing_pairs=pair_count)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path(__file__).with_name('observer.json'));a=p.parse_args()
    out=dict(levels=[control(E) for E in (1,2,3)]);out['checks']=checks
    out['scope']='finite controls for the all-E group-generation and observer lower-bound proof'
    report=json.dumps(out,indent=2,sort_keys=True)+'\n';a.output.write_text(report);print(report,end='')
