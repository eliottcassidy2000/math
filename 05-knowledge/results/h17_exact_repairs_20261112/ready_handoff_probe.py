"""Finite H29 retained-star handoff and hostile block enlargement."""
from pathlib import Path
from itertools import combinations
import importlib.util
import json

HERE=Path(__file__).resolve().parent
p=HERE.parent/'eight_nine_twentyseven_20261112/constructors/escape.py'
spec=importlib.util.spec_from_file_location('handoff_graphs',p)
g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)


def need(ok,message):
    if not ok:raise RuntimeError(message)


def inflate(b):
    base=(22,8,2,21,6); blocks=[];n=0
    for size in (1,1,b,1,1):
        blocks.append(list(range(n,n+size)));n+=size
    rows=[0]*n
    for i,A in enumerate(blocks):
        for a in A:
            for z in A:
                if a<z:rows[a]|=1<<z
        for j,B in enumerate(blocks):
            if base[i]>>j&1:
                for a in A:
                    for z in B:rows[a]|=1<<z
    return tuple(rows),blocks


def delete(t,vertex):
    names=[i for i in range(len(t)) if i!=vertex]
    return tuple(sum(1<<j for j,b in enumerate(names) if t[a]>>b&1) for a in names),names


def main():
    source=(22,8,2,21,6)
    old,r=g.child(source,0,True);old,r=g.child(old,r)
    need(old==(54,8,2,21,6,94,31) and r==5,'named H29 source')
    need(g.root_state(old,r)==(29,9,9,29),'raw root state')
    edges=[(i,j) for i in range(7) for j in range(i+1,7) if r not in (i,j)]
    attempts=0;hits=[]
    for radius in (1,2):
        for choices in combinations(edges,radius):
            t=old
            for i,j in choices:t=g.flip(t,i,j)
            z=g.root_state(t,r);attempts+=1
            if z[0]==29 and z[2]-z[1]==1:hits.append(choices)
    need(attempts==120 and not hits,'radius two boundary')
    repaired=old
    for edge in ((0,2),(1,3),(4,6)):repaired=g.flip(repaired,*edge)
    need(repaired==(50,0,3,23,70,94,15),'three-arc repair')
    need(repaired[r]==old[r] and all((repaired[i]>>r&1)==(old[i]>>r&1) for i in range(7)),'root star retained')
    need(g.root_state(repaired,r)==(29,8,9,30),'unit-ready exact target')
    small,names=delete(repaired,1);smallroot=names.index(r)
    need(g.strong(small) and g.root_state(small,smallroot)==(29,8,9,30),'neutral sink removal')
    nxt,nr=g.child(small,smallroot,True);nxt,nr=g.child(nxt,nr)
    need(g.root_state(nxt,nr)==(89,29,29,89),'one more exact target')
    trials=[]
    for b in range(1,7):
        t,blocks=inflate(b);initial=g.root_state(t,0)
        raw,root=g.child(t,0,True);raw,root=g.child(raw,root)
        modified=raw
        for j in blocks[2]:modified=g.flip(modified,0,j)
        modified=g.flip(modified,1,blocks[3][0])
        modified=g.flip(modified,blocks[4][0],len(raw)-1)
        z=g.root_state(raw,root);w=g.root_state(modified,root)
        need(w[2]-w[1]==1 and modified[1]==0,'enlarged readiness and sink')
        trials.append(dict(b=b,initial=initial,raw=z,modified=w,target_preserved=z[0]==w[0]))
    need(trials[0]['target_preserved'] and trials[1]['raw'][0]==53 and trials[1]['modified'][0]==59,'first failed enlargement')
    receipt=dict(status='FINITE-EXACT: fixed-source local minimum and six declared enlargement controls',
                 source_rows=old,root=r,nonroot_radius_one_two_attempts=attempts,ready_hits=hits,
                 repair_pairs=[[0,2],[1,3],[4,6]],repaired_rows=repaired,repaired_state=[29,8,9,30],
                 deleted_physical_vertex=1,six_vertex_rows=small,six_vertex_name_map=names,
                 next_rows=nxt,next_root=nr,next_state=[89,29,29,89],enlargements=trials)
    (HERE/'ready_handoff_probe.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','minimal_arc_edits':3,'next_count':89,'inflation_hostile':[53,59]}))


if __name__=='__main__':main()
