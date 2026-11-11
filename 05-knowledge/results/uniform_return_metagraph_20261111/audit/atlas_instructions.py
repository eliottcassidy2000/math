"""Independent named-graph check of all six sink/apex atlas instructions."""
from pathlib import Path
import json
import retained_reset as audit
from ready_rebase import frame


def source(k,tail):
    chain=tuple('P'+str(i) for i in range(k))+tuple('T'+str(i) for i in range(len(tail)))
    names=('r','a','b')+chain
    edges={('a','r'),('b','r'),('a','b')}
    edges.update(('r',x) for x in chain)
    edges.update((x,y) for i,x in enumerate(chain) for y in chain[i+1:])
    for x,letter in zip(chain,(3,)*k+tail):
        for j,port in enumerate(('a','b')):
            edges.add((x,port) if letter>>j&1 else (port,x))
    return names,edges


def main():
    recipes=[((0,2),1,2,('b',),9,-1,2,9,10),
             ((0,1),0,1,('a','b'),11,1,2,13,8),
             ((2,2,1),0,2,('T1','x'),34,5,2,41,24),
             ((2,2,1),1,2,('T0','x'),34,5,2,41,24),
             ((2,3,1),1,2,('T0','x'),34,5,2,41,24),
             ((2,1,2),0,1,('T2','a','b'),30,5,4,33,14)]
    results=[]
    for tail,drop,sinks,cut,A,B,offset,q,b in recipes:
        for k in range(1,7):
            old,old_edges=source(k,tail)
            before=frame(old,old_edges,'r')
            need_drop='T'+str(drop)
            names=tuple(v for v in old if v!=need_drop)
            edges={(u,v) for u,v in old_edges if u!=need_drop and v!=need_drop}
            kept=names; kept_edges=edges.copy()
            for new in ('x','y')[:sinks]:
                edges.update((v,new) for v in names); names+=(new,)
            edges.update(('z',v) if v in cut else (v,'z') for v in names)
            names+=('z',)
            after=frame(names,edges,'r')
            audit.need(before[0]==A*2**k+B and before[4]-before[0]==4,'source ray')
            audit.need(after[0]==3*before[0]+offset,'exact branch')
            audit.need(after[2]==after[5]==q and after[4]-after[0]==b,'root response')
            audit.need(audit.strong(audit.masks(names,edges)),'strong target')
            audit.need({(u,v) for u,v in edges if u in kept and v in kept}==kept_edges,'persistent old arcs')
            results.append(dict(tail=tail,drop=drop,k=k,H=after[0],root_frame=after))
    report=dict(status='PASS', universe='Six exact sink/apex instructions, k1..6; independent persistent named graph gauge',
                graph_cases=len(results),records=results)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('PASS: all36 independent atlas graphs, both ternary offsets and every retained old arc')


if __name__=='__main__':main()
