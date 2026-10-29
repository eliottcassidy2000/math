"""Coarsest level-preserving equitable quotient of the folded-path DAG."""
from collections import Counter,defaultdict
from pathlib import Path
from itertools import product
import sys,json
import probes as p
core=p.core;need=p.need;r=p.r
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'session90_even_graphs_odd_descent_20261029/oddness'
sys.path.insert(0,str(OLD))
import check_centralizer_path_dp as old


def minimize(rows,tau,selectors=()):
    root=next(i for i in range(len(tau)) if tau[i]==i)
    pairs=[(i,tau[i]) for i in range(len(tau)) if i<tau[i]];n=len(pairs)
    vertices=[v for pair in pairs for v in pair];owner={v:i for i,vs in enumerate(pairs) for v in vs}
    if not n:return dict(count=1,states=0,cells=0)
    states=[(S,v) for S in range(1,1<<n) for v in vertices if S>>owner[v]&1]
    successors={s:[(s[0]|1<<owner[w],w) for w in vertices if not s[0]>>owner[w]&1 and core.edge(rows,s[1],w)] for s in states}
    def terminal(s):return int(s[0]==(1<<n)-1 and core.edge(rows,s[1],root))
    def obs(s):return s[0].bit_count(),terminal(s),tuple(f(s) for f in selectors)
    colours={};table={}
    for s in states:colours[s]=table.setdefault(obs(s),len(table))
    history=[len(table)]
    while True:
        table={};new={}
        for s in states:
            sig=colours[s],tuple(sorted(Counter(colours[t] for t in successors[s]).items()))
            new[s]=table.setdefault(sig,len(table))
        if all(new[s]==colours[s] for s in states):break
        stable=len(table)==len(set(colours.values()));colours=new
        if stable:break
        history.append(len(table))
    groups=defaultdict(list)
    for s,c in colours.items():groups[c].append(s)
    out={c:Counter(colours[t] for t in successors[ss[0]]) for c,ss in groups.items()}
    for c,ss in groups.items():
        need(all(Counter(colours[t] for t in successors[s])==out[c] and obs(s)==obs(ss[0]) for s in ss),'equitable cell')
    # Independent bottom-up signatures on this acyclic graph certify the same partition.
    bottom={};intern={}
    for s in sorted(states,key=lambda s:-s[0].bit_count()):
        sig=obs(s),tuple(sorted(Counter(bottom[t] for t in successors[s]).items()))
        bottom[s]=intern.setdefault(sig,len(intern))
    correspondence={(colours[s],bottom[s]) for s in states}
    need(len(correspondence)==len(groups)==len(intern),'bottom-up minimal partition equality')
    mass=Counter()
    for c,ss in groups.items():
        if ss[0][0].bit_count()==1:mass[c]=len(ss)
    for c in sorted(groups,key=lambda c:groups[c][0][0].bit_count()):
        for d,m in out[c].items():mass[d]+=mass[c]*m
    labelled=Counter({s:1 for s in states if s[0].bit_count()==1})
    for s in sorted(states,key=lambda s:s[0].bit_count()):
        for t in successors[s]:labelled[t]+=labelled[s]
    need(all(mass[c]==sum(labelled[s] for s in ss) for c,ss in groups.items()),'every equitable cell has its full labelled mass')
    count=sum(mass[c] for c,ss in groups.items() if terminal(ss[0]))
    need(count==old.folded_paths(rows,tau),'independent labelled folded recurrence')
    return dict(count=count,states=len(states),cells=len(groups),refinement_counts=history,
                cell_sizes=dict(sorted(Counter(map(len,groups.values())).items())),
                layers=dict(sorted(Counter(ss[0][0].bit_count() for ss in groups.values()).items())),
                weighted_entries=sum(len(row) for row in out.values()),
                classes=[dict(states=ss,out=sorted(out[c].items()),terminal=terminal(ss[0])) for c,ss in sorted(groups.items())])


def main():
    rows,_=core.substitute(p.cycle3(),(p.cycle3(),)*3)
    tau=tuple(3*((-(v//3))%3)+(-(v%3))%3 for v in range(9))
    G=r.oracle.isomorphisms(rows,rows,True)
    C=[g for g in G if all(g[tau[i]]==tau[g[i]] for i in range(9))]
    baseline=old.compressed(rows,tau,C)
    plain=minimize(rows,tau)
    one=minimize(rows,tau,(lambda s:s[1]==3,))
    vertex=minimize(rows,tau,(lambda s:s[1],))
    full=minimize(rows,tau,(lambda s:s[1],lambda s:s[0]))
    cases=0;totals=Counter()
    for n in range(4):
        for x in range(1<<(n*n)):
            t,a=old.lift(x,n);z=minimize(t,a);g=old.normalized_centralizer(t,n);orbit=old.compressed(t,a,g)
            need(z['cells']<=orbit['orbits'],'minimum refines no more than group orbits')
            cases+=1;totals.update(labelled=z['states'],equitable=z['cells'],orbit=orbit['orbits'])
    result=dict(status='FINITE-EXACT',baseline=baseline,minimal=plain,selected_vertex=one,
                all_last_vertices=vertex,full_state=full,small_lifts=dict(cases=cases,totals=dict(totals)))
    (HERE/'folded-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(baseline=baseline,minimal={k:v for k,v in plain.items() if k!='classes'},
                         selector_cells=[one['cells'],vertex['cells'],full['cells']],small_lifts=result['small_lifts'])));return result


if __name__=='__main__':main()
