"""Standalone finite controls for cell-separator localization of original repair.

No mathematical producer is imported. The whole-instance reference enumerates
triple hitting choices; the localized route uses Tarjan blocks, two-state
messages, and exact bounded weighted graph-capacity DP. This implementation is
an exact finite verifier, not an implementation timing claim about the theorem.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path
import sys

sys.stdout.reconfigure(newline="\n")
GATES = Counter()

def check(condition, name):
    GATES[name] += 1
    if not condition:
        raise RuntimeError('always-active gate failed: '+name)

SLOPES = (1, -1, 2)

def lines(points, slopes=SLOPES):
    by = defaultdict(list)
    for i, (r,c) in enumerate(points):
        for s in slopes:
            by[s,c-s*r].append(i)
    keys = sorted(k for k,v in by.items() if len(v)>2)
    return keys, [tuple(by[k]) for k in keys]

def exact(ls, deleted=(), retained=()):
    dm=sum(1<<p for p in deleted)
    rm=sum(1<<p for p in retained)
    if dm & rm:
        raise ValueError('contradictory forced state')
    clauses=set()
    for line in ls:
        for tr in combinations(line,3):
            bits=sum(1<<p for p in tr)
            if not bits & dm:
                clauses.add(bits & ~rm)
    @lru_cache(None)
    def solve(cs):
        if not cs: return 0
        if cs[0]==0: return 10**9
        first=min(cs,key=int.bit_count)
        choices=[1<<p for p in range(first.bit_length()) if first>>p&1]
        return 1+min(solve(tuple(t for t in cs if not t&b)) for b in choices)
    return len(deleted)+solve(tuple(sorted(clauses)))

def decompose(ls):
    active=sorted({p for l in ls for p in l})
    n=len(ls)
    pointnode={p:n+i for i,p in enumerate(active)}
    nodepoint={v:k for k,v in pointnode.items()}
    adj=defaultdict(list)
    for a,line in enumerate(ls):
        for p in line:
            b=pointnode[p]
            adj[a].append(b); adj[b].append(a)
    times={}; low={}; stack=[]; blocks=[]
    def dfs(a,parent):
        times[a]=low[a]=len(times)
        for b in adj[a]:
            if b==parent: continue
            if b not in times:
                stack.append((a,b)); dfs(b,a); low[a]=min(low[a],low[b])
                if low[b]>=times[a]:
                    block=[]
                    while True:
                        e=stack.pop();block.append(e)
                        if e==(a,b): break
                    blocks.append(block)
            elif times[b]<times[a]:
                stack.append((a,b));low[a]=min(low[a],times[b])
    for a in sorted(adj):
        if a not in times: dfs(a,None)
    parent=list(range(len(blocks)))
    def find(a):
        while parent[a]!=a:
            parent[a]=parent[parent[a]]; a=parent[a]
        return a
    memberships=defaultdict(list)
    for i,edges in enumerate(blocks):
        for node in {v for e in edges for v in e}: memberships[node].append(i)
    for node,bs in memberships.items():
        if node<n:
            for b in bs[1:]: parent[find(b)]=find(bs[0])
    merged=defaultdict(set)
    for i,edges in enumerate(blocks):
        merged[find(i)].update(v for e in edges for v in e)
    bags=[]
    for nodes in merged.values():
        bags.append({'lines':tuple(sorted(v for v in nodes if v<n)),
                     'cells':tuple(sorted(nodepoint[v] for v in nodes if v>=n))})
    bags.sort(key=lambda b:b['lines'])
    occ=defaultdict(list)
    for b,bag in enumerate(bags):
        for p in bag['cells']: occ[p].append(b)
    separators={p:tuple(bs) for p,bs in occ.items() if len(bs)>1}
    return bags,separators

def localized(ls, first_root=None):
    bags,seps=decompose(ls)
    line_sets=list(map(set,ls))
    calls=0; negative=[]; local_h=[]; fixed_high_parent=0
    for b,bag in enumerate(bags):
        local_h.append(sum(sum(p in line_sets[a] for a in bag['lines'])>=3 for p in bag['cells']))
    @lru_cache(None)
    def cell_msg(p,parent_b,s):
        return 1-s+sum(bag_msg(b,p,s) for b in seps[p] if b!=parent_b)
    @lru_cache(None)
    def bag_msg(b,parent_p,parent_s):
        nonlocal calls, fixed_high_parent
        bag=bags[b]; aa=bag['lines']
        if parent_p is not None and sum(parent_p in line_sets[a] for a in aa)>=3:
            fixed_high_parent+=1
        cap=tuple(2-int(parent_p in line_sets[a]) for a in aa) if parent_s else (2,)*len(aa)
        ordinary=[]; high=[]; baseline=0
        for p in bag['cells']:
            if p==parent_p: continue
            c0,c1=(cell_msg(p,b,0),cell_msg(p,b,1)) if p in seps else (1,0)
            weight=c0-c1;baseline+=c0
            if weight<0: negative.append((b,p,weight))
            inc=tuple(i for i,a in enumerate(aa) if p in line_sets[a])
            (high if len(inc)>=3 else ordinary).append((p,inc,weight))
        best=10**9
        for states in product((0,1),repeat=len(high)):
            residual=list(cap); reward=0
            for s,(_,inc,w) in zip(states,high):
                if s:
                    reward+=w
                    for i in inc: residual[i]-=1
            if min(residual,default=0)<0: continue
            calls+=1
            @lru_cache(None)
            def matching(j,capacity):
                if j==len(ordinary):return 0
                _,inc,w=ordinary[j]
                value=matching(j+1,capacity)
                if all(capacity[i]>0 for i in inc):
                    nxt=list(capacity)
                    for i in inc:nxt[i]-=1
                    value=max(value,w+matching(j+1,tuple(nxt)))
                return value
            best=min(best,baseline-reward-matching(0,tuple(residual)))
        return best
    seen=set();roots=[]
    def visit(b):
        if b in seen:return
        seen.add(b)
        for p in bags[b]['cells']:
            for child in seps.get(p,()): visit(child)
    order=list(range(len(bags)))
    if first_root is not None:
        order.remove(first_root);order.insert(0,first_root)
    for b in order:
        if b not in seen:roots.append(b);visit(b)
    result=sum(bag_msg(b,None,0) for b in roots)
    return {'tau':result,'bags':len(bags),'local_h':local_h,'calls':calls,
            'bound':2*sum(2**h for h in local_h),'negative':negative,'separators':seps,
            'fixed_high_parent':fixed_high_parent}

def fill_private(centres, extra_shared=(), slopes=SLOPES):
    selected={(s,c-s*r) for r,c in centres for s in slopes}
    points=list(centres)+list(extra_shared)
    usedrows={r for r,c in points};usedcols={c for r,c in points}
    used={s:{c-s*r for r,c in points} for s in slopes}
    for slope,label in sorted(selected):
        need=3-sum(c-slope*r==label for r,c in points)
        if need<0:raise ValueError('selected line already too long')
        for _ in range(need):
            for r in range(100000):
                c=slope*r+label
                if r in usedrows or c in usedcols:continue
                if any((s,c-s*r) in selected or c-s*r in used[s] for s in slopes if s!=slope):continue
                points.append((r,c));usedrows.add(r);usedcols.add(c)
                for s in slopes:used[s].add(c-s*r)
                break
            else:raise RuntimeError('private-point search limit')
    rr=min(r for r,c in points); cc=min(c for r,c in points)
    return [(r-rr,c-cc) for r,c in points]

def chain(k):
    return fill_private([(2*i,2*i+2*(i//2)) for i in range(k)])

def synergy(slopes=SLOPES):
    return fill_private([(0,0),(6,0)],[(3,3),(4,-4),(-6,-12)],slopes)

def adjacency(ls):
    graph=defaultdict(set)
    for a,line in enumerate(ls):
        for p in line:
            graph['l',a].add(('p',p));graph['p',p].add(('l',a))
    return graph

def components(graph, omitted=frozenset()):
    remaining=set(graph)-set(omitted);result=[]
    while remaining:
        seed=min(remaining);seen={seed};todo=[seed];remaining.remove(seed)
        while todo:
            a=todo.pop()
            for b in graph[a]:
                if b in remaining:
                    remaining.remove(b);seen.add(b);todo.append(b)
        result.append(seen)
    return result

def validate(ls, reroot=False, max_directions=None):
    check(all(len(line)>=3 and len(set(line))==len(line) for line in ls),'original overfull line types')
    bagdata,seps=decompose(ls);graph=adjacency(ls)
    memberships=Counter(a for bag in bagdata for a in bag['lines'])
    check(memberships==Counter(range(len(ls))),'complete unique line ownership')
    for bag in bagdata:
        check(set(bag['cells'])=={p for a in bag['lines'] for p in ls[a]},'bag contains exactly owned-line cells')
    pointbags=defaultdict(set)
    quotient=defaultdict(set)
    for b,bag in enumerate(bagdata):
        quotient['b',b]
        for p in bag['cells']:pointbags[p].add(b)
    for p,bb in pointbags.items():
        check((len(bb)>1)==(p in seps),'separator multiplicity')
        if len(bb)>1:
            check(bb==set(seps[p]),'separator membership')
            for b in bb:
                quotient['b',b].add(('p',p));quotient['p',p].add(('b',b))
    qcomponents=components(quotient)
    edgecount=sum(map(len,quotient.values()))//2
    check(edgecount==len(quotient)-len(qcomponents),'bag-cell quotient is forest')
    check(len(qcomponents)==len(components(graph)),'component count preserved')
    base_components=len(components(graph))
    articulations={p for kind,p in graph if kind=='p' and len(components(graph,{('p',p)}))>base_components}
    check(articulations==set(seps),'separator cells equal brute-force articulations')
    inc=Counter(p for l in ls for p in l)
    if max_directions is not None:
        check(max(inc.values(),default=0)<=max_directions,'geometric direction incidence bound')
    if max_directions==3:
        local=Counter(p for bag in bagdata for p in bag['cells'] if sum(p in ls[a] for a in bag['lines'])>=3)
        expected=Counter({p:1 for p,d in inc.items() if d==3 and p not in articulations})
        check(local==expected,'three-direction nonarticulation corollary')
    reference=exact(ls);result=localized(ls)
    check(result['tau']==reference,'localized equals independent triple hitting')
    check(result['calls']<=result['bound'],'declared graph oracle bound')
    if reroot:
        for root in range(len(bagdata)):
            alternate=localized(ls,root)
            check(alternate['tau']==reference,'arbitrary bag root same optimum')
            check(alternate['calls']<=alternate['bound'],'arbitrary root call bound')
    return result,inc

def boards(n):
    residual=[2]*n;rows=[]
    def visit(r):
        if r==n:
            if all(v==0 for v in residual):yield tuple(rows)
            return
        for a,b in combinations(range(n),2):
            if not residual[a] or not residual[b]:continue
            residual[a]-=1;residual[b]-=1;rows.append((a,b))
            if max(residual,default=0)<=n-r-1:yield from visit(r+1)
            rows.pop();residual[a]+=1;residual[b]+=1
    yield from visit(0)

def row_points(rows):
    return [(r,c) for r,pair in enumerate(rows) for c in pair]

def abstract_bank(v):
    """All shared 2/3-cells, no repeated line pair, then private triple fill."""
    candidates=list(combinations(range(v),2))+list(combinations(range(v),3))
    for bits in range(1<<len(candidates)):
        chosen=[c for i,c in enumerate(candidates) if bits>>i&1]
        degree=Counter(a for c in chosen for a in c)
        pairs=Counter(pair for c in chosen for pair in combinations(c,2))
        if max(degree.values(),default=0)>3 or max(pairs.values(),default=0)>1:continue
        ls=[[] for _ in range(v)]
        for p,c in enumerate(chosen):
            for a in c:ls[a].append(p)
        nxt=len(chosen)
        for a in range(v):
            while len(ls[a])<3:ls[a].append(nxt);nxt+=1
        yield list(map(tuple,ls))

def naive_quotient_cycle(ls):
    graph=adjacency(ls);initial=len(components(graph))
    arts={node for node in graph if node[0]=='p' and len(components(graph,{node}))>initial}
    parts=components(graph,arts);quotient=defaultdict(set)
    for j,part in enumerate(parts):
        quotient['c',j]
        for p in arts:
            if any(x in graph[p] for x in part):
                quotient['c',j].add(p);quotient[p].add(('c',j))
    edges=sum(map(len,quotient.values()))//2
    return {'articulation_cells':sorted(p[1] for p in arts),'vertices':len(quotient),
            'edges':edges,'components':len(components(quotient)),
            'cycle_rank':edges-len(quotient)+len(components(quotient))}

def compact(result):
    return {key:result[key] for key in ('tau','bags','local_h','calls','bound','fixed_high_parent')} | {
        'negative_rewards':result['negative'],'separators':result['separators']}

def main():
    certificate={'scope':'Finite exact original-cell checks; no probabilistic or extremal asymptotic claim.',
                 'producer_routes':['whole-instance triple hitting','block-cut plus weighted capacity DP']}
    censuses={};total=0
    for n in range(2,6):
        count=0;high=0;hist=Counter();calls=0;neg=0
        for rows in boards(n):
            points=row_points(rows);_,ls=lines(points)
            check(len(points)==2*n and len(set(points))==2*n,'simple complete board')
            check(set(Counter(c for r,c in points).values())=={2},'column degree two')
            result,inc=validate(ls,max_directions=3)
            count+=1;calls+=result['calls'];neg+=bool(result['negative'])
            high+=any(d>=3 for d in inc.values())
            hist[max(result['local_h'],default=0)]+=1
        check(count=={2:1,3:6,4:90,5:2040}[n],'complete two-regular universe count')
        censuses[n]={'boards':count,'with_global_exception':high,'max_local_h_histogram':dict(hist),
                     'total_oracle_calls':calls,'negative_reward_boards':neg}
        total+=count
    check(total==2137,'complete board total')
    certificate['complete_boards']=censuses
    abstract={}
    for v in range(5):
        count=0;local_hist=Counter()
        for ls in abstract_bank(v):
            result,_=validate(ls,reroot=True)
            local_hist[max(result['local_h'],default=0)]+=1;count+=1
        check(count=={0:1,1:1,2:2,3:9,4:96}[v],'complete abstract linear triple universe')
        abstract[v]={'systems':count,'max_local_h_histogram':dict(local_hist)}
    certificate['abstract_systems']=abstract
    control_rows={
        'triangle':((0,3),(1,2),(0,4),(3,4),(1,2)),
        'five_cycle':((0,3),(1,2),(0,4),(1,4),(3,5),(2,5)),
        'concurrent_triple':((2,4),(0,1),(2,3),(3,4),(0,1)),
        'triangle_pendant':((0,1),(0,2),(3,4),(1,3),(2,4)),
        'larger_line':((0,1),(0,1),(2,3),(2,4),(3,4)),
        'first_saturated_core':((0,1),(0,1),(2,3),(4,5),(3,5),(2,4))}
    controls={}
    for name,rows in control_rows.items():
        points=row_points(rows);modes={}
        for slopes in (SLOPES,(1,-1,2,-2)):
            labels,ls=lines(points,slopes);result,_=validate(ls,reroot=True,max_directions=len(slopes))
            modes[str(slopes)]={'lines':list(zip(labels,ls)),**compact(result)}
        controls[name]={'rows':rows,'modes':modes}
    check(controls['triangle']['modes'][str(SLOPES)]['tau']==2,'inherited triangle integral control')
    check(controls['five_cycle']['modes'][str(SLOPES)]['tau']==3,'inherited odd five-cycle control')
    check(controls['concurrent_triple']['modes'][str(SLOPES)]['tau']==1,'literal concurrency control')
    check(controls['triangle_pendant']['modes'][str(SLOPES)]['tau']==2,'attachment parity control')
    core=controls['first_saturated_core']['modes'][str(SLOPES)]
    check(core['tau']==4 and core['bags']==1 and core['local_h']==[1],'dimension-six genuine local core')
    check(all(d['max_local_h_histogram']=={0:d['boards']} for d in censuses.values()),'no smaller saturated local core')
    for j,rows in enumerate(boards(6),1):
        _,ls=lines(row_points(rows));result=localized(ls)
        if any(result['local_h']):
            check(j==14 and rows==control_rows['first_saturated_core'],'bounded lexicographic core discovery')
            break
    else:raise RuntimeError('missing expected early dimension-six core')
    certificate['saturated_core_discovery']={'dimension':6,'boards_examined':j,'exhaustive_dimension_six':False,
        'smaller_dimensions_complete':True,'row_pairs':rows}
    certificate['inherited_controls']=controls
    family=[]
    for k in range(1,9):
        points=chain(k);labels,ls=lines(points);result,inc=validate(ls,reroot=True,max_directions=3)
        graph=adjacency(ls)
        check(len(points)==4*k+3 and len(ls)==2*k+1,'chain geometric cardinalities')
        check(all(len(l)==3 for l in ls),'chain all intended triples')
        check(len(graph)==6*k+4 and sum(map(len,graph.values()))//2==6*k+3,'chain tree cardinalities')
        check(len(components(graph))==1,'chain connected incidence')
        check({p for p,d in inc.items() if d>=3}==set(range(k)),'chain exact exceptional centres')
        check(not any(result['local_h']),'chain no local exceptional cells')
        check(result['tau']==k,'chain exact unbounded-family formula')
        check(result['calls']==4*k+1 and result['bound']==4*k+2,'chain linear graph-call count')
        check(len({r for r,c in points})==len(points)==len({c for r,c in points}),'chain all occupied rows and columns one')
        all_labels=Counter((s,c-s*r) for r,c in points for s in SLOPES)
        check(set(all_labels.values())<={1,3},'chain nonselected lines singleton')
        if k==1:check(any(w==-1 for b,p,w in result['negative']),'negative continuation reward control')
        family.append({'k':k,'points':points,'lines':list(zip(labels,ls)),**compact(result)})
    certificate['connected_chain_family']=family
    points=synergy();labels,ls=lines(points);result,inc=validate(ls,reroot=True,max_directions=3)
    state_cost={''.join(map(str,states)):exact(ls,deleted=[i for i,s in enumerate(states) if s],
                  retained=[i for i,s in enumerate(states) if not s]) for states in product((0,1),repeat=2)}
    check(state_cost=={'00':3,'01':4,'10':4,'11':2},'inseparable two-cell synergy table')
    check(result['bags']==1 and result['local_h']==[2],'theta core retains both exceptions')
    check(len(points)==11 and len(ls)==6 and {p for p,d in inc.items() if d==3}=={0,1},'theta original geometry')
    certificate['geometric_joint_choice']={'points':points,'lines':list(zip(labels,ls)),
                                          'forced_deletion_cost':state_cost,**compact(result)}
    points=fill_private([(0,0),(6,6),(2,-2)])
    labels,ls=lines(points);result,_=validate(ls,reroot=True,max_directions=3)
    naive=naive_quotient_cycle(ls)
    check(naive['cycle_rank']==1 and naive['articulation_cells']==[0,1,2],'naive articulation quotient is cyclic')
    check(result['bags']==4 and not any(result['local_h']) and result['tau']==3,'whole cyclic bag with cell arms control')
    certificate['naive_component_hostile']={'points':points,'lines':list(zip(labels,ls)),
                                           'naive_quotient':naive,**compact(result)}
    points=synergy((1,-1,2,-2));labels,ls=lines(points,(1,-1,2,-2))
    result,inc=validate(ls,reroot=True,max_directions=4)
    check({p for p,d in inc.items() if d>=3}=={0,1}==set(result['separators']),'four-direction exceptional articulations')
    check(max(result['local_h'])==2,'four-direction articulations still locally exceptional')
    fixed=0
    for b in range(result['bags']):fixed+=localized(ls,b)['fixed_high_parent']
    check(fixed>0,'high-incidence fixed parent state exercised')
    certificate['four_direction_boundary']={'points':points,'lines':list(zip(labels,ls)),
                                           'all_root_fixed_high_states':fixed,**compact(result)}
    certificate['gate_counts']=dict(sorted(GATES.items()))
    certificate['total_gates']=sum(GATES.values())
    stem=Path(__file__).stem
    directory=Path(__file__).resolve().parent
    if directory.name=='04-computation':directory=directory.parent/'05-knowledge'/'results'
    target=directory/(stem+'_certificate.json')
    target.write_bytes((json.dumps(certificate,sort_keys=True,indent=2)+'\n').encode('utf-8'))
    print('Continuing13 cell-separator repair: PASS')
    print('Complete simple two-regular boards n=2..5: '+str(total))
    print('Abstract complete linear triple-incidence systems v=0..4: '+str(sum(d['systems'] for d in abstract.values())))
    print('Smallest saturated-board local-core dimension: 6; bounded discovery examined 14 boards')
    print('Connected geometric family k=1..8: global H=k, local h=0, tau=k')
    print('Geometric joint-choice forced-deletion costs: 00=3, 01=4, 10=4, 11=2')
    print('Naive articulation-removal quotient cycle rank: 1; correct quotient: forest')
    print('Negative continuation reward and four-direction high-incidence parent controls: PASS')
    print('Always-active exact gates: '+str(sum(GATES.values())))
    print('Certificate: '+target.name)

if __name__=='__main__':main()
