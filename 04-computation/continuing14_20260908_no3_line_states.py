"""Exact unmerged-block repair with binary cell and capacity-two line states.

Standalone: no mathematical producer import. NetworkX's integer-weight blossom
is the graph consumer; whole-instance triple hitting and literal subset controls
independently test the compiled predicate and costs on the declared finite bank.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path
import sys
import networkx as nx

sys.stdout.reconfigure(newline='\n')
GATES=Counter()
INF=10**18
SLOPES=(1,-1,2)

def check(value,label):
    GATES[label]+=1
    if not value:raise RuntimeError('always-active gate failed: '+label)

def lines(points,slopes=SLOPES):
    by=defaultdict(list)
    for p,(r,c) in enumerate(points):
        for a in slopes:by[a,c-a*r].append(p)
    labels=sorted(a for a in by if len(by[a])>2)
    return labels,[tuple(by[a]) for a in labels]

def exact(ls):
    clauses=tuple(sorted({sum(1<<p for p in triple) for line in ls for triple in combinations(line,3)}))
    @lru_cache(None)
    def solve(cs):
        if not cs:return 0
        first=min(cs,key=int.bit_count)
        return 1+min(solve(tuple(t for t in cs if not t&(1<<p))) for p in range(first.bit_length()) if first>>p&1)
    return solve(clauses)

def raw_blocks(ls):
    adj=defaultdict(set)
    for a,line in enumerate(ls):
        for p in line:adj['L',a].add(('P',p));adj['P',p].add(('L',a))
    clock={};low={};stack=[];edge_blocks=[]
    def dfs(u,parent):
        clock[u]=low[u]=len(clock)
        for v in sorted(adj[u]):
            if v==parent:continue
            if v not in clock:
                stack.append((u,v));dfs(v,u);low[u]=min(low[u],low[v])
                if low[v]>=clock[u]:
                    piece=[]
                    while True:
                        edge=stack.pop();piece.append(edge)
                        if edge==(u,v):break
                    edge_blocks.append(piece)
            elif clock[v]<clock[u]:
                stack.append((u,v));low[u]=min(low[u],clock[v])
    for u in sorted(adj):
        if u not in clock:dfs(u,None)
    blocks=[]
    for edges in edge_blocks:
        incidence=defaultdict(list)
        for u,v in edges:
            a,p=(u[1],v[1]) if u[0]=='L' else (v[1],u[1])
            incidence[p].append(a)
        blocks.append({'lines':tuple(sorted({a for aa in incidence.values() for a in aa})),
                       'cells':tuple(sorted(incidence)),
                       'inc':{p:tuple(sorted(aa)) for p,aa in incidence.items()}})
    blocks.sort(key=lambda b:(b['lines'],b['cells']))
    occurrence=defaultdict(list)
    for j,block in enumerate(blocks):
        for a in block['lines']:occurrence['L',a].append(j)
        for p in block['cells']:occurrence['P',p].append(j)
    return blocks,{v:tuple(bb) for v,bb in occurrence.items() if len(bb)>1}

def weighted_kernel(cells,profiles,exact_line=None,exact_count=None):
    """Max sum(cell rewards)+sum(line degree rewards), with capacities<=2.

    cells are (original_id, incident_line_tuple, reward). Profiles map each
    line to its rewards for degrees0..capacity. Optional one exact degree is
    enforced lexicographically below mandatory original cell-port coverage.
    """
    if exact_line is not None:
        if exact_count<0 or exact_count>=len(profiles[exact_line]):return None
        profiles=dict(profiles)
        profiles[exact_line]=profiles[exact_line][:exact_count+1]
    graph=nx.Graph();slots={};charge={};constant=0
    for a,rewards in profiles.items():
        capacity=len(rewards)-1
        check(capacity in (0,1,2),'kernel line capacity type')
        slots[a]=[('slot',a,i) for i in range(capacity)]
        graph.add_nodes_from(slots[a])
        if capacity==0:constant+=rewards[0];charge[a]=0
        elif capacity==1:
            constant+=rewards[0];charge[a]=rewards[1]-rewards[0]
        else:
            r0,r1,r2=rewards
            shift=min(r1-r2,r0-r1)
            leaf_reward=r1-r2-shift
            pair_reward=r0-r2-2*shift
            check(pair_reward>=leaf_reward>=0,'arbitrary degree-profile gadget signs')
            charge[a]=-shift;constant+=r2+2*shift
            leaf=('leaf',a)
            for slot in slots[a]:graph.add_edge(slot,leaf,base=leaf_reward,ports=0)
            graph.add_edge(*slots[a],base=pair_reward,ports=0)
    ports=[]
    for p,inc,reward in cells:
        check(len(inc) in (1,2) and len(set(inc))==len(inc),'kernel original ordinary cell type')
        aa=('port',p,0);bb=('port',p,1);ports.extend((aa,bb))
        graph.add_edge(aa,bb,base=0,ports=2)
        for j,a in enumerate(inc):
            port=aa if j==0 else bb
            for slot in slots[a]:
                graph.add_edge(port,slot,base=charge[a]+(reward if j==0 else 0),ports=1,
                               selected_line=a,original_cell=p)
        if len(inc)==1:
            graph.add_edge(bb,('private',p),base=0,ports=1)
    exact_bonus=1+sum(abs(d['base']) for u,v,d in graph.edges(data=True))
    for u,v,d in graph.edges(data=True):
        d['pre']=d['base']+(exact_bonus if exact_line is not None and d.get('selected_line')==exact_line else 0)
    coverage_bonus=1+sum(abs(d['pre']) for u,v,d in graph.edges(data=True))
    for u,v,d in graph.edges(data=True):d['weight']=d['pre']+coverage_bonus*d['ports']
    matching=nx.max_weight_matching(graph,weight='weight')
    partner={}
    for u,v in matching:partner[u]=v;partner[v]=u
    check(all(p in partner for p in ports),'all original cell ports covered')
    chosen=[];degrees={a:0 for a in profiles}
    for p,inc,reward in cells:
        aa=('port',p,0);bb=('port',p,1)
        if partner[aa]!=bb:
            check(partner[aa][0]=='slot' and partner[bb][0] in ('slot','private'),'binary original-cell gadget')
            chosen.append(p)
            for a in inc:degrees[a]+=1
    score=sum(w for p,inc,w in cells if p in chosen)+sum(profiles[a][degrees[a]] for a in profiles)
    matching_value=sum(graph[u][v]['weight'] for u,v in matching)
    recovered=matching_value-len(ports)*coverage_bonus+constant
    if exact_line is not None:recovered-=exact_bonus*degrees[exact_line]
    check(score==recovered,'literal vertex costs equal matching gadget objective')
    if exact_line is not None and degrees[exact_line]!=exact_count:
        return None
    return score,tuple(chosen),degrees

def literal_kernel(cells,profiles,exact_line=None,exact_count=None):
    best=None
    for states in product((0,1),repeat=len(cells)):
        degrees={a:0 for a in profiles};reward=0
        for state,(p,inc,w) in zip(states,cells):
            if state:
                reward+=w
                for a in inc:degrees[a]+=1
        if any(degrees[a]>=len(profiles[a]) for a in profiles):continue
        if exact_line is not None and degrees[exact_line]!=exact_count:continue
        value=reward+sum(profiles[a][degrees[a]] for a in profiles)
        best=value if best is None else max(best,value)
    return best

def solve_blocks(ls,first_root=None,check_kernels=False):
    blocks,separators=raw_blocks(ls)
    calls=0;line_messages={};cell_messages={}
    local_h=[sum(len(aa)>=3 for aa in b['inc'].values()) for b in blocks]
    @lru_cache(None)
    def cell_message(p,parent):
        value=tuple(1-state+sum(block_message(b,('P',p),state) for b in separators['P',p] if b!=parent)
                    for state in (0,1))
        cell_messages[p,parent]=value
        return value
    @lru_cache(None)
    def line_message(a,parent):
        cost={0:0}
        for b in separators['L',a]:
            if b==parent:continue
            nxt={}
            for used,value in cost.items():
                for retained in range(3-used):
                    child=block_message(b,('L',a),retained)
                    if child>=INF:continue
                    total=used+retained
                    nxt[total]=min(nxt.get(total,INF),value+child)
            cost=nxt
        value=tuple(min(v for used,v in cost.items() if used<=2-parent_used) for parent_used in range(3))
        line_messages[a,parent]=value
        return value
    @lru_cache(None)
    def block_message(b,parent,parent_state):
        nonlocal calls
        block=blocks[b]
        fixed_parent_cell=parent[1] if parent is not None and parent[0]=='P' else None
        parent_line=parent[1] if parent is not None and parent[0]=='L' else None
        high=[];ordinary=[];base=0
        for p in block['cells']:
            if p==fixed_parent_cell:continue
            c0,c1=cell_message(p,b) if ('P',p) in separators else (1,0)
            base+=c0;item=(p,block['inc'][p],c0-c1)
            (high if len(item[1])>=3 else ordinary).append(item)
        costs={a:line_message(a,b) if ('L',a) in separators and a!=parent_line else (0,0,0)
               for a in block['lines']}
        best=INF
        for states in product((0,1),repeat=len(high)):
            fixed={a:0 for a in block['lines']};high_reward=0
            if fixed_parent_cell is not None and parent_state:
                for a in block['inc'][fixed_parent_cell]:fixed[a]+=1
            for state,(p,inc,w) in zip(states,high):
                if state:
                    high_reward+=w
                    for a in inc:fixed[a]+=1
            capacities={a:(parent_state if a==parent_line else 2)-fixed[a] for a in block['lines']}
            if min(capacities.values(),default=0)<0:continue
            profiles={a:tuple(-costs[a][fixed[a]+j] for j in range(capacities[a]+1)) for a in block['lines']}
            calls+=1
            exact_count=capacities[parent_line] if parent_line is not None else None
            result=weighted_kernel(ordinary,profiles,parent_line,exact_count)
            if check_kernels:
                brute=literal_kernel(ordinary,profiles,parent_line,exact_count)
                check((None if result is None else result[0])==brute,'compiled block kernel equals literal subsets')
            if result is not None:best=min(best,base-high_reward-result[0])
        return best
    seen=set();roots=[]
    def visit(b):
        if b in seen:return
        seen.add(b)
        block=blocks[b]
        for node in [('L',a) for a in block['lines']]+[('P',p) for p in block['cells']]:
            for child in separators.get(node,()):visit(child)
    order=list(range(len(blocks)))
    if first_root is not None:order.remove(first_root);order.insert(0,first_root)
    for b in order:
        if b not in seen:roots.append(b);visit(b)
    value=sum(block_message(b,None,0) for b in roots)
    return {'tau':value,'blocks':len(blocks),'local_h':local_h,'calls':calls,
            'bound':3*sum(2**h for h in local_h),
            'line_messages':[[a,b,list(v)] for (a,b),v in sorted(line_messages.items())],
            'cell_messages':[[p,b,list(v)] for (p,b),v in sorted(cell_messages.items())]}

def theta_chain(k):
    base=[(6,12),(12,12),(9,15),(10,8),(0,0),(8,10),(7,17),(11,11),(13,19),(15,18),(14,28)]
    return [(1000**j*(r+100),1000**j*(c+100)) for j in range(k) for r,c in base]

def row_points(rows):return [(r,c) for r,pair in enumerate(rows) for c in pair]

def gadget_bank():
    for rewards in product(range(-2,3),repeat=3):
        for cells in ([(0,(0,),1),(1,(0,),-1),(2,(0,),2)],
                      [(0,(0,1),1),(1,(0,1),-1),(2,(1,),2)]):
            profiles={0:rewards}
            if any(1 in inc for p,inc,w in cells):profiles[1]=(0,-2,1)
            for count in (None,0,1,2):
                trial=weighted_kernel(cells,profiles,None if count is None else 0,count)
                brute=literal_kernel(cells,profiles,None if count is None else 0,count)
                check((None if trial is None else trial[0])==brute,'profile bank exact matching optimum')
    naive=nx.Graph();naive.add_edge('p0','p1',weight=0);naive.add_edge('p0','slot',weight=10)
    match=nx.max_weight_matching(naive,weight='weight')
    check(sum(naive[u][v]['weight'] for u,v in match)==10,'unprotected half-cell hostile')
    protected=weighted_kernel([(0,(0,1),10)],{0:(0,0),1:(0,)})
    check(protected[0]==0 and protected[1]==(),'coverage bonus repairs half-cell hostile')
    return {'reward_profiles':125,'cell_banks':2,'exact_degree_modes':4,'comparisons':1000,
            'unprotected_half_cell_reward':10,'true_original_cell_reward':0}

def connected_parts(graph,removed=frozenset()):
    remaining=set(graph)-set(removed);out=[]
    while remaining:
        seed=min(remaining);remaining.remove(seed);part={seed};todo=[seed]
        while todo:
            u=todo.pop()
            for v in graph[u]:
                if v in remaining:remaining.remove(v);part.add(v);todo.append(v)
        out.append(part)
    return out

def validate(ls,reroot=False):
    check(all(len(line)>=3 and len(set(line))==len(line) for line in ls),'original line universe types')
    blocks,separators=raw_blocks(ls)
    graph=defaultdict(set);owned=Counter();quotient=defaultdict(set)
    for a,line in enumerate(ls):
        for p in line:graph['L',a].add(('P',p));graph['P',p].add(('L',a))
    for b,block in enumerate(blocks):
        quotient['B',b]
        local=defaultdict(set)
        for p,aa in block['inc'].items():
            for a in aa:
                owned[a,p]+=1;local['L',a].add(('P',p));local['P',p].add(('L',a))
        edge_count=sum(map(len,local.values()))//2
        check(len(connected_parts(local))==1,'ordinary block connected')
        if edge_count>1:
            for v in local:check(len(connected_parts(local,{v}))==1,'ordinary block has no articulation')
            beta=edge_count-len(local)+1
            high=sum(len(aa)>=3 for aa in block['inc'].values())
            check(high<=2*beta-2,'ordinary-block cycle-excess bound')
        else:check(sum(len(aa)>=3 for aa in block['inc'].values())==0,'bridge has no local hyperedge')
        for v in [('L',a) for a in block['lines']]+[('P',p) for p in block['cells']]:
            if v in separators:quotient['B',b].add(v);quotient[v].add(('B',b))
    check(owned==Counter((a,p) for a,line in enumerate(ls) for p in line),'every original incidence owned once')
    parts=connected_parts(graph)
    arts={v for v in graph if len(connected_parts(graph,{v}))>len(parts)}
    check(arts==set(separators),'both separator types equal deletion-based articulation test')
    qparts=connected_parts(quotient)
    check(sum(map(len,quotient.values()))//2==len(quotient)-len(qparts),'unmerged block-cut forest')
    check(len(qparts)==len(parts),'original components preserved')
    expected=exact(ls)
    result=solve_blocks(ls,check_kernels=True)
    check(result['tau']==expected,'full recurrence equals independent original triple hitting')
    check(result['calls']<=result['bound'],'unmerged matching-call bound')
    if reroot:
        for b in range(len(blocks)):
            other=solve_blocks(ls,b,check_kernels=True)
            check(other['tau']==expected,'every tested root preserves full optimum')
            check(other['calls']<=other['bound'],'every tested root respects graph-call bound')
    return result,blocks,separators

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

def abstract_bank(v):
    shared=list(combinations(range(v),2))+list(combinations(range(v),3))
    for bits in range(1<<len(shared)):
        chosen=[c for i,c in enumerate(shared) if bits>>i&1]
        degree=Counter(a for c in chosen for a in c)
        pair_degree=Counter(p for c in chosen for p in combinations(c,2))
        if max(degree.values(),default=0)>3 or max(pair_degree.values(),default=0)>1:continue
        ls=[[] for a in range(v)]
        for p,inc in enumerate(chosen):
            for a in inc:ls[a].append(p)
        nxt=len(chosen)
        for a in range(v):
            while len(ls[a])<3:ls[a].append(nxt);nxt+=1
        yield list(map(tuple,ls))

def compact(result):return result

def main():
    certificate={'scope':'Exact finite original-cell repair; no random-board frequency or extremal constant.',
                 'matching_consumer':'NetworkX '+nx.__version__+' integer-weight blossom',
                 'gadget_bank':gadget_bank()}
    census={}
    for n in range(2,6):
        count=0;calls=0;cell_separators=0;line_separators=0
        for rows in boards(n):
            points=row_points(rows);_,ls=lines(points)
            check(len(points)==2*n and len(set(points))==2*n,'simple two-regular board cells')
            check(set(Counter(c for r,c in points).values())=={2},'complete column degree two')
            result,bb,sep=validate(ls)
            calls+=result['calls'];count+=1
            cell_separators+=any(v[0]=='P' for v in sep)
            line_separators+=any(v[0]=='L' for v in sep)
        check(count=={2:1,3:6,4:90,5:2040}[n],'complete saturated universe count')
        census[n]={'boards':count,'matching_calls':calls,'boards_with_cell_separator':cell_separators,
                   'boards_with_line_separator':line_separators}
    certificate['complete_boards']=census
    abstract={}
    for v in range(5):
        count=0
        for ls in abstract_bank(v):validate(ls,reroot=True);count+=1
        check(count=={0:1,1:1,2:2,3:9,4:96}[v],'complete abstract linear triple universe count')
        abstract[v]=count
    certificate['abstract_systems']=abstract
    controls={
        'triangle':((0,3),(1,2),(0,4),(3,4),(1,2)),
        'five_cycle':((0,3),(1,2),(0,4),(1,4),(3,5),(2,5)),
        'concurrent_triple':((2,4),(0,1),(2,3),(3,4),(0,1)),
        'triangle_pendant':((0,1),(0,2),(3,4),(1,3),(2,4)),
        'larger_line':((0,1),(0,1),(2,3),(2,4),(3,4)),
        'saturated_core':((0,1),(0,1),(2,3),(4,5),(3,5),(2,4))}
    output_controls={}
    for name,rows in controls.items():
        modes={}
        for slopes in (SLOPES,(1,-1,2,-2)):
            points=row_points(rows);labels,ls=lines(points,slopes);result,_,_=validate(ls,reroot=True)
            modes[str(slopes)]={'lines':list(zip(labels,ls)),**compact(result)}
        output_controls[name]={'rows':rows,'modes':modes}
    check(output_controls['concurrent_triple']['modes'][str(SLOPES)]['tau']==1,'inherited true hyperedge control')
    check(output_controls['saturated_core']['modes'][str(SLOPES)]['tau']==4,'inherited saturated inseparable core')
    certificate['inherited_controls']=output_controls
    family=[]
    for k in range(1,6):
        points=theta_chain(k);labels,ls=lines(points);result,bb,seps=validate(ls,reroot=k<=2)
        a=labels.index((1,0));inc=Counter(p for line in ls for p in line)
        check(len(points)==11*k and len(ls)==5*k+1,'connected theta family sizes')
        check(len(ls[a])==3*k and all(len(line)==3 for j,line in enumerate(ls) if j!=a),'one original shared line and remaining triples')
        check(not any(v[0]=='P' for v in seps),'old cell separators cannot split theta family')
        check(result['blocks']==7*k and Counter(result['local_h'])==Counter({0:6*k,2:k}),'unmerged theta block parameters')
        check(result['tau']==4*k-2 and result['bound']==30*k,'exact theta repair and linear oracle bound')
        check(max(Counter(r for r,c in points).values())<=2 and max(Counter(c for r,c in points).values())<=2,'theta occupied rows and columns at most two')
        high=sorted(p for p,d in inc.items() if d>=3)
        check(len(high)==2*k,'old merged exceptional count')
        valid=0
        for states in product((0,1),repeat=len(high)):
            selected={p for p,s in zip(high,states) if s}
            valid+=all(len(selected.intersection(line))<=2 for line in ls)
        formula=2**k*(1+k+k*(k-1)//2)
        check(valid==formula,'old valid-branch count after capacity filter')
        family.append({'k':k,'points':points,'lines':list(zip(labels,ls)),
                       'old_merged_h':len(high),'old_valid_branches':valid,**compact(result)})
    certificate['connected_theta_family']=family
    labels,ls=lines(theta_chain(1));shared=ls[labels.index((1,0))];other=[l for l in ls if l!=shared]
    phi={r:INF for r in range(3)}
    for bits in range(1<<11):
        kept={p for p in range(11) if bits>>p&1};r=len(kept.intersection(shared))
        if r<3 and all(len(kept.intersection(line))<=2 for line in other):phi[r]=min(phi[r],11-len(kept))
    check(phi=={0:4,1:3,2:2},'literal shared-line contribution cost table')
    certificate['theta_contribution_costs']=phi
    double=[(100,100),(101,102),(102,104),(103,103),(104,102),
            (100000,100000),(101000,102000),(102000,104000),(103000,103000),(104000,102000)]
    labels,ls=lines(double);result,_,_=validate(ls,reroot=True)
    common=labels.index((1,0))
    check(any(a==common and costs==[1,2,2] for a,b,costs in result['line_messages']),'actual nonconvex line continuation')
    check(result['tau']==3 and not any(result['local_h']),'double-triangle graph-only exact repair')
    certificate['nonconvex_geometric_control']={'points':double,'lines':list(zip(labels,ls)),**compact(result)}
    certificate['gate_counts']=dict(sorted(GATES.items()));certificate['total_gates']=sum(GATES.values())
    directory=Path(__file__).resolve().parent
    if directory.name=='04-computation':directory=directory.parent/'05-knowledge'/'results'
    target=directory/(Path(__file__).stem+'_certificate.json')
    target.write_bytes((json.dumps(certificate,sort_keys=True,indent=2)+'\n').encode())
    print('Continuing14 unmerged line/cell states: PASS')
    print('Complete simple two-regular boards n=2..5: 2137')
    print('Complete abstract linear triple-incidence systems v=0..4: 109')
    print('Arbitrary degree-cost matching compiler: 1000 exact comparisons plus coverage hostile')
    print('Connected geometric theta family k=1..5: tau=4k-2; old branches exponential, new bound30k')
    print('Actual nonconvex line continuation [1,2,2], both separator types, and exact parent degrees: PASS')
    print('Always-active exact gates: '+str(sum(GATES.values())))
    print('Certificate: '+target.name)

if __name__=='__main__':main()
