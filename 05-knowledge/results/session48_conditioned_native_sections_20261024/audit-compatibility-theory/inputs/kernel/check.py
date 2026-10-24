#!/usr/bin/env python3
"""Independent finite controls for native equality holonomy and fibre bridges."""
from collections import deque
from itertools import combinations, permutations, product
from pathlib import Path
import argparse, hashlib, json

HERE = Path(__file__).resolve().parent
CHECKS = 0

def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)

# The compiler below never uses path records; literal comparisons below do.
def compile_voltage(n, edges, pins=()):
    """edge (i,j,d) means phase[j]-phase[i]=d modulo 3."""
    adj = [[] for _ in range(n)]
    for i,j,d in edges:
        adj[i].append((j,d % 3)); adj[j].append((i,(-d) % 3))
    component = [-1]*n; potential = [None]*n; members=[]
    for root in range(n):
        if component[root] >= 0: continue
        c=len(members); members.append([]); component[root]=c; potential[root]=0
        todo=[root]
        while todo:
            i=todo.pop(); members[c].append(i)
            for j,d in adj[i]:
                expected=(potential[i]+d)%3
                if component[j] < 0:
                    component[j]=c; potential[j]=expected; todo.append(j)
                elif potential[j] != expected:
                    return None
    fixed={}
    for i,k in pins:
        c=component[i]; t=(k-potential[i])%3
        if c in fixed and fixed[c]!=t: return None
        fixed[c]=t
    return dict(component=component,potential=potential,members=members,fixed=fixed)

def compiled_states(f):
    if f is None: return []
    domains=[(f['fixed'][c],) if c in f['fixed'] else range(3) for c in range(len(f['members']))]
    return sorted(tuple((roots[f['component'][i]]+p)%3 for i,p in enumerate(f['potential'])) for roots in product(*domains))

P=((3,0,6,7,1,5,4,2),(3,1,6,7,0,5,4,2))
G=(4,0,2,3,1,5,6,7)

def source(r):
    core=[[0]*8 for _ in range(8)]
    for z,(i,j) in enumerate(combinations(range(8),2)):
        core[i][j]=(29558785>>(27-z))&1; core[j][i]=1-core[i][j]
    q=8*r
    return tuple(tuple(0 if i==j else core[i%8][j%8] if i<q and j<q and i//8==j//8 else int(i<j) if i<q and j<q else int(i<q) if (i<q)!=(j<q) else int((j-i)%3==1) for j in range(q+3)) for i in range(q+3))

def path(r,b,k):
    return tuple(8*j+x for j in range(r) for x in P[(b>>j)&1])+tuple(8*r+(k+t)%3 for t in range(3))

def partner(a,p):
    pos={u:i for i,u in enumerate(p)}
    return tuple(tuple(0 if u==v else a[u][v] if abs(pos[u]-pos[v])==1 else a[v][u] for v in range(len(a))) for u in range(len(a)))

def word(a,p):
    return tuple(a[p[i]][p[j]] for i in range(len(p)) for j in range(i+2,len(p)))

def literal_states(n,edges,pins=(),r=None):
    # Realizes d by equality of slot i clock q+d and slot j clock q.
    r = max(1,(n-1).bit_length()) if r is None else r
    pp=[[path(r,i,k) for k in range(3)] for i in range(n)]
    return [ks for ks in product(range(3),repeat=n) if all(pp[i][ks[i]][8*r+d%3]==pp[j][ks[j]][8*r] for i,j,d in edges) and all(pp[i][ks[i]][8*r]==8*r+k%3 for i,k in pins)]

def native_checks():
    rows=[]
    for r in range(1,5):
        a=source(r); n=len(a); rev=tuple(zip(*a))
        pp={(b,k):path(r,b,k) for b in range(1<<r) for k in range(3)}
        partners={z:partner(a,p) for z,p in pp.items()}; ww={z:word(a,p) for z,p in pp.items()}
        family=set(ww.values()); reverse={word(rev,p[::-1]) for p in pp.values()}
        reservoir=family|reverse|{tuple(1-z for z in w) for w in family|reverse}
        need(len(family)==1<<r,'prefix keys must be distinct source resources')
        need(len(reservoir)==4*(1<<r),'normalized K/converse tile reservoir')
        for b in range(1<<r):
            need(len({ww[b,k] for k in range(3)})==1,'three lifts consume one unit')
            prefix=pp[b,0][:-3]
            completions={prefix+tail for tail in permutations(range(8*r,n)) if all(a[u][v] for u,v in zip(prefix+tail,(prefix+tail)[1:]))}
            need(completions=={pp[b,k] for k in range(3)},'complete prefix-pinned native menu')
        for (b,k),p in pp.items():
            need(sorted(p)==list(range(n)) and all(a[u][v] for u,v in zip(p,p[1:])),'actual Hamiltonian record')
            need(partner(rev,p[::-1])==tuple(zip(*partners[b,k])),'literal converse partner square')
            need(all(p[::-1][2-t]==p[8*r+t] for t in range(3)),'all indexed final clocks reflected correctly')
            for delta in (1,2):
                perm=tuple(u if u<8*r else 8*r+(u-8*r+delta)%3 for u in range(n))
                target=(b,(k+delta)%3)
                need(tuple(perm[u] for u in p)==pp[target],'literal native kernel owner transport')
                need(all(a[u][v]==a[perm[u]][perm[v]] for u in range(n) for v in range(n)),'actual source automorphism')
                need(all(partners[b,k][u][v]==partners[target][perm[u]][perm[v]] for u in range(n) for v in range(n)),'actual kernel partner gluing')
                need(sum(x!=y for x,y in zip(p,pp[target]))==3,'exact kernel clock work')
            for j in range(r):
                g=G if not (b>>j)&1 else tuple(G.index(u) for u in range(8))
                perm=tuple(8*j+g[u%8] if u//8==j else u for u in range(n)); target=(b^(1<<j),k)
                need(all(partners[b,k][u][v]==partners[target][perm[u]][perm[v]] for u in range(n) for v in range(n)),'endpoint-marked prefix partner gluing')
        rows.append(dict(prefix_blocks=r,source_order=n,full_records=len(pp),normalized_tiles=1<<r,reservoir_words=len(reservoir)))
    # Authenticate inherited rigidity separately, without an isomorphism library.
    a=source(1); seed=tuple(row[:8] for row in a[:8]); autos=[]
    deg=[sum(row) for row in seed]
    for p in permutations(range(8)):
        if all(deg[i]==deg[p[i]] for i in range(8)) and all(seed[i][j]==seed[p[i]][p[j]] for i in range(8) for j in range(8)):
            autos.append(p)
    need(autos==[tuple(range(8))],'rigid seed has no hidden resource identifications')
    # Arbitrary prefix/final indexed tests, not only the standard realizing clocks.
    r=2; allp=[[path(r,i,k) for k in range(3)] for i in range(4)]; tests=0
    for i,j in combinations(range(4),2):
        for x,y in product(range(19),repeat=2):
            actual={(u,v) for u,v in product(range(3),repeat=2) if allp[i][u][x]==allp[j][v][y]}
            if x<16 and y<16:
                expected=set(product(range(3),repeat=2)) if allp[i][0][x]==allp[j][0][y] else set()
            elif (x<16)!=(y<16): expected=set()
            else: expected={(u,v) for u,v in product(range(3),repeat=2) if (v-u)%3==(x-y)%3}
            need(actual==expected,'indexed-clock equality compiler'); tests+=1
    pin_tests=0
    for i,x,owner in product(range(4),range(19),range(19)):
        actual={k for k in range(3) if allp[i][k][x]==owner}
        if x<16: expected=set(range(3)) if allp[i][0][x]==owner else set()
        elif owner<16: expected=set()
        else: expected={(owner-x)%3}
        need(actual==expected,'arbitrary indexed owner pin compiler');pin_tests+=1
    for start in range(8):
        seen={start};todo=[start]
        while todo:
            i=todo.pop()
            for j in range(8):
                if seed[i][j] and j not in seen:seen.add(j);todo.append(j)
        need(len(seen)==8,'seed strong component is genuine')
    return dict(records=rows,seed_automorphisms=len(autos),indexed_clock_pair_tests=tests,indexed_owner_pin_tests=pin_tests)

def graph_catalogue(n):
    pairs=list(combinations(range(n),2))
    for code in product(range(4),repeat=len(pairs)):
        yield [(i,j,d-1) for (i,j),d in zip(pairs,code) if d]

def classification_checks():
    rows=[]; catalogues={}
    for n in range(1,5):
        graph_cases=pin_cases=frustrated=0; unique={}
        for edges in graph_catalogue(n):
            f=compile_voltage(n,edges); actual=literal_states(n,edges)
            need(compiled_states(f)==actual,'graph holonomy iff literal owner feasibility')
            graph_cases+=1; frustrated+=f is None
            if f is not None: unique[tuple(actual)]=f
            if n<=3:
                for pv in product(range(4),repeat=n):
                    pins=[(i,v-1) for i,v in enumerate(pv) if v]
                    need(compiled_states(compile_voltage(n,edges,pins))==literal_states(n,edges,pins),'pin consistency iff literal owner feasibility')
                    pin_cases+=1
        catalogues[n]=unique
        rows.append(dict(slots=n,voltage_graphs=graph_cases,pinned_cases=pin_cases,frustrated_graphs=frustrated,distinct_nonempty_unpinned_fibres=len(unique)))
    # Parallel contradictory edges and self-loop constraints are also allowed.
    for edges in [[(0,1,0),(0,1,1)],[(0,0,1)],[(0,0,0)],[(0,1,2),(1,0,1)]]:
        need(compiled_states(compile_voltage(2,edges))==literal_states(2,edges),'parallel and loop relation control')
    return rows,catalogues

def overlay(f,g):
    # Independent fibres give a left and right component-root vertex per component.
    offset=len(f['members']); n=offset+len(g['members'])
    edges=[(f['component'][i],offset+g['component'][i],(f['potential'][i]-g['potential'][i])%3) for i in range(len(f['component']))]
    pins=list(f['fixed'].items())+[(offset+c,t) for c,t in g['fixed'].items()]
    return n,edges,pins

def frustration_distance(f,g,weights):
    n,edges,pins=overlay(f,g); best=None; witness=None
    for mask in range(1<<len(edges)):
        cost=sum(w for i,w in enumerate(weights) if mask>>i&1)
        if best is not None and cost>=best: continue
        keep=[e for i,e in enumerate(edges) if not(mask>>i&1)]
        model=compile_voltage(n,keep,pins)
        if model is not None:
            best=cost; witness=compiled_states(model)[0]
    if best is None: raise RuntimeError('own fibre pins cannot conflict across disjoint root sets')
    left=len(f['members'])
    x=tuple((witness[f['component'][i]]+f['potential'][i])%3 for i in range(len(weights)))
    y=tuple((witness[left+g['component'][i]]+g['potential'][i])%3 for i in range(len(weights)))
    return best,x,y

def boundary_checks(catalogues):
    rows=[]
    for n in range(1,5):
        cases=0; positive=0; maxdelta=0
        entries=list(catalogues[n].items())
        for (aa,f),(bb,g) in product(entries,repeat=2):
            weights=tuple(1+(i%2) for i in range(n))
            delta,x,y=frustration_distance(f,g,weights)
            brute=min(sum(w for w,u,v in zip(weights,a,b) if u!=v) for a in aa for b in bb)
            need(delta==brute,'overlay frustration index equals full literal Hamming minimum')
            need(x in aa and y in bb,'overlay minimizer reconstructs actual fibre records')
            need(sum(w for w,u,v in zip(weights,x,y) if u!=v)==delta,'boundary witness attains optimum')
            need((delta==0)==bool(set(aa)&set(bb)),'positive weights detect actual fibre intersection')
            cases+=1; positive+=delta>0; maxdelta=max(maxdelta,delta)
        rows.append(dict(slots=n,fibre_pairs=cases,positive_distance_pairs=positive,max_weighted_distance=maxdelta,weights=list(1+i%2 for i in range(n))))
    # Pinned fibre controls include a conflict along a tree, not only cycle frustration.
    fs=[]
    for n in (2,3):
        seen={}
        for edges in graph_catalogue(n):
            for pins in [[],[(0,0)],[(0,1)],[(0,0),(n-1,2)]]:
                f=compile_voltage(n,edges,pins)
                if f is not None: seen[tuple(compiled_states(f))]=f
        count=0
        for (aa,f),(bb,g) in product(seen.items(),repeat=2):
            weights=(1,)*n;delta,x,y=frustration_distance(f,g,weights)
            need(delta==min(sum(u!=v for u,v in zip(a,b)) for a in aa for b in bb),'pinned overlay distance including path pin incompatibility')
            need(x in aa and y in bb,'pinned boundary witness');count+=1
        fs.append(dict(slots=n,distinct_fibres=len(seen),pairs=count))
    return dict(unpinned=rows,pinned=fs)

def components(states,budget,cost):
    todo=set(states); out=[]
    while todo:
        root=min(todo); todo.remove(root); found={root}; queue=[root]
        while queue:
            x=queue.pop()
            new={y for y in todo if cost(x,y)<=budget}
            todo-=new; found|=new; queue.extend(new)
        out.append(found)
    return sorted(len(c) for c in out)

def compact_native_bridge():
    r=3;n=4
    ea=[(0,1,0),(2,3,0)];eb=[(0,1,1),(2,3,1)]
    f=compile_voltage(n,ea);g=compile_voltage(n,eb)
    aa=literal_states(n,ea,r=r);bb=literal_states(n,eb,r=r)
    need(aa==compiled_states(f) and bb==compiled_states(g),'compact native chart completeness')
    delta,x,y=frustration_distance(f,g,(1,)*n)
    need(delta==2 and len(aa)==len(bb)==9,'compact two-parallel-cycle boundary')
    states=[(0,)+a for a in aa]+[(1,)+b for b in bb]
    counts={str(b):components(states,b,lambda x,y:sum(a!=c for a,c in zip(x,y))) for b in (2,3)}
    need(counts=={'2':[9,9],'3':[18]},'compact sharp budget component counts')
    a=source(r);pp={(i,k):path(r,i,k) for i in range(5) for k in range(3)}
    ww={z:word(a,p) for z,p in pp.items()}
    for state in states:
        phases=state[1:]+(state[0],)
        need(len({ww[i,k] for i,k in enumerate(phases)})==5,'compact five actual resources remain injective')
    need(sum(sum(u!=v for u,v in zip(pp[i,x[i]],pp[i,y[i]])) for i in range(n))==6,'compact old boundary clock work')
    need(sum(u!=v for u,v in zip(pp[4,0],pp[4,1]))==3,'compact helper clock work')
    return dict(source_order=27,distinct_resources=5,total_states=18,fibre_sizes=[9,9],sharp_slot_budget=3,sharp_clock_budget=9,component_sizes_by_slot_budget=counts,witness_left=list(x),witness_right=list(y))

def native_hostiles():
    # All three edge constraints are individually 3-by-3 perfect matchings.
    edges=[(0,1,0),(1,2,0),(2,0,1)]
    edge_counts=[len(literal_states(3,[e])) for e in edges]
    need(edge_counts==[9,9,9] and not literal_states(3,edges),'locally total matching does not glue around a frustrated cycle')
    # Eight old resources + one actual helper resource, all in S_4.
    r=4;n=8
    ea=[(s,s+1,0) for s in (0,2,4,6)]
    eb=[(s,s+2,0) for s in (0,4)]+[(s+1,s+3,1) for s in (0,4)]
    f=compile_voltage(n,ea);g=compile_voltage(n,eb)
    aa=literal_states(n,ea,r=r);bb=literal_states(n,eb,r=r)
    need(aa==compiled_states(f) and bb==compiled_states(g),'native conditional equality charts are complete')
    need(len(aa)==len(bb)==81 and not(set(aa)&set(bb)),'two neutral fibres are nonempty and disjoint')
    need(all({a[i] for a in aa}=={a[i] for a in bb}=={0,1,2} for i in range(n)),'every old slot has full local phase coverage in both charts')
    delta,x,y=frustration_distance(f,g,(1,)*n)
    need(delta==2,'two disjoint frustrated overlay cycles force two changed old slots')
    states=[(0,)+a for a in aa]+[(1,)+b for b in bb]
    counts={}
    for budget in (1,2,3):
        counts[str(budget)]=components(states,budget,lambda x,y:sum(a!=b for a,b in zip(x,y)))
    need(counts['2']==[81,81] and counts['3']==[162],'sharp native helper crossing budget is three')
    # Compare full physical records and reservoir capacity for every arrival.
    allpaths={(i,k):path(r,i,k) for i in range(9) for k in range(3)};a=source(r)
    allwords={(i,k):word(a,p) for (i,k),p in allpaths.items()}
    for state in states:
        phases=state[1:]+(state[0],)
        need(len({allwords[i,k] for i,k in enumerate(phases)})==9,'native helper and old slots use nine distinct units, not twenty-seven')
    oldcost=sum(sum(a!=b for a,b in zip(allpaths[i,x[i]],allpaths[i,y[i]])) for i in range(8))
    helpercost=sum(a!=b for a,b in zip(allpaths[8,0],allpaths[8,1]))
    need(oldcost==6 and helpercost==3,'sharp bridge costs nine actual owner-clock occurrences')
    return dict(frustrated_triangle=dict(source_order=19,distinct_resources=3,accepted_states=0,single_edge_accepted_states=edge_counts,each_edge_two_slot_matches=3),conditional_helper=dict(source_order=35,distinct_resources=9,helper_record_phases=[0,1],old_fibre_sizes=[len(aa),len(bb)],neutral_component_size=2,boundary_old_changed_slots=delta,sharp_total_changed_slot_budget=3,sharp_total_clock_budget=9,component_sizes_by_slot_budget=counts,witness_old_left=list(x),witness_old_right=list(y),guard_scope='explicit helper-conditioned higher-arity equality predicate; not inferred from pairwise owner equalities'))

def distance_enumerator_checks(catalogues):
    cases=budgets=routes=0
    for n,cat in catalogues.items():
        for state_tuple,f in cat.items():
            states=list(state_tuple); reference=states[0]
            polynomial={(0,0):1}
            free=[c for c in range(len(f['members'])) if c not in f['fixed']]
            for c in free:
                m=len(f['members'][c]); new={}
                for (x,y),value in polynomial.items():
                    new[x,y]=new.get((x,y),0)+value
                    new[x+m,y+3*m]=new.get((x+m,y+3*m),0)+2*value
                polynomial=new
            observed={}
            for state in states:
                dist=sum(a!=b for a,b in zip(reference,state))
                observed[dist,3*dist]=observed.get((dist,3*dist),0)+1
            need(polynomial==observed,'all-size joint packet distance enumerator finite controls');cases+=1
            for budget in range(n+1):
                actual=components(states,budget,lambda x,y:sum(a!=b for a,b in zip(x,y)))
                heavy=sum(len(f['members'][c])>budget for c in free)
                need(actual==[3**(len(free)-heavy)]*(3**heavy),'exact heavy-phase budget component count');budgets+=1
            # For every endpoint pair, direct component flips must never return
            # to a departed full record and never change the resource key.
            for start,target in product(states,repeat=2):
                current=list(start);visited=[{k} for k in current]
                for c in free:
                    if current[f['members'][c][0]]==target[f['members'][c][0]]:continue
                    for i in f['members'][c]:
                        need(target[i] not in visited[i],'direct kernel route never returns to a departed record')
                        visited[i].add(target[i]);current[i]=target[i]
                    need(tuple(current) in state_tuple,'direct component route remains accepted')
                need(tuple(current)==target,'direct no-return route reaches exact record endpoint');routes+=1
    return dict(fibres=cases,budget_graphs=budgets,direct_no_record_return_routes=routes,resource_semantics='one fixed prefix key continuously occupied per slot; no release or reacquisition')

def quotient_checks():
    # Abstract exact additive sidecar: independent from native zero-charge controls.
    # Exhaustive small component charge maps, with both neutral and charged factors.
    cases=0
    for n in range(1,4):
        weights=tuple(1+i%2 for i in range(n))
        for codes in product(range(4),repeat=n):
            charge=[(0,0,0) if code==0 else tuple((code*k)%4 for k in range(3)) for code in codes]
            for target in range(4):
                states=[x for x in product(range(3),repeat=n) if sum(charge[i][k] for i,k in enumerate(x))%4==target]
                for budget in range(1,5):
                    visible=[i for i in range(n) if len(set(charge[i]))>1 or weights[i]>budget]
                    quotient=sorted({tuple(x[i] for i in visible) for x in states})
                    full=components(states,budget,lambda x,y:sum(w for w,a,b in zip(weights,x,y) if a!=b))
                    qcomp=components(quotient,budget,lambda x,y:sum(weights[i] for i,a,b in zip(visible,x,y) if a!=b))
                    factor=3**(n-len(visible))
                    need(full==[x*factor for x in qcomp],'full additive fibre packet graph vs exact neutral quotient')
                    cases+=1
    return dict(cases=cases,charge_group='Z/4',components='1..3 independent C3 factors, complete charge maps from the stated four-function family',budgets=[1,2,3,4])

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    classified,catalogues=classification_checks()
    result=dict(status='FINITE-EXACT PASS; all-size claims proved in RELATIONAL_KERNEL_HOLONOMY.md',native=native_checks(),classification=classified,boundary=boundary_checks(catalogues),distance_enumerators=distance_enumerator_checks(catalogues),native_controls=native_hostiles(),compact_native_bridge=compact_native_bridge(),charge_quotient=quotient_checks())
    result['source_sha256_lf']={p.name:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in sorted((HERE/'inputs').iterdir()) if p.is_file() and p.suffix in ('.md','.txt')}
    result['checks_before_receipt']=CHECKS
    out=HERE/'receipt.json'
    if args.write:out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:need(result==json.loads(out.read_text()),'frozen receipt exact replay')
    print('PASS',CHECKS,'checks;', 'wrote' if args.write else 'replayed','receipt.json')

if __name__=='__main__': main()
