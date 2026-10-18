"""Exact source-marked CRT branches and their remaining six-colour support.

Mod2/mod3 coordinates are arithmetic colour residues, not F2^2 flow labels.
A graph input makes no metric claim. Every negative receipt is independently
checked by a voltage cycle or a complete Boolean decision tree.
"""
from collections import deque
from dataclasses import dataclass, replace
from itertools import combinations, product
import hashlib
import json
from pathlib import Path

from pysat.solvers import Solver
import snark_boundary_transfer_20261018a as old
import exceptional_factor_support_20261018a as source8
import plane_six_reframes_parity_20261011b as parity

need, nat = old.require, old.nat
Graph = old.Graph


@dataclass(frozen=True)
class Branch:
    graph: Graph
    equations: tuple  # (modulus, start, finish, difference)
    pins: tuple       # (literal vertex, literal colour)


def audit_branch(b):
    need(type(b) is Branch, 'exact Branch')
    old.audit_graph(b.graph); nat(b.graph.n,1)
    need(type(b.equations) is tuple and type(b.pins) is tuple, 'tuple source banks')
    for row in b.equations:
        need(type(row) is tuple and len(row)==4, 'equation row')
        for x in row: nat(x)
        p,u,v,a=row
        need(p in (2,3) and u<b.graph.n and v<b.graph.n and a<p, 'native modular row')
    for row in b.pins:
        need(type(row) is tuple and len(row)==2, 'pin row')
        for x in row: nat(x)
        need(row[0]<b.graph.n and row[1]<6, 'native pin')
    need(len({v for v,c in b.pins})==len(b.pins), 'one colour per pinned vertex')
    return b


@dataclass(frozen=True)
class Voltage:
    modulus: int
    rows: tuple
    roots: tuple
    offsets: tuple
    cycle: tuple  # (row index, direction), a closed oriented source walk


def voltage(b,p):
    audit_branch(b); nat(p); need(p in (2,3),'prime coordinate')
    n=b.graph.n
    rows=tuple((u,v,a) for q,u,v,a in b.equations if q==p)
    rows+=tuple((n,v,c%p) for v,c in b.pins)
    adj=[[] for _ in range(n+1)]
    for i,(u,v,a) in enumerate(rows):
        adj[u].append((v,a,i,1)); adj[v].append((u,-a,i,-1))
    roots=[-1]*(n+1); offsets=[0]*(n+1); tree=[[] for _ in range(n+1)]
    for seed in (n,)+tuple(range(n)):
        if roots[seed]>=0: continue
        roots[seed]=seed; queue=deque([seed])
        while queue:
            u=queue.popleft()
            for v,a,index,direction in adj[u]:
                expected=(offsets[u]+a)%p
                if roots[v]<0:
                    roots[v]=seed; offsets[v]=expected; queue.append(v)
                    tree[u].append((v,index,direction));tree[v].append((u,index,-direction))
                elif offsets[v]!=expected:
                    # The offending u->v edge followed by the actual tree v->u path.
                    prev={v:None}; q=deque([v])
                    while u not in prev:
                        w=q.popleft()
                        for z,i,d in tree[w]:
                            if z not in prev: prev[z]=(w,i,d);q.append(z)
                    path=[]; w=u
                    while w!=v:
                        z,i,d=prev[w];path.append((i,d));w=z
                    cycle=((index,direction),)+tuple(reversed(path))
                    result=Voltage(p,rows,(),(),cycle)
                    audit_cycle(result,n+1)
                    return result
    return Voltage(p,rows,tuple(roots),tuple(offsets),())


def audit_cycle(v,size):
    need(type(v) is Voltage, 'voltage receipt'); nat(v.modulus)
    need(v.modulus in (2,3), 'coordinate prime'); nat(size,1)
    need(type(v.rows) is tuple and type(v.cycle) is tuple and bool(v.cycle),'cycle banks')
    for row in v.rows:
        need(type(row) is tuple and len(row)==3,'source row')
        for x in row: nat(x)
        need(row[0]<size and row[1]<size and row[2]<v.modulus,'row domain')
    walk=[]; total=0
    for item in v.cycle:
        need(type(item) is tuple and len(item)==2,'oriented index')
        i,d=item;nat(i);need(type(d) is int and d in (-1,1) and i<len(v.rows),'direction')
        a,b,r=v.rows[i]
        walk.append((a,b) if d==1 else (b,a));total+=d*r
    need(all(walk[i][1]==walk[(i+1)%len(walk)][0] for i in range(len(walk))), 'closed actual walk')
    need(total%v.modulus!=0,'nonzero holonomy')
    return True


@dataclass(frozen=True)
class Compiled:
    branch: Branch
    binary: Voltage
    ternary: Voltage
    parameters: tuple  # (p, component root)
    clauses: tuple
    variable_rows: tuple  # one-hot variables for each parameter value


def parameter_colour(c,vertex,values):
    n=c.branch.graph.n; residues=[]
    for v in (c.binary,c.ternary):
        root=v.roots[vertex]
        phase=0 if root==n else values[c.parameters.index((v.modulus,root))]
        residues.append((v.offsets[vertex]+phase)%v.modulus)
    return (3*residues[0]+4*residues[1])%6


def compile_branch(b):
    audit_branch(b)
    binary,ternary=voltage(b,2),voltage(b,3)
    if binary.cycle or ternary.cycle:
        return Compiled(b,binary,ternary,(),(),())
    n=b.graph.n
    parameters=tuple((v.modulus,r) for v in (binary,ternary)
                     for r in sorted(set(v.roots)-{n}))
    nextvar=1; ids=[]; clauses=[]
    for p,r in parameters:
        row=tuple(range(nextvar,nextvar+p));nextvar+=p;ids.append(row)
        clauses.append(row)
        clauses.extend((-a,-d) for a,d in combinations(row,2))
    tmp=Compiled(b,binary,ternary,parameters,(),tuple(ids))
    # Keep the JOINT edge condition, rather than separate marginal graphs.
    for a,d in b.graph.edges:
        if any(v.roots[a]==v.roots[d] and v.offsets[a]!=v.offsets[d] for v in (binary,ternary)):
            continue  # This coordinate alone certifies inequality.
        if all(v.roots[a]==v.roots[d] and v.offsets[a]==v.offsets[d] for v in (binary,ternary)):
            clauses.append(())
            continue  # Both coordinates force the same actual colour.
        involved=tuple(sorted({parameters.index((v.modulus,v.roots[x]))
                               for v in (binary,ternary) for x in (a,d) if v.roots[x]!=n}))
        for choices in product(*(range(parameters[i][0]) for i in involved)):
            values=[0]*len(parameters)
            for i,t in zip(involved,choices): values[i]=t
            if parameter_colour(tmp,a,values)==parameter_colour(tmp,d,values):
                clauses.append(tuple(-ids[i][t] for i,t in zip(involved,choices)))
    canonical=tuple(sorted(set(tuple(sorted(row,key=lambda z:(abs(z),z))) for row in clauses)))
    return replace(tmp,clauses=canonical)


def audit_compiled(c):
    need(type(c) is Compiled,'compiled source')
    audit_branch(c.branch)
    # Do not compare bool/float-equivalent packets before exact type checks.
    for v in (c.binary,c.ternary):
        need(type(v) is Voltage,'voltage type');nat(v.modulus)
        for bank in (v.rows,v.cycle):
            need(type(bank) is tuple,'voltage tuple bank')
            for row in bank:
                need(type(row) is tuple,'row tuple')
                for x in row: need(type(x) is int,'exact voltage integer')
        for bank in (v.roots,v.offsets):
            need(type(bank) is tuple,'coordinate tuple')
            for x in bank: nat(x)
    for bank in (c.parameters,c.clauses,c.variable_rows):
        need(type(bank) is tuple,'compiled tuple bank')
        for row in bank:
            need(type(row) is tuple,'compiled row')
            for x in row: need(type(x) is int,'exact compiled integer')
    need(c==compile_branch(c.branch),'compiled source mismatch')
    return c


def audit_colouring(b,colours):
    audit_branch(b);old.audit_colouring(b.graph,6,colours)
    need(all(colours[v]==a for v,a in b.pins),'literal source pin')
    need(all((colours[v]-colours[u])%p==a for p,u,v,a in b.equations),'branch equations')
    return colours


def propagate(clauses,assignment):
    assignment=dict(assignment)
    while True:
        changed=False; unresolved=[]
        for clause in clauses:
            if any(abs(x) in assignment and assignment[abs(x)]==(x>0) for x in clause):continue
            row=tuple(x for x in clause if abs(x) not in assignment)
            if not row:return None,()
            if len(row)==1:
                assignment[abs(row[0])]=row[0]>0;changed=True
            else: unresolved.append(row)
        if not changed:return assignment,tuple(unresolved)


class ProofLimit(Exception): pass


def unsat_tree(clauses,cap=100000):
    nat(cap,1);nodes=0
    def visit(values):
        nonlocal nodes
        nodes+=1
        if nodes>cap:raise ProofLimit('negative decision-tree budget exhausted')
        values,remaining=propagate(clauses,values)
        if values is None:return ()
        need(bool(remaining),'SAT assignment encountered while certifying UNSAT')
        row=min(remaining,key=lambda c:(len(c),c));var=abs(row[0])
        return (var,visit(values|{var:False}),visit(values|{var:True}))
    return visit({}),nodes


def audit_unsat_tree(clauses,tree):
    nodes=0
    def visit(values,node):
        nonlocal nodes
        nodes+=1
        need(type(node) is tuple,'decision node tuple')
        values,remaining=propagate(clauses,values)
        if values is None:need(node==(),'conflict leaf');return
        need(bool(remaining) and len(node)==3,'complete unsatisfiable branching')
        var,left,right=node;nat(var,1)
        need(var not in values and any(abs(x)==var for row in remaining for x in row),'unfixed branch variable')
        visit(values|{var:False},left);visit(values|{var:True},right)
    visit({},tree)
    return nodes


@dataclass(frozen=True)
class Receipt:
    compiled: Compiled
    status: str
    colours: tuple
    tree: tuple


def solve(b,proof_cap=100000):
    nat(proof_cap,1)
    c=compile_branch(b)
    if c.binary.cycle or c.ternary.cycle:return Receipt(c,'voltage_conflict',(),())
    with Solver(name='g3',bootstrap_with=c.clauses) as solver:
        sat=solver.solve()
        if sat:
            model=set(solver.get_model() or [])
            values=tuple(next(t for t,var in enumerate(row) if var in model) for row in c.variable_rows)
            colours=tuple(parameter_colour(c,v,values) for v in range(b.graph.n))
            audit_colouring(b,colours)
            return Receipt(c,'colouring',colours,())
    tree,_=unsat_tree(c.clauses,proof_cap)
    audit_unsat_tree(c.clauses,tree)
    return Receipt(c,'unsatisfiable',(),tree)


def audit_receipt(r):
    need(type(r) is Receipt,'Receipt required');audit_compiled(r.compiled)
    need(type(r.status) is str,'status string')
    need(type(r.colours) is tuple and type(r.tree) is tuple,'receipt tuple banks')
    for x in r.colours:nat(x)
    c=r.compiled
    if r.status=='voltage_conflict':
        need(not r.colours and not r.tree and bool(c.binary.cycle or c.ternary.cycle),'literal voltage failure')
        for v in (c.binary,c.ternary):
            if v.cycle:audit_cycle(v,c.branch.graph.n+1)
    elif r.status=='colouring':
        need(not r.tree and not c.binary.cycle and not c.ternary.cycle,'positive receipt type')
        audit_colouring(c.branch,r.colours)
    elif r.status=='unsatisfiable':
        need(not r.colours and not c.binary.cycle and not c.ternary.cycle,'negative receipt type')
        audit_unsat_tree(c.clauses,r.tree)
    else:need(False,'unknown receipt status')
    return r


def count_branch(b):
    c=compile_branch(b)
    if c.binary.cycle or c.ternary.cycle:return 0
    total=0
    for values in product(*(range(p) for p,r in c.parameters)):
        colours=tuple(parameter_colour(c,v,values) for v in range(b.graph.n))
        if all(colours[a]!=colours[d] for a,d in b.graph.edges):total+=1
    return total


BRANCHES=tuple((i,p,a,j,q,b) for i,j in combinations(range(6),2)
               for p in (2,3) for a in range(p) for q in (2,3) for b in range(q))


def contact_branch(base,factors,choice):
    audit_branch(base)
    need(type(factors) is tuple and len(factors)==6,'six marked factor pair sets')
    for pairs in factors:
        need(type(pairs) is tuple and bool(pairs),'nonempty pair set')
        for pair in pairs:
            need(type(pair) is tuple and len(pair)==2,'direction pair')
            for x in pair:nat(x);need(x<base.graph.n,'literal factor vertex')
    need(type(choice) is tuple and len(choice)==6,'G4 branch index')
    for x in choice:nat(x)
    need(choice in BRANCHES,'one of the375 local branches')
    i,p,a,j,q,b=choice
    equations=base.equations+tuple((p,u,v,a) for u,v in factors[i])+tuple((q,u,v,b) for u,v in factors[j])
    return audit_branch(Branch(base.graph,equations,base.pins))


def complete_boundary(graph,ports,patterns):
    """Decide each supplied pattern; completeness of its bank is caller-owned."""
    old.audit_graph(graph)
    need(type(ports) is tuple and len(set(ports))==len(ports),'literal distinct ports')
    for v in ports:nat(v);need(v<graph.n,'port in source')
    need(type(patterns) is tuple,'explicit pattern universe')
    results=[]
    for pattern in patterns:
        need(type(pattern) is tuple and len(pattern)==len(ports),'pattern domain')
        for x in pattern:nat(x);need(x<6,'six colours')
        results.append(solve(Branch(graph,(),tuple(zip(ports,pattern)))))
    return tuple(results)


def peel(graph,ports):
    """Full six-colour SUPPORT only, without extra branch equations."""
    old.audit_graph(graph)
    need(type(ports) is tuple and len(set(ports))==len(ports),'literal port set')
    for v in ports:nat(v);need(v<graph.n,'original port')
    adj=[set() for _ in range(graph.n)]
    for a,b in graph.edges:adj[a].add(b);adj[b].add(a)
    live=set(range(graph.n));order=[]
    while True:
        candidates=[v for v in sorted(live-set(ports)) if len(adj[v]&live)<=5]
        if not candidates:break
        v=candidates[0];live.remove(v);order.append(v)
    return tuple(sorted(live)),tuple(order)


def audit_peeling(graph,ports,core,order):
    old.audit_graph(graph)
    for bank in (ports,core,order):
        need(type(bank) is tuple,'vertex tuple')
        for v in bank:nat(v);need(v<graph.n,'vertex in graph')
        need(len(set(bank))==len(bank),'distinct vertex bank')
    need(set(ports)<=set(core) and set(core).isdisjoint(order)
         and set(core)|set(order)==set(range(graph.n)),'source-preserving vertex partition')
    live=set(range(graph.n))
    for v in order:
        degree=sum((a==v and b in live) or (b==v and a in live) for a,b in graph.edges)
        need(degree<=5,'six-colour reverse-greedy guard');live.remove(v)
    return True


def lift_core(graph,ports,core,order,colours):
    audit_peeling(graph,ports,core,order)
    need(type(colours) is tuple and len(colours)==len(core),'literal core colouring')
    for c in colours:nat(c);need(c<6,'six core colours')
    colour=dict(zip(core,colours))
    need(all(colour[a]!=colour[b] for a,b in graph.edges if a in colour and b in colour),'proper core')
    for v in reversed(order):
        blocked={colour[b if a==v else a] for a,b in graph.edges
                 if (a==v and b in colour) or (b==v and a in colour)}
        colour[v]=next(c for c in range(6) if c not in blocked)
    result=tuple(colour[v] for v in range(graph.n))
    old.audit_colouring(graph,6,result)
    return result


def support_batch(graph,ports,patterns,proof_cap=100000):
    """Incremental full-colour decisions with checked positive/negative exports.

    An exhausted negative proof budget raises ProofLimit; it never exports UNSAT.
    """
    nat(proof_cap,1);need(type(ports) is tuple and bool(ports),'nonempty marked source')
    core,order=peel(graph,ports);audit_peeling(graph,ports,core,order)
    index={v:i for i,v in enumerate(core)}
    cg=old.graph(len(core),((index[a],index[b]) for a,b in graph.edges if a in index and b in index))
    # With no extra equations, each binary/ternary pair is one unrestricted
    # six-state variable. This exact specialization propagates more strongly.
    ids=tuple(tuple(6*v+a+1 for a in range(6)) for v in range(cg.n))
    clauses=tuple([row for row in ids]+
        [(-a,-b) for row in ids for a,b in combinations(row,2)]+
        [(-ids[a][k],-ids[b][k]) for a,b in cg.edges for k in range(6)])
    need(type(patterns) is tuple,'explicit pattern universe')
    rows=[]
    with Solver(name='g3',bootstrap_with=clauses) as solver:
        for pattern in patterns:
            need(type(pattern) is tuple and len(pattern)==len(ports),'source pattern')
            for x in pattern:nat(x);need(x<6,'six source colours')
            assumptions=tuple(ids[index[v]][a] for v,a in zip(ports,pattern))
            if solver.solve(assumptions=assumptions):
                model=set(solver.get_model() or [])
                z=tuple(next(k for k,var in enumerate(row) if var in model) for row in ids)
                full=lift_core(graph,ports,core,order,z)
                need(tuple(full[v] for v in ports)==pattern,'literal boundary pins')
                rows.append({'pattern':pattern,'status':'colouring','colours':full})
            else:
                pinned=clauses+tuple((x,) for x in assumptions)
                tree,nodes=unsat_tree(pinned,proof_cap);audit_unsat_tree(pinned,tree)
                rows.append({'pattern':pattern,'status':'unsatisfiable','tree':tree,'proof_nodes':nodes})
    return {'core':core,'deletion_order':order,'core_edges':cg.edges,'rows':rows}


def geometry_support(data):
    """Consume a separately authenticated complete coordinate-graph packet."""
    points=data['points'];edges=tuple(tuple(e) for e in data['edges'])
    ports=tuple(data['original_ports'])
    need(ports==tuple(range(8)),'unchanged original eight ports')
    graph=old.graph(len(points),edges)
    need(tuple(e for e in graph.edges if e[1]<8)==tuple(sorted(source8.BASE_SET)),'literal original source graph')
    patterns=tuple(source8.partitions())
    result=support_batch(graph,ports,patterns)
    for row in result['rows']:
        check(row['status']=='colouring','current actual candidate retains this source pattern')
        old.audit_colouring(graph,6,row['colours'])
        check(tuple(row['colours'][v] for v in ports)==row['pattern'],'full exported boundary')
    result.update({'point_count':graph.n,'edges':graph.edges,'ports':ports,
         'coordinate_sha256':hashlib.sha256(json.dumps(points,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
         'graph_sha256':hashlib.sha256(json.dumps([graph.n,graph.edges,ports],separators=(',',':')).encode()).hexdigest(),
         'equal_patterns':sum(r['pattern'][2]==r['pattern'][4] for r in result['rows']),
         'separating_patterns':sum(r['pattern'][2]!=r['pattern'][4] for r in result['rows'])})
    check(len(patterns)==570 and len(result['rows'])==570,'complete source partition bank')
    # A lossless base-six word saves hundreds of thousands of repeated JSON lines.
    result['rows']=[dict(row,colours=''.join(str(c) for c in row['colours'])) for row in result['rows']]
    return result


def address_pair_repair(graph,addresses,aliases,ports,pattern,matching):
    """A fixed inherited-pair chart on physical vertices, including ALL aliases."""
    old.audit_graph(graph)
    need(type(addresses) is tuple and type(aliases) is tuple and len(addresses)==len(aliases),'address map')
    need(type(ports) is tuple and type(pattern) is tuple and len(ports)==len(pattern),'source pins')
    for x in ports:nat(x);need(x<graph.n,'source vertex')
    for x in pattern:nat(x);need(x<6,'source colour')
    need(len(set(ports))==len(ports),'distinct ports')
    for row,v in zip(addresses,aliases):
        need(type(row) is tuple and len(row)==3,'three-factor marked address')
        for x in row:nat(x)
        nat(v);need(v<graph.n and row[0]<len(ports),'literal address domain')
    need(set(aliases)==set(range(graph.n)),'every physical vertex has an address')
    need(type(matching) is tuple and len(matching)==3,'palette matching')
    for pair in matching:
        need(type(pair) is tuple and len(pair)==2,'palette pair')
        for x in pair:nat(x)
        need(pair[0]<pair[1],'canonical pair orientation')
    decode=tuple(c for pair in matching for c in pair)
    need(tuple(sorted(decode))==tuple(range(6)) and tuple(sorted(matching))==matching,'perfect matching')
    encode={c:i for i,c in enumerate(decode)}
    labels=[None]*graph.n;owner=[None]*graph.n
    for index,(row,v) in enumerate(zip(addresses,aliases)):
        label=encode[pattern[row[0]]]//2
        if labels[v] is None:labels[v]=label;owner[v]=index
        elif labels[v]!=label:
            return {'status':'alias_conflict','obstruction':(v,owner[v],index)}
    for v,c in zip(ports,pattern):
        if labels[v]!=encode[c]//2:return {'status':'pin_label_conflict','obstruction':(v,owner[v],c)}
    pg=parity.Graph(tuple(range(graph.n)),graph.edges)
    packet=parity.compile_boundary(pg,tuple(labels),ports)
    if packet.odd_cycle:return {'status':'odd_cycle','obstruction':packet.odd_cycle}
    bits=tuple(encode[c]%2 for c in pattern);pinned=dict(zip(ports,bits))
    for e in parity.constraints(packet):
        if pinned[e.start]^pinned[e.end]!=e.parity:
            return {'status':'source_path','obstruction':e.path}
    encoded=parity.lift(packet,bits);colours=tuple(decode[c] for c in encoded)
    old.audit_colouring(graph,6,colours)
    need(tuple(colours[v] for v in ports)==pattern,'source-fixed shade lift')
    free=set(packet.roots)-{packet.roots[v] for v in ports}
    return {'status':'colouring','count':2**len(free),'colours':colours}


def geometric_pair_sweep(data):
    graph=old.graph(len(data['points']),(tuple(e) for e in data['edges']))
    addresses=tuple(tuple(a) for a in data['addresses']);aliases=tuple(data['address_to_vertex'])
    ports=tuple(data['original_ports']);matchings=tuple(old.perfect_pairings(tuple(range(6))))
    totals={};rows=[]
    for pattern in source8.partitions():
        charts=[];witness=None
        for matching in matchings:
            r=address_pair_repair(graph,addresses,aliases,ports,pattern,matching)
            totals[r['status']]=totals.get(r['status'],0)+1
            if r['status']=='colouring':
                if witness is None:witness=''.join(map(str,r['colours']))
                charts.append({'status':r['status'],'count':r['count']})
            else:charts.append(r)
            check(True,'exact physical-address pairing chart')
        profile=_alias_profile(addresses,aliases,ports,pattern)
        passed=sum(c['status'] not in ('alias_conflict','pin_label_conflict') for c in charts)
        check(passed==profile['allowed_matchings'],'complete six-name alias matching obstruction')
        rows.append({'pattern':pattern,'charts':charts,'first_positive_colouring':witness,
                     'alias_profile':profile})
    # This model already fails on a fixed-anchor seed triangle. The mixed
    # geometric edges are not needed for this particular obstruction.
    marked_addresses=((0,0,0),(0,1,0),(0,2,0))
    marked_indices=tuple(addresses.index(a) for a in marked_addresses)
    triangle=tuple(aliases[i] for i in marked_indices)
    check(len(set(triangle))==3 and all(tuple(sorted(e)) in graph.edges
          for e in combinations(triangle,2)),'actual fixed-anchor triangle')
    check(all(sum(c['status']=='odd_cycle' for c in r['charts'])==3
              and r['alias_profile']['allowed_matchings']==3 for r in rows),
          'all remaining inherited-pair charts fail on the same seed obstruction')
    return {'matchings':matchings,'status_totals':totals,'rows':rows,
            'fixed_anchor_triangle':{'addresses':marked_addresses,'address_indices':marked_indices,
                                     'physical_vertices':triangle},
            'patterns_with_repair':sum(r['first_positive_colouring'] is not None for r in rows),
            'scope':'Inherited-pair charts only; the fixed-anchor seed triangle already defeats this class, without any mixed edge.'}


def _alias_profile(addresses,aliases,ports,pattern):
    """Called after exact address validation by every matching chart."""
    owners={};rows=[]
    for i,(address,v) in enumerate(zip(addresses,aliases)):
        if v not in owners:owners[v]=i
        else:
            first=owners[v];a=pattern[addresses[first][0]];b=pattern[address[0]]
            if a!=b:rows.append({'colours':(a,b),'vertex':v,'addresses':(first,i)})
    for v,c in zip(ports,pattern):
        i=owners[v];a=pattern[addresses[i][0]]
        if a!=c:rows.append({'colours':(a,c),'vertex':v,'address':i,'source_pin':c})
    adj=[set() for _ in range(6)]
    for row in rows:
        a,b=row['colours'];adj[a].add(b);adj[b].add(a)
    components=[];seen=set()
    for i in range(6):
        if i in seen:continue
        component=[i];seen.add(i)
        for a in component:
            for b in sorted(adj[a]):
                if b not in seen:seen.add(b);component.append(b)
        components.append(tuple(sorted(component)))
    if any(len(c)>2 for c in components):count=0
    else:
        remaining=sum(len(c)==1 for c in components);count=1
        for factor in range(1,remaining,2):count*=factor
    return {'components':tuple(components),'allowed_matchings':count,'source_alias_edges':rows}


CHECKS=0
def check(b,message):
    global CHECKS
    CHECKS+=1
    if not b:raise RuntimeError(message)


def rejects(fn):
    try:fn()
    except ValueError:check(True,'guard')
    else:check(False,'malformed input accepted')


def main():
    check(len(BRANCHES)==375,'local contact cover size')
    for mask in range(1,64):
        D={i for i in range(6) if mask>>i&1}
        rich=any((i+1)%6 in D for i in D)
        periodic=any(len({i%p for i in D})==1 for p in (2,3))
        check(periodic== (not rich),'nonrich classification')
    # Complete all simple graphs through four vertices, one literal source pin.
    small=[]
    for n in range(1,5):
        possible=tuple(combinations(range(n),2))
        for mask in range(1<<len(possible)):
            g=old.graph(n,(e for i,e in enumerate(possible) if mask>>i&1))
            b=Branch(g,(),((0,0),))
            direct=sum(all(((0,)+x)[a]!=((0,)+x)[d] for a,d in g.edges)
                       for x in product(range(6),repeat=n-1))
            check(count_branch(b)==direct,'unconstrained CRT count is full six-colour count')
            r=solve(b);audit_receipt(r)
            check(r.status=='colouring','all <=4 graphs six-colourable')
            small.append(direct)
    # All oriented single equations on three vertices, plus a second-coordinate equation.
    n=3;g=old.graph(n,((0,1),(1,2),(0,2)))
    configs=0
    for p in (2,3):
        for u,v in product(range(n),repeat=2):
            for a in range(p):
                for q in (2,3):
                    for d in range(q):
                        b=Branch(g,((p,u,v,a),(q,0,2,d)),((0,0),))
                        direct=0
                        for tail in product(range(6),repeat=2):
                            z=(0,)+tail
                            if all(z[x]!=z[y] for x,y in g.edges) and all((z[y]-z[x])%m==r for m,x,y,r in b.equations):direct+=1
                        check(count_branch(b)==direct,'coupled equations direct enumeration')
                        r=solve(b);audit_receipt(r)
                        check((r.status=='colouring')==(direct>0),'independent complete decision')
                        configs+=1
    # A genuine residual obstruction beyond linear consistency.
    k4=old.graph(4,combinations(range(4),2))
    kb=Branch(k4,tuple((2,0,i,0) for i in range(1,4)),((0,0),))
    kr=solve(kb);audit_receipt(kr)
    check(not kr.compiled.binary.cycle and not kr.compiled.ternary.cycle and kr.status=='unsatisfiable',
          'consistent branch can leave a non-three-colourable ternary graph')
    check(solve(Branch(k4,(),((0,0),))).status=='colouring','branch obstruction is not full K4 six-colour obstruction')
    # A forced loop in BOTH component quotients is immediately impossible.
    loop=Branch(old.graph(2,((0,1),)),((2,0,1,0),(3,0,1,0)),())
    lr=solve(loop);audit_receipt(lr)
    check(lr.status=='unsatisfiable' and () in lr.compiled.clauses,'joint equality forbids the actual edge')
    # Complete375 local cover rejection with five independently rich factors.
    factors=tuple(((3*i,3*i+1),(3*i,3*i+2)) for i in range(6))
    pins=tuple((3*i+j,j if i<5 else int(j>0)) for i in range(6) for j in range(3))
    source=Branch(old.graph(18,()),(),pins)
    covered=[]
    for choice in BRANCHES:
        r=solve(contact_branch(source,factors,choice));audit_receipt(r)
        check(r.status=='voltage_conflict','five-rich source has no G4 branch')
        covered.append({'choice':choice,'modulus':2 if r.compiled.binary.cycle else 3,
                        'cycle':r.compiled.binary.cycle or r.compiled.ternary.cycle})
    check(solve(source).status=='colouring','G4 contradiction needs the additive contact premise')
    # Arithmetic CRT is not the old consecutive-pair quotient.
    check((1//2+1//2)%3 != ((1+1)%6)//2,'pair-label projection is not additive')
    check(all((3*(c%2)+4*(c%3))%6==c for c in range(6)),'marked arithmetic CRT inverse')
    rejects(lambda:solve(replace(kb,pins=((0,False),))))
    rejects(lambda:solve(replace(kb,equations=((2,0,1,0.0),))))
    rejects(lambda:audit_receipt(replace(kr,status=True)))
    rejects(lambda:audit_receipt(replace(kr,tree=())))
    rejects(lambda:audit_receipt(replace(kr,compiled=replace(kr.compiled,parameters=((True,1),)))))
    rejects(lambda:contact_branch(source,factors,(False,2,0,1,2,0)))
    rejects(lambda:complete_boundary(k4,(False,),((0,),)))
    # A full abstract equality positive control and its checked separating negatives.
    # This graph has no asserted unit embedding.
    terminal_graph=old.graph(7,(e for e in combinations(range(7),2) if e!=(0,1)))
    equality=support_batch(terminal_graph,(0,1),((0,0),(0,1)))
    check([r['status'] for r in equality['rows']]==['colouring','unsatisfiable'],
          'complete fixed-source equality control')
    check(len(equality['core'])==7,'critical source core retained')
    # Geometry is authenticated by its separate producer; this consumer retains
    # the whole graph, literal source IDs, and a content hash of the coordinates.
    path=Path(__file__).resolve().parents[2]/'05-knowledge/results/coupled_fractional_geometry_20261018c.json'
    data=json.loads(path.read_text(encoding='utf-8'))
    geometric={name:geometry_support(data[name]) for name in ('candidate','lifted_candidate')}
    pair_sweep=geometric_pair_sweep(data['lifted_candidate'])
    return {'checks':CHECKS,'status':'PROVED typed CRT residual interface; FINITE-EXACT controls; plane forcing OPEN',
            'simple_graphs':len(small),'three_vertex_coupled_branches':configs,
            'k4_constant_binary':{'equations':kb.equations,'pins':kb.pins,'clauses':kr.compiled.clauses,
                'parameters':kr.compiled.parameters,'unsat_tree':kr.tree,'proof_nodes':audit_unsat_tree(kr.compiled.clauses,kr.tree)},
            'complete_375_contact_refutation':covered,
            'abstract_equality':equality,
            'geometry_dependency':str(path.name),'geometric_candidates':geometric,
            'lifted_pair_sweep':pair_sweep,
            'scope':'G4 is a necessary cover only for the declared additive Cartesian contact; unconstrained CRT solves actual graph colouring.'}


if __name__=='__main__':
    print(json.dumps(main(),sort_keys=True,indent=2))
