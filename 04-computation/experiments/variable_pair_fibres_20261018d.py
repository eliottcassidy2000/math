"""Exact variable ternary labels with pinned bipartite shade fibres.

All learned cuts carry an actual source walk. A bounded unfinished search
returns pending, never UNSAT. No geometric assertion follows from a graph.
"""
from collections import deque
from dataclasses import dataclass, replace
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
from pysat.solvers import Solver
import snark_boundary_transfer_20261018a as old
import coupled_support_branches_20261018c as prior
import exceptional_factor_support_20261018a as source8

need, nat = old.require, old.nat


@dataclass(frozen=True)
class Source:
    graph: old.Graph
    ports: tuple


@dataclass(frozen=True)
class Cut:
    kind: str
    walk: tuple


@dataclass(frozen=True)
class Fibre:
    labels: tuple
    cuts: tuple
    colours: tuple
    free_components: int


@dataclass(frozen=True)
class Result:
    source: Source
    pattern: tuple
    status: str
    cuts: tuple
    labels: tuple
    colours: tuple
    tree: tuple
    iterations: int
    free_components: int


def audit_source(s):
    need(type(s) is Source,'exact Source');old.audit_graph(s.graph)
    need(type(s.ports) is tuple,'literal port tuple')
    for v in s.ports:nat(v);need(v<s.graph.n,'port in graph')
    need(len(set(s.ports))==len(s.ports),'distinct ports')
    return s


def audit_pattern(s,pattern):
    audit_source(s);need(type(pattern) is tuple and len(pattern)==len(s.ports),'source pattern')
    for c in pattern:nat(c);need(c<6,'six source colours')
    return pattern


def path_between(u,v,parents):
    left=[];x=u
    while x!=-1:left.append(x);x=parents[x]
    right=[];x=v
    while x not in left:right.append(x);x=parents[x]
    return tuple(left[:left.index(x)+1]+list(reversed(right)))


def canonical_cycle(walk):
    row=walk[:-1]
    variants=[r[i:]+r[:i] for r in (row,tuple(reversed(row))) for i in range(len(row))]
    best=min(variants)
    return best+(best[0],)


def audit_cut(s,c):
    audit_source(s)
    return _audit_cut(s,c,frozenset(s.graph.edges))


def _audit_cut(s,c,edges):
    """Internal fast path after exact source validation; never equality-cached."""
    need(type(c) is Cut and type(c.kind) is str,'exact Cut')
    need(type(c.walk) is tuple,'source walk tuple')
    for v in c.walk:nat(v);need(v<s.graph.n,'walk vertex')
    need(all(tuple(sorted(e)) in edges for e in zip(c.walk,c.walk[1:])),
         'every walk edge is an actual source edge')
    if c.kind=='cycle':
        need(len(c.walk)>=4 and c.walk[0]==c.walk[-1] and len(c.walk)%2==0,'odd closed cycle')
        need(len(set(c.walk[:-1]))==len(c.walk)-1,'simple cycle')
        need(c.walk==canonical_cycle(c.walk),'canonical cycle orientation')
    elif c.kind=='path':
        need(len(c.walk)>=2 and len(set(c.walk))==len(c.walk),'simple nonempty path')
        need(c.walk[0] in s.ports and c.walk[-1] in s.ports,'retained pinned endpoints')
        need(c.walk<=tuple(reversed(c.walk)),'canonical path orientation')
    else:need(False,'unknown cut kind')
    return c


def shade_fibre(s,pattern,labels):
    audit_pattern(s,pattern)
    need(type(labels) is tuple and len(labels)==s.graph.n,'one label per physical vertex')
    for x in labels:nat(x);need(x<3,'ternary labels')
    need(all(labels[v]==c%3 for v,c in zip(s.ports,pattern)),'literal source ternary labels')
    adj=[[] for _ in labels]
    for a,b in s.graph.edges:
        if labels[a]==labels[b]:adj[a].append(b);adj[b].append(a)
    roots=[-1]*len(labels);bits=[0]*len(labels);parents=[-1]*len(labels);cycles=set()
    for start in range(len(labels)):
        if roots[start]>=0:continue
        roots[start]=start;q=deque([start])
        while q:
            u=q.popleft()
            for v in adj[u]:
                if roots[v]<0:
                    roots[v]=start;bits[v]=bits[u]^1;parents[v]=u;q.append(v)
                elif bits[v]==bits[u]:
                    cycles.add(Cut('cycle',canonical_cycle(path_between(u,v,parents)+(u,))))
    if cycles:
        cuts=tuple(sorted(cycles,key=lambda c:(len(c.walk),c.walk)))
        edges=frozenset(s.graph.edges)
        for c in cuts:_audit_cut(s,c,edges)
        return Fibre(labels,cuts,(),0)
    first={};flips={};paths=set()
    for v,c in zip(s.ports,pattern):
        r=roots[v];required=(c%2)^bits[v]
        if r not in first:first[r]=v;flips[r]=required
        elif flips[r]!=required:
            path=path_between(first[r],v,parents);path=min(path,tuple(reversed(path)))
            paths.add(Cut('path',path))
    if paths:
        cuts=tuple(sorted(paths,key=lambda c:(len(c.walk),c.walk)))
        edges=frozenset(s.graph.edges)
        for c in cuts:_audit_cut(s,c,edges)
        return Fibre(labels,cuts,(),0)
    shades=tuple(bits[v]^flips.get(roots[v],0) for v in range(len(labels)))
    colours=tuple((3*b+4*t)%6 for b,t in zip(shades,labels))
    old.audit_colouring(s.graph,6,colours)
    need(tuple(colours[v] for v in s.ports)==pattern,'literal pinned lift')
    return Fibre(labels,(),colours,len(set(roots)-set(first)))


def ids(s):
    return tuple(tuple(3*v+t+1 for t in range(3)) for v in range(s.graph.n))


def shade_ids(s):
    return {v:3*s.graph.n+i+1 for i,v in enumerate(s.ports)}


def base_clauses(s):
    audit_source(s);rows=ids(s)
    return tuple(list(rows)+[(-a,-b) for row in rows for a,b in combinations(row,2)])


def cut_clauses(s,c):
    audit_cut(s,c)
    return _cut_clauses(s,c)


def _cut_clauses(s,c):
    rows=ids(s);clauses=[]
    if c.kind=='cycle':
        return tuple(tuple(-rows[v][t] for v in c.walk[:-1]) for t in range(3))
    u,v=c.walk[0],c.walk[-1];bits=shade_ids(s);parity=(len(c.walk)-1)%2
    for t in range(3):
        labels=tuple(-rows[w][t] for w in c.walk)
        for a,b in product(range(2),repeat=2):
            if a^b!=parity:
                clauses.append(labels+((-bits[u] if a else bits[u]),(-bits[v] if b else bits[v])))
    return tuple(clauses)


def assumptions(s,pattern):
    audit_pattern(s,pattern);rows=ids(s);bits=shade_ids(s)
    return tuple(x for v,c in zip(s.ports,pattern)
                 for x in (rows[v][c%3],bits[v] if c%2 else -bits[v]))


def audit_result(r):
    need(type(r) is Result,'exact Result');audit_pattern(r.source,r.pattern)
    need(type(r.status) is str and type(r.cuts) is tuple,'typed result banks')
    edges=frozenset(r.source.graph.edges)
    for c in r.cuts:_audit_cut(r.source,c,edges)
    need(len(set(r.cuts))==len(r.cuts),'distinct retained cuts')
    for bank in (r.labels,r.colours):
        need(type(bank) is tuple,'tuple colouring bank')
        for x in bank:nat(x)
    need(type(r.tree) is tuple,'tree tuple');nat(r.iterations);nat(r.free_components)
    if r.status=='colouring':
        need(not r.tree,'positive receipt has no negative tree')
        f=shade_fibre(r.source,r.pattern,r.labels)
        need(not f.cuts and f.colours==r.colours and f.free_components==r.free_components,'exact shade fibre')
    elif r.status=='unsatisfiable':
        need(not r.labels and not r.colours and r.free_components==0,'negative fields')
        clauses=base_clauses(r.source)+tuple(row for c in r.cuts for row in _cut_clauses(r.source,c))
        clauses+=tuple((x,) for x in assumptions(r.source,r.pattern))
        prior.audit_unsat_tree(clauses,r.tree)
    elif r.status=='pending':
        need(not r.labels and not r.colours and not r.tree and r.free_components==0,'pending is no conclusion')
    else:need(False,'unknown result')
    return r


class Engine:
    """Source-fixed lazy label solver; cycle/path cuts survive new source colours."""
    def __init__(self,source,initial_cuts=()):
        audit_source(source);need(type(initial_cuts) is tuple,'initial cuts tuple')
        self.source=source;self.cuts=[];self.seen=set();self.edges=frozenset(source.graph.edges)
        self.solver=Solver(name='g3',bootstrap_with=base_clauses(source))
        for c in initial_cuts:self.add(c)
    def close(self):self.solver.delete()
    def add(self,c):
        _audit_cut(self.source,c,self.edges)
        if c in self.seen:return False
        for row in _cut_clauses(self.source,c):self.solver.add_clause(row)
        self.seen.add(c);self.cuts.append(c);return True
    def solve(self,pattern,max_rounds=10000,proof_cap=100000):
        s=self.source;a=assumptions(s,pattern);nat(max_rounds);nat(proof_cap,1)
        for iteration in range(max_rounds):
            if not self.solver.solve(assumptions=a):
                clauses=base_clauses(s)+tuple(row for c in self.cuts for row in _cut_clauses(s,c))
                clauses+=tuple((x,) for x in a)
                try:tree,_=prior.unsat_tree(clauses,proof_cap)
                except prior.ProofLimit:
                    return Result(s,pattern,'pending',tuple(self.cuts),(),(),(),iteration,0)
                r=Result(s,pattern,'unsatisfiable',tuple(self.cuts),(),(),tree,iteration,0)
                return audit_result(r)
            model=set(self.solver.get_model() or [])
            labels=tuple(next(t for t,var in enumerate(row) if var in model) for row in ids(s))
            f=shade_fibre(s,pattern,labels)
            if not f.cuts:
                return Result(s,pattern,'colouring',tuple(self.cuts),labels,f.colours,(),iteration+1,f.free_components)
            need(any(c not in self.seen for c in f.cuts),'each failed label assignment is eliminated')
            for c in f.cuts:self.add(c)
        return Result(s,pattern,'pending',tuple(self.cuts),(),(),(),max_rounds,0)


def solve_fibres(graph,ports,pattern,max_rounds=10000,initial_cuts=(),proof_cap=100000):
    e=Engine(Source(graph,ports),initial_cuts)
    try:return e.solve(pattern,max_rounds,proof_cap)
    finally:e.close()


def encode_result(r):
    audit_result(r)
    return {'vertices':r.source.graph.n,'edges':r.source.graph.edges,'ports':r.source.ports,
            'pattern':r.pattern,'status':r.status,
            'cuts':[{'kind':c.kind,'walk':c.walk} for c in r.cuts],
            'labels':''.join(map(str,r.labels)),'colours':''.join(map(str,r.colours)),
            'tree':r.tree,'iterations':r.iterations,'free_components':r.free_components}


CHECKS=0
def check(b,message):
    global CHECKS
    CHECKS+=1
    if not b:raise RuntimeError(message)


def rejects(fn):
    try:fn()
    except ValueError:check(True,'exact guard')
    else:check(False,'malformed source accepted')


def small_controls():
    cases=0
    for n in range(1,5):
        possible=tuple(combinations(range(n),2))
        for mask in range(1<<len(possible)):
            g=old.graph(n,(e for i,e in enumerate(possible) if mask>>i&1))
            ports=tuple(range(min(2,n)));pattern=(0,3)[:len(ports)];s=Source(g,ports)
            total=0
            for tail in product(range(3),repeat=n-len(ports)):
                labels=tuple(c%3 for c in pattern)+tail;f=shade_fibre(s,pattern,labels)
                direct=0
                for bits in product(range(2),repeat=n-len(ports)):
                    shades=tuple(c%2 for c in pattern)+bits
                    colours=tuple((3*b+4*t)%6 for b,t in zip(shades,labels))
                    direct+=all(colours[a]!=colours[b] for a,b in g.edges)
                check(direct==(0 if f.cuts else 2**f.free_components),'complete fixed-label fibre count')
                for c in f.cuts:
                    audit_cut(s,c)
                    # Full original six-colour solutions obey every learned cut.
                    for free in product(range(6),repeat=n-len(ports)):
                        z=pattern+free
                        if not all(z[a]!=z[b] for a,b in g.edges):continue
                        values={3*v+t+1:z[v]%3==t for v in range(n) for t in range(3)}
                        values.update({v:z[p]%2==1 for p,v in shade_ids(s).items()})
                        check(all(any(values[abs(lit)]==(lit>0) for lit in row) for row in cut_clauses(s,c)),
                              'learned cut preserves every full source colouring')
                total+=direct;cases+=1
            result=solve_fibres(g,ports,pattern);audit_result(result)
            check((result.status=='colouring')==(total>0),'full lazy decision versus enumeration')
    # A fixed cycle basis misses a monochromatic nonbasis triangle in K4.
    g=old.graph(4,combinations(range(4),2));s=Source(g,())
    labels=(1,0,0,0);basis=((0,1,2,0),(0,1,3,0),(0,2,3,0))
    check(all(len({labels[v] for v in C})>1 for C in basis),'fixed basis passes')
    f=shade_fibre(s,(),labels)
    check(Cut('cycle',(1,2,3,1)) in f.cuts,'nonbasis actual odd triangle rejected')
    # Same ternary pins, incompatible versus compatible endpoint shade.
    p=old.graph(3,((0,1),(1,2)));s=Source(p,(0,2))
    bad=shade_fibre(s,(0,3),(0,0,0));good=shade_fibre(s,(0,0),(0,0,0))
    check(bad.cuts==(Cut('path',(0,1,2)),) and not good.cuts,'pinned path parity is essential')
    engine=Engine(s,bad.cuts)
    try:
        check(engine.solve((0,0)).status=='colouring','conditional path cut survives another source pattern')
        check(engine.solve((0,3)).status=='colouring','variable middle label repairs bad path')
    finally:engine.close()
    k7=old.graph(7,(e for e in combinations(range(7),2) if e!=(0,1)))
    positive=solve_fibres(k7,(0,1),(0,0));negative=solve_fibres(k7,(0,1),(0,1))
    audit_result(positive);audit_result(negative)
    check(positive.status=='colouring' and negative.status=='unsatisfiable','abstract equality exact control')
    pending=solve_fibres(k7,(0,1),(0,1),max_rounds=0);audit_result(pending)
    check(pending.status=='pending','zero budget is not UNSAT')
    rejects(lambda:shade_fibre(Source(p,(0,2)),(0,3),(False,0,0)))
    rejects(lambda:audit_cut(Source(p,(0,2)),Cut('path',(0,2))))
    rejects(lambda:audit_result(replace(negative,tree=())))
    rejects(lambda:solve_fibres(k7,(False,1),(0,1)))
    rejects(lambda:solve_fibres(k7,(0,1),(0,1),max_rounds=1.0))
    return {'fibre_cases':cases,'fixed_basis_hostile':{'labels':labels,'basis':basis,'missed_cycle':f.cuts[0].walk},
            'pinned_path':{'walk':bad.cuts[0].walk,'bad_pattern':(0,3),'good_pattern':(0,0)},
            'abstract_equal':encode_result(positive),'abstract_separating':encode_result(negative)}


def _graph(data):
    return old.graph(len(data['points']),(tuple(e) for e in data['edges']))


def _hash(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def _prepare_core(graph,ports):
    core,order=prior.peel(graph,ports);prior.audit_peeling(graph,ports,core,order)
    index={v:i for i,v in enumerate(core)}
    reduced=old.graph(len(core),((index[a],index[b]) for a,b in graph.edges if a in index and b in index))
    adj=[[] for _ in range(graph.n)]
    for a,b in graph.edges:adj[a].append(b);adj[b].append(a)
    return core,order,index,reduced,adj


def _lift_prechecked(graph,core,order,adj,colours):
    """Private loop after _prepare_core authenticated the complete deletion order."""
    values=[-1]*graph.n
    for v,c in zip(core,colours):values[v]=c
    for v in reversed(order):
        blocked={values[w] for w in adj[v] if values[w]>=0}
        values[v]=next(c for c in range(6) if c not in blocked)
    result=tuple(values);old.audit_colouring(graph,6,result)
    return result


def geometry_controls(progress=None):
    """Actual graph control; exact metric authentication belongs to its producer."""
    root=Path(__file__).resolve().parents[2]/'05-knowledge/results'
    previous=json.loads((root/'coupled_fractional_geometry_20261018c.json').read_text(encoding='utf-8'))['lifted_candidate']
    saved_path=root/'coupled_support_branches_20261018c.json'
    saved_text=saved_path.read_text(encoding='utf-8')
    saved=json.loads(saved_text)['geometric_candidates']['lifted_candidate']
    graph=_graph(previous);ports=tuple(range(8));patterns=tuple(source8.partitions())
    need(tuple(tuple(r['pattern']) for r in saved['rows'])==patterns,'complete unchanged old source bank')
    need(saved['point_count']==graph.n and tuple(map(tuple,saved['edges']))==graph.edges,'old word source binding')
    core,order,index,reduced,adj=_prepare_core(graph,ports)
    engine=Engine(Source(reduced,tuple(index[v] for v in ports)))
    selected=(0,next(i for i,p in enumerate(patterns) if p[2]==p[4]),len(patterns)-1)
    old_rows=[]
    try:
        for i in selected:
            pattern=patterns[i];before=len(engine.cuts);r=engine.solve(pattern,max_rounds=2000)
            check(r.status=='colouring','selected old variable-label witness')
            full=_lift_prechecked(graph,core,order,adj,r.colours)
            f=shade_fibre(Source(graph,ports),pattern,tuple(c%3 for c in full))
            check(not f.cuts,'full lifted label fibre')
            full=f.colours
            different=[j for j,(address,v) in enumerate(zip(previous['addresses'],previous['address_to_vertex']))
                       if full[v]%3!=pattern[address[0]]%3]
            physical={previous['address_to_vertex'][j] for j in different}
            old_rows.append({'index':i,'pattern':pattern,'iterations':r.iterations,
                'new_cuts':len(engine.cuts)-before,'colours':''.join(map(str,full)),
                'free_components':f.free_components,'migrated_addresses':len(different),
                'migrated_physical_vertices':len(physical)})
        old_cuts=tuple(Cut(c.kind,tuple(core[v] for v in c.walk)) for c in engine.cuts)
        for c in old_cuts:audit_cut(Source(graph,ports),c)
    finally:engine.close()
    path=root/'surviving_repair_geometry_20261018d.json'
    data=json.loads(path.read_text(encoding='utf-8'))['candidate'];new=_graph(data)
    need(new.n>=graph.n and data['points'][:graph.n]==previous['points'],'literal old-coordinate inclusion')
    need(tuple(e for e in new.edges if e[1]<graph.n)==graph.edges,'full old induced graph retained')
    fixed=tuple(range(graph.n));nc,no,ni,ng,na=_prepare_core(new,fixed)
    source=Source(ng,tuple(ni[v] for v in fixed));engine=Engine(source)
    rows=[]
    try:
        for i,old_row in enumerate(saved['rows']):
            word=old_row['colours'];need(type(word) is str and len(word)==graph.n and set(word)<=set('012345'),'old literal colour word')
            pinned=tuple(map(int,word));old.audit_colouring(graph,6,pinned)
            need(tuple(pinned[v] for v in ports)==patterns[i],'old boundary source authentication')
            before=len(engine.cuts);r=engine.solve(pinned,max_rounds=2000)
            if r.status!='colouring':
                rows.append({'index':i,'pattern':patterns[i],'status':r.status,'receipt':encode_result(r)})
                break
            full=_lift_prechecked(new,nc,no,na,r.colours)
            check(full[:graph.n]==pinned,'every old physical colour is unchanged')
            f=shade_fibre(Source(new,fixed),pinned,tuple(c%3 for c in full))
            check(not f.cuts and f.colours[:graph.n]==pinned,'exact full new fibre with all old pins')
            full=f.colours
            rows.append({'index':i,'pattern':patterns[i],'status':'colouring','iterations':r.iterations,
                         'new_cuts':len(engine.cuts)-before,'free_components':f.free_components,
                         'colours':''.join(map(str,full))})
            if progress and (i+1)%50==0:progress(i+1,len(engine.cuts))
        full_cuts=tuple(Cut(c.kind,tuple(nc[v] for v in c.walk)) for c in engine.cuts)
        fs=Source(new,fixed);edge_set=frozenset(new.edges)
        for c in full_cuts:_audit_cut(fs,c,edge_set)
    finally:engine.close()
    complete=len(rows)==len(patterns) and all(r['status']=='colouring' for r in rows)
    return {'old_variable_labels':{'point_count':graph.n,'source_hash':_hash([graph.n,graph.edges,ports]),
                'selected':old_rows,'cuts':[{'kind':c.kind,'walk':c.walk} for c in old_cuts]},
            'literal_old_extensions':{'status':'all_saved_words_extend' if complete else 'partial',
                'old_word_bank_lf_sha256':hashlib.sha256(saved_text.replace('\r\n','\n').encode()).hexdigest(),
                'geometry_dependency':path.name,'coordinate_sha256':_hash(data['points']),
                'graph_sha256':_hash([new.n,new.edges,fixed]),'point_count':new.n,'edges':new.edges,
                'fixed_vertices':fixed,'core_vertices':nc,'core_edges':ng.edges,'deletion_order':no,
                'rows':rows,'cuts':[{'kind':c.kind,'walk':c.walk} for c in full_cuts],
                'cycle_cuts':sum(c.kind=='cycle' for c in full_cuts),
                'pinned_path_cuts':sum(c.kind=='path' for c in full_cuts),
                'scope':'Each row retains all 714 old physical colours, not merely the original eight-port pattern.'}}


if __name__=='__main__':
    controls=small_controls();geometry=geometry_controls()
    print(json.dumps({'checks':CHECKS,'controls':controls,'geometry':geometry},sort_keys=True,indent=2))
