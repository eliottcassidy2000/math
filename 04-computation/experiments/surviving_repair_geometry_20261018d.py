"""Bounded source-fixed true-cubic insertions into the frozen 714-point graph."""
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations
import hashlib
import json
import sys
import time

from coupled_fractional_geometry_20261018c import Search, relative_core, independent_edges
from cubic_boundary_search_20261018a import encode_element
from escaped_ports_20261017d import Field


ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'05-knowledge/results/coupled_fractional_geometry_20261018c.json'


def need(ok,message):
    if not ok: raise ValueError(message)


def decode_element(row,n=48):
    need(type(n) is int and n>=1,'positive exact field dimension')
    need(type(row) is list,'literal sparse element')
    out=[Q(0)]*n; last=-1
    for term in row:
        need(type(term) is list and len(term)==3 and all(type(v) is int for v in term),
             'exact sparse coefficient')
        i,a,b=term
        need(last<i<n and a!=0 and b>0,'ordered nonzero coefficient with positive denominator')
        out[i]=Q(a,b); last=i
    return tuple(out)


class Engine:
    def __init__(self):
        self.s=Search();self.ctx=self.s.ctx;self.f=self.s.f
        self.payload=json.loads(OLD.read_text(encoding='utf8'))
        self.old=self.payload['lifted_candidate']
        self.points=tuple(decode_element(x,self.f.n) for x in self.old['points'])
        self.edges=tuple(tuple(e) for e in self.old['edges'])
        self.seed=tuple(decode_element(x,self.f.n) for x in self.old['spindle'])
        self.shifts=tuple(decode_element(x,self.f.n) for x in self.old['fractional_shifts'])
        self.q=decode_element(self.old['phase'],self.f.n)
        self.gens=tuple(decode_element(x,self.f.n) for x in self.old['fractional_generators'])
        self.u=self.ctx.cubic_unit
        need(self.points[:8]==self.ctx.ports,'literal old source')
        need(len(self.points)==714 and len(self.edges)==2914,'fixed inherited geometry')
        need(self.ctx.norm(self.u)==self.f.one,'unit cubic rotation')
        need(any(x for x,b in zip(self.u,self.f.basis) if b[2]),'genuine cubic coefficient')
        self.R=tuple(dict.fromkeys(self.f.add(c,a) for c in self.ctx.ports for a in self.seed))
        self.k=Field(self.f.polys[:2])
        self.kbar=self.k.substitution((self.k.power(self.k.gen(0),14),self.k.scale(self.k.gen(1),-1)))

    def coefficients(self,x):
        self.ctx.element(x)
        return tuple(tuple(x[self.f.index[(a,b,c)]] for a,b in self.k.basis) for c in range(3))

    def embed_k(self,x):
        need(type(x) is tuple and len(x)==self.k.n and all(type(v) is Q for v in x),'exact K coefficient')
        out=list(self.f.zero)
        for (a,b),v in zip(self.k.basis,x):out[self.f.index[(a,b,0)]]=v
        return tuple(out)

    def k_neighbour(self,x):
        """Exact complete K-unit-neighbour reader for x in F minus K."""
        self.ctx.element(x);k=self.k
        t=self.coefficients(x);bar=self.coefficients(self.ctx.bar(x));N=self.coefficients(self.ctx.norm(x))
        need(t[1]!=k.zero or t[2]!=k.zero,'exterior point required')
        det=k.sub(k.mul(bar[1],t[2]),k.mul(bar[2],t[1]))
        # A singular two-row system makes bar(x)=lambda*x+mu over K.
        # A K-neighbour would then make x satisfy a nonzero quadratic,
        # impossible for an exterior element of this cubic extension.
        if det==k.zero:return None
        inv=k.inv(det)
        a=k.mul(k.sub(k.mul(N[1],t[2]),k.mul(N[2],t[1])),inv)
        ab=k.mul(k.sub(k.mul(bar[1],N[2]),k.mul(bar[2],N[1])),inv)
        if k.apply(a,self.kbar)!=ab:return None
        out=self.embed_k(a)
        return out if self.ctx.norm(self.f.sub(x,out))==self.f.one else None

    def rebased_shifts(self):
        f=self.f;z=f.mul(self.u,self.q)
        return tuple(dict.fromkeys((f.zero,)+tuple(f.scale(f.sub(self.transform(z,s),z),2) for s in self.gens)))

    def transform(self,z,s):
        self.ctx.element(z);self.ctx.element(s)
        f=self.f
        need(self.ctx.norm(z)==f.one and self.ctx.norm(s)==f.one,'unit fractional inputs')
        den=f.sub(f.scalar(2),f.mul(self.ctx.bar(s),z))
        need(den!=f.zero,'no fractional pole')
        out=f.mul(f.sub(f.scale(z,2),s),f.inv(den))
        need(self.ctx.norm(out)==f.one,'unit fractional output')
        return out

    def recipe(self,name):
        f=self.f;u=self.u
        if name=='coupled-seed-phase':
            return self.recipe('rotated-seed')+self.recipe('cubic-rebased-phase')
        if name=='rotated-seed':
            return tuple(f.add(c,f.add(f.mul(u,a),t))
                         for c in self.ctx.ports for a in self.seed for t in self.shifts)
        z=f.mul(u,self.q)
        if name=='cubic-phase':
            phases=tuple(dict.fromkeys((z,)+tuple(self.transform(z,s) for s in self.gens)))
            shifts=tuple(f.scale(f.sub(x,self.q),2) for x in phases)
        elif name=='cubic-rebased-phase':
            shifts=self.rebased_shifts()
        else: raise ValueError('unknown bounded recipe')
        return tuple(f.add(r,t) for r in self.R for t in shifts)

    def projected_recipe(self,name):
        """The proposed old-point target at every literal new address."""
        f=self.f
        if name=='coupled-seed-phase':
            return self.projected_recipe('rotated-seed')+self.projected_recipe('cubic-rebased-phase')
        if name=='rotated-seed':
            return tuple(f.add(c,f.add(a,t)) for c in self.ctx.ports for a in self.seed for t in self.shifts)
        need(name in ('cubic-phase','cubic-rebased-phase'),'declared projection recipe')
        raw=self.recipe(name); layers=len(raw)//len(self.R)
        need(len(raw)==layers*len(self.R),'complete rectangular address list')
        return tuple(r for r in self.R for _ in range(layers))

    def build(self,name):
        raw=self.recipe(name)
        pts=list(self.points);index={x:i for i,x in enumerate(pts)};aliases=[]
        for x in raw:
            if x not in index:index[x]=len(pts);pts.append(x)
            aliases.append(index[x])
        self.ctx.geo.pl=[p for p in self.s.good_places
                         if all(c.denominator%p[0] for x in pts for c in x)]
        need(self.ctx.geo.pl,'authenticated denominator-safe split place')
        edges=tuple(self.ctx.geo.unit_edges(pts)); n0=len(self.points)
        need(tuple(e for e in edges if e[1]<n0)==self.edges,'exact old induced inclusion')
        exterior,peel=relative_core(len(pts),edges)
        relative,oldpeel=relative_core(len(pts),edges,tuple(range(n0)))
        core=set(exterior)|set(range(8));new=set(range(n0,len(pts)))
        cross=tuple(e for e in edges if e[0]<n0<=e[1])
        retained=tuple(e for e in cross if set(e)<=core)
        degree=Counter(v for e in cross for v in e if v>=n0)
        retract=list(range(n0))+[None]*(len(pts)-n0)
        for v,y in zip(aliases,self.projected_recipe(name)):
            need(y in index and index[y]<n0,'literal old image of every new address')
            target=index[y]
            need(retract[v] is None or retract[v]==target,'all quotient aliases preserve projected image')
            retract[v]=target
        retract=tuple(retract)
        audit_retraction(len(pts),edges,n0,self.edges,retract)
        summary={'recipe':name,'raw_new_addresses':len(raw),'vertices':len(pts),'edges':len(edges),
                 'new_vertices':len(new),'new_cubic_vertices':sum(any(c for c,b in zip(x,self.f.basis) if b[2]) for x in pts[n0:]),
                 'core_vertices':len(core),'new_core_vertices':len(core&new),
                 'old_fixed_new_core':len(relative),'old_new_edges':len(cross),
                 'old_new_edges_retained':len(retained),'new_old_degree_histogram':dict(sorted(Counter(degree[v] for v in new).items())),
                 'sieve_primes':[p[0] for p in self.ctx.geo.pl]}
        return {'summary':summary,'points':tuple(pts),'edges':edges,'old_inclusion':tuple(range(n0)),
                'new_address_to_vertex':tuple(aliases),'relative_core':exterior,'deletion_order':peel,
                'old_fixed_new_core':relative,'old_fixed_deletion_order':oldpeel,'old_new_edges':cross,
                'retained_old_new_edges':retained,'retraction':retract}


def graph(n,edges):
    need(type(n) is int and n>=1,'positive exact vertex count')
    need(type(edges) is tuple,'literal edge tuple')
    need(all(type(e) is tuple and len(e)==2 and all(type(x) is int for x in e)
             and 0<=e[0]<e[1]<n for e in edges),'canonical typed endpoints')
    need(tuple(sorted(set(edges)))==edges,'canonical distinct edges')


def audit_retraction(n,edges,old_n,old_edges,retraction):
    graph(n,edges);graph(old_n,old_edges)
    need(old_n<=n and tuple(e for e in edges if e[1]<old_n)==old_edges,'literal induced old graph')
    need(type(retraction) is tuple and len(retraction)==n
         and all(type(v) is int and 0<=v<old_n for v in retraction),'exact map into the old graph')
    need(retraction[:old_n]==tuple(range(old_n)),'every old vertex fixed')
    E=set(old_edges)
    need(all(tuple(sorted((retraction[a],retraction[b]))) in E for a,b in edges),
         'every actual edge preserved without collapse')
    return True


def extend_colouring(packet,old_edges,colours):
    """Authenticate an old colouring and pull it back along a checked retraction."""
    need(type(packet) is dict,'literal packet')
    n=len(packet['points']);old_n=len(packet['old_inclusion'])
    need(type(packet['old_inclusion']) is tuple
         and all(type(v) is int for v in packet['old_inclusion'])
         and packet['old_inclusion']==tuple(range(old_n)),'literal old index inclusion')
    audit_retraction(n,packet['edges'],old_n,old_edges,packet['retraction'])
    need(type(colours) is tuple and len(colours)==old_n
         and all(type(c) is int and 0<=c<6 for c in colours),'exact six-colour source word')
    need(all(colours[a]!=colours[b] for a,b in old_edges),'proper old colouring')
    out=tuple(colours[v] for v in packet['retraction'])
    need(out[:old_n]==colours and all(out[a]!=out[b] for a,b in packet['edges']),
         'literal source preserved and all full edges checked')
    return out


def probe():
    e=Engine();start=time.monotonic()
    for name in ('rotated-seed','cubic-phase','cubic-rebased-phase','coupled-seed-phase'):
        p=e.build(name)
        print(json.dumps(p['summary'],sort_keys=True),round(time.monotonic()-start,3),flush=True)


def export_packet(p):
    out=dict(p);out['points']=[encode_element(x) for x in p['points']]
    raw=json.dumps({'points':out['points'],'edges':out['edges']},sort_keys=True,separators=(',',':')).encode()
    out['point_edge_sha256']=hashlib.sha256(raw).hexdigest()
    return out


def verify():
    e=Engine();f=e.f;checks=0
    def check(ok,msg):
        nonlocal checks
        checks+=1;need(ok,msg)
    def rejects(call):
        try:call()
        except (ValueError,KeyError,TypeError):check(True,'invalid exact source rejected')
        else:check(False,'invalid exact source accepted')
    check(e.ctx.norm(e.u)==f.one,'unit cubic rotation')
    w=f.gen(2)
    check(f.sub(f.mul(e.ctx.geo.j,
                      f.mul(f.add(f.one,e.u),f.inv(f.sub(f.one,e.u)))),f.scalar(3))==w,
          'recover cubic generator from unit rotation')
    check(e.k_neighbour(w) is None,'real exterior has no K-unit neighbour: singular determinant')
    for a in (f.zero,f.one,e.ctx.ports[2]):
        for v in (e.u,f.power(e.u,2),f.mul(e.u,e.q)):
            check(e.k_neighbour(f.add(a,v))==a,'planted unique K-neighbour')
    rejects(lambda:e.k_neighbour(f.one))
    rejects(lambda:e.k_neighbour(tuple(float(v) for v in e.u)))
    rejects(lambda:decode_element([[False,1,1]]))
    rejects(lambda:decode_element([[0,1.0,1]]))
    rejects(lambda:decode_element([[0,1,0]]))
    shifts=e.rebased_shifts();catalogue=[]
    check(len(shifts)==7,'six distinct genuine-cubic branches')
    for i,j in combinations(range(7),2):
        diff=f.sub(shifts[j],shifts[i]);a=e.k_neighbour(diff)
        expected=f.scale(e.gens[j-1],-1) if i==0 else None
        check(a==expected,'complete star K-contact catalogue')
        catalogue.append({'layers':[i,j],'K_neighbour':None if a is None else encode_element(a)})
    packets={}
    old_edge_set=set(e.edges)
    for name in ('rotated-seed','cubic-phase','cubic-rebased-phase','coupled-seed-phase'):
        p=e.build(name);packets[name]=p
        check(p['points'][:714]==e.points,'literal 714-point inclusion')
        check(p['summary']['new_cubic_vertices']==p['summary']['new_vertices'],'every new point is outside K')
        for a,b in p['edges']:
            check(tuple(sorted((p['retraction'][a],p['retraction'][b]))) in old_edge_set,
                  'actual edge maps to an old edge')
        check(all(v<=1 for v in map(int,p['summary']['new_old_degree_histogram'])),
              'each genuine cubic point has at most one old neighbour')
    p=packets['cubic-rebased-phase']
    check(independent_edges(e.s,p['points'])==p['edges'],'independent one-place complete incidence')
    check((len(p['points']),len(p['edges']),len(p['old_fixed_new_core']))==(1482,6118,491),
          'frozen main candidate dimensions')
    # Compare the six entire new branches, not only selected contacts.
    index={x:i for i,x in enumerate(p['points'])}
    branches=tuple(tuple(index[f.add(r,t)] for r in e.R) for t in shifts[1:])
    check(len(set(v for branch in branches for v in branch))==768,'complete disjoint new branch addresses')
    expected_cross=set()
    old_index={x:i for i,x in enumerate(e.points)}
    for s,t,branch in zip(e.gens,shifts[1:],branches):
        check(e.ctx.norm(f.add(t,s))==f.one,'universal fractional contact identity')
        for r,v in zip(e.R,branch):
            other=f.sub(r,s)
            if other in old_index:expected_cross.add((old_index[other],v))
    check(tuple(sorted(expected_cross))==p['old_new_edges'],'all and only translated old contacts')
    same_branch={v:i for i,b in enumerate(branches) for v in b}
    check(all(same_branch[a]==same_branch[b] for a,b in p['edges'] if a>=714),
          'no new edges between distinct fractional branches')
    union=packets['coupled-seed-phase'];ui={x:i for i,x in enumerate(union['points'])}
    ra={ui[x] for x in e.recipe('rotated-seed')}-set(range(714))
    rb={ui[x] for x in e.recipe('cubic-rebased-phase')}-set(range(714))
    check(not ra&rb,'two exterior families have no aliases')
    check(not any((a in ra and b in rb) or (b in ra and a in rb) for a,b in union['edges']),
          'last bounded union adds no cross-family edge')
    # Saved physical words are independent inputs, not inferred from a successful chart.
    support_path=ROOT/'05-knowledge/results/coupled_support_branches_20261018c.json'
    bank=json.loads(support_path.read_text(encoding='utf8'))['geometric_candidates']['lifted_candidate']['rows']
    check(len(bank)==570,'complete inherited source-pattern witness bank')
    witness_hash=hashlib.sha256();first=None
    for row in bank:
        colours=tuple(int(c) for c in row['colours'])
        out=extend_colouring(p,e.edges,colours)
        check(out[:8]==tuple(row['pattern']),'source-pattern and full-old-word preservation')
        word=''.join(map(str,out));witness_hash.update(word.encode());witness_hash.update(b'\n')
        if first is None:first={'old_colours':row['colours'],'full_colours':word}
    rejects(lambda:extend_colouring(p,e.edges,tuple(False for _ in e.points)))
    bad=dict(p);bad['retraction']=(False,)+p['retraction'][1:]
    rejects(lambda:extend_colouring(bad,e.edges,tuple(int(c) for c in bank[0]['colours'])))
    bad=dict(p);bad['old_inclusion']=(False,)+p['old_inclusion'][1:]
    rejects(lambda:extend_colouring(bad,e.edges,tuple(int(c) for c in bank[0]['colours'])))
    bad=dict(p);bad['retraction']=tuple(range(714))+(0,)*(len(p['points'])-714)
    rejects(lambda:extend_colouring(bad,e.edges,tuple(int(c) for c in bank[0]['colours'])))
    rejects(lambda:audit_retraction(3,((0,1),(1,2)),2,((0,1),),(0,1,True)))
    check(audit_retraction(3,((0,1),(1,2)),2,((0,1),),(0,1,0)),'positive path retraction control')
    square=(f.zero,f.one,e.u,f.add(f.one,e.u))
    square_edges=tuple((i,j) for i,j in combinations(range(4),2)
                       if e.ctx.norm(f.sub(square[j],square[i]))==f.one)
    check(square_edges==((0,1),(0,2),(1,3),(2,3)),'planted cubic perfect-matching square')
    rejects(lambda:audit_retraction(4,square_edges,2,((0,1),),(0,1,0,1)))
    square_colours=(0,1,1,2)
    check(all(square_colours[a]!=square_colours[b] for a,b in square_edges),
          'failure of forgetting the layer is not failure of six-colouring')
    return {'status':'PROVED scoped cubic contact/retraction mechanism; FINITE-EXACT complete candidate geometry',
            'checks':checks,
            'field_polynomials':e.f.polys,'basis':e.f.basis,
            'old_dependency':str(OLD.relative_to(ROOT)).replace('\\','/'),
            'old_point_edge_sha256':e.old['point_edge_sha256'],
            'cubic_unit':encode_element(e.u),'phase':encode_element(e.q),
            'rebased_shifts':[encode_element(t) for t in shifts],
            'contact_catalogue':catalogue,'candidate':export_packet(p),
            'probe_summaries':[x['summary'] for x in packets.values()],
            'old_word_pullback':{'words':570,'all_full_word_sha256':witness_hash.hexdigest(),'first':first},
            'nonretraction_control':{'points':[encode_element(x) for x in square],
                                    'edges':square_edges,'colours':square_colours},
            'scope':'Every old colouring pulls back. No source pattern is eliminated. Marked terminals (2,4) remain nonunit.'}


if __name__=='__main__':
    if '--probe' in sys.argv:probe()
    else:print(json.dumps(verify(),sort_keys=True,indent=2))
