"""Exact bounded coupled fractional candidates; source ports remain literal."""
from fractions import Fraction as Q
from itertools import combinations
from collections import Counter
import hashlib
import json
import sys
import time
import numpy as np

from cubic_boundary_search_20261018a import Context, encode_element
from escaped_ports_20261017d import Field
from shared_phase_geometry_20261018b import PLACES


def need(ok, message):
    if not ok:
        raise ValueError(message)


def relative_core(n, edges, fixed=tuple(range(8))):
    need(type(n) is int and n>=1,'positive exact vertex count')
    need(type(edges) is tuple and type(fixed) is tuple,'literal edge and fixed-port tuples')
    for e in edges:
        need(type(e) is tuple and len(e)==2 and all(type(x) is int for x in e)
             and 0<=e[0]<e[1]<n,'canonical typed edge')
    need(tuple(sorted(set(edges)))==edges,'complete canonical edge list')
    need(all(type(x) is int and 0<=x<n for x in fixed) and len(set(fixed))==len(fixed),
         'literal distinct fixed vertices')
    adj = [set() for _ in range(n)]
    for a,b in edges:
        adj[a].add(b); adj[b].add(a)
    live = set(range(n)); order=[]; fixed=set(fixed)
    while True:
        batch = [v for v in sorted(live-fixed) if len(adj[v]&live)<=5]
        if not batch:
            break
        for v in batch:
            live.remove(v); order.append(v)
    return tuple(sorted(live-fixed)), tuple(order)


class Search:
    def __init__(self):
        self.ctx=Context(); self.f=self.ctx.f
        self.good_places=[]
        for place in PLACES:
            p=place.p
            rs=(place.z,place.j,place.w)
            for poly,r in zip(self.f.polys,rs):
                need(sum(c*pow(r,i,p) for i,c in enumerate(poly))%p==0,'source place relation')
            rsc=(pow(place.z,-1,p),-place.j%p,place.w)
            def bs(roots):
                return np.array([pow(roots[0],a,p)*pow(roots[1],b,p)*pow(roots[2],c,p)%p
                                 for a,b,c in self.f.basis],dtype=object)
            self.good_places.append((p,bs(rs),bs(rsc)))
        self.h=Field([self.f.polys[0]])
        h=self.h; z=h.gen(0)
        self.bar_basis=h.substitution((h.power(z,14),))
        delta=h.power(z,5)
        alpha=h.power(z,3)
        self.eta=h.add(h.one,delta)
        root5=h.add(h.one,h.scale(h.add(alpha,h.power(alpha,4)),2))
        root15=h.mul(root5,h.add(h.one,h.scale(delta,2)))
        self.q=h.scale(h.mul(delta,h.sub(root15,h.one)),Q(1,4))
        self.qbar=self.bar(self.q)
        need(h.mul(self.q,self.qbar)==h.one,'unit starting phase')

    def bar(self,x):
        need(type(x) is tuple and len(x)==self.h.n and all(type(v) is Q for v in x),
             'exact cyclotomic element')
        return self.h.apply(x,self.bar_basis)

    def embed(self,x):
        need(type(x) is tuple and len(x)==self.h.n and all(type(v) is Q for v in x),
             'exact cyclotomic element')
        out=list(self.f.zero)
        for a,c in enumerate(x):
            out[self.f.index[(a,0,0)]]=c
        return tuple(out)

    def transform(self,x,s):
        self.bar(x); sb=self.bar(s)
        h=self.h
        need(h.mul(x,self.bar(x))==h.one and h.mul(s,sb)==h.one,'unit source and direction')
        # T_s(z)=(2z-s)/(2-bar(s)z), preserving the unit circle.
        den=h.sub(h.scalar(2),h.mul(sb,x))
        need(den!=h.zero,'no pole')
        y=h.mul(h.sub(h.scale(x,2),s),h.inv(den))
        need(h.mul(y,self.bar(y))==h.one,'unit phase after fractional move')
        return y

    def orbit(self, depth, alphabet='triangular'):
        need(type(depth) is int and 0<=depth<=3,'bounded exact word depth')
        need(type(alphabet) is str,'named alphabet')
        h=self.h
        if alphabet=='triangular':
            generators=tuple(h.power(self.eta,i) for i in range(6))
        elif alphabet=='pentagonal':
            alpha=h.power(h.gen(0),3)
            generators=tuple(h.scale(h.power(alpha,i),s) for s in (1,-1) for i in range(5))
        elif alphabet=='both':
            alpha=h.power(h.gen(0),3)
            generators=tuple(dict.fromkeys(tuple(h.power(self.eta,i) for i in range(6))+
                              tuple(h.scale(h.power(alpha,i),s) for s in (1,-1) for i in range(5))))
        else:
            raise ValueError('unknown alphabet')
        vals=[self.q]; words=[()]; seen={self.q:0}; front=[0]
        for _ in range(depth):
            new=[]
            for i in front:
                for j,s in enumerate(generators):
                    y=self.transform(vals[i],s)
                    if y not in seen:
                        seen[y]=len(vals); new.append(len(vals)); vals.append(y); words.append(words[i]+(j,))
            front=new
        shifts=tuple(self.embed(h.scale(h.sub(v,self.q),2)) for v in vals)
        return shifts,words,generators

    def build(self, shifts, label):
        need(type(shifts) is tuple and shifts and shifts[0]==self.f.zero,'source-owned zero shift first')
        for x in shifts:self.ctx.element(x)
        f=self.f
        points=list(self.ctx.ports); lookup={x:i for i,x in enumerate(points)}; raw=[]
        for layer,shift in enumerate(shifts):
            for a,c in enumerate(self.ctx.ports):
                x=f.add(c,shift)
                if x not in lookup:
                    lookup[x]=len(points); points.append(x)
                raw.append(lookup[x])
        self.ctx.geo.pl=[row for row in self.good_places
                         if all(v.denominator%row[0] for x in points for v in x)]
        need(self.ctx.geo.pl,'at least one authenticated good reduction')
        edges=tuple(self.ctx.geo.unit_edges(points))
        core,order=relative_core(len(points),edges)
        return {'label':label,'points':tuple(points),'edges':edges,'raw_to_vertex':raw,
                'layers':len(shifts),'relative_core':core,'deletion_order':order,
                'sieve_primes':tuple(row[0] for row in self.ctx.geo.pl)}


def probe():
    s=Search(); start=time.monotonic()
    for alphabet,depth in [('triangular',1),('triangular',2),('triangular',3),('both',1),('both',2)]:
        shifts,words,generators=s.orbit(depth,alphabet)
        p=s.build(shifts,f'{alphabet}-{depth}')
        print(p['label'],len(shifts),len(p['points']),len(p['edges']),len(p['relative_core']),
              round(time.monotonic()-start,2),flush=True)
        if p['relative_core']:
            return p
    shifts,words,generators=s.orbit(1,'triangular')
    for power in (2,3):
        if power==2:
            sums=tuple(dict.fromkeys(s.f.add(a,b) for a in shifts for b in shifts))
        else:
            sums=tuple(dict.fromkeys(s.f.add(a,b) for a in sums for b in shifts))
        p=s.build(sums,f'sum-{power}')
        print(p['label'],len(sums),len(p['points']),len(p['edges']),len(p['relative_core']),
              round(time.monotonic()-start,2),flush=True)
        if p['relative_core']:
            return p
    return None


def build_candidate(s, lift=False):
    need(type(lift) is bool,'literal triangle-lift flag')
    f=s.f
    shifts,words,gens=s.orbit(1,'triangular')
    M=tuple(s.ctx.geo.M)
    triangle=(f.zero,f.one,s.embed(s.eta))
    if lift:
        M=tuple(dict.fromkeys(f.add(a,b) for a in M for b in triangle))
    banks=[('M+triangle' if lift else 'M',M)]
    for label,seed in banks:
        sums=tuple(dict.fromkeys(f.add(a,b) for a in seed for b in shifts))
        p=s.build(sums,label+'+fractional')
        if p['relative_core']:
            p['phase']=s.embed(s.q)
            p['fractional_shifts']=shifts
            p['fractional_words']=words
            p['fractional_generators']=tuple(s.embed(g) for g in gens)
            phase_values=tuple(s.h.add(s.q,s.h.scale(tuple(x[s.f.index[(a,0,0)]] for a in range(s.h.n)),Q(1,2)))
                               for x in shifts)
            p['fractional_move_to_layer']=tuple(phase_values.index(s.transform(s.q,g)) for g in gens)
            p['spindle']=M
            p['sum_shifts']=sums
            lookup={x:i for i,x in enumerate(p['points'])}
            p['addresses']=tuple((a,m,t) for a in range(8) for m in range(len(M)) for t in range(len(shifts)))
            p['address_to_vertex']=tuple(lookup[f.add(s.ctx.ports[a],f.add(M[m],shifts[t]))]
                                          for a,m,t in p['addresses'])
            factors=(s.ctx.ports,M,shifts)
            factor_edges=tuple(tuple(s.ctx.geo.unit_edges(factor)) for factor in factors)
            p['factor_edges']=factor_edges
            owner={address:v for address,v in zip(p['addresses'],p['address_to_vertex'])}
            cartesian=set()
            for address,v in owner.items():
                for k,ee in enumerate(factor_edges):
                    for aa,bb in ee:
                        if address[k] not in (aa,bb):
                            continue
                        other=list(address);other[k]=bb if address[k]==aa else aa
                        w=owner[tuple(other)]
                        if v!=w:
                            cartesian.add(tuple(sorted((v,w))))
            p['cartesian_edges']=tuple(sorted(cartesian))
            p['noncartesian_edges']=tuple(e for e in p['edges'] if e not in cartesian)
            keep=set(p['relative_core'])|set(range(8))
            p['core_vertices']=tuple(sorted(keep))
            p['core_edges']=tuple(e for e in p['edges'] if set(e)<=keep)
            p['retained_noncartesian_edges']=tuple(e for e in p['noncartesian_edges'] if set(e)<=keep)
            p['triangle_lift']=lift
            return p
    raise ValueError('declared candidate unexpectedly has empty core')


def export_packet(p):
    out=dict(p)
    for key in ('points','fractional_shifts','fractional_generators','spindle','sum_shifts'):
        out[key]=[encode_element(x) for x in p[key]]
    out['phase']=encode_element(p['phase'])
    out['original_ports']=list(range(8))
    return out


def independent_edges(s, points):
    """Complete alternative one-place rejection, then exact surviving norms."""
    from golden_vsquare_unit_census_20261017 import red
    row=next(row for row in reversed(s.good_places)
             if all(v.denominator%row[0] for x in points for v in x))
    prime,basis,conjugate=row
    rr=[red(x,prime,basis) for x in points]
    bb=[red(x,prime,conjugate) for x in points]
    out=[]
    for a,b in combinations(range(len(points)),2):
        if (rr[b]-rr[a])*(bb[b]-bb[a])%prime==1:
            if s.ctx.norm(s.f.sub(points[b],points[a]))==s.f.one:
                out.append((a,b))
    return tuple(out)


def verify():
    s=Search(); checks=0
    def check(ok,message):
        nonlocal checks
        checks+=1
        need(ok,message)
    def rejects(call):
        try:call()
        except ValueError:check(True,'invalid type/source rejected')
        else:check(False,'invalid type/source accepted')
    summaries=[]
    for alphabet,depth in [('triangular',1),('triangular',2),('triangular',3),('both',1),('both',2)]:
        shifts,_,_=s.orbit(depth,alphabet)
        p=s.build(shifts,f'{alphabet}-{depth}')
        check(not p['relative_core'],'word ball peels completely')
        summaries.append({'alphabet':alphabet,'depth':depth,'layers':len(shifts),
                          'vertices':len(p['points']),'edges':len(p['edges']),
                          'relative_core_size':len(p['relative_core'])})
    packets=[]
    for lift in (False,True):
        p=build_candidate(s,lift);points=p['points'];edges=p['edges'];n=len(points)
        expected=(714,2914,464,1903,180,94) if lift else (316,1004,56,190,45,0)
        check((n,len(edges),len(p['relative_core']),len(p['core_edges']),
               len(p['noncartesian_edges']),len(p['retained_noncartesian_edges']))==expected,
              'declared complete candidate universe')
        check(points[:8]==s.ctx.ports and len(set(points))==n,'source identity and exact distinctness')
        for x in points:
            check(all(not coeff or basis[2]==0 for coeff,basis in zip(x,s.f.basis)),
                  'actual coefficients lie in Q(zeta15,j), no cubic coordinate')
        check(edges==independent_edges(s,points),'independent complete single-place incidence')
        check(set(p['cartesian_edges'])<=set(edges),'every declared factor edge is actual')
        check(set(edges)==set(p['cartesian_edges'])|set(p['noncartesian_edges']),
              'exact factor/extra edge decomposition')
        for address,v in zip(p['addresses'],p['address_to_vertex']):
            a,m,t=address
            check(points[v]==s.f.add(s.ctx.ports[a],s.f.add(p['spindle'][m],p['fractional_shifts'][t])),
                  'every retained source address authenticates')
        for i,layer in enumerate(p['fractional_move_to_layer']):
            generator=tuple(p['fractional_generators'][i][s.f.index[(a,0,0)]] for a in range(s.h.n))
            expected_shift=s.embed(s.h.scale(s.h.sub(s.transform(s.q,generator),s.q),2))
            check(p['fractional_shifts'][layer]==expected_shift,'all six raw fractional moves retained')
        adjacent=[set() for _ in range(n)]
        for a,b in edges:
            adjacent[a].add(b);adjacent[b].add(a)
            check(s.ctx.norm(s.f.sub(points[b],points[a]))==s.f.one,'all reported edges exact')
        live=set(range(n))
        for v in p['deletion_order']:
            check(v>=8 and v in live and len(adjacent[v]&live)<=5,'each deletion valid')
            live.remove(v)
        check(live==set(p['core_vertices']),'exact complete relative core')
        for v in p['relative_core']:
            check(len(adjacent[v]&live)>=6,'no remaining exterior deletion')
        p['core_exterior_degree_histogram']=dict(sorted(Counter(len(adjacent[v]&live) for v in p['relative_core']).items()))
        p['point_edge_sha256']=hashlib.sha256(json.dumps([p['addresses'],p['address_to_vertex'],edges],
                                                        separators=(',',':')).encode()).hexdigest()
        packets.append(p)
    # Explicit old-source terminal remains a nonunit pentagon diagonal.
    f=s.f; phi=f.add(f.one,f.add(s.ctx.geo.alpha,f.power(s.ctx.geo.alpha,4)))
    check(s.ctx.norm(f.sub(s.ctx.ports[4],s.ctx.ports[2]))==f.mul(phi,phi)
          and f.mul(phi,phi)!=f.one,'terminal2/4 remains NONUNIT')
    first,second=packets
    second_lookup={x:i for i,x in enumerate(second['points'])}
    inclusion=tuple(second_lookup[x] for x in first['points'])
    check(inclusion[:8]==tuple(range(8)),'source-fixed triangle inclusion')
    check(all(tuple(sorted((inclusion[a],inclusion[b]))) in set(second['edges'])
              for a,b in first['edges']),'all old induced edges retained')
    # A second order of valid peeling reaches the same core.
    for p in packets:
        n=len(p['points']);adj=[set() for _ in range(n)]
        for a,b in p['edges']:adj[a].add(b);adj[b].add(a)
        live=set(range(n))
        while True:
            removable=[v for v in sorted(live,reverse=True) if v>=8 and len(adj[v]&live)<=5]
            if not removable:break
            live.remove(removable[0])
        check(live==set(p['core_vertices']),'peeling order independence')
    # Abstract positive/negative controls, not geometric K7 claims.
    check(relative_core(7,tuple(combinations(range(7),2)),(0,))[0]==tuple(range(1,7)),
          'nonempty degree6 abstract control')
    check(not relative_core(7,tuple((i,i+1) for i in range(6)),(0,))[0],
          'path peels completely')
    for bad in (True,1.0):
        rejects(lambda bad=bad:s.orbit(bad))
        rejects(lambda bad=bad:relative_core(8,((0,bad),)))
        rejects(lambda bad=bad:s.embed((Q(1),bad)+(Q(0),)*6))
    rejects(lambda:build_candidate(s,1))
    rejects(lambda:s.transform(s.h.zero,s.h.one))
    rejects(lambda:relative_core(8,((0,1),(0,1))))
    return {'status':'FINITE-EXACT complete coupled configurations; equality forcing OPEN',
            'checks':checks,'field_polynomials':s.f.polys,'basis':s.f.basis,
            'actual_declared_subfield':'Q(zeta15,sqrt(-11)); every w coefficient is zero',
            'word_ball_controls':summaries,'candidate':export_packet(first),
            'lifted_candidate':export_packet(second),'inclusion_in_lift':inclusion,
            'terminal':{'vertices':(2,4),'squared_distance':encode_element(f.mul(phi,phi)),
                        'unit':False}}


if __name__=='__main__':
    if '--probe' in sys.argv:
        probe()
    else:
        print(json.dumps(verify(),indent=2,sort_keys=True))
