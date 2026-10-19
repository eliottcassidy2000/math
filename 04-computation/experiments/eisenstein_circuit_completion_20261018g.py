#!/usr/bin/env python3
"""Exact circuit completion of supplied native Nyquist packets.

Consumes the complete source/route receipt of native_cancellation_105; no
orbit search, frontier enlargement, or new source supply. Default writes
this packet's result; --verify reproduces and compares the saved receipt.
"""
from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from pathlib import Path
import argparse
import json

import native_cancellation_105_20261018g as native

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent / '05-knowledge/results'
SOURCE = RESULTS / 'native_cancellation_105_20261018g.json'
RESULT = RESULTS / 'eisenstein_circuit_completion_20261018g.json'
SCALES = (4096, 5997, 11002, 16384)
CHECKS = 0

# Counterclockwise in a+b*(-1+i*sqrt(3))/2. No floating angle decisions.
DIRECTIONS = ((1,0),(2,1),(1,1),(1,2),(0,1),(-1,1),
              (-1,0),(-2,-1),(-1,-1),(-1,-2),(0,-1),(1,-1))
FUNCTIONALS = ((1,4),(-1,5),(-4,5),(-5,4),(-5,1),(-4,-1),
               (-1,-4),(1,-5),(4,-5),(5,-4),(5,-1),(4,1))
ALPHABET = tuple(sorted(DIRECTIONS))


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def integer(x, minimum=0):
    if type(x) is not int or x < minimum:
        raise ValueError('exact integer outside domain')
    return x


def vector(v):
    if type(v) is not tuple or len(v) != 2 or any(type(x) is not int for x in v):
        raise ValueError('exact integer coordinate tuple required')
    return v


def norm(v):
    a,b = vector(v)
    return a*a-a*b+b*b


def det(v,w):
    a,b = vector(v)
    c,d = vector(w)
    return a*d-b*c


def dot(f,v):
    a,b = vector(f)
    c,d = vector(v)
    return a*c+b*d


def total(vectors):
    a=b=0
    for v in vectors:
        x,y=vector(v)
        a+=x; b+=y
    return a,b


def digest(obj):
    return sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def enumerate_circuits():
    ans=[]
    for v,w in combinations(ALPHABET,2):
        if total((v,w)) == (0,0):
            ans.append(((v,w),(1,1)))
    for vs in combinations(ALPHABET,3):
        cs=(det(vs[1],vs[2]),det(vs[2],vs[0]),det(vs[0],vs[1]))
        if max(cs)<0:
            cs=tuple(-c for c in cs)
        if min(cs)<=0:
            continue
        g=gcd(gcd(cs[0],cs[1]),cs[2])
        ans.append((vs,tuple(c//g for c in cs)))
    return tuple(sorted(ans, key=lambda r:(sum(r[1]),len(r[0]),r[0],r[1])))


CIRCUITS=enumerate_circuits()


def audit_counts(counts):
    if type(counts) is not dict:
        raise ValueError('exact charge multiplicity dictionary required')
    for v,c in counts.items():
        vector(v); integer(c)
        if v not in ALPHABET:
            raise ValueError('charge outside the twelve-type alphabet')
    return counts


def available(counts,circuit):
    vs,cs=circuit
    return all(counts.get(v,0)>=c for v,c in zip(vs,cs))


def completed_bound(counts):
    """A source-agnostic algebra certificate; native authentication is separate."""
    audit_counts(counts)
    if any(available(counts,c) for c in CIRCUITS):
        raise ValueError('available primitive circuit remains')
    heavy={v for v,c in counts.items() if c>=3}
    choices=[i for i in range(12) if heavy <= {DIRECTIONS[(i+j)%12] for j in range(6)}]
    if not choices:
        raise RuntimeError('heavy half-plane theorem failed')
    i=choices[0]
    arc={DIRECTIONS[(i+j)%12] for j in range(6)}
    f=FUNCTIONALS[i]
    V=total(tuple(c*x for x in v) for v,c in counts.items())
    M=sum(c*norm(v) for v,c in counts.items())
    positive_defect=sum(c*max(norm(v)-dot(f,v),0) for v,c in counts.items())
    for v,c in counts.items():
        check(v in arc or c<=2, 'outside-arc charges are light')
    check(positive_defect<=80, 'actual light-defect bill')
    check(M-dot(f,V)<=positive_defect, 'linear majorant retains signed light terms')
    check(dot(f,V)**2<=28*norm(V), 'exact dual Cauchy bound')
    check(max(M-80,0)**2<=28*norm(V), 'residual mass bound')
    return {'arc':i,'functional':list(f),'vector':list(V),'source_norm':M,
            'actual_positive_light_defect':positive_defect,
            'linear_gap':M-dot(f,V),'quadratic_left':max(M-80,0)**2,
            'quadratic_right':28*norm(V),
            'counts':[[*v,counts.get(v,0)] for v in ALPHABET]}


def authenticate_source(data):
    """Replay supplied finite DAG edges and weights; never discover a route."""
    ranks={int(line.split()[0]):int(line.split()[1]) for line in native.FRONTIER.read_text().splitlines()
           if line.strip() and not line.startswith('#')}
    V=set(ranks)
    check(len(V)==627 and ranks[1]==0, 'fixed inherited frontier')
    check(data['frontier']==sorted(V), 'exact frontier source binding')
    check(data['frontier_ranks']==[[s,ranks[s]] for s in sorted(V)], 'frontier ranks retained')
    for s in V:
        native.primitive(s)
        y,dt=native.inherited.edge(s)
        check(y in V and (s==1 or ranks[y]<ranks[s]), 'finite strict frontier ROOT rank')
    check(data['frontier_sha256']==digest(sorted(V)), 'frontier hash')
    producer=sha256((HERE/'root_clock_channels_20261018o.py').read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    check(data['producer_sha256']==producer, 'unchanged actual source producer')
    edges={}
    for row in data['route_dag']:
        if type(row) is not list or len(row)!=6 or any(type(x) is not int for x in row):
            raise ValueError('literal integer route row required')
        s,y,dt,u,t,d=row
        native.primitive(s); native.primitive(y)
        integer(dt,1); integer(t,1); integer(d,1)
        if s in edges or s in V or u not in V:
            raise ValueError('duplicate, premature or invalid first-arrival route')
        check((y,dt)==native.inherited.edge(s), 'actual canonical primitive edge')
        edges[s]=(y,dt,u,t,d)
    check(data['route_dag_sha256']==digest(data['route_dag']), 'source DAG digest')
    for s,(y,dt,u,t,d) in edges.items():
        if y in V:
            tail=(y,0,0)
        elif y in edges:
            tail=edges[y][2:]
        else:
            raise ValueError('missing route suffix')
        check((u,t,d)==(tail[0],tail[1]+dt,tail[2]+1), 'first-arrival depth and clock induction')
    check([x['statistics']['N'] for x in data['scales']]==list(SCALES), 'fixed four-scale source universe')
    for item in data['scales']:
        n=item['statistics']['N']; rows=item['source_rows']
        expected=[s for s in range(1,4*n,2) if s%8!=5]
        if type(rows) is not list or len(rows)!=len(expected):
            raise ValueError('incomplete physical source universe')
        for s,row in zip(expected,rows):
            if type(row) is not list or len(row)!=9 or any(type(x) is not int for x in row):
                raise ValueError('literal physical source row required')
            u,t,d=(s,0,0) if s in V else edges[s][2:]
            a,b=native.charge(s,n)
            check((a,b)==native.independent_charge(s,n), 'independent actual shell charge')
            check(row==[s,s%15,s%35,u,t,d,a,b,norm((a,b))], 'source/clock/charge/register authentication')
        check(item['statistics']['source_rows_sha256']==digest(rows), 'complete source rows digest')
        check(item['statistics']['packets_sha256']==digest(item['packets']), 'inherited packet digest')
    return data


def packet_receipt(ids, rows, depth_horizon, clock_horizon):
    """Consume already-authenticated source rows at both complete-arrival cuts."""
    integer(depth_horizon); integer(clock_horizon)
    if type(ids) is not list or not 2<=len(ids)<=6:
        raise ValueError('two to six named physical sources required')
    if any(type(s) is not int for s in ids) or len(set(ids))!=len(ids):
        raise ValueError('distinct exact source IDs required')
    if any(s not in rows for s in ids):
        raise ValueError('source absent from authenticated row')
    selected=[rows[s] for s in ids]
    if len({r[3] for r in selected})!=1 or any(r[8]==0 for r in selected):
        raise ValueError('one actual terminal and nonzero charges required')
    if max(r[5] for r in selected)>depth_horizon or max(r[4] for r in selected)>clock_horizon:
        raise ValueError('packet is not fully arrived at the declared cut')
    vectors=[native.signed(r) for r in selected]
    if total(vectors)!=(0,0):
        raise ValueError('nonzero packet')
    mass=sum(r[8] for r in selected)
    debt=sum(native.cross(v,w) for v,w in combinations(vectors,2))
    if debt!=-mass:
        raise ValueError('unpaid internal cross debt')
    return {'sources':ids,'terminal':selected[0][3],'source_norm':mass,'cross_debt':debt,
            'depth_horizon':depth_horizon,'clock_horizon':clock_horizon,
            'source_clocks':[[r[0],r[4],r[5]] for r in selected]}


def consume(ids,rows,spent):
    if spent.intersection(ids):
        raise ValueError('physical source spent twice')
    rec=packet_receipt(ids,rows,max(rows[s][5] for s in ids),max(rows[s][4] for s in ids))
    spent.update(ids)
    return rec


def complete_scale(item):
    n=item['statistics']['N']; source_rows=item['source_rows']
    rows={r[0]:r for r in source_rows}; spent=set(); inherited_mass=0
    for ids in item['packets']:
        native.packet_data(ids,rows)
        inherited_mass+=consume(ids,rows,spent)['source_norm']
    check(inherited_mass==item['statistics']['removed_source_norm'], 'inherited extraction retained')
    old_spent=set(spent)
    buckets=defaultdict(lambda:defaultdict(deque))
    for r in source_rows:
        if r[8] and r[0] not in spent:
            buckets[r[3]][native.signed(r)].append(r[0])
    extra=[]; bounds=[]
    for u in sorted({r[3] for r in source_rows}):
        q=buckets[u]
        for index,circuit in enumerate(CIRCUITS):
            vs,cs=circuit
            while all(len(q[v])>=c for v,c in zip(vs,cs)):
                ids=sorted(q[v].popleft() for v,c in zip(vs,cs) for _ in range(c))
                rec=consume(ids,rows,spent); rec['circuit']=index
                extra.append(rec)
        counts={v:len(q[v]) for v in ALPHABET}
        cert=completed_bound(counts); cert['terminal']=u
        bounds.append(cert)
    before=native.terminal_vectors(source_rows)
    after=native.terminal_vectors(source_rows,spent)
    check(before==after, 'every actual terminal vector unchanged')
    extra_mass=sum(r['source_norm'] for r in extra)
    mass=sum(r[8] for r in source_rows)
    residual=sum(r[8] for r in source_rows if r[0] not in spent)
    check(mass-inherited_mass-extra_mass==residual, 'source norm never spent twice')
    check(sum(r['source_norm'] for r in bounds)==residual, 'all residual ports retained')
    check(all(tuple(r['vector'])==before.get(r['terminal'],(0,0)) for r in bounds), 'certificate vector is actual terminal sum')
    energy=sum(norm(v) for v in before.values())
    check(energy==item['statistics']['terminal_energy'], 'original Nyquist energy retained')
    # A spent packet contributes zero cross terms against every other packet
    # or residual vector at its terminal, by bilinearity; none is merged across ports.
    stats={'N':n,'inherited_paid_norm':inherited_mass,'inherited_residual_norm':mass-inherited_mass,
           'extra_packets':len(extra),'extra_sources':len(spent-old_spent),'extra_paid_norm':extra_mass,
           'extra_cross_debt':-extra_mass,'residual_source_norm':residual,'unchanged_terminal_energy':energy,
           'packet_size_histogram':dict(sorted(Counter(len(r['sources']) for r in extra).items())),
           'actual_terminal_count':len(bounds),'max_actual_light_defect':max(r['actual_positive_light_defect'] for r in bounds),
           'source_rows_sha256':digest(source_rows),'extra_packets_sha256':digest(extra)}
    return {'statistics':stats,'extra_packets':extra,'terminal_bounds':bounds}


def reject(fn,why):
    try:
        fn()
    except (ValueError,RuntimeError):
        check(True,why)
    else:
        check(False,why)


def controls():
    check(len(set(DIRECTIONS))==12 and set(DIRECTIONS)==set(native.ALPHABET), 'exact twelve signed charges')
    check(len(CIRCUITS)==46, 'six pair and forty three-support circuits')
    check(Counter(len(vs) for vs,cs in CIRCUITS)=={2:6,3:40}, 'circuit support census')
    check(Counter(sum(cs) for vs,cs in CIRCUITS)=={2:6,3:10,4:12,5:6,6:12}, 'circuit lengths')
    for vs,cs in CIRCUITS:
        check(max(cs)<=3 and gcd(*cs)==1, 'primitive bounded circuit coefficients')
        check(total(tuple(c*x for x in v) for v,c in zip(vs,cs))==(0,0), 'exact positive circuit')
        check(all(det(v,w)!=0 for v,w in combinations(vs,2)) if len(vs)==3 else True, 'minimal three-support dependence')
    certificates=[]
    for i,f in enumerate(FUNCTIONALS):
        arc=[DIRECTIONS[(i+j)%12] for j in range(6)]
        check(all(dot(f,v)>=norm(v) for v in arc), 'six-arc linear majorant')
        check(4*(f[0]*f[0]+f[0]*f[1]+f[1]*f[1])==3*28, 'dual norm exactly twenty-eight')
        check(2*sum(max(norm(v)-dot(f,v),0) for v in ALPHABET)==80, 'universal light bill exactly eighty')
        certificates.append({'arc':i,'directions':arc,'functional':f,'dual_norm_squared':28,'light_bill':80})
    # Exhaustive support-only alternative to the planar positive-dependence argument.
    for mask in range(1<<12):
        support={DIRECTIONS[i] for i in range(12) if mask>>i&1}
        has_circuit=any(set(vs)<=support for vs,cs in CIRCUITS)
        has_arc=any(support<={DIRECTIONS[(i+j)%12] for j in range(6)} for i in range(12))
        check(has_circuit!=has_arc, 'all 4096 support subsets: circuit or strict half-plane')
    # Independent small multiplicity tests on every support of at most three types.
    tested=0
    for size in range(4):
        for vs in combinations(ALPHABET,size):
            for cs in product(range(1,5),repeat=size):
                counts=dict(zip(vs,cs))
                if not any(available(counts,c) for c in CIRCUITS):
                    completed_bound(counts); tested+=1
    hostile={(1,0):1,(0,1):2,(-1,-2):1}
    check(total(tuple(c*x for x in v) for v,c in hostile.items())==(0,0), 'four-member zero hostile')
    check(not any(sum(cs)<=3 and available(hostile,(vs,cs)) for vs,cs in CIRCUITS), 'no pair or triple consumes the four-member hostile')
    check(any(sum(cs)==4 and available(hostile,(vs,cs)) for vs,cs in CIRCUITS), 'four-member circuit resolves hostile')
    reject(lambda:completed_bound(hostile), 'incomplete residual certificate rejected')
    # Minimal real supports need not fit an integer multiset's capacities.
    integer_hostile={(-2,-1):1,(0,-1):1,(1,0):1,(1,2):1}
    check(total(integer_hostile)==(0,0), 'four-support integer zero sum')
    check(not any(available(integer_hostile,c) for c in CIRCUITS), 'primitive support circuits miss this capacity-limited zero sum')
    bound=completed_bound(integer_hostile)
    check(bound['source_norm']==8 and bound['quadratic_right']==0, 'positive light allowance cannot be discarded')
    reject(lambda:completed_bound({(True,0):3}), 'Boolean charge rejected')
    reject(lambda:completed_bound({(1,0):3.0}), 'float multiplicity rejected')
    check(completed_bound({})['source_norm']==0, 'empty residual positive control')
    # Exact tight homogeneous coefficient: five copies of direction 0 and three
    # of direction 150 have Q(V)=7 and source mass14, so 14^2=28*7.
    sharp={(1,0):5,(-1,1):3}
    check(sum(c*norm(v) for v,c in sharp.items())**2==28*norm(total(tuple(c*x for x in v) for v,c in sharp.items())), 'sharp homogeneous coefficient')
    return certificates,{'multiplicity_cases':tested,'four_member_hostile':[[*v,c] for v,c in sorted(hostile.items())],
                         'integer_capacity_hostile':[[*v,c] for v,c in sorted(integer_hostile.items())],
                         'sharp_homogeneous_control':[[*v,c] for v,c in sorted(sharp.items())]}


def run():
    global CHECKS
    CHECKS=0
    certificates,small=controls()
    raw=SOURCE.read_bytes().replace(b'\r\n',b'\n')
    data=authenticate_source(json.loads(raw))
    output=[complete_scale(item) for item in data['scales']]
    for item in output:
        print(json.dumps(item['statistics'],sort_keys=True),flush=True)
    first=next(item for item in output if item['extra_packets'])
    rec=first['extra_packets'][0]
    item=next(x for x in data['scales'] if x['statistics']['N']==first['statistics']['N'])
    rows={r[0]:r for r in item['source_rows']}; ids=rec['sources']
    reject(lambda:packet_receipt(ids,rows,rec['depth_horizon']-1,rec['clock_horizon']), 'premature P-depth deletion rejected')
    reject(lambda:packet_receipt(ids,rows,rec['depth_horizon'],rec['clock_horizon']-1), 'premature integer-clock cut rejected')
    reject(lambda:packet_receipt([float(ids[0])]+ids[1:],rows,rec['depth_horizon'],rec['clock_horizon']), 'float source alias rejected')
    reject(lambda:packet_receipt([ids[0]]*len(ids),rows,rec['depth_horizon'],rec['clock_horizon']), 'repeated source rejected')
    reject(lambda:consume(ids,rows,set(ids)), 'existing packet cannot be paid again')
    reject(lambda:packet_receipt(ids,rows,True,rec['clock_horizon']), 'Boolean horizon rejected')
    return {'status':'PROVED finite-alphabet compression bound; FINITE-EXACT native residual completion; uniform energy OPEN',
            'checks':CHECKS,'source_receipt_lf_sha256':sha256(raw).hexdigest(),
            'source_route_dag_sha256':data['route_dag_sha256'],
            'circuits':[{'vectors':vs,'coefficients':cs,'members':sum(cs)} for vs,cs in CIRCUITS],
            'linear_certificates':certificates,'controls':small,'scales':output}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--verify',action='store_true'); args=ap.parse_args()
    payload=run(); text=json.dumps(payload,sort_keys=True,separators=(',',':'))+'\n'
    if args.verify:
        if RESULT.read_text(encoding='utf-8')!=text:
            raise RuntimeError('saved circuit receipt differs')
    else:
        RESULT.write_text(text,encoding='utf-8',newline='\n')
    print('checks',CHECKS,'LF_SHA256',sha256(text.encode()).hexdigest())


if __name__=='__main__':
    main()
