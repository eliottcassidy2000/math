"""Weighted destination-coupled tail gate; consumes frozen prefix states only.

No native trajectory, unseen suffix, or new source census is evaluated.
All arithmetic in the new gate and its optimizer is rational/exact.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path
import argparse
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
GEOMETRY_SHA = 'a6fa97c41fa8d364f174addcaa231da19a1d2a93406a6eb2146ce7a142c74dbf'
PREFIX_SHA = '7731a9a9ddf8e5f0478b8d6455e53c6130870f2004353db2e0fefa345bb509c0'
FRONTIER_SHA = '391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f'
PULSE_SHA = 'c0be85b63ccb3c00d7263b7180c63a88f3c4fcc6d5c4225728578a189c90a4ef'
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def read_pinned(path, expected):
    raw = path.read_bytes().replace(b'\r\n', b'\n')
    if sha256(raw).hexdigest() != expected:
        raise ValueError('frozen dependency hash mismatch: '+str(path))
    return raw


path = HERE.parent/'tail_geometry_20261021c/check.py'
read_pinned(path, GEOMETRY_SHA)
spec = importlib.util.spec_from_file_location('weighted_tail_inherited_geometry', path)
G = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = G
spec.loader.exec_module(G)
State = G.State


def rational(x):
    if type(x) not in (int, Fraction):
        raise ValueError('literal rational required')
    return Fraction(x)


def unit_interval(x):
    x = rational(x)
    if not 0 <= x <= 1:
        raise ValueError('alpha outside [0,1]')
    return x


def pair_fraction(x):
    x = Fraction(x)
    return [x.numerator, x.denominator]


def g(v):
    return G.add(v[:2], v[2:])


def local_lines(state, v, w, degree):
    """Affine lines (intercept,slope), with exact marked maximizing choices.

    These maximize degree times one pair factor. Each unary has weight
    1/degree; its interaction has weight one. The coefficient of alpha is
    retained, rather than separately maximizing its endpoints.
    """
    G.audit_state(state); G.vec(v); G.vec(w); G.integer(degree, 1)
    fv = [(t, G.cross(c, v)) for t, c in state.ports]
    fw = [(t, G.cross(c, w)) for t, c in state.ports]
    tv = sorted([(abs(f), t, f) for t, f in fv], reverse=True)[:2]
    tw = sorted([(abs(f), t, f) for t, f in fw], reverse=True)[:2]
    dv, dw = G.delta(v), G.delta(w)
    qv, qw = G.Q(g(v)), G.Q(g(w))
    k, kg = G.cross(v, w), G.beta(g(v), g(w))
    lines = [(0, 0, (None, None))]
    for i, (d, q, top) in enumerate(((dv, qv, tv), (dw, qw, tw))):
        field, terminal, signed = top[0]
        choice = [None, None]; choice[i] = (terminal, 1 if signed >= 0 else -1)
        lines.append((d+field, -q, tuple(choice)))
    for sv, sw in product((-1, 1), repeat=2):
        field, terminal = max((sv*a+sw*b, t) for (t, a), (_, b) in zip(fv, fw))
        lines.append((dv+dw+field+degree*sv*sw*k,
                      -qv-qw-degree*sv*sw*kg,
                      ((terminal, sv), (terminal, sw))))
    choices = [(a+b, t, u, fa, fb) for a, t, fa in tv for b, u, fb in tw if t != u]
    if choices:
        field, t, u, fa, fb = max(choices)
        lines.append((dv+dw+field, -qv-qw,
                      ((t, 1 if fa >= 0 else -1), (u, 1 if fb >= 0 else -1))))
    return tuple(lines)


def local_value(state, v, w, degree, alpha, choices):
    """Independent literal factor evaluator, used by bounded controls."""
    G.audit_state(state); G.vec(v); G.vec(w); G.integer(degree, 1)
    alpha = unit_interval(alpha)
    if type(choices) is not tuple or len(choices) != 2:
        raise ValueError('two marked choices')
    ports = dict(state.ports); value = Fraction(0)
    for b, choice in zip((v, w), choices):
        if choice is None:
            continue
        if type(choice) is not tuple or len(choice) != 2:
            raise ValueError('marked port/sign')
        t, sign = choice
        G.integer(t)
        if type(sign) is not int or sign not in (-1, 1) or t not in ports:
            raise ValueError('illegal literal port/sign')
        value += G.delta(b)-alpha*G.Q(g(b))+sign*G.cross(ports[t], b)
    if choices[0] is not None and choices[1] is not None and choices[0][0] == choices[1][0]:
        value += degree*choices[0][1]*choices[1][1]*(G.cross(v, w)-alpha*G.beta(g(v), g(w)))
    return value


def upper_hull(lines):
    """Upper affine envelope on the real line; start=None denotes -infinity."""
    by_slope = {}
    for b, m, _ in lines:
        by_slope[m] = max(b, by_slope.get(m, b))
    hull = []
    for m, b in sorted(by_slope.items()):
        start = None
        while hull:
            old_b, old_m, old_start = hull[-1]
            start = Fraction(old_b-b, m-old_m)
            if old_start is None or start > old_start:
                break
            hull.pop()
        if not hull:
            start = None
        hull.append((b, m, start))
    return tuple(hull)


def compile_factors(state):
    G.audit_state(state)
    n = len(state.pending)
    if n < 2:
        raise ValueError('pair-factor compiler requires at least two pending nodes')
    bank = Counter(v for _, v in state.pending)
    factors = []
    for v, w in combinations_with_replacement(sorted(bank), 2):
        count = bank[v]*(bank[v]-1)//2 if v == w else bank[v]*bank[w]
        if count:
            lines = local_lines(state, v, w, n-1)
            factors.append((v, w, count, lines, upper_hull(lines)))
    return tuple(factors)


def _raw_from_factors(state, alpha, factors):
    G.audit_state(state); alpha = unit_interval(alpha)
    base = G.physical(v for _, v in state.ports)
    n = len(state.pending)
    if n == 0:
        return Fraction(base)
    if n == 1:
        b = state.pending[0][1]
        return Fraction(base+max(0, G.delta(b)+G.support(state, b)))
    D = G.zonotope(tuple(g(v) for _, v in state.pending))[0]
    residual = sum(count*max(b+alpha*m for b, m, _ in lines) for _, _, count, lines, _ in factors)
    return base+alpha*D+residual/(n-1)


def raw_bound(state, alpha):
    G.audit_state(state); alpha=unit_interval(alpha)
    factors=compile_factors(state) if len(state.pending)>=2 else ()
    return _raw_from_factors(state,alpha,factors)


def optimize(state):
    """Exact optimum of this one-parameter relaxation on [0,1], not of routes."""
    G.audit_state(state)
    n = len(state.pending)
    if n < 2:
        value = raw_bound(state, 0)
        return {'state_sha256': G.state_hash(state), 'alpha': [0, 1],
                'rational_bound': pair_fraction(value), 'integer_bound': value.numerator//value.denominator,
                'pending_nodes': n, 'pair_factors': 0, 'type_pair_factors': 0,
                'interior_breakpoints': 0, 'subgradient_numerators': [0, 0],
                'subgradient_denominator': 1}
    factors = compile_factors(state)
    D = G.zonotope(tuple(g(v) for _, v in state.pending))[0]
    slope = (n-1)*D; events = defaultdict(int)
    for _, _, count, _, hull in factors:
        current = hull[0]
        for row in hull[1:]:
            if row[2] <= 0:
                current = row
        slope += count*current[1]
        for old, new in zip(hull, hull[1:]):
            if 0 < new[2] < 1:
                events[new[2]] += count*(new[1]-old[1])
    alpha = Fraction(0)
    if slope < 0:
        alpha = Fraction(1)
        for x, jump in sorted(events.items()):
            slope += jump
            if slope >= 0:
                alpha = x
                break
    left = right = (n-1)*D
    for _, _, count, lines, _ in factors:
        value = max(b+alpha*m for b, m, _ in lines)
        slopes = [m for b, m, _ in lines if b+alpha*m == value]
        left += count*min(slopes); right += count*max(slopes)
    if not ((alpha == 0 and right >= 0) or (alpha == 1 and left <= 0) or (0 < alpha < 1 and left <= 0 <= right)):
        raise ValueError('convex optimality certificate failed')
    value = _raw_from_factors(state, alpha, factors)
    return {'state_sha256': G.state_hash(state), 'alpha': pair_fraction(alpha),
            'rational_bound': pair_fraction(value), 'integer_bound': value.numerator//value.denominator,
            'pending_nodes': n, 'pair_factors': n*(n-1)//2, 'type_pair_factors': len(factors),
            'interior_breakpoints': len(events), 'subgradient_numerators': [left, right],
            'subgradient_denominator': n-1}


def audit_certificate(state, certificate):
    G.audit_state(state)
    if type(certificate) is not dict:
        raise ValueError('source-bound dictionary required')
    expected = optimize(state)
    if G.enc(certificate) != G.enc(expected):
        raise ValueError('weighted certificate/source mismatch')
    return expected['integer_bound']


def allocation_bound(state, alpha, weights):
    """A general feasible unary allocation, with each pair interaction once.

    Rows (u,v,weight_at_u,weight_at_v), u<v, include EVERY unordered node pair.
    This is an evaluator, not an assertion that these weights solve the LP.
    """
    G.audit_state(state); alpha=unit_interval(alpha)
    if type(weights) is not tuple or len(state.pending)<2:
        raise ValueError('tuple allocation and at least two pending nodes')
    bank=dict(state.pending); wanted=set(combinations(sorted(bank),2)); found=set()
    mass={u:Fraction(0) for u in bank}; total=Fraction(0)
    ports=dict(state.ports)
    for row in weights:
        if type(row) is not tuple or len(row)!=4:
            raise ValueError('literal four-entry allocation row')
        u,v,a,b=row; G.integer(u); G.integer(v); a=rational(a); b=rational(b)
        if (u,v) not in wanted or (u,v) in found or a<0 or b<0:
            raise ValueError('pair address/positive allocation')
        found.add((u,v)); mass[u]+=a; mass[v]+=b
        x,y=bank[u],bank[v]
        dx,dy=G.delta(x)-alpha*G.Q(g(x)),G.delta(y)-alpha*G.Q(g(y))
        k=G.cross(x,y)-alpha*G.beta(g(x),g(y))
        fx=[(t,G.cross(c,x)) for t,c in state.ports]
        fy=[(t,G.cross(c,y)) for t,c in state.ports]
        tx=sorted(((abs(f),t) for t,f in fx),reverse=True)[:2]
        ty=sorted(((abs(f),t) for t,f in fy),reverse=True)[:2]
        best=max(0,a*(dx+tx[0][0]),b*(dy+ty[0][0]))
        for (t,f),(_,h) in zip(fx,fy):
            for sx,sy in product((-1,1),repeat=2):
                best=max(best,a*(dx+sx*f)+b*(dy+sy*h)+sx*sy*k)
        for f,t in tx:
            for h,z in ty:
                if t!=z:
                    best=max(best,a*(dx+f)+b*(dy+h))
        total+=best
    if found!=wanted or any(v!=1 for v in mass.values()):
        raise ValueError('every unary must be allocated exactly once')
    return G.physical(ports.values())+alpha*G.zonotope(tuple(g(v) for v in bank.values()))[0]+total


def rejects(call):
    try:
        call()
    except (ValueError, TypeError):
        need(True, 'hostile rejected')
    else:
        need(False, 'hostile accepted')


def controls():
    # Literal local factors: exhaustive ports/signs versus the eight-line decoder.
    vectors = ((0,0,1,0),(0,0,0,1),(0,0,-1,1),(1,0,-1,0),(-1,0,-1,-1),(-1,1,1,-1))
    states = [State(((0,(0,0,0,0)),), ()),
              State(((0,(0,0,1,0)),(1,(0,0,0,1))), ())]
    for state in states:
        options = (None,)+tuple((t,s) for t,_ in state.ports for s in (-1,1))
        for v,w in product(vectors, repeat=2):
            for degree in (1,2,5):
                lines = local_lines(state,v,w,degree)
                hull = upper_hull(lines)
                for alpha in (0,Fraction(1,3),Fraction(2,3),1):
                    expected = max(local_value(state,v,w,degree,alpha,c) for c in product(options,repeat=2))
                    need(expected == max(b+alpha*m for b,m,_ in lines), 'complete pair-factor maximum')
                    need(expected == max(b+alpha*m for b,m,_ in hull), 'convex envelope exact')
                    for b,m,c in lines:
                        need(local_value(state,v,w,degree,alpha,c) == b+alpha*m, 'marked maximizing-line witness')
    # Full small linked/mixed assignments. This checks optional nonarrival, not only sign vertices.
    small = vectors[:4]
    for blocks in product(small, repeat=3):
        state = State(states[1].ports, tuple((10+i,v) for i,v in enumerate(blocks)))
        options = (None,)+tuple((t,s) for t,_ in state.ports for s in (-1,1))
        actual = max(G.assignment(state,c) for c in product(options,repeat=3))
        cert = optimize(state); alpha = Fraction(*cert['alpha'])
        need(actual <= cert['integer_bound'], 'whole mixed/optional assignment dominated')
        for x in (0,Fraction(1,4),Fraction(1,2),Fraction(3,4),1):
            need(raw_bound(state,alpha) <= raw_bound(state,x), 'exact optimum beats fixed rational trials')
        audit_certificate(state,cert)
    # Empty/singleton and source-bound literal types; negative delta cannot force arrival.
    negative = State(((0,(0,0,0,0)),), ((1,(-5,0,8,0)),))
    need(optimize(negative)['integer_bound'] == 0, 'nonarrival retains zero over negative active increment')
    need(optimize(states[1])['integer_bound'] == G.physical(v for _,v in states[1].ports), 'empty exact')
    rejects(lambda: raw_bound(negative,True)); rejects(lambda: raw_bound(negative,0.5))
    rejects(lambda: raw_bound(negative,Fraction(-1,2)))
    rejects(lambda: raw_bound(negative,Fraction(3,2)))
    rejects(lambda: local_lines(states[1],vectors[0],vectors[1],True))
    bad = optimize(negative); bad['integer_bound'] = 0.0
    rejects(lambda: audit_certificate(negative,bad))
    rejects(lambda: audit_certificate(State(negative.ports,((2,negative.pending[0][1]),)),optimize(negative)))
    # Same scalar supports/coupling do not retain the literal terminal coincidence.
    v,w=(0,0,1,0),(0,0,0,1)
    same=State(((0,(0,0,3,3)),(1,(0,0,0,0))),((10,v),(11,w)))
    split=State(((0,(0,0,1,-1)),(1,(0,0,-1,1))),((10,v),(11,w)))
    need(all(G.support(s,v)==G.support(s,w)==3 for s in (same,split)), 'identical marginal field supports')
    local=[max(b for b,_,_ in local_lines(s,v,w,1)) for s in (same,split)]
    need(local==[7,9], 'literal joint terminal fields differ despite same marginal supports')
    # Fixed-alpha adaptive LP positive control. A full actual maximizing choice is
    # simultaneously a lower bound for every sound upper relaxation, proving LP optimality here.
    adaptive=State(states[1].ports,((10,(0,0,1,0)),(11,(1,0,-1,0)),(12,(1,0,-1,0))))
    weights=((10,11,1,0),(10,12,0,0),(11,12,1,1))
    options=(None,)+tuple((t,s) for t,_ in adaptive.ports for s in (-1,1))
    witness=max(product(options,repeat=3),key=lambda c:G.assignment(adaptive,c))
    exact=G.assignment(adaptive,witness); allocated=allocation_bound(adaptive,0,weights)
    need(exact==allocated==5 and optimize(adaptive)['integer_bound']==6, 'adaptive allocation exact LP optimum5 versus uniform6')
    uniform=tuple((u,v,Fraction(1,2),Fraction(1,2)) for u,v in combinations((10,11,12),2))
    for alpha in (0,Fraction(1,3),1):
        need(allocation_bound(adaptive,alpha,uniform)==raw_bound(adaptive,alpha),'general evaluator reproduces uniform allocation')
    rejects(lambda: allocation_bound(adaptive,0,weights[:-1]))
    rejects(lambda: allocation_bound(adaptive,0,((10,11,True,0),)+weights[1:]))
    rejects(lambda: allocation_bound(adaptive,0,((10,11,Fraction(1,2),0),)+weights[1:]))
    rejects(lambda: allocation_bound(adaptive,0,weights+(weights[0],)))
    return {'optional_negative': optimize(negative),
            'destination_controls': {'marginal_supports':[3,3], 'same_port_local':local[0], 'split_field_local':local[1]},
            'adaptive_control':{'state_sha256':G.state_hash(adaptive),'weights':weights,'alpha':[0,1],
                                'uniform_bound':6,'adaptive_bound':5,'actual_maximizing_choices':witness,'actual_value':exact}}


def native():
    prefix = json.loads(read_pinned(HERE.parent/'future_tail_blocks_20261021b/receipt.json',PREFIX_SHA))
    frontier = json.loads(read_pinned(HERE.parent/'session18_native_annulus/frontier.json',FRONTIER_SHA))['frontier']
    results=[]
    for row in prefix['cases']:
        N,D=row['N'],row['D']; G.integer(N,1); G.integer(D)
        ports={r[0]:tuple(r[1:5]) for r in row['arrived_ports']}
        need(set(ports)<=set(frontier),'literal frontier ports')
        ports.update({t:(0,0,0,0) for t in frontier if t not in ports})
        state=State(tuple(sorted(ports.items())),tuple((r[0],tuple(r[1:5])) for r in row['pending']))
        witnesses={w['node']:w for w in row['witnesses']}
        for r in row['pending']:
            node,_,_,_,_,count,lo,hi,s,clock=r; w=witnesses[node]
            need(count>0 and lo<=clock<=hi and s<256*N,'carrier/count/clock range')
            need(w['source']==s and w['P_depth']==D and w['integer_clock']==clock and not w['arrived']
                 and sum(w['edge_clocks'])==clock,'full actual primitive prefix clock')
            need(w['mod15']==s%15 and w['mod35']==s%35 and w['mod15']%5==w['mod35']%5,
                 'joint source registers retained')
        cert=optimize(state); audit_certificate(state,cert)
        angular=G.projected_certificate(state)['integer_bound']
        need(all(not any(v[:2]) for _,v in state.pending),'these finite states have pure-new pending aggregates')
        need(raw_bound(state,1)==angular,'pure-new endpoint is the previous angular gate')
        need(cert['integer_bound']<=angular,'optimized weighted gate retains angular endpoint')
        need((cert['integer_bound']<=153*N)==(N>=4096),'positive controls and N7 hostile')
        vectors=[v for _,v in state.pending]; alpha=Fraction(*cert['alpha']); n=len(vectors)
        Dz=G.zonotope(tuple(g(v) for v in vectors))[0]
        L=sum(G.support(state,v) for v in vectors)
        gram_slack=sum(G.Q(g(v)) for v in vectors)+sum(abs(G.cross(v,w)) for v,w in combinations(vectors,2))-Dz
        residual=raw_bound(state,alpha)-G.physical(ports.values())-alpha*Dz
        loss=(1-alpha)*(gram_slack+Dz)+L-residual
        independent_loss=Fraction(0)
        for v,w,count,lines,_ in compile_factors(state):
            baseline=(1-alpha)*(G.Q(g(v))+G.Q(g(w))+(n-1)*abs(G.cross(v,w)))+G.support(state,v)+G.support(state,w)
            local_loss=baseline-max(b+alpha*m for b,m,_ in lines)
            need(local_loss>=0,'each local deficit is lost relaxation slack, not negative physical credit')
            independent_loss+=count*local_loss/(n-1)
        need(independent_loss==loss,'source-typed pair multiplicities conserve local deficits')
        gain=loss-(1-alpha)*gram_slack
        need(G.physical(ports.values())+Dz+L-gain==Fraction(*cert['rational_bound']), 'exact joint-destination deficit accounting')
        results.append({'N':N,'D':D,'budget':153*N,'source_rows_sha256':row['source_rows_sha256'],
                        'prefix_edges_sha256':row['prefix_edges_sha256'],'U':row['envelope'][-1],
                        'angular':angular,'weighted':cert,'gain_over_angular':angular-cert['integer_bound'],
                        'destination_deficit':{'D':Dz,'L':L,'unsigned_gram_slack':gram_slack,
                                               'local_factor_loss':pair_fraction(loss),'net_gain':pair_fraction(gain)},
                        'all_earlier_cuts_pass':row['early_cuts_pass']})
    return results


def pulse_native():
    """Independent typed consumer of supplied finite prefix paths, not discovery."""
    folder=HERE.parent/'lagged_tail_pulses_20261022'
    manifest=json.loads(read_pinned(folder/'receipt.json',PULSE_SHA))
    paths=json.loads(read_pinned(folder/'paths.json',manifest['paths_sha256_lf']))
    exports=json.loads((folder/'relaxations.json').read_bytes())
    prefix=json.loads(read_pinned(HERE.parent/'future_tail_blocks_20261021b/receipt.json',PREFIX_SHA))
    frontier=json.loads(read_pinned(HERE.parent/'session18_native_annulus/frontier.json',FRONTIER_SHA))['frontier']
    V=set(frontier); answer=[]
    need(len(paths)==len(exports)==len(prefix['cases'])==len(manifest['cases']),'complete pulse case bank')
    for old,proof,export,meta in zip(prefix['cases'],paths,exports,manifest['cases']):
        ports={r[0]:tuple(r[1:5]) for r in old['arrived_ports']}
        state=State(tuple((t,ports.get(t,(0,0,0,0))) for t in sorted(V)),
                    tuple((r[0],tuple(r[1:5])) for r in old['pending']))
        need(proof['input_state_sha256']==G.state_hash(state)==meta['input_state_sha256'],'pulse original fixed source')
        bank=dict(state.pending); used=set(); removed=set(); early=[]
        for pulse in proof['pulses']:
            u,v=pulse['late_node'],pulse['early_node']; G.integer(u,1); G.integer(v,1)
            G.integer(pulse['native_delay'],1); G.integer(pulse['integer_clock'],1)
            need(u!=v and u in bank and v in bank and u not in used and v not in used,'disjoint two-fibre pulse owners')
            used.update((u,v)); removed.add(u); early.append(v)
            x=u; clock=0
            need(type(pulse['path']) is list and len(pulse['path'])==pulse['native_delay'],'retained literal delay')
            for edge in pulse['path']:
                need(type(edge) is list and len(edge)==3,'path triple shape')
                source,target,a=edge
                for z in edge:G.integer(z,1)
                need(source==x and source not in V and target not in V,'prefix avoids first-arrival boundary')
                need(source%2==1 and source%8!=5,'primitive native source')
                cost=2 if source%8==1 else 1
                y=(3*source+1)//(2**cost)
                while y%8==5:
                    y=(y-1)//4; cost+=2
                need((target,a)==(y,cost),'independent exact primitive edge and clock')
                x=target; clock+=a
            need(x==v and clock==pulse['integer_clock'],'complete path endpoint and clock')
            need(all(bank[v][j]+(-1)**clock*bank[u][j]==0 for j in range(4)),'opposite full linked channels')
        relaxed=State(state.ports,tuple((u,b) for u,b in state.pending if u not in removed))
        expected={'N':old['N'],'D':old['D'],'kind':'optional pulse relaxation, NOT synchronous native state',
                  'ports':relaxed.ports,'pending':relaxed.pending,'pulse_early_nodes':early,
                  'relaxation_sha256':G.state_hash(relaxed)}
        need(G.enc(expected)==G.enc(export),'source-bound typed pulse export, no synchronous reinterpretation')
        need(meta['relaxation_sha256']==G.state_hash(relaxed),'pulse manifest binding')
        cert=optimize(relaxed); audit_certificate(relaxed,cert)
        angular=G.projected_certificate(relaxed)['integer_bound']
        need(angular==meta['pulse_angular'] and cert['integer_bound']<=angular,'composed pulse and weighted gates')
        answer.append({'N':old['N'],'D':old['D'],'kind':export['kind'],'pulses':len(early),
                       'original_state_sha256':G.state_hash(state),'relaxation_sha256':G.state_hash(relaxed),
                       'original_source_rows_sha256':old['source_rows_sha256'],
                       'original_prefix_edges_sha256':old['prefix_edges_sha256'],
                       'pulse_angular':angular,'weighted':cert,'gain_over_pulse_angular':angular-cert['integer_bound']})
    return answer


def run():
    ctrl=controls(); cases=native(); pulses=pulse_native()
    return {'status':'PROVED weighted linked-channel relaxation; FINITE-EXACT frozen prefixes; uniform stopping OPEN',
            'geometry_code_sha256_lf':GEOMETRY_SHA,'prefix_sha256_lf':PREFIX_SHA,'frontier_sha256_lf':FRONTIER_SHA,
            'pulse_receipt_sha256_lf':PULSE_SHA,'controls':ctrl,'cases':cases,'pulse_cases':pulses,'checks':CHECKS}


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true'); args=parser.parse_args()
    result=run(); raw=G.enc(result); target=HERE/'receipt.json'
    if args.write:
        target.write_bytes(raw)
    else:
        need(target.read_bytes().replace(b'\r\n',b'\n')==raw,'saved exact receipt comparison')
    print(json.dumps({'checks':CHECKS,'sha256_lf':sha256(raw).hexdigest(),
                      'cases':[[r['N'],r['D'],r['weighted']['alpha'],r['weighted']['integer_bound'],r['gain_over_angular']]
                               for r in result['cases']],
                      'pulse_cases':[[r['N'],r['weighted']['alpha'],r['weighted']['integer_bound'],r['gain_over_pulse_angular']]
                                     for r in result['pulse_cases']]}))
