"""Source-owned charge rules and complete eleventh-root port relations.

Run from the repository root with Python3, and repeat with -O.
Dependencies: the two explicitly imported maintained exact engines, SymPy,
and PySAT. No orbit discovery is called by the production dispatcher.
"""
from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
from math import lcm
import hashlib
import importlib.util
import json
from pysat.solvers import Solver

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0

def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)

def need(ok, label):
    if not ok:
        raise ValueError(label)

def engine(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT/'04-computation/experiments'/filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

C = engine('charged_family', 'collatz_charged_terminal_family_20261012a.py')
G = engine('cyclotomic_geometry', 'charge_ports_eleven_20261011.py')

def E(n):
    return (n+1).bit_length()-1

def ceil_div(a, b):
    return -((-a)//b)

@dataclass(frozen=True)
class ExitRule:
    head: tuple
    final_min: int
    drop: int
    potential_power: int = 1


def rule_data(rule):
    """Compile an all-height affine inequality, never a sampled inequality."""
    need(type(rule) is ExitRule, 'exact ExitRule')
    C.word(rule.head)
    need(type(rule.final_min) is int and rule.final_min >= 1, 'positive final lower bound')
    need(type(rule.drop) is int and rule.drop >= 1, 'positive bit drop')
    need(type(rule.potential_power) is int and rule.potential_power in (1, 2), 'rank power 1 or 2')
    p, q, b = C.carrier(rule.head)
    P, B, M = 3*p, 3*b+q, q << rule.final_min
    gap = M-(1 << rule.drop)*P
    need(gap > 0, 'strict all-height slope gap')
    cut = max(3, (1 << (rule.drop+1))-1,
              ceil_div((1 << rule.drop)*B+((1 << rule.drop)-1)*M, gap))
    for j in range(1, len(rule.head)+1):
        pp, qq, bb = C.carrier(rule.head[:j])
        cut = max(cut, ceil_div(3*qq-bb, pp))  # every head state is >=3
    residue = -B*pow(P, -1, M) % M
    first = residue+max(0, ceil_div(cut-residue, M))*M
    fee = sum(a % 2 for a in rule.head)+1
    power = rule.potential_power
    need(fee <= 3*(E(first)**power-(E(first)-rule.drop)**power), 'all-height charge invoice')
    return dict(P=P, Q=q, B=B, modulus=M, residue=residue,
                first=first, cut=cut, charge_cap=fee, drop=rule.drop,
                head=list(rule.head), final_min=rule.final_min, power=power)


def apply_rule(rule, n):
    C.odd(n)
    data = rule_data(rule)
    if n < data['first'] or n % data['modulus'] != data['residue']:
        return None
    y = C.replay(n, rule.head)
    a = C.valuation(3*y+1)
    h = (3*y+1) >> a
    delta = sum(x % 2 for x in rule.head)+a % 2
    need(a >= rule.final_min and (h+1)*(1 << rule.drop) <= n+1, 'native paid exit')
    need(delta <= 3*(E(n)**rule.potential_power-E(h)**rule.potential_power), 'exact source invoice')
    return h, rule.head+(a,), delta

RUNS = [ExitRule((1, 2), 3, 1), ExitRule((1, 2, 2), 3, 1),
        ExitRule((1, 2, 2, 2), 2, 1)]
R127 = ExitRule((1,)*6+(2, 4, 3, 3), 3, 3)
Q127 = ExitRule((1,)*6+(2, 4, 3), 3, 2, 2)
R135 = ExitRule((1, 1, 2), 4, 1)
Q159 = ExitRule((1,)*4+(2,)+(1,)*4+(4, 2, 2), 4, 1, 2)


def inverse_tables():
    tables = {}
    count = 0
    for length in range(1, 9):
        p = 3**length
        for w in C.compositions(length, (p-1).bit_length()-1):
            _, q, b = C.carrier(w)
            k = 0
            while q*4**k < p:
                c, d = q*4**k, q*(4**k-1)//3-b
                residue = -d*pow(c, -1, p) % p
                tables.setdefault((length, residue), []).append((c, d, -sum(a % 2 for a in w)))
                count += 1
                k += 1
    check(count == 1045, 'complete inherited inverse template universe')
    return tables


def grammar(limit, rules=(), power=1):
    """Proof DAG from supplied seeds and checked smaller dependencies only."""
    C.natural(limit)
    need(type(power) is int and power in (1, 2), 'exact rank power')
    need(all(r.potential_power <= power for r in rules), 'compatible rule budget')
    raw = json.loads((ROOT/'05-knowledge/results/charge_ports_eleven_20261011.json').read_text())
    words = {int(n): tuple(w) for n, w in raw['controllers']['bank_words'].items()}
    need(set(words) == set(range(1, 30, 2)), 'exact inherited seed bank')
    ledger = {}
    for n, w in words.items():
        need(C.replay(n, w) == 1, 'supplied strict seed')
        ledger[n] = (None, None, 'seed', sum(a % 2 for a in w))
    tables = inverse_tables()
    compiled = [(r, rule_data(r)) for r in rules]
    for n in range(31, limit+1, 2):
        choices = []
        if n % 4 == 1:
            a = C.valuation(3*n+1)
            choices.append(((3*n+1) >> a, a % 2, 'forward'))
        if n % 8 == 5:
            choices.append(((n-1)//4, 0, 'sibling'))
        k = C.valuation(n+1)
        if k >= 2:
            t = (n+1) >> k
            r = C.valuation(3**k*t-1)
            if r >= 2:
                choices.append(((n-1)//2, 1, 'reset'))
            else:
                y = (3**k*t-1)//2
                a = C.valuation(3*y+1)
                y = (3*y+1) >> a
                if a == 6:
                    choices.append(((n-1)//2, -1, 'sporadic26'))
                if a == 2 and C.valuation(3*y+1) == 10:
                    choices.append(((n-1)//2, -1, 'sporadic2210'))
        for length in range(1, 9):
            p = 3**length
            for c, d, charge in tables.get((length, n % p), ()):
                h = (c*n+d)//p
                if 0 < h < n:
                    choices.append((h, charge, 'inverse'))
        if n % 1024 == 43:
            choices.append(((81*n+101)//512, 1, '43-cylinder'))
        if n % 1024 == 39:
            choices.append(((243*n+251)//512, 3, '39-cylinder'))
        for function, label in [(C.coarse75, '75-coarse'), (C.coarse123, '123-coarse')]:
            receipt = function(n)
            if receipt:
                h, _, charge = receipt
                choices.append((h, charge, label))
        for index, (rule, data) in enumerate(compiled):
            if n >= data['first'] and n % data['modulus'] == data['residue']:
                h, _, charge = apply_rule(rule, n)
                choices.append((h, charge, 'new-'+str(index)))
        for h, delta, _ in choices:
            check(0 < h < n and h % 2 == 1, 'strict smaller odd dependency')
            check(delta <= 3*(E(n)**power-E(h)**power), 'policy invoice')
        paid = [r for r in choices if r[0] in ledger]
        if paid:
            h, delta, name = min(paid, key=lambda r: (ledger[r[0]][3]+r[1], r))
            ledger[n] = (h, delta, name, ledger[h][3]+delta)
    return ledger


def control_orbit(n):
    """CONTROL ONLY: bounded discovery, never called by grammar/apply_rule."""
    word, states = [], [n]
    while n != 1 and len(word) < 10000:
        a = C.valuation(3*n+1)
        n = (3*n+1) >> a
        word.append(a)
        states.append(n)
    check(n == 1, 'declared control orbit cap')
    return tuple(word), states


def collatz_controls():
    rules = RUNS+[R127, Q127, R135, Q159]
    for rule in rules:
        d = rule_data(rule)
        for t in list(range(129))+[2**100+17]:
            n = d['first']+d['modulus']*t
            h, w, delta = apply_rule(rule, n)
            check(C.replay(n, w) == h, 'literal native word')
            check(delta == sum(a % 2 for a in w), 'literal odd charge')
            check((h+1)*(1 << rule.drop) <= n+1, 'symbolically compiled inequality control')
        check(apply_rule(rule, d['first']+2) is None, 'wrong-source rejection')
    # Independent fixed rational carriers and coarse guards, not orbit sampling.
    check(rule_data(R127)['residue'] == 127 and rule_data(R127)['modulus'] == 2**21,
          '127 linear cylinder')
    check(apply_rule(R127, 127) == (11, (1,)*6+(2,4,3,3,3), 9), '127 paid seed route')
    check(apply_rule(Q127, 127)[0] == 29, 'quadratic rule stops at first band')
    check([(rule_data(r)['residue'], rule_data(r)['modulus']) for r in RUNS]
          == [(11,64),(43,256),(171,512)], 'folded even-run guards')
    for r in range(1, 65):
        for t in (1,7,13,61):
            n = (2**(2*r+3)*t+1)//3
            a = 2+C.valuation(1+3**(r+1)*t)
            h = (1+3**(r+1)*t) >> (a-2)
            check(C.replay(n, (1,)+(2,)*r+(a,)) == h, 'unbounded-run formula control')
            check(2*(h+1) <= n+1 and 1+a%2 <= 2, 'run invoice')
            check(any(apply_rule(rule, n) for rule in RUNS), 'bounded receiver covers run family')
    # Small exhaustive compiler universe, with all-height proof separate.
    compiler_count = 0
    for length in range(4):
        for w in product(range(1,5), repeat=length):
            for a in range(1,8):
                for drop in (1,2,3):
                    r = ExitRule(w, a, drop, 2)
                    try:
                        d = rule_data(r)
                    except ValueError:
                        continue
                    for t in (0,1,7):
                        n = d['first']+d['modulus']*t
                        h, ww, fee = apply_rule(r, n)
                        check(C.replay(n, ww) == h and fee <= 3*(E(n)**2-E(h)**2), 'general compiler control')
                    compiler_count += 1
    for bad in [ExitRule((True,),3,1), ExitRule((1,),3.0,1), ExitRule((1,),3,True),
                ExitRule((1,),3,1,True), ExitRule((1,)*6,1,3), ExitRule((1,)*6+(2,4,3),3,2)]:
        try:
            rule_data(bad)
        except ValueError:
            check(True, 'bad or unpaid rule rejected')
        else:
            check(False, 'bad rule accepted')
    # All-depth shield uses a finite union of authenticated complete paths.
    pool = set()
    seed_words = {}
    for n in range(1,127,2):
        w, states = control_orbit(n)
        seed_words[str(n)] = list(w)
        pool.update(states)
    siblings = []
    x = 127
    while x <= max(pool):
        siblings.append(x)
        x = 4*x+1
    check((len(pool),max(pool),siblings) == (104,3077,[127,509,2037]), 'exact inverse shield universe')
    check(not pool.intersection(siblings), 'all-depth smaller-source sibling exclusion')
    # Least rate-three hostile, independent bounded replay of every smaller odd.
    hostile = None
    for n in range(1,704,2):
        w, _ = control_orbit(n)
        if sum(a % 2 for a in w) > 3*E(n)+16 and hostile is None:
            hostile = (n, sum(a % 2 for a in w), E(n))
    check(hostile == (703,44,9), 'least linear-rate-three hostile in declared universe')
    policies = [('incoming',(),1),('127',(R127,),1),
                ('folded',tuple(RUNS+[R127,R135]),1),
                ('quadratic',tuple(RUNS+[R127,Q127,R135,Q159]),2)]
    results = {}
    for label, extra, power in policies:
        ledger = grammar(8191,extra,power)
        if label == 'incoming':
            check(ledger == C.finite_grammar(8191,True,True), 'independent inherited dispatcher comparison')
        for n, row in ledger.items():
            w, _ = control_orbit(n)
            check(sum(a % 2 for a in w) == row[3], 'control-only replay of compiled exact charge')
            budget = 3*E(n)+16 if power == 1 else 3*(E(n)**2-1)
            check(row[3] <= budget, 'global telescoped source charge')
        pending = [n for n in range(1,8192,2) if n not in ledger]
        results[label] = {'accepted':len(ledger),'first_pending':pending[:12],
                          'rank_power':power,'ledger':[[n,*ledger[n]] for n in sorted(ledger)]}
    check([results[k]['accepted'] for k in ('incoming','127','folded')] == [1278,1362,1650], 'matched policy gains')
    w2047, states2047 = control_orbit(2047)
    j = states2047.index(127)
    check((j,sum(w2047[:j]),sum(a % 2 for a in w2047[:j])) == (41,69,27), '2047 actual marked prefix')
    return {'rule_plans':[rule_data(r) for r in rules], 'compiler_plans_tested':compiler_count,
            'shield127':{'pool':sorted(pool),'smaller_source_words':seed_words,'siblings_below_pool_max':siblings},
            'least_linear_rate3_hostile':list(hostile),'hostile703_word':list(control_orbit(703)[0]),
            '2047_to_127':list(w2047[:j]),'policies':results}


def relative_core(n, edges, ports, k):
    need(type(n) is int and n >= 0 and type(k) is int and k >= 1, 'graph size/palette')
    need(len(set(ports)) == len(ports) and all(type(x) is int and 0 <= x < n for x in ports), 'distinct marked ports')
    adj = [set() for _ in range(n)]
    for u,v in edges:
        need(type(u) is int and type(v) is int and 0 <= u < v < n, 'simple ordered graph edge')
        adj[u].add(v);adj[v].add(u)
    frozen=set(ports);live=set(range(n));degrees=list(map(len,adj))
    queue=deque(v for v in range(n) if v not in frozen and degrees[v] < k)
    order=[]
    while queue:
        v=queue.popleft()
        if v not in live:
            continue
        need(v not in frozen and len(adj[v]&live) < k, 'actual paid elimination')
        live.remove(v);order.append(v)
        for u in sorted(adj[v]&live):
            degrees[u]-=1
            if u not in frozen and degrees[u] < k:
                queue.append(u)
    return sorted(live),order,adj


def partitions(size, k, edges, prefix=()):
    if len(prefix) == size:
        yield prefix
        return
    for c in range(min(k, max(prefix,default=-1)+2)):
        if any(v == len(prefix) and prefix[u] == c for u,v in edges):
            continue
        yield from partitions(size,k,edges,prefix+(c,))


def boundary_packet(points, edges, palette=6):
    n=len(points);ports=list(range(8))
    live,order,adj=relative_core(n,edges,ports,palette)
    ids={v:i for i,v in enumerate(live)}
    pe=[(u,v) for u,v in edges if u < 8 and v < 8]
    core_edges=[(ids[u],ids[v]) for u,v in edges if u in ids and v in ids]
    clauses=[[palette*i+c+1 for c in range(palette)] for i in range(len(live))]
    clauses += [[-(palette*i+c+1),-(palette*i+d+1)] for i in range(len(live))
                for c in range(palette) for d in range(c)]
    clauses += [[-(palette*u+c+1),-(palette*v+c+1)] for u,v in core_edges for c in range(palette)]
    counts=Counter();lifts=[]
    with Solver(name='g3',bootstrap_with=clauses) as solver:
        for p in partitions(8,palette,pe):
            check(solver.solve(assumptions=[palette*ids[v]+c+1 for v,c in enumerate(p)]),
                  'every bare boundary partition has a lift')
            model=set(solver.get_model())
            colour={v:next(c for c in range(palette) if palette*i+c+1 in model) for v,i in ids.items()}
            check(tuple(colour[v] for v in ports) == p, 'original labelled ports fixed')
            core_string=''.join(str(colour[v]) for v in live)
            for v in reversed(order):
                colour[v]=next(c for c in range(palette) if all(colour.get(u) != c for u in adj[v]))
            for u,v in edges:
                check(colour[u] != colour[v], 'complete geometric edge in restored witness')
            counts[max(p)+1]+=1
            lifts.append([list(p),core_string])
    check(dict(counts) == {3:30,4:180,5:247,6:113}, 'complete canonical partition universe')
    return {'ports':ports,'port_edges':pe,'core_vertices':live,'elimination_order':order,
            'partition_counts':dict(counts),'all_570_core_lifts':lifts}


def graph_record(points, edges, den, survivors):
    coords=[[int(x*den) for x in p] for p in points]
    return {'vertices':len(points),'edges':len(edges),'denominator':den,'filter_survivors':survivors,
            'coordinate_hash':hashlib.sha256(json.dumps(coords).encode()).hexdigest(),
            'edge_hash':hashlib.sha256(json.dumps(edges).encode()).hexdigest()}


def geometry_controls():
    # Reconstruct the previous context from its saved exact generator, not scratch files.
    old=json.loads((ROOT/'05-knowledge/results/charge_ports_eleven_20261011.json').read_text())
    u=tuple(Q(x) for x in old['ports']['mixed_generator'])
    directions=[G.powers[11*j] for j in range(5)]+[G.powers[5*j] for j in range(11)]
    directions += [u,G.mul(u,G.conj(G.a)),G.mul(u,G.conj(G.g))]
    directions=list(dict.fromkeys(directions+[G.sc(v,-1) for v in directions]))
    D=lcm(*(Q(x).denominator for p in directions+G.ports for x in p))
    dd=[tuple(int(x*D) for x in p) for p in directions]
    ball={G.O};shell={G.O}
    for _ in range(3):
        shell={G.add(x,y) for x,y in product(shell,dd)}-ball;ball|=shell
    points=list(dict.fromkeys(G.ports+[tuple(Q(x,D) for x in p) for p in sorted(ball)]))
    edges,surv,den=G.graph(points)
    rec=graph_record(points,edges,den,surv)
    check(rec['coordinate_hash'] == old['ports']['graph']['coordinate_hash'], 'same original mixed context')
    check(rec['edge_hash'] == old['ports']['graph']['edge_hash'], 'all original geometric edges')
    old_boundary=boundary_packet(points,edges)
    check(len(old_boundary['core_vertices']) == 874, 'relative-six-core reduction')
    # New denominator-three connector and a real eleventh-root assembly.
    v=G.sc(G.add(G.sc(G.I,-5),G.b),Q(1,6))
    check(G.mul(v,G.conj(v)) == G.I, 'new unit exact norm')
    check(G.add(G.add(G.sc(G.I,2),v),G.conj(v)) == G.sc(G.I,Q(1,3)), 'four-edge thirds connector')
    path=[G.O,G.I,G.sc(G.I,2),G.add(G.sc(G.I,2),v),G.sc(G.I,Q(1,3))]
    path_edges,_,_=G.graph(path)
    check(path_edges == [(0,1),(1,2),(2,3),(3,4)], 'complete induced connector P5')
    check(G.add(G.sc(G.I,5),G.sc(G.add(v,G.conj(v)),3)) == G.O, 'non-additive-colour relation')
    # Complete unit lattice A/3+Bv: norm=(2A-5B)^2+11B^2 over36.
    lattice_directions=[]
    c5=[0,1,1,2,0]  # proper colours for edges a -> a+3 modulo5
    for A in range(-12,13):
        for B in range(-4,5):
            actual=(2*A-5*B)**2+11*B*B == 36
            declared=(A,B) in [(3,0),(-3,0),(0,1),(0,-1),(5,1),(-5,-1)]
            check(actual == declared,'full lattice norm equation, exterior bound |B|<=1')
            if actual:
                lattice_directions.append((A,B))
    for A,B in product(range(15),range(3)):
        for a,b in lattice_directions:
            check((c5[A%5]+B)%3 != (c5[(A+a)%5]+B+b)%3, 'nonlinear all-lattice three-colouring')
    # A sparse modular factor gives the same simple K unit rotated out of K.
    gamma=G.g
    numerator=G.O
    for j,c in [(7,1),(5,-2),(4,1),(2,1),(1,1),(0,1)]:
        numerator=G.add(numerator,G.sc(G.powers[5*j%55],c))
    u3=G.sc(numerator,Q(1,3))
    check(u3 == G.mul(G.powers[25],G.conj(v)), 'degree-five factor simplifies to rotated K unit')
    check(G.mul(u3,G.conj(u3)) == G.I, 'eleventh-root denominator-three unit')
    x=G.z
    f1=x**5-x**3+x**2-x-1
    f2=x**5+x**4-x**3+x**2-1
    phi11=G.s.cyclotomic_poly(11,x)
    check(all(int(c)%3 == 0 for c in G.s.Poly(f1*f2-phi11,x).all_coeffs()),
          'complete modular factorization at three')
    check(G.s.resultant(phi11,f1,x) == 3**5, 'exact cyclotomic factor norm')
    fx=G.O
    for j,c in [(5,1),(3,-1),(2,1),(1,-1),(0,-1)]:
        fx=G.add(fx,G.sc(G.powers[5*j%55],c))
    check(G.mul(u3,G.conj(fx)) == fx, 'actual factor-to-unit ratio')
    patch=[G.add(G.sc(G.I,Q(a,3)),G.sc(v,b)) for a in range(-4,5) for b in range(-1,2)]
    rotations=[G.powers[5*j] for j in range(11)]
    fresh=list(dict.fromkeys(G.ports+[G.add(p,G.mul(r,x)) for p in G.ports for r in rotations for x in patch]))
    fresh_edges,surv2,den2=G.graph(fresh)
    fresh_record=graph_record(fresh,fresh_edges,den2,surv2)
    check((len(fresh),len(fresh_edges)) == (2278,6405), 'new complete eleventh-root thirds context')
    colouring=G.proper_colour(len(fresh),fresh_edges,3)
    fresh_boundary=boundary_packet(fresh,fresh_edges)
    check(fresh_boundary['core_vertices'] == list(range(8)), 'new context peels entirely to original ports')
    # Their actual geometric union has extra cross-edges. Authenticate those
    # too; a support intersection at the eight ports would miss them.
    joined=list(dict.fromkeys(points+fresh))
    index={p:i for i,p in enumerate(joined)}
    inherited={tuple(sorted((index[points[i]],index[points[j]]))) for i,j in edges}
    inherited |= {tuple(sorted((index[fresh[i]],index[fresh[j]]))) for i,j in fresh_edges}
    joined_edges,surv3,den3=G.graph(joined)
    joined_set=set(joined_edges)
    check(inherited <= joined_set, 'both original complete contexts retained')
    new_edges=sorted(joined_set-inherited)
    check((len(joined),len(joined_edges),len(new_edges)) == (10598,30506,88),
          'complete changing-generator union')
    live3,order3,_=relative_core(len(joined),joined_edges,list(range(8)),6)
    check(live3 == old_boundary['core_vertices'], 'joined context retains exactly the old marked core')
    core_set=set(live3)
    check({(i,j) for i,j in joined_edges if i in core_set and j in core_set}
          == {(i,j) for i,j in edges if i in core_set and j in core_set},
          'old core has exactly the same complete edge relation')
    joined_colouring=G.proper_colour(len(joined),joined_edges,3)
    joined_record=graph_record(joined,joined_edges,den3,surv3)
    joined_record.update(new_cross_edges=new_edges,core_vertices=live3,
                         elimination_order=order3,three_colouring=joined_colouring,
                         shared_vertices=len(points)+len(fresh)-len(joined))
    # Independent unfiltered geometry on the connector, ports and one small patch.
    small=list(dict.fromkeys(G.ports+path+patch))
    se,_,_=G.graph(small);se=set(se)
    for i,j in combinations(range(len(small)),2):
        difference=G.add(small[i],G.sc(small[j],-1))
        check(((i,j) in se) == (G.mul(difference,G.conj(difference)) == G.I), 'unfiltered complete norm control')
    # A hostile core control: all six distinct precolours on star leaves block the hub.
    star=[(i,6) for i in range(6)]
    core,_,_=relative_core(7,star,list(range(6)),6)
    check(core == list(range(7)), 'essential degree-six obstruction is retained')
    core,_,_=relative_core(6,[(i,5) for i in range(5)],list(range(5)),6)
    check(core == list(range(5)), 'degree-five hub always extends')
    return {'old_context':rec,'old_boundary':old_boundary,
            'thirds_context':fresh_record,'thirds_boundary':fresh_boundary,'thirds_three_colouring':colouring,
            'joined_context':joined_record,
            'v_coefficients':[str(x) for x in v],'u3_coefficients':[str(x) for x in u3],
            'lattice_periodic_colouring':c5,'lattice_directions':lattice_directions}


def number_controls():
    check(2**11-1 == 23*89 and 127 == 10**2+3**3, 'literal requested identities')
    for p in (23,89):
        check(pow(2,11,p) == 1 and all(pow(2,k,p) != 1 for k in range(1,11)), 'exact order eleven')
    units=[G.unit(2,j) for j in range(1,11)]
    total=G.O
    for u in units:
        check(G.mul(u,G.conj(u)) == G.I, 'Mersenne norm direction')
        total=G.add(total,u)
    check(total == G.sc(G.I,Q(14849,2047)), 'complete eleventh-root trace')
    check(G.add(G.sc(G.I,1799),G.sc(total,-248)) == G.sc(G.I,Q(1,2047)), 'explicit denominator connector')
    # The attractive 127 rotation belongs to the excluded Eisenstein extension.
    check(73**2+3*60**2 == 127**2, '127 Eisenstein unit norm')
    check(G.s.totient(165) == 2*G.s.totient(55), 'sqrt(-3) is not in the required field')
    return {'Mersenne_factorization':[23,89], 'unit_trace':'14849/2047',
            'connector_identity':'1799 - 248 * sum_(j=1)^10 u_(2,j) = 1/2047',
            'connector_unit_steps':1799+2480,
            '127_rotation':'(73+60*sqrt(-3))/127; outside Q(zeta55)',
            'connection_scope':'exact arithmetic and denominator identities; no dynamical conjugacy'}


if __name__ == '__main__':
    result={'collatz':collatz_controls(),'geometry':geometry_controls(),'numbers':number_controls()}
    result['checks']=CHECKS
    result['inherited_engine_checks']=G.checks
    print(json.dumps(result,indent=2,sort_keys=True))
