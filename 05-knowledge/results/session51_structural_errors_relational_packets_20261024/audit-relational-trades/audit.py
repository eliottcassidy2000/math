#!/usr/bin/env python3
"""Independent Session51 audit: raw native records, DSU filtration, monotone DAG.

No candidate or inherited mathematical checker is imported. Frozen authentic
native-menu and filtration scripts are separately replayed by reproduce.py.
"""
from pathlib import Path
from collections import defaultdict, Counter, deque
from itertools import combinations, product
from hashlib import sha256
from math import comb
import argparse, csv, json

P = Path(__file__).resolve().parent
COUNT = 0
U, V = 1965, 2234
FULL = (1 << 17) - 1

def require(value, why):
    global COUNT
    COUNT += 1
    if not value:
        raise RuntimeError(why)

def read(name):
    return json.loads((P / 'native/inputs' / name).read_text())

def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()

def add(a, b):
    return ((a[0] + b[0]) % 128,) + tuple(x ^ y for x, y in zip(a[1:], b[1:]))

ZERO = (0, 0, 0, 0, 0)

def total(values):
    out = ZERO
    for value in values:
        out = add(out, value)
    return out

def encode(value):
    return 16 * value[0] + sum(x << i for i, x in enumerate(value[1:]))

def source_pins():
    pins = json.loads((P / 'SOURCEPINS.json').read_text())
    names = {f.relative_to(P).as_posix() for d in ('native', 'sources')
             for f in (P / d).rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.suffix != '.pyc'}
    require(names == set(pins), 'exact source closure')
    for name, row in pins.items():
        b = (P / name).read_bytes()
        require(len(b) == row['bytes'] and sha256(b).hexdigest() == row['sha256'], 'source hash ' + name)

def native():
    domain = read('native_packet_decomposition_20261023/native/domains.json')['sources']
    supply = read('compatible_charged_helpers_20261023/supply/selected_supply.json')
    repair = read('guarded_component_completion_20261024/sections/native_2125_receipt.json')
    host = read('guarded_component_completion_20261024/native/xor_tree_receipt.json')
    electric = read('metagraph_order8_electrical_20261023/receipt.json')['order8']
    X = electric['higher_layers'][-1]['lifts'][0]
    Y = electric['radical_lifts_mod4']
    potentials, moduli = [X] + Y, [256] + [4] * len(Y)
    labels = {s: (X[s] % 128,) + tuple(f[s] % 2 for f in Y[1:5]) for s in range(len(X))}
    records = {r['source']: sorted(r['exact_all_shared_clocks_target_domain'], key=lambda t: t['target']) for r in domain}
    core = sorted(records)
    for row in supply['helpers']:
        records[row['source']] = [row['selected']['before'], row['selected']['after']]
    records[2125] = repair['records']
    prices = {s: sum(a != b for a, b in zip(rs[0]['path'], rs[1]['path'])) for s, rs in records.items()}
    base = sorted(s for s in core if prices[s] == 2 and s not in {176, 1425, 2071, 1457})
    triple = [1457, 1841, 3015]
    groups = [[s] for s in base] + [triple]
    require(base == [1394,1477,1480,1549,1638,1691,1921,1974,2176,2185,2581,2837,2905,2945,3000,3369], 'authentic source inventory')
    require(len(records) == 82 and [prices[s] for s in triple] == [2,3,3], 'authentic eighty-two menus and triple clock prices')
    frames = {368: list(range(8)), 776: [8,9,4,10,11,12,13,3], 3139: [14,15,16,17,18,3,4,19], 2125: list(range(20,28))}
    next_owner = 28
    for s in sorted(set(records) - set(frames)):
        frames[s] = list(range(next_owner, next_owner + 8)); next_owner += 8
    require(next_owner == 652 and frames == {int(s): f for s, f in host['source_frames'].items()}, 'independent full source-frame reconstruction')
    global_paths = {s: [[frames[s][v] for v in r['path']] for r in rs] for s, rs in records.items()}
    require(set(frames[1457]).isdisjoint(frames[1394]), 'new relation does not identify private owners')
    table = {(194,162), (191,159)}
    truth = []
    for a, b in product((0,1), repeat=2):
        pair = (global_paths[1457][a][1], global_paths[1394][b][2])
        reflected = (global_paths[1457][a][::-1][6], global_paths[1394][b][::-1][5])
        require(pair == reflected, 'literal same relation on converse-reflected ports')
        require((pair in table) == (a == b), 'actual owner table equivalent to branch bit relation')
        require(pair[0] != pair[1], 'allowed relation is not literal equality')
        truth.append({'bits': [a,b], 'owners': pair, 'accepted': pair in table})
    for s, rs in records.items():
        require([r['target'] for r in rs] == [U,V], 'original full target records')
        require(prices[s] == sum(a != b for a,b in zip(*global_paths[s])), 'clocks are literal per-slot owner disagreements')
    for s in base + triple:
        require(all(2*(f[s]-f[U]) % m == 2*(f[s]-f[V]) % m for f,m in zip(potentials,moduli)), 'all source rows unchanged under native paired switch')
    # Whole reservations are read from literal masks and reconstructed reflections.
    chords = [(i,j) for i in range(8) for j in range(i+2,8)]
    def reflected(mask):
        return sum(((mask >> chords.index((7-j,7-i))) & 1) << k for k,(i,j) in enumerate(chords))
    all_words = set(); words = {}
    for s, rs in records.items():
        for bit, r in enumerate(rs):
            m = r['mask']; rr = reflected(m); full = (1 << 21) - 1
            unit = {m, rr, full ^ m, full ^ rr}
            require(len(unit) == 4 and not unit & all_words, 'disjoint full alternatives, including within-source choices')
            all_words |= unit; words[s,bit] = unit
    require(len(all_words) == 656, 'unchanged 164 reserved units,656 closed words')
    return {'base':base, 'triple':triple, 'groups':groups, 'labels':labels, 'prices':prices,
            'potentials':potentials, 'moduli':moduli, 'records':records, 'global_paths':global_paths,
            'words':words, 'truth':truth, 'table':table, 'core':core,
            'free':sorted(set(records)-set(core)-{368,776,3139})}

def cost(mask):
    b, t = (mask & 65535).bit_count(), mask >> 16
    return b + 3*t, 2*b + 8*t

def states(n):
    weights = [len(g) for g in n['groups']]
    labels = [total(n['labels'][s] for s in g) for g in n['groups']]
    sums = [tuple(sum(f[s] for s in g) % m for f,m in zip(n['potentials'],n['moduli'])) for g in n['groups']]
    complete_sum = tuple(sum(f[s] for g in n['groups'] for s in g) % m for f,m in zip(n['potentials'],n['moduli']))
    grades = [0] * (FULL+1); charges = [ZERO] * (FULL+1); potentials = [(0,)*6] * (FULL+1)
    original = defaultdict(list); restricted = defaultdict(list); signatures = {}
    for mask in range(FULL+1):
        if mask:
            low = mask & -mask; j = low.bit_length()-1; previous = mask ^ low
            grades[mask] = grades[previous] + weights[j]
            charges[mask] = add(charges[previous], labels[j])
            potentials[mask] = tuple((a+b) % m for a,b,m in zip(potentials[previous],sums[j],n['moduli']))
        key = grades[mask], encode(charges[mask])
        rows = tuple((2*((19-grades[mask])*f[U] - (tt-p)) % m,
                      2*(grades[mask]*f[V]-p) % m)
                     for f,m,tt,p in zip(n['potentials'],n['moduli'],complete_sum,potentials[mask]))
        if key in signatures:
            require(rows == signatures[key], 'literal U/V target diagonals and paired factor2 on every grade/charge fibre')
        else:
            signatures[key] = rows
        original[key].append(mask)
        owner_pair = (n['global_paths'][1457][mask >> 16][1], n['global_paths'][1394][mask & 1][2])
        if owner_pair in n['table']:
            restricted[key].append(mask)
    require(sum(map(len,restricted.values())) == 65536, 'all accepted native record assignments')
    return dict(sorted(original.items())), dict(sorted(restricted.items())), grades, charges

def ordinary_dsu(v, matrix):
    """All-pair minimax by edge filtration, not Floyd or a spanning-tree solver."""
    m = len(v); clusters = {j:{j} for j in range(m)}; owner = list(range(m)); answer = [[(0,0)]*m for _ in v]
    for price,i,j in sorted((matrix[i][j],i,j) for i,j in combinations(range(m),2)):
        a,b = owner[i],owner[j]
        if a == b: continue
        A,B = clusters[a],clusters.pop(b)
        for u in A:
            for w in B:
                answer[u][w] = answer[w][u] = price
        for w in B: owner[w] = a
        A.update(B)
    return answer

def strict_dag(v, matrix):
    """A root-relative subset DAG directly realizes strict non-reacquisition."""
    answer = []
    for root,x in enumerate(v):
        diff = [x ^ y for y in v]
        order = sorted(range(len(v)), key=lambda j: diff[j].bit_count())
        best = [(99,999)] * len(v); best[root] = (0,0)
        for j in order:
            if j == root: continue
            best[j] = min(max(best[i],matrix[i][j]) for i in order if i != j and diff[i] & diff[j] == diff[i])
        answer.append(best)
    return answer

def intervals(v):
    return {(i,j): [z for z in v if (z ^ v[i]) & (v[i] ^ v[j]) == (z ^ v[i])]
            for i,j in combinations(range(len(v)),2)}

def route(v, x, y, budget):
    previous = {x:None}; queue = deque([x])
    while queue:
        a = queue.popleft()
        for b in v:
            if b not in previous and all(p <= q for p,q in zip(cost(a ^ b),budget)):
                previous[b] = a; queue.append(b)
    if y not in previous:
        return None, sorted(previous)
    path = [y]
    while previous[path[-1]] is not None:
        path.append(previous[path[-1]])
    return path[::-1], sorted(previous)

def full_route(n, path, require_pin):
    fixed = set(n['free'][:2]) | {776,3139}; arrivals=[]; first_rows=None
    for mask in path:
        active = fixed | {s for j,g in enumerate(n['groups']) if mask >> j & 1 for s in g}
        used=set(); rows={s:[0]*len(n['potentials']) for s in list(n['records'])+[U,V]}
        for s,rs in n['records'].items():
            bit=int(s in active); r=rs[bit]; unit=n['words'][s,bit]
            require(not unit & used, 'literal arrival cannot duplicate a normalized reserved word');used |= unit
            for j,f in enumerate(n['potentials']):
                delta=2*(f[s]-f[r['target']]); rows[s][j]+=delta;rows[r['target']][j]-=delta
        if first_rows is None: first_rows=rows
        else: require(all((rows[s][j]-first_rows[s][j]) % m == 0 for s in rows for j,m in enumerate(n['moduli'])), 'all original physical source/target rows restored at each literal arrival')
        require(len(used)==328, 'exactly82 occupied native inverse/converse units at arrival')
        pair=(n['global_paths'][1457][int(1457 in active)][1],n['global_paths'][1394][int(1394 in active)][2])
        if require_pin: require(pair in n['table'], 'new actual relation holds at every arrival')
        arrivals.append({'mask':mask,'core_V':sorted(active & set(n['core'])),'helpers_V':sorted(active-set(n['core'])),'owner_pair':pair,'relation_accepted':pair in n['table']})
    changes=Counter();events=[]
    for a,b in zip(path,path[1:]):
        changed=[s for j,g in enumerate(n['groups']) if (a^b)>>j&1 for s in g]
        actual=(len(changed),sum(n['prices'][s] for s in changed))
        require(actual==cost(a^b), 'event cost is exact literal changed-owner-clock sum')
        changes.update(changed);events.append({'sources':sorted(changed),'cost':actual})
    return {'arrivals':arrivals,'events':events,'repeated_source_changes':{s:c for s,c in changes.items() if c>1}}

def classify(original, accepted, grades, charges):
    stats={k:{'states':0,'charges':0,'pairs':0,'primitive_pairs':0,'ordinary_clocks':0,'strict_clocks':0} for k in range(20)}
    hist={'ordinary':Counter(),'strict':Counter(),'primitive':Counter()}; trades={}; all_edges=set(); pair_hash=sha256(); fibre_rows=[]; lowers={}
    for (k,q),v in accepted.items():
        matrix=[[cost(x^y) for y in v] for x in v]; ordinary=ordinary_dsu(v,matrix);strict=strict_dag(v,matrix);ints=intervals(v)
        primitive=[ij for ij,ins in ints.items() if len(ins)==2]
        ob=max((ordinary[i][j] for i,j in ints),default=(0,0));sb=max((strict[i][j] for i,j in ints),default=(0,0))
        require(sb==max((matrix[i][j] for i,j in primitive),default=(0,0)), 'strict DAG agrees with independent primitive-uniform criterion')
        r=stats[k];r['states']+=len(v);r['charges']+=1;r['pairs']+=len(ints);r['primitive_pairs']+=len(primitive)
        r['ordinary_clocks']=max(r['ordinary_clocks'],ob[1]);r['strict_clocks']=max(r['strict_clocks'],sb[1]);fibre_rows.append([k,q,len(v),ob[1],sb[1],len(primitive)])
        if ob[1] and (k not in lowers or ob[1]>lowers[k]['cost'][1]):
            i,j=next(ij for ij in ints if ordinary[ij[0]][ij[1]]==ob);x,y=v[i],v[j]
            path,_=route(v,x,y,ob);_,slotcut=route(v,x,y,(ob[0]-1,999));_,clockcut=route(v,x,y,(99,ob[1]-1))
            require(y not in slotcut and y not in clockcut, 'sharp ordinary lower cut in each budget coordinate')
            lowers[k]={'charge':q,'complete_fibre':v,'endpoints':[x,y],'cost':ob,'route':path,'below_slot_cut':slotcut,'below_clock_cut':clockcut}
        for i,j in ints:
            require(ordinary[i][j]<=strict[i][j]<=matrix[i][j], 'ordinary/strict/direct containment')
            require(strict[i][j]==strict[j][i], 'strict reversal isometry')
            all_edges.add(matrix[i][j]);hist['ordinary'][ordinary[i][j][1]]+=1;hist['strict'][strict[i][j][1]]+=1
            pair_hash.update(canonical([k,q,v[i],v[j],ordinary[i][j],strict[i][j]]))
        for i,j in primitive:
            x,y=v[i],v[j];a,b=sorted((x & ~y,y & ~x));hist['primitive'][matrix[i][j][1]]+=1
            if (a,b) not in trades:
                prior=[z for z in original[grades[a],encode(charges[a])] if z not in (a,b) and (z^a)&(a^b)==(z^a)]
                trades[a,b]={'masks':[a,b],'cost':matrix[i][j],'original_interval_interior':prior,'originally_primitive':not prior}
    edges=sorted(all_edges)
    require(all(a[0]<=b[0] and a[1]<=b[1] for a,b in zip(edges,edges[1:])), 'actual paired event costs form one coordinatewise chain')
    require(all(s%2==0 and c in (2*s,2*s+2) for s,c in edges), 'exact even changed-slot/native clock relation')
    # Direct subset-pair generation independently classifies normalized trades.
    # For disjoint sides, an interior is a proper grade/charge-balanced exchange.
    independently_primitive=set()
    for (k,q),v in accepted.items():
        for a,b in combinations(v,2):
            if a & b: continue
            # Same-fibre submasks at full grade cannot detect smaller exchanges;
            # instead enumerate all proper submasks of each side and meet by margin/charge.
            left=defaultdict(list);u=a
            while True:
                left[grades[u],encode(charges[u])].append(u)
                if u==0:break
                u=(u-1)&a
            proper=False;w=b
            while True:
                for u in left.get((grades[w],encode(charges[w])),[]):
                    if (u,w) not in ((0,0),(a,b)):
                        # Intermediate must obey the physical shared-port relation.
                        z=(a ^ u) | w
                        if (z&1)==z>>16:proper=True;break
                if proper or w==0:break
                w=(w-1)&b
            if not proper: independently_primitive.add((a,b))
    require(independently_primitive==set(trades), 'independent signed balanced-subtrade classifier agrees with every normalized interval primitive')
    norm=Counter(r['cost'][1] for r in trades.values());new=Counter(r['cost'][1] for r in trades.values() if not r['originally_primitive'])
    completion=Counter()
    for (a,b),row in trades.items():
        changed=a^b; t=changed>>16; singleton=(changed&65535).bit_count()-t
        require((a&1)==a>>16 and (b&1)==b>>16, 'stripping common selected records preserves new relation')
        for j in range(16-singleton):
            for z in range(2-t):completion[grades[a]+j+4*z,row['cost'][1]]+=comb(15-singleton,j)
    actual=Counter()
    for (k,q),v in accepted.items():
        for (i,j),ins in intervals(v).items():
            if len(ins)==2:actual[k,cost(v[i]^v[j])[1]]+=1
    require(completion==actual, 'entire primitive-context generating function by grade and clocks')
    critical={k:c for (k,p),c in actual.items() if p==34}
    return {'summary':stats,'histograms':{k:dict(sorted(v.items())) for k,v in hist.items()},'trades':list(trades.values()),
            'normalized_histogram':dict(sorted(norm.items())),'newly_primitive_histogram':dict(sorted(new.items())),
            'pair_digest':pair_hash.hexdigest(),'fibre_summaries':fibre_rows,'lower_witnesses':lowers,'event_cost_chain':edges,
            'critical_by_grade':critical,'primitive_context_counts':[[k,c,v] for (k,c),v in sorted(actual.items())]}

def main():
    source_pins();n=native();old,accepted,grades,charges=states(n)
    result=classify(old,accepted,grades,charges)
    # Independently regenerate charge holes by a tuple-valued generating function.
    dp={(0,ZERO):1}
    for s in n['base'][1:]:
        nxt=Counter(dp)
        for (k,q),count in dp.items():nxt[k+1,add(q,n['labels'][s])]+=count
        dp=dict(nxt)
    block_charge=total(n['labels'][s] for s in [1394]+n['triple']);combined=Counter();branch=[]
    for bit in (0,1):
        support=[set() for _ in range(20)]
        for (k,q),count in dp.items():
            grade=k+4*bit;qq=encode(add(q,block_charge) if bit else q);combined[grade,qq]+=count;support[grade].add(qq)
        for k,vals in enumerate(support):
            branch.append({'branch':bit,'grade':k,'support_count':len(vals),'holes':sorted(set(range(2048))-vals),'support_sha256':sha256(sum(1<<q for q in vals).to_bytes(256,'little')).hexdigest()})
    require(combined==Counter({key:len(v) for key,v in accepted.items()}), 'four-block polynomial recovers every fibre multiplicity')
    holes={str(k):sorted(set(range(2048))-{q for g,q in accepted if g==k}) for k in range(20)}
    for k in range(20):require(set(holes[str(k)])==set(branch[k]['holes']) & set(branch[20+k]['holes']), 'conditioned hole intersection, including out-of-range branch emptiness')
    require(encode(block_charge)==1009, 'authentic four-block charge')
    # New primitive hostile, retaining its full pre-pin interval and an old strict route.
    x,y,z=53028,77843,67124;key=(grades[x],encode(charges[x]));v=old[key];newv=accepted[key]
    interval=[w for w in v if (x^w)&(x^y)==(x^w)]
    newinterval=[w for w in newv if (x^w)&(x^y)==(x^w)]
    require(key==(8,165) and interval==[40483,53028,63747,67124,77843,90388], 'complete original grade8q165 interval')
    require(newinterval==[x,y] and all((w&1)!=(w>>16) for w in interval if w not in (x,y)), 'new actual pin deletes all four and only interval-interior records')
    require(cost(x^z)==(8,18) and cost(z^y)==(8,16) and not (x^z)&(z^y), 'old8/18 then8/16 route is strict with disjoint changed support')
    original_matrix=[[cost(a^b) for b in v] for a in v];new_matrix=[[cost(a^b) for b in newv] for a in newv]
    oldordinary=ordinary_dsu(v,original_matrix)[v.index(x)][v.index(y)];oldstrict=strict_dag(v,original_matrix)[v.index(x)][v.index(y)]
    newordinary=ordinary_dsu(newv,new_matrix)[newv.index(x)][newv.index(y)];newstrict=strict_dag(newv,new_matrix)[newv.index(x)][newv.index(y)]
    require(oldordinary==oldstrict==(8,18) and newstrict==(16,34), 'hostile exact old strict/ordinary and new strict optimum')
    newpath,_=route(newv,x,y,newordinary);_,newslotcut=route(newv,x,y,(newordinary[0]-1,999));_,newclockcut=route(newv,x,y,(99,newordinary[1]-1))
    require(newordinary==(10,22) and y not in newslotcut and y not in newclockcut, 'hostile new ordinary10/22 independently certified by route and both lower cuts')
    hostile={'grade':8,'charge':165,'endpoints':[x,y],'old_complete_fibre':v,'new_complete_fibre':newv,'old_interval':interval,'new_interval':newinterval,
             'old_ordinary':oldordinary,'old_strict':oldstrict,'new_ordinary':newordinary,'new_strict':newstrict,
             'old_literal_route':full_route(n,[x,z,y],False),'new_literal_route':full_route(n,[x,y],True),
             'new_ordinary_literal_route':full_route(n,newpath,True),'new_ordinary_below_slot_cut':newslotcut,'new_ordinary_below_clock_cut':newclockcut}
    # Every original source is separately tied to the triple: exact finite census.
    link_census=[]
    for j,s in enumerate(n['base']):
        hist={k:{'charges':0,'ordinary_clocks':0,'strict_clocks':0} for k in range(20)}
        for (k,q),v0 in old.items():
            v=[x for x in v0 if ((x>>j)&1)==(x>>16)]
            if not v:continue
            matrix=[[cost(x^y) for y in v] for x in v];ordinary=ordinary_dsu(v,matrix);ins=intervals(v)
            r=hist[k];r['charges']+=1;r['ordinary_clocks']=max(r['ordinary_clocks'],max(map(max,ordinary))[1])
            r['strict_clocks']=max(r['strict_clocks'],max((cost(v[i]^v[l])[1] for (i,l),iv in ins.items() if len(iv)==2),default=0))
        link_census.append({'source':s,'grades':hist})
    # Compare only AFTER independent derivation. Candidate code is never imported.
    candidate=json.loads((P/'sources/candidate-certificate.json.txt').read_text())
    normalized=lambda value:json.loads(json.dumps(value,sort_keys=True))
    comparisons={'summary':result['summary'],'pair_digest':result['pair_digest'],'ordinary_histogram':result['histograms']['ordinary'],
                 'strict_histogram':result['histograms']['strict'],'primitive_endpoint_histogram':result['histograms']['primitive'],
                 'normalized_primitive_histogram':result['normalized_histogram'],'newly_primitive_histogram':result['newly_primitive_histogram'],
                 'critical_by_grade':result['critical_by_grade'],'exact_holes':holes,'branch_conditioned_holes':branch,'all_link_census':link_census,
                 'fibre_summaries':result['fibre_summaries']}
    for key,value in comparisons.items():require(normalized(value)==candidate[key], 'independent agreement '+key)
    trade_keys=('masks','cost','original_interval_interior','originally_primitive')
    candidate_trades=[{key:r[key] for key in trade_keys} for r in candidate['primitive_catalogue']]
    require(sorted(normalized(result['trades']),key=lambda r:r['masks'])==sorted(candidate_trades,key=lambda r:r['masks']), 'every normalized trade, price and deleted original interval exactly matches')
    require(candidate['new_primitive_witness']['old_strict_route']==[x,z,y] and candidate['new_primitive_witness']['old_strict_optimum']==list(oldstrict), 'candidate hostile annotation agrees with independent full-route audit')
    require((sum(map(len,accepted.values())),len(accepted),sum(len(v)*(len(v)-1)//2 for v in accepted.values()),max(map(len,accepted.values())))==(65536,21696,105112,11), 'complete exact domain census')
    require(len(result['trades'])==762 and sum(result['newly_primitive_histogram'].values())==34, '762 normalized primitives,34 new')
    require(result['normalized_histogram'][34]==100 and result['critical_by_grade']=={8:100,9:300,10:300,11:100}, 'exact critical polynomial100 z8(1+z)^3')
    result.update({'status':'PASS finite-exact original nineteen-slot receiver, additional actual two-port relation, frozen exterior',
        'state_count':65536,'fibre_count':21696,'endpoint_pairs':105112,'largest_fibre':11,'relation_truth':n['truth'],
        'four_block_charge':encode(block_charge),'four_block_clocks':sum(n['prices'][s] for s in [1394]+n['triple']),
        'branch_conditioned_holes':branch,'exact_holes':holes,'all_link_census':link_census,'new_primitive_hostile':hostile,'checks':COUNT})
    return normalized(result)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=main();target=P/'AUDIT_RESULT.json'
    if args.write:
        target.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
        with (P/'FIBRE_SUMMARY.tsv').open('w') as f:
            f.write('grade\tcharge\tstates\tordinary_clocks\tstrict_clocks\tprimitive_pairs\n')
            for row in result['fibre_summaries']:f.write('\t'.join(map(str,row))+'\n')
    else:require(result==json.loads(target.read_text()), 'complete independently frozen audit result')
    print(json.dumps({'status':'PASS','states':result['state_count'],'fibres':result['fibre_count'],'pairs':result['endpoint_pairs'],
        'normalized_primitives':len(result['trades']),'new_primitives':sum(result['newly_primitive_histogram'].values()),
        'pair_digest':result['pair_digest'],'checks':result['checks']},sort_keys=True))
