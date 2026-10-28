"""Literal-permutation audit, independent of the C++ refinement canonicalizer."""
from collections import Counter, defaultdict
from functools import lru_cache, reduce
from itertools import permutations, combinations
from math import gcd
from pathlib import Path
import gzip, json, importlib.util

HERE = Path(__file__).resolve().parent
CHECKS = 0

def need(b, message):
    global CHECKS
    CHECKS += 1
    if not b:
        raise RuntimeError(message)

def rows(n, x):
    r = [0]*n
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            forward = j == i+1 or (x >> k & 1)
            if j != i+1:
                k += 1
            a,b = (i,j) if forward else (j,i)
            r[a] |= 1 << b
    return tuple(r)

def relabel(r, p):
    return tuple(sum(((r[v] >> w) & 1) << j for j,w in enumerate(p)) for v in p)

def encode(r):
    return sum(((r[i] >> j) & 1) << k for k,(i,j) in
               enumerate(combinations(range(len(r)), 2)))

def freecode(r):
    return sum(((r[i] >> j) & 1) << k for k,(i,j) in enumerate(
        (i,j) for i in range(len(r)) for j in range(i+2,len(r))))

def delete(r, j):
    return relabel(r, tuple(i for i in range(len(r)) if i != j))

def converse(r):
    full = (1 << len(r))-1
    return tuple(v ^ full ^ (1 << i) for i,v in enumerate(r))

@lru_cache(None)
def canon(r):
    return min(encode(relabel(r,p)) for p in permutations(range(len(r))))

def pointed(r, root):
    return min((encode(relabel(r,p)),p.index(root)) for p in permutations(range(len(r))))

def intern(values):
    ids = {}
    return [ids.setdefault(x,len(ids)) for x in values]

def same_partition(a,b):
    forward,backward = {},{}
    for x,y in zip(a,b):
        need(forward.setdefault(x,y)==y and backward.setdefault(y,x)==x, 'partition equivalence')

def probability(counter):
    d = reduce(gcd, counter.values())
    return tuple((k,v//d) for k,v in sorted(counter.items()))

def literal_carriers():
    result = []
    for n in range(2,6):
        saved = json.loads(gzip.decompress((HERE/f'closed_flags{n}.json.gz').read_bytes()))
        count = n << ((n-1)*(n-2)//2)
        q, observed, transitions = [], [], []
        for s in range(count):
            x,j = divmod(s,n)
            r = rows(n,x)
            q.append(pointed(r,j))
            observed.append((canon(r),canon(delete(r,j))))
            reset = Counter()
            for p in permutations(range(n)):
                if all(r[p[k]] >> p[k+1] & 1 for k in range(n-1)):
                    reset[freecode(relabel(r,p))*n + p.index(j)] += 1
            need(bool(reset), 'Hamilton path exists')
            transitions.append(reset)
        same_partition(q, saved['pointed_id_by_state'])
        groups = defaultdict(list)
        for s,t in enumerate(q):
            groups[t].append(s)
        for s,reset in enumerate(transitions):
            need(set(reset)==set(groups[q[s]]), 'literal reset fibre equals pointed type')
            need(len(set(reset.values()))==1, 'rooted automorphism multiplicity is uniform')
        colors = intern(observed)
        full = count//n-1
        tau = [((s//n)^full)*n+s%n for s in range(count)]
        trajectory = [len(set(colors))]
        while True:
            profiles = []
            for transition in transitions:
                counter = Counter()
                for t,multiplicity in transition.items():
                    counter[colors[t]] += multiplicity
                profiles.append(probability(counter))
            new = intern((colors[s],colors[tau[s]],profiles[s]) for s in range(count))
            if new==colors:
                break
            colors=new
            trajectory.append(len(set(colors)))
        same_partition(colors, saved['closed_color_by_state'])
        result.append({'order':n,'literal_states':count,'refinement_sizes':trajectory})
    return result

def forest_identity():
    cases = 0
    for n in range(2,7):
        total = 1 << ((n-1)*(n-2)//2)
        for x in range(total):
            r,partner = rows(n,x),rows(n,x^(total-1))
            for j in range(n):
                keep = [i for i in range(n) if i!=j]
                actual = delete(partner,j)
                predicted = list(converse(delete(r,j)))
                for a in range(n-2):
                    if keep[a+1]==keep[a]+1:
                        predicted[a] ^= 1 << (a+1)
                        predicted[a+1] ^= 1 << a
                need(tuple(predicted)==actual, 'delete native complement equals ordered forest complement')
                for k in keep:
                    unchanged = abs(k-j)==1
                    need(((r[j]>>k)&1)^((partner[j]>>k)&1)==(not unchanged), 'root mask neighbour exception')
                cases += 1
    # The skipped predecessor/successor pair is NOT a restored path edge.
    r=rows(3,1)
    need(delete(rows(3,0),1)!=delete(r,1),'minimal interior deletion flips the gap edge')
    return {'cases':cases,'minimal_interior_witness':{'order':3,'freecode':1,'root':1,
             'parent_rows':delete(r,1),'partner_parent_rows':delete(rows(3,0),1)}}

def native_self_converse_boundary():
    result={}
    for n in range(3,7):
        if n%2 in result:
            continue
        total=1 << ((n-1)*(n-2)//2)
        for x in range(total):
            r,partner=rows(n,x),rows(n,x^(total-1))
            if canon(r)==canon(converse(r)) and canon(partner)!=canon(converse(partner)):
                result[n%2]={'order':n,'path':list(range(n)),'freecode':x,'rows':r,'partner_rows':partner,
                    'ordinary_codes':[canon(r),canon(partner),canon(converse(partner))]}
                break
    need(len(result)==2,'odd and even native closure hostiles')
    return result

def score_roots():
    cases=0
    for n in range(1,4):
        edges=list(combinations(range(n),2));e=len(edges)
        for odd in (False,True):
            for x in range(1 << (2*e+(n if odd else 0))):
                N=2*n+odd;r=[0]*N;twice_b=[0]*n
                def arc(a,b):
                    r[a] |= 1 << b
                for k,(i,j) in enumerate(edges):
                    a=(x>>k)&1;g=(x>>(e+k))&1
                    arc(i,j) if a else arc(j,i)
                    arc(n+j,n+i) if a else arc(n+i,n+j)
                    arc(i,n+j) if g else arc(n+j,i)
                    arc(j,n+i) if g else arc(n+i,j)
                    twice_b[i] += (2*a-1)+(2*g-1)
                    twice_b[j] += -(2*a-1)+(2*g-1)
                for i in range(n):
                    s=(x>>(2*e+i))&1 if odd else 1
                    arc(i,n+i) if s else arc(n+i,i)
                    if odd:
                        arc(2*n,i);arc(n+i,2*n)
                        twice_b[i] += 2*s-1
                    twice_score=2*r[i].bit_count()-(N-1)
                    need(twice_score==twice_b[i]+(-1 if odd else 1),'B/D root score with gauge shift')
                    need(r[i].bit_count()+r[n+i].bit_count()==N-1,'anti-paired score reflection')
                cases+=1
    return {'labelled_inputs':cases,'odd_shift':'score_plus = B_score - 1/2',
            'even_shift':'score_plus = D_score + 1/2'}

def hostile():
    source=HERE.parent/'relation_tower_frobenius_ports_20261028/check_root_switch.py'
    spec=importlib.util.spec_from_file_location('inherited_root_switch',source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    result=module.hostile()
    need(result['automorphism_orders']==[81,27],'inherited independent group search')
    need(result['labelled_cyclic_triples']==67,'full labelled triple control')
    result['source']='../relation_tower_frobenius_ports_20261028/check_root_switch.py'
    return result

def witness():
    n=4
    saved=json.loads(gzip.decompress((HERE/'closed_flags4.json.gz').read_bytes()))
    meta=saved['pointed_types'];q=saved['pointed_id_by_state']
    def output(s):
        return meta[q[s]][2],meta[q[s]][0]
    a,b=0,4
    tau=lambda s: ((s//n)^7)*n+s%n
    need(output(a)==output(b),'same initial deletion observations')
    need(output(tau(a))!=output(tau(b)),'one complement distinguishes')
    return {'states':[a,b],'rows':[rows(n,a//n),rows(n,b//n)],'roots':[a%n,b%n],
            'initial_output':output(a),'after_complement':[output(tau(a)),output(tau(b))]}

def main():
    report={'status':'PASS','audit':'literal permutations and literal Hamilton resets',
            'carriers':literal_carriers(),'deletion_forest':forest_identity(),
            'self_converse_native_nonclosure':native_self_converse_boundary(),
            'child_parent_insufficiency':witness(),'root_system_scores':score_roots(),
            'hostile_81_27':hostile()}
    report['checks']=CHECKS
    (HERE/'structure_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report))

if __name__=='__main__':
    main()
