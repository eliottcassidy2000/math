#!/usr/bin/env python3
"""Independent referee controls. No producer imports; explicit exceptions under -O."""
from collections import Counter
from itertools import combinations, product
from math import comb
from pathlib import Path
import hashlib
import json
import subprocess
import argparse

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repo',type=Path,help='optional repository for pinned path-to-Git-object verification')
parser.add_argument('--write-constructor',action='store_true')
args=parser.parse_args()

BASE = Path(__file__).resolve().parent
SOURCE = BASE / 'input'
counts = Counter()

def verify(condition, label, witness=None):
    counts[label] += 1
    if not condition:
        raise RuntimeError((label, witness))

def digest(data):
    return hashlib.sha256(data).hexdigest()

def v2(x):
    if x == 0:
        raise ValueError('valuation of zero excluded')
    return (abs(x) & -abs(x)).bit_length() - 1

def signed_edge(n, b):
    if b not in (-1, 1) or n % 2 != 1:
        raise ValueError('signed odd source required')
    a = v2(3*n+b)
    return (3*n+b)//(2**a), a

def cycles(p):
    verify(sorted(p) == list(range(len(p))), 'permutation_type')
    used = set()
    out = []
    for start in range(len(p)):
        if start in used:
            continue
        row = []
        x = start
        while x not in used:
            used.add(x)
            row.append(x)
            x = p[x]
        verify(x == start, 'cycle_closure')
        out.append(row)
    return out

def partitions(n, minimum=1):
    if n == 0:
        yield ()
    for k in range(minimum, n+1):
        for tail in partitions(n-k, k):
            yield (k,) + tail

def canonical_permutation(lengths):
    p = []
    for length in lengths:
        start = len(p)
        p.extend(start+(j+1)%length for j in range(length))
    return p

def pair_data(p):
    edges = list(combinations(range(len(p)), 2))
    index = {edge: i for i, edge in enumerate(edges)}
    action = [index[tuple(sorted((p[u], p[v])))] for u, v in edges]
    flips = [int(p[u] > p[v]) for u, v in edges]
    return edges, action, flips

def defect_bits(bits, action, flips):
    return sum(bits[i] ^ bits[j] ^ twist
               for i, (j, twist) in enumerate(zip(action, flips)))

def minimum_orientation(p):
    edges, action, flips = pair_data(p)
    bits = [None]*len(edges)
    orbits = cycles(action)
    for orbit in orbits:
        bits[orbit[0]] = 0
        for i, j in zip(orbit, orbit[1:]):
            bits[j] = bits[i] ^ flips[i]
    negative = sum(sum(flips[i] for i in row) % 2 for row in orbits)
    verify(defect_bits(bits, action, flips) == negative, 'construct_minimum')
    return edges, bits, orbits, negative

def histogram_formula(action, flips):
    poly = [1]
    for row in cycles(action):
        parity = sum(flips[i] for i in row) % 2
        factor = [2*comb(len(row), k) if k%2 == parity else 0
                  for k in range(len(row)+1)]
        out = [0]*(len(poly)+len(factor)-1)
        for i, a in enumerate(poly):
            for j, b in enumerate(factor):
                out[i+j] += a*b
        poly = out
    return {i: value for i, value in enumerate(poly) if value}

def tournament_reader(edges, bits):
    table = {edge: bits[i] for i, edge in enumerate(edges)}
    def edge(u, v):
        return table[min(u,v), max(u,v)] ^ int(u > v)
    return edge

def direct_defect(p, edge):
    return sum(edge(u,v) != edge(p[u],p[v])
               for u,v in combinations(range(len(p)),2))

def antipodes(p):
    cs = cycles(p)
    verify(all(len(row)%2 == 0 for row in cs), 'even_cycle_domain')
    blocks = [tuple(sorted((row[i], row[i+len(row)//2])))
              for row in cs for i in range(len(row)//2)]
    where = {v:i for i, block in enumerate(blocks) for v in block}
    pi = [where[p[block[0]]] for block in blocks]
    verify(all(where[p[u]] == where[p[v]] for u,v in blocks), 'invariant_partition')
    return cs, blocks, where, pi

def substitute(blocks, where, quotient, inside):
    def edge(u,v):
        i,j = where[u],where[v]
        if i == j:
            return inside[i] ^ int(u > v)
        return quotient(i,j)
    return edge

# Offline default: authenticate archived bytes and their Git blob object hashes.
# Optional --repo rechecks each pinned repository path, as the initial audit did.
manifest = json.loads((SOURCE/'MANIFEST.json').read_text())
for row in manifest['files']:
    data = (SOURCE/row['path']).read_bytes()
    verify(len(data) == row['bytes'] and digest(data) == row['sha256'], 'packet_hash')
source_manifest = json.loads((SOURCE/'SOURCE_MANIFEST.json').read_text())
pin = source_manifest['pin']
for row in source_manifest['sources']:
    data = (SOURCE/row['local']).read_bytes()
    git_blob_hash=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    verify(len(data)==row['bytes'] and digest(data)==row['sha256'] and git_blob_hash==row['git_blob'],
           'archived_source_blob_identity')
    if args.repo:
        obj = subprocess.check_output(['git','-C',str(args.repo),'rev-parse',pin+':'+row['path']]).decode().strip()
        gitdata = subprocess.check_output(['git','-C',str(args.repo),'cat-file','blob',obj])
        verify(obj == row['git_blob'] and data == gitdata, 'optional_pinned_repository_identity')
verify((SOURCE/'normal.json').read_bytes() == (SOURCE/'optimized.json').read_bytes(), 'producer_sealed_modes_equal')

# Complete signed-source boundary probes, not just CRT-selected positive sources.
for n in range(-5001,5002,2):
    for b in (-1,1):
        q,a = signed_edge(n,b)
        verify(q%2 == 1 and q%3 != 0 and (q>0)==(n>0), 'signed_boundary')
        verify(signed_edge(-n,-b) == (-q,a), 'negation_forward_conjugacy')
        for j in (0,1,2,3,7,16):
            carry = b*(2**j-(-1)**j)//3
            nn,bb = 2**j*n+carry, (-1)**j*b
            verify(signed_edge(nn,bb) == (q,a+j), 'ordered_integer_carry')

# Guard iff for every height in a complete binary period, plus negative copies.
negative_carries = 0
for m in range(3,36,2):
    for r,b,a in product(range(m),(-1,1),range(1,7)):
        denominator = 2**a
        s = ((3*r+b)*pow(denominator,-1,m))%m
        c = (3*r+b-denominator*s)//m
        negative_carries += int(c < 0)
        verify(3*r+b == denominator*s+m*c, 'signed_carry_integer')
        for t in range(-2**(a+1),2**(a+1)):
            n = r+m*t
            numerator = 3*n+b
            guard = (3*t+c-(1-s)*denominator)%(2*denominator)==0
            exact = n%2==1 and v2(numerator)==a
            verify(guard == exact, 'guard_iff_all_heights', (m,r,b,a,t))
            if exact:
                y = numerator//denominator
                verify(y == s+m*((3*t+c)//denominator) and y%2==1, 'signed_height_identity')
                bad = n+m*denominator
                verify(v2(3*bad+b)>a and bad%m==r and bad%denominator==n%denominator,
                       'missing_top_bit_hostile')
verify(negative_carries > 0, 'negative_carry_coverage')
verify(signed_edge(1,1)==(1,2) and signed_edge(211,1)==(317,1), 'one_211_hostile')
verify(signed_edge(3,1)[0]==signed_edge(7,-1)[0]==5 and signed_edge(5,1)[0]!=signed_edge(5,-1)[0], 'false_sheet_continuation')

# Full finite signed-port actions, checked against an independent doubling carrier.
ports = [(n,b) for n in range(1,210,2) for b in (-1,1)]
port_index = {port:i for i,port in enumerate(ports)}
p = [port_index[((2*n+b)%210,-b)] for n,b in ports]
cs = cycles(p)
verify(Counter(map(len,cs)) == {6:7,12:14}, 'full_cycle_type')
zs = [(3*n+b)%315 for n,b in ports]
verify(sorted(zs) == [z for z in range(315) if z%3], 'phi_bijection')
for i,(n,b) in enumerate(ports):
    verify(zs[p[i]] == 2*zs[i]%315, 'phi_intertwines')
    # Inverse Phi: recover b, n mod105, then select the odd lift.
    z = zs[i]
    bb = 1 if z%3==1 else -1
    r = ((z-bb)//3)%105
    nn = r if r%2 else r+105
    verify((nn,bb)==(n,b), 'phi_inverse')
for row in cs:
    verify(len(row)==(6 if zs[row[0]]%5==0 else 12), 'five_cycle_gate')
    n,b = ports[row[0]]
    ell = len(row)
    lifted = 2**ell*n+b*(2**ell-(-1)**ell)//3
    verify(lifted != n and (lifted-n)%210 == 0, 'residue_cycle_is_not_integer_cycle')
factor_data = {}
for modulus,odd in ((14,True),(15,False)):
    states=[(r,b) for r in range(modulus) if not odd or r%2 for b in (-1,1)]
    index={x:i for i,x in enumerate(states)}
    action=[index[((2*r+b)%modulus,-b)] for r,b in states]
    factor_data[modulus]=(states,action)
    expected={2:1,6:2} if modulus==14 else {6:1,12:2}
    verify(Counter(map(len,cycles(action)))==expected,'factor_cycle_type')
matched = mismatched = 0
for (r,b),(s,c) in product(factor_data[14][0],factor_data[15][0]):
    n = (r+14*((s-r)*pow(14,-1,15)%15))%210
    if b==c:
        matched += 1
        verify((n,b) in port_index,'common_sheet_CRT')
        dest = ((2*n+b)%210,-b)
        verify((dest[0]%14,dest[1])==((2*r+b)%14,-b) and
               (dest[0]%15,dest[1])==((2*s+c)%15,-c),'CRT_map_commutes')
    else:
        mismatched += 1
        verify(not any((n,d) in port_index and d==b and d==c for d in (-1,1)), 'counterfeit_gluing_rejected')
        verify(-b != -c, 'counterfeit_mismatch_persists')
verify((matched,mismatched)==(210,210),'fibre_product_sizes')

# Ternary native orientation: all pairs, including modular wrap boundaries.
def ternary_unit_digit(x):
    if x == 0:
        raise ValueError('distinct vertices required')
    while x%3 == 0:
        x //= 3
    return x%3

for k in range(6):
    modulus=3**(k+1)
    units=[x for x in range(modulus) if x%3]
    ii={x:i for i,x in enumerate(units)}
    pp=[ii[(2*x)%modulus] for x in units]
    verify(len(cycles(pp))==1,'ternary_single_cycle')
    for x,y in combinations(units,2):
        verify(ternary_unit_digit((2*y)%modulus-(2*x)%modulus)==3-ternary_unit_digit(y-x),
               'native_ternary_anti_action')

# Exhaust every tournament for every permutation cycle type through six vertices.
unrestricted_cases = []
for n in range(7):
    for lengths in partitions(n):
        sig=canonical_permutation(lengths)
        edges,action,flips=pair_data(sig)
        histogram=Counter()
        for mask in range(1<<len(edges)):
            bits=[(mask>>i)&1 for i in range(len(edges))]
            defect=defect_bits(bits,action,flips)
            histogram[defect]+=1
            if n<=4:
                verify(defect==direct_defect(sig,tournament_reader(edges,bits)), 'direct_small_defect')
        verify(dict(histogram)==histogram_formula(action,flips),'full_defect_enumerator',lengths)
        optimum=sum(length%2==0 for length in lengths)
        verify(min(histogram)==optimum,'unrestricted_minimum',lengths)
        ee,bb,orbits,neg=minimum_orientation(sig)
        verify(neg==optimum,'negative_orbit_classification',lengths)
        unrestricted_cases.append({'cycles':lengths,'orientations':sum(histogram.values()),'minimum':min(histogram)})

# Full module-constrained orientation universes for every even type through eight.
module_cases=[]
for n in range(0,9,2):
    for halves in partitions(n//2):
        lengths=tuple(2*x for x in halves)
        sig=canonical_permutation(lengths)
        cc,blocks,where,pi=antipodes(sig)
        qedges,qa,qf=pair_data(pi)
        histogram=Counter()
        for mask in range(1<<(len(qedges)+len(blocks))):
            qbits=[(mask>>i)&1 for i in range(len(qedges))]
            inside=[(mask>>(len(qedges)+i))&1 for i in range(len(blocks))]
            qe=tournament_reader(qedges,qbits)
            edge=substitute(blocks,where,qe,inside)
            external=direct_defect(pi,qe)
            internal=sum(edge(u,v)!=edge(sig[u],sig[v]) for u,v in blocks)
            total=direct_defect(sig,edge)
            verify(total==4*external+internal,'module_weighted_identity')
            verify(internal>=len(cc) and internal%2==len(cc)%2,'internal_holonomy_parity')
            histogram[total]+=1
        target=len(cc)+4*sum(d%4==0 for d in lengths)
        verify(min(histogram)==target,'module_minimum',lengths)
        module_cases.append({'cycles':lengths,'orientations':sum(histogram.values()),'minimum':min(histogram)})

# Explicit actual 210-vertex minimizing tournaments, not merely formulas.
edges,unrestricted_bits,pair_orbits,negative=minimum_orientation(p)
unrestricted_edge=tournament_reader(edges,unrestricted_bits)
verify(direct_defect(p,unrestricted_edge)==21,'actual_unrestricted_210')
cc,blocks,where,pi=antipodes(p)
qedges,qbits,qorbits,qnegative=minimum_orientation(pi)
qe=tournament_reader(qedges,qbits)
inside=[]
# Use the already constructed unrestricted orientation only on internal pairs;
# these form one negative orbit per original vertex cycle.
for u,v in blocks:
    inside.append(unrestricted_edge(u,v))
module_edge=substitute(blocks,where,qe,inside)
internal=sum(module_edge(u,v)!=module_edge(p[u],p[v]) for u,v in blocks)
for u,v in blocks:
    for w in range(210):
        if w not in (u,v):
            verify(module_edge(u,w)==module_edge(v,w),'actual_homogeneous_modules')
verify(Counter(map(len,cycles(pi)))=={3:7,6:14},'actual_quotient_cycle_type')
verify((direct_defect(pi,qe),internal,direct_defect(p,module_edge))==(14,21,77),'actual_restricted_210')
module_bits=[module_edge(u,v) for u,v in edges]
constructor={'vertex_order':ports,'permutation':p,'blocks':blocks,'quotient_permutation':pi,
             'orientation_convention':'bit 1 means edge[i][0] -> edge[i][1]; edges are lexicographic unordered pairs',
             'unrestricted_orientation_bits':''.join(map(str,unrestricted_bits)),
             'module_orientation_bits':''.join(map(str,module_bits)),
             'quotient_orientation_bits':''.join(map(str,qbits)),
             'unrestricted_defect':21,'module_defect':77,'quotient_defect':14,'internal_defect':21}
constructor_bytes=(json.dumps(constructor,indent=2,sort_keys=True)+'\n').encode()
# Emit certificate separately when called with this explicit flag.
if args.write_constructor:
    (BASE/'constructor-210.json').write_bytes(constructor_bytes)

print(json.dumps({'status':'PASS','source_pin':pin,'checks':dict(sorted(counts.items())),
                 'total_checks':sum(counts.values()),'negative_carries':negative_carries,
                 'unrestricted_cases':unrestricted_cases,'module_cases':module_cases,
                 'unrestricted_orientations':sum(x['orientations'] for x in unrestricted_cases),
                 'module_orientations':sum(x['orientations'] for x in module_cases),
                 'actual_210':{'cycle_type':{'6':7,'12':14},'unrestricted':21,'module':77,
                               'internal':21,'weighted_quotient':56,'blocks':105,
                               'certificate_sha256':digest(constructor_bytes)},
                 'counterfeit_glue':{'matching_sheets':matched,'mismatched_sheets':mismatched,
                                     'witness':{'mod14_port':[1,1],'mod15_port':[1,-1]}}},
                indent=2,sort_keys=True))
