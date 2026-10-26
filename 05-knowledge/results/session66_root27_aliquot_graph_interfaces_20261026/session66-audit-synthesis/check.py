#!/usr/bin/env python3
"""Independent, small cross-carrier controls; imports no producer code.

Run from this directory with --repo ../session66-repo for Git-tree checks.
Without --repo, archived local snapshot SHA-256/Git-blob IDs are still checked;
the root's repository-only source references require the optional checkout.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
from math import isqrt
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ap = argparse.ArgumentParser()
ap.add_argument('--repo', type=Path)
args = ap.parse_args()
checks = Counter()


def require(ok, tag):
    checks[tag] += 1
    if not ok:
        raise RuntimeError(tag)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def blob(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()


def v2(n):
    n = Fraction(n)
    require(n != 0, 'valuation_nonzero')
    a, b = abs(n.numerator), n.denominator
    return (a & -a).bit_length() - (b & -b).bit_length()


def constants(word):
    L = len(word)
    C = sum(int(c) * 2**i * 3**(L-1-i) for i, c in enumerate(word))
    return L, 3**L, 2**L, C, 3**L-2**L


def native(z, q):
    c = z % 2
    if q < 1 or (q == 1 and c == 0):
        raise ValueError('blocked native macro')
    return (3*z+c)//2, q+3*c-1, c


z, q = 7, 1
states = [(z, q)]
for _ in range(110):
    z, q, _ = native(z, q)
    states.append((z, q))
require(states[:5] == [(7,1),(11,3),(17,5),(26,7),(39,6)], 'root_sections')
require(states[86] == (10952205219944621,47), 'actual_macro86')
P, Q = constants('100'), constants('101000')
lp, ap_, bp, cp, dp = P
lq, aq, bq, cq, dq = Q
delta = cq*dp-cp*dq
z, q = states[86]
y = (ap_*z+cp)//bp
require((v2(dp*z+cp), v2(dp*y+cp), v2(dq*y+cq)) == (5,2,9), 'actual_refill_5_2_9')
require(delta == 684 and v2(delta) == 2, 'signed_defect')
require(bp*dp*(dq*y+cq) == ap_*dq*(dp*z+cp)+bp*delta, 'cocycle')
rP, rQ = -Fraction(cp,dp), -Fraction(cq,dq)
require(rP != rQ, 'finite_bank_hostile_two_centers')
require(Fraction(2,3)**lp < 1 and Fraction(2,3)**lq < 1, 'real_contraction')
try:
    native(4, 1)
except ValueError:
    checks['blocked_symbol_rejected'] += 1
else:
    raise RuntimeError('blocked_symbol_executed')
require(native(4,2)[:2] == (6,1), 'q2_not_q1_lifetime')
require(v2(4-rP) == v2(4-rQ) == 0, 'q2_center_bank_hostile')


def cycles(perm):
    todo, result = set(perm), []
    while todo:
        start = min(todo)
        current, cycle = start, []
        while current in todo:
            todo.remove(current)
            cycle.append(current)
            current = perm[current]
        require(current == start, 'permutation_cycles')
        result.append(cycle)
    return result


ports = [(n,b) for n in range(1,210,2) for b in (-1,1)]
lift = {(n,b):((2*n+b)%210,-b) for n,b in ports}
ctype = Counter(map(len,cycles(lift)))
require(ctype == {6:7,12:14}, 'signed_port_cycle_census')
phi = {(n,b):(3*n+b)%315 for n,b in ports}
require(len(set(phi.values())) == 210, 'port_phi_bijection')
require(all(phi[lift[p]] == 2*phi[p]%315 for p in ports), 'port_doubling_conjugacy')


def negative_pair_orbits(perm):
    todo = set(itertools.combinations(sorted(perm),2))
    count = 0
    while todo:
        a,b = min(todo)
        start = (a,b)
        while True:
            todo.remove(tuple(sorted((a,b))))
            a,b = perm[a],perm[b]
            if tuple(sorted((a,b))) == start:
                count += (a,b) == start[::-1]
                break
    return count


require(negative_pair_orbits(lift) == 21, 'signed_port_minimum_defect')
port_cycles = cycles(lift)
blocks = [frozenset((cyc[i],cyc[i+len(cyc)//2])) for cyc in port_cycles for i in range(len(cyc)//2)]
idx = {b:i for i,b in enumerate(blocks)}
quotient = {i:idx[frozenset(lift[x] for x in b)] for i,b in enumerate(blocks)}
qtype = Counter(map(len,cycles(quotient)))
require(qtype == {3:7,6:14}, 'clone_quotient_cycles')
qdef = negative_pair_orbits(quotient)
require(qdef == 14 and 21+4*qdef == 77, 'clone_weighted_defect_77')


def order2(m):
    r, d = 2 % m, 1
    while r != 1:
        r = r*2 % m
        d += 1
    return d


def prime_factors(m):
    factors=[]
    p=3
    while p*p <= m:
        if m % p == 0:
            factors.append(p)
            while m % p == 0:
                m //= p
        p += 2
    if m > 1:
        factors.append(m)
    return factors


for m in range(3,1000,2):
    d=order2(m)
    actual=(d % 2 == 0 and pow(2,d//2,m)==m-1)
    two_parts={v2(order2(p)) for p in prime_factors(m)}
    predicted=(len(two_parts)==1 and next(iter(two_parts)) > 0)
    require(actual==predicted,'odd_table_shift_classification')
require(pow(2,1,3)==2 and pow(2,2,5)==4 and all(pow(2,t,15)!=14 for t in range(4)), 'mod3_mod5_not_mod15')


def s(n):
    return sum(d for d in range(1,n) if n%d == 0)


divisors = [d for d in range(1,211) if 210%d == 0]
interior = divisors[1:-1]
iota = {d:210//d for d in interior}
require((len(divisors),len(divisors[:-1]),len(interior)) == (16,15,14), 'divisor_carrier_sizes')
require(210//1 not in divisors[:-1], 'proper_divisor_boundary_not_stable')
require(negative_pair_orbits(iota) == 7, 'divisor_minimum_defect_7')
require(sum(iota[a]>iota[b] for a,b in itertools.combinations(interior,2)) == 91, 'divisor_actual_increasing_defect_91')
pairs = [(d,210//d) for d in interior if d<210//d]
require(min(pairs,key=lambda p:p[1]-p[0]) == (14,15), 'unique_closest_pair')
require(sum((b-a)==1 for a,b in pairs) == 1, 'closest_unique')
require(sum(divisors[:-1]) == 366, 'divisor_weighted_evaluation')
require(s(4)==3 and s(9)==4, 'prime_relabelling_loses_weights')
require(s(2)==s(3)==1 and s(4)!=s(6), 'fiber_not_monoid_congruence')
require(4 == (s(4)-1)**2, 'inclusive_fiber_boundary')


def independent_evaluation(n, edges):
    result=0
    for mask in range(1<<n):
        if all(not(mask>>a&1 and mask>>b&1) for a,b in edges):
            result += 2**mask.bit_count()
    return result


require(independent_evaluation(3,[(0,1),(0,2),(1,2)])==7, 'arbitrary_graph_K3_value7')
require(independent_evaluation(4,[(0,1),(1,2),(2,3)])==21, 'arbitrary_graph_P4_value21')
triangle={(0,1),(1,2),(2,0)}
sigma={0:1,1:0,2:2}
triangle_image={(sigma[u],sigma[v]) for u,v in triangle}
def port_edges(arcs):
    return {frozenset((('t',u),('h',v))) for u,v in arcs}
tdef=len(triangle ^ triangle_image)//2
gdiff=len(port_edges(triangle) ^ port_edges(triangle_image))
require((tdef,gdiff)==(3,6) and gdiff==2*tdef,'graph_difference_twice_tournament_defect')

# Exact small completion of the cited large Goldbach computation, no FFT.
bound=19562
prime=bytearray(b'\1')*(bound+1)
prime[:2]=b'\0\0'
for p in range(2,isqrt(bound)+1):
    if prime[p]:
        prime[p*p:bound+1:p]=b'\0'*len(range(p*p,bound+1,p))
primes=[p for p in range(2,bound+1) if prime[p]]
witnesses=[]
for e in range(8,bound+1,2):
    pair=next(((p,e-p) for p in primes if 2*p<e and prime[e-p]),None)
    require(pair is not None,'distinct_goldbach_small_completion')
    witnesses.append((e,*pair))

# Verify source hashes, recorded Git blob IDs, and optionally committed bytes.
source_records=[]
specs=[('session66-root27','SOURCE_PINS.json','sources','snapshot','repo_path'),
       ('session66-modular','SOURCE_MANIFEST.json','sources','local','path'),
       ('session66-graph-pairs','SOURCE_PINS.json','files','local_path','source_path')]
for package,name,key,localkey,pathkey in specs:
    manifest=json.loads((BASE/package/name).read_text())
    for entry in manifest[key]:
        data=(BASE/package/entry[localkey]).read_bytes()
        require(sha(data)==entry['sha256'],'archived_source_sha256')
        require(blob(data)==entry['git_blob'],'archived_source_git_blob')
        commit=entry.get('commit',entry.get('source_commit',manifest.get('pin')))
        if args.repo:
            committed=subprocess.check_output(['git','-C',str(args.repo),'show',f'{commit}:{entry[pathkey]}'])
            require(committed==data,'committed_source_bytes')
        source_records.append((package,entry[pathkey],entry['sha256']))
if args.repo:
    manifest=json.loads((BASE/'session66-root/SOURCES.json').read_text())
    for entry in manifest['repository']:
        data=subprocess.check_output(['git','-C',str(args.repo),'show',entry['commit']+':'+entry['path']])
        require(sha(data)==entry['sha256'],'root_repository_source_sha256')
        require(blob(data)==entry['git_blob'],'root_repository_source_git_blob')
        source_records.append(('session66-root',entry['path'],entry['sha256']))

output={
 'status':'PASS',
 'scope':'Small independently implemented cross-carrier and boundary controls; universal claims require the audited proofs.',
 'checks':dict(sorted(checks.items())),
 'signed_port_cycle_type':dict(sorted(ctype.items())),
 'clone_quotient_cycle_type':dict(sorted(qtype.items())),
 'different_defect_observables':{'signed_ports_unrestricted':21,'signed_ports_antipodal_modules':77,'divisor_interior_unrestricted':7,'divisor_increasing_orientation':91},
 'goldbach_small_completion':{'even_range':[8,bound],'count':len(witnesses),'witness_sha256':sha(json.dumps(witnesses,separators=(',',':')).encode())},
 'git_commit_bytes_checked':bool(args.repo),
 'source_records':source_records,
}
print(json.dumps(output,indent=2,sort_keys=True))
