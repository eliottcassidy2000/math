#!/usr/bin/env python3
"""Independent exact audit of the frozen divisor/norm-ball packet; stdlib only.

No producer functions are imported. Independent checks use doubled real coordinates,
finite-ring quotient groups, factor-pair form enumeration, and Gram matrices.
The producer is additionally rerun in normal/-O modes to check its receipts.
"""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import permutations, product
import json
from math import gcd, isqrt
from pathlib import Path
import subprocess
import sys

EXPECTED_MANIFEST = '7f00140b9b85d0b670231f8481f9cc78ec5f615975a60dde36a7e423f9072265'
checks = Counter()
def require(ok, category, detail):
    checks[category] += 1
    if not ok:
        raise RuntimeError(f'{category}: {detail}')

def doubled_chart(k, a, b, c):
    # Coordinates (2 Re(z), 2 Im(z)/sqrt(3)), independently of producer basis.
    return 2*a-k*(b+c), k*(b-c)

def norm4(z):
    u, v = z
    return u*u+3*v*v

def rotate(z):
    u,v=z
    return Fraction(-u-3*v,2), Fraction(u-v,2)

def direct_ball(k):
    # Exhaustive enlarged rectangle. u and v have matching parity;
    # v is a multiple of k. Norm4 bound guarantees |u|<=2k, |v|<=2k.
    return {(u,v) for u in range(-2*k,2*k+1)
            for v in range(-2*k,2*k+1,k)
            if (u+v)%2==0 and 0<norm4((u,v))<=4*k*k}

def phi(n):
    return sum(gcd(x,n)==1 for x in range(n))

def finite_ring_classes(f):
    # E/(f), in the Eisenstein basis. A pair is a unit iff its norm is coprime to f.
    units={(x,t) for x in range(f) for t in range(f)
           if gcd(x*x-x*t+t*t,f)==1}
    scalars=[s for s in range(f) if gcd(s,f)==1]
    roots=[(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]
    H={(s*x%f,s*t%f) for s in scalars for x,t in roots}
    require(H<=units,'finite_ring','scalar-unit image is contained in units')
    require(len(H)==3*phi(f),'finite_ring','unit index exactly three')
    require(len(units)%len(H)==0,'finite_ring','quotient order integral')
    bad=next(((x,t) for x,t in sorted(units)
              if ((x*x-t*t)%f,(2*x*t-t*t)%f) not in H),None)
    return len(units)//len(H),bad,len(units),len(H)

def factor_pair_forms(d):
    # Discriminant -4d. Enumerate b=B/2 and AC=b²+d, rather than A,B.
    result=[]
    for b in range(isqrt(d//3)+1):
        n=b*b+d
        for A in range(max(1,2*b),isqrt(n)+1):
            if n%A:
                continue
            C=n//A
            if gcd(gcd(A,2*b),C)!=1:
                continue
            result.append((A,2*b,C))
            if b and 2*b<A<C:
                result.append((A,-2*b,C))
    return sorted(result)

def self_inverse(t):
    A,B,C=t
    return B==0 or abs(B)==A or A==C

def gram(z,w):
    return z[0]*w[0]+3*z[1]*w[1]

parser=argparse.ArgumentParser()
parser.add_argument('--packet',required=True,type=Path)
parser.add_argument('--repository',type=Path)
args=parser.parse_args()
packet=args.packet.resolve()
manifest_bytes=(packet/'PUBLICATION_MANIFEST.json').read_bytes()
require(sha256(manifest_bytes).hexdigest()==EXPECTED_MANIFEST,'integrity','frozen manifest hash')
manifest=json.loads(manifest_bytes)
for item in manifest['files']:
    data=(packet/Path(item['source']).name).read_bytes()
    require(sha256(data).hexdigest()==item['sha256'],'integrity',item['target'])
for line in (packet/'SHA256SUMS').read_text().splitlines():
    expected,name=line.split(maxsplit=1)
    require(sha256((packet/name).read_bytes()).hexdigest()==expected,'integrity',name)

normal=subprocess.check_output([sys.executable,'-B','verify_norm_ball.py'],cwd=packet)
optimized=subprocess.check_output([sys.executable,'-B','-O','verify_norm_ball.py'],cwd=packet)
require(normal==optimized,'producer_replay','normal and optimized bytes')
require(normal==(packet/'verification.json').read_bytes(),'producer_replay','normal frozen receipt')
require(optimized==(packet/'verification.optimized.json').read_bytes(),'producer_replay','optimized frozen receipt')
receipt=json.loads(normal)
require(receipt['status']=='PASS' and receipt['live_checks']==72384,'producer_replay','reported check total')

isometry_census=[]
for k in range(1,257):
    decoded={doubled_chart(k,a,b,c):(a,b,c)
             for a in range(k+1) for b,c in product(range(2),repeat=2)
             if (a,b,c) not in ((0,0,0),(k,1,1))}
    ball=direct_ball(k)
    require(len(decoded)==4*k+2,'ball','injectivity')
    require(set(decoded)==ball,'ball','exhaustive target')
    require(sum(u for u,v in ball)==sum(v for u,v in ball)==0,'ball','centroid zero')
    core={z for z in ball if rotate(z) in ball}
    require(len(core)==6,'rotor','six admissible points')
    require(core=={z for z in ball if z[1]%k==0 and (z[0]+z[1])%(2*k)==0},'rotor','conductor intersection')
    require(all(norm4(z)==4*k*k for z in core),'rotor','outer unit shell')
    require({rotate(z) for z in core}==core,'rotor','core is rotation closed')
    require((core==ball)==(k==1),'rotor','whole-carrier closure boundary')
    for (u,v),(a,b,c) in decoded.items():
        require(doubled_chart(k,k-a,1-b,1-c)==(-u,-v),'reflection','complement')
        require(doubled_chart(k,a,c,b)==(u,-v),'reflection','prime exchange')
        require(norm4((u,v))%4==0,'ball','integral field norm')
    require(sum(v==0 for u,v in ball)==2*k,'reflection','Q fixed count')
    require(sum(u==0 for u,v in ball)==(2 if k%2==0 else 0),'reflection','CQ fixed count')
    require(not any(u==v==0 for u,v in ball),'reflection','C has no fixed point')
    if k>=2:
        require({z for z in ball if norm4(z)==4}=={(-2,0),(2,0)},'isometry','unique shortest pair')
    if k<=32:
        # A norm isometry is determined by images of 1 and k*omega.
        first=[z for z in ball if norm4(z)==4]
        second=[z for z in ball if norm4(z)==4*k*k]
        matrices=[]
        for z,w in product(first,second):
            if gram(z,w)==-2*k:
                # Source determinant in doubled real coordinates is 2k.
                require(abs(z[0]*w[1]-z[1]*w[0])==2*k,'isometry','lattice basis determinant')
                matrices.append((z,w))
        require(len(matrices)==(12 if k==1 else 4),'isometry','full Gram-isometry count')
        if k<=8:
            isometry_census.append({'k':k,'isometries':len(matrices)})

# Divisor sets built by trial division of actual integers, not exponent generation.
for p,q,r in ((2,3,5),(3,2,7),(5,3,2)):
    for k in range(1,9):
        N=p**k*q*r
        divisors=set()
        for d in range(1,isqrt(N)+1):
            if N%d==0:
                divisors.update((d,N//d))
        divisors-={1,N}
        points=set()
        for d in divisors:
            value=d
            abc=[]
            for prime in (p,q,r):
                e=0
                while value%prime==0:
                    value//=prime
                    e+=1
                abc.append(e)
            require(value==1,'integer_divisors','all prime factors decoded')
            points.add(doubled_chart(k,*abc))
        require(points==direct_ball(k),'integer_divisors','numerical-divisor image equals ball')

# All conductor elements in enlarged boxes, also for negative coefficients and zero.
for k in range(1,65):
    for x,y in product(range(-2*k,2*k+1),range(-4,5)):
        u,v=2*x-k*y,k*y
        ru,rv=rotate((u,v))
        in_order=rv.denominator==1 and rv%k==0 and (ru+rv)%2==0
        require(in_order==(x%k==0),'whole_lattice','one-step integrality gate')

ring_rows=[]
for f in range(2,193):
    h,bad,units,H=finite_ring_classes(f)
    if f%2==0:
        d=3*(f//2)**2
        fs=factor_pair_forms(d)
        require(len(fs)==h,'class_group','factor-pair form count versus finite-ring quotient')
        require(all(self_inverse(t) for t in fs)==(bad is None),'class_group','ambiguity versus direct group squares')
        require((bad is None)==(f in (2,4,8)),'class_group','even-family exponent-two classification')
        if f & (f-1):
            require(h%3==0,'class_group','every odd conductor factor yields factor three')
        if f in (2,4,6,8,16,32,64,128,192):
            ring_rows.append({'conductor':f,'unit_count':units,'scalar_root_subgroup_size':H,
                              'class_number':h,'square_obstruction':bad,'self_inverse_classes':sum(map(self_inverse,fs))})
        if f in (2,4,8,16):
            producer_row=next(row for row in receipt['even_idoneal_rows'] if row['conductor']==f)
            require(fs==[tuple(t) for t in producer_row['forms']],'class_group','complete literal producer list')

# Big integers plus direct symbolic derivation in the written audit.
for exponent in range(4,1025):
    f=2**exponent
    A,B,C=(13,8,16) if exponent==4 else (16,8,1+3*4**(exponent-3))
    require(B*B-4*A*C==-3*f*f,'dyadic','exact discriminant')
    require(gcd(gcd(A,B),C)==1,'dyadic','primitive witness')
    require(0<B<A<C,'dyadic','strictly interior reduced witness')

# Labelled reflection obstruction: every involutive vertex permutation fixes an edge coordinate.
for n in range(2,9):
    for perm in permutations(range(n)):
        if all(perm[perm[i]]==i for i in range(n)):
            require(any({perm[i],perm[j]}=={i,j} for i in range(n) for j in range(i+1,n)),
                    'labelled_obstruction','a fixed unordered pair forbids complement reversal')

native={(1,0),(2,1),(1,1),(1,2),(0,1),(-1,1),(-1,0),(-2,-1),(-1,-1),(-1,-2),(0,-1),(1,-1)}
native_norms=Counter(x*x-x*t+t*t for x,t in native)
ball_norms=Counter(norm4(z)//4 for z in direct_ball(2))
require(native_norms=={1:6,3:6},'nonidentification','native two-shell inventory')
require(ball_norms=={1:2,3:2,4:6},'nonidentification','divisor three-shell inventory')
require(len(ball_norms)>len(native_norms),'nonidentification','no nonzero common scalar norm matching')
require(doubled_chart(2,2,0,0)!=(2,0),'nonidentification','p squared does not map multiplicatively')

repository_sources_checked=0
if args.repository:
    for source in json.loads((packet/'SOURCE_PROVENANCE.json').read_text())['repository_sources']:
        data=subprocess.check_output(['git','-C',str(args.repository),'show',source['repository_commit']+':'+source['path']])
        require(sha256(data).hexdigest()==source['sha256'],'provenance',source['path'])
        repository_sources_checked+=1

# Recheck the producer freeze after all computation.
require(sha256((packet/'PUBLICATION_MANIFEST.json').read_bytes()).hexdigest()==EXPECTED_MANIFEST,'integrity','manifest unchanged')
for item in manifest['files']:
    require(sha256((packet/Path(item['source']).name).read_bytes()).hexdigest()==item['sha256'],'integrity','payload unchanged')

result={'status':'PASS','live_checks':sum(checks.values()),'checks_by_category':dict(sorted(checks.items())),
        'reviewed_producer_manifest_sha256':EXPECTED_MANIFEST,
        'producer_checks_per_mode':receipt['live_checks'],'producer_modes_byte_identical':True,
        'scope':{'balls_and_reflections_k':[1,256],'gram_isometry_k':[1,32],
                 'direct_divisors_k':[1,8],'prime_assignments':3,'whole_lattice_k':[1,64],
                 'finite_ring_quotients_f':[2,192],'even_reduced_forms_f':[2,192],
                 'dyadic_witness_exponents':[4,1024],'involutive_vertex_permutations_n':[2,8],
                 'repository_sources_checked':repository_sources_checked},
        'finite_ring_group_rows':ring_rows,'isometry_rows':isometry_census,
        'uniform_claim_basis':'Written mathematical proof, not extrapolation from these finite controls.',
        'nonclaims':['No global 65-idoneal completeness proof.','No native source, clock, or cancellation identity.',
                     'No divisor/ring multiplicativity or divisor/ideal-class bijection.']}
print(json.dumps(result,indent=2,sort_keys=True))
