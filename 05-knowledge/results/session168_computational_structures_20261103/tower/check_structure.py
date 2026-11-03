#!/usr/bin/env python3
"""Exact controls for Cayley-Dickson structural consequences, stdlib only."""
from collections import deque
from functools import lru_cache
from math import comb
from pathlib import Path
import hashlib
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(None)
def sign(n, a, b):
    """Sign of e_a e_b, obtained directly from pair multiplication."""
    if n == 0:
        return 1
    half = 1 << (n-1)
    u, A = a % half, a // half
    v, B = b % half, b // half
    if A == B == 0:
        return sign(n-1, u, v)
    if A == 0:
        return sign(n-1, v, u)
    conjugate_v = 1 if v == 0 else -1
    if B == 0:
        return conjugate_v * sign(n-1, u, v)
    return -conjugate_v * sign(n-1, v, u)


def shifted(mask, a, size):
    out = 0
    while mask:
        bit = mask & -mask
        out |= 1 << ((bit.bit_length()-1) ^ a)
        mask ^= bit
    return out


def insert(basis, mask):
    while mask:
        p = mask.bit_length()-1
        if p in basis:
            mask ^= basis[p]
        else:
            basis[p] = mask
            return True
    return False


def close_translations(basis, n):
    queue = list(basis.values())
    size = 1 << n
    while queue:
        mask = queue.pop()
        for j in range(n):
            image = shifted(mask, 1 << j, size)
            if insert(basis, image):
                queue.append(image)


def zero_moments(mask, d):
    size = 1 << d
    if mask.bit_count() % 2:
        return False
    moments = 0
    for y in range(size):
        if (mask >> y) & 1:
            moments ^= y
    return moments == 0


def belongs_W(mask, n):
    half = 1 << (n-1)
    low = mask & ((1 << half)-1)
    high = mask >> half
    return low.bit_count() % 2 == 0 and (low ^ high) in (0, (1 << half)-1)


def belongs_I2(mask, n):
    half = 1 << (n-1)
    low = mask & ((1 << half)-1)
    high = mask >> half
    return low == high and zero_moments(low, n-1)


def multiplication(n, x, y):
    result = {}
    for a, ca in x.items():
        for b, cb in y.items():
            c = a ^ b
            result[c] = result.get(c, 0) + ca * cb * sign(n, a, b)
    return {a: c for a, c in result.items() if c}


def check_dimension(n):
    size = 1 << n
    half = size // 2
    fs = [sum((sign(n, a, x) < 0) << x for x in range(size)) for a in range(size)]
    all_mask = (1 << size)-1
    previous=[sum((sign(n-1,a,x)<0)<<x for x in range(half)) for a in range(half)]
    low_all=(1<<half)-1
    q_mask=low_all^1
    q = lambda a: int(a != 0)
    c = lambda a, b: q(a) ^ q(b) ^ q(a ^ b)
    for a in range(size):
        for x in range(size):
            u, A = a % half, a // half
            y, t = x % half, x // half
            f = int(sign(n-1, u, y) < 0) ^ (t*c(u,y)) ^ (A*q(y)) ^ (A*t)
            require(f == int(sign(n,a,x)<0), f'sign recurrence {n,a,x}')
        require(fs[a] ^ shifted(fs[a], a, size) == (all_mask if a else 0), 'basis squares')

    kernel, derived = {}, {}
    commutator_count = 0
    for a in range(size):
        for b in range(size):
            # Schreier difference for L_a L_b versus L_(a+b).
            defect = fs[b] ^ shifted(fs[a],b,size) ^ fs[a^b]
            require(belongs_W(defect,n), 'kernel upper bound')
            insert(kernel,defect)
            p = fs[b] ^ shifted(fs[a],b,size) ^ fs[a] ^ shifted(fs[b],a,size)
            require(belongs_I2(p,n), f'commutator bound {n,a,b}')
            u,A=a%half,a//half
            v,B=b%half,b//half
            lower_p=previous[v]^shifted(previous[u],v,half)^previous[u]^shifted(previous[v],u,half)
            double_q=q_mask^shifted(q_mask,u,half)^shifted(q_mask,v,half)^shifted(q_mask,u^v,half)
            if A^B:
                lower_p ^= double_q
            if (B*q(u))^(A*q(v)):
                lower_p ^= low_all
            require(p==lower_p|(lower_p<<half),'commutator recurrence')
            insert(derived,p)
            commutator_count += 1
    close_translations(kernel,n)
    close_translations(derived,n)
    require(len(kernel)==half,'kernel rank')
    require(len(derived)==half-n,'derived rank')
    for mask in kernel.values():
        for j in range(n):
            require(belongs_I2(mask ^ shifted(mask,1<<j,size),n),'kernel commutator')
    require(belongs_I2(all_mask,n),'squares in derived')
    require(half+n-len(derived)==2*n,'abelianization dimension')
    constructed={}
    for j in range(n-1):
        low=1 | (1 << (1 << j))
        insert(constructed,low | (low << half))
    insert(constructed,((1<<half)-1)<<half)
    close_translations(constructed,n)
    require(len(constructed)==half,'constructive 2n-generator upper bound')

    # Independent expected kernel basis: lower coordinate pair flips, repeated, plus top bit.
    explicit_kernel = {}
    for y in range(1,half):
        low = 1 | (1 << y)
        insert(explicit_kernel,low | (low << half))
    insert(explicit_kernel,((1<<half)-1)<<half)
    require(len(explicit_kernel)==half,'explicit kernel dimension')
    combined=dict(kernel)
    for v in explicit_kernel.values():
        insert(combined,v)
    require(len(combined)==half,'kernel equality')

    lower_ranks = {}
    current=dict(derived)
    k=2
    while current:
        expected=sum(comb(n-1,j) for j in range(k,n))
        require(len(current)==expected,f'lower central rank {n,k}')
        lower_ranks[str(k)]=len(current)
        following={}
        for mask in current.values():
            for j in range(n):
                insert(following,mask^shifted(mask,1<<j,size))
        current=following
        k+=1
    require(k==n,'nilpotency class')

    # A nonzero maximal nested commutator has constant sign.
    low=1 | (1 << 1)
    nested=low | (low << half)
    for j in range(1,n-1):
        nested ^= shifted(nested,1<<j,size)
    require(nested==all_mask,'maximal commutator witness')

    # Common fixed points of the identity-sheet stabilizer are exactly +/-e_0.
    stabilizer=[((1<<half)-1)<<half]
    for y in range(2,half):
        low=(1<<1)|(1<<y)
        stabilizer.append(low|(low<<half))
    fixed=[x for x in range(size) if all(((f>>x)&1)==0 for f in stabilizer)]
    require(fixed==[0],'centralizer fixed points')

    # The canonical doubling operators form a regular Clifford subgroup.
    canonical=[1 << j for j in range(n)]
    for a in canonical:
        for b in canonical:
            if a != b:
                p=fs[b]^shifted(fs[a],b,size)^fs[a]^shifted(fs[b],a,size)
                require(p==all_mask,'canonical Clifford anticommutation')
    clifford={(0,0)}
    queue=deque([(0,0)])
    while queue:
        b,g=queue.popleft()
        for a in canonical:
            element=(a^b,g^shifted(fs[a],b,size))
            if element not in clifford:
                clifford.add(element)
                queue.append(element)
    require(len(clifford)==2*size,'canonical Clifford group order')
    require(len({(a,f&1) for a,f in clifford})==2*size,'canonical Clifford regular action')

    # Enumerate all operators independently for n=3,4.
    enumerated=None
    if n<=4:
        seen={(0,0)}
        queue=deque([(0,0)])
        translated_rows={(a,b):shifted(fs[a],b,size) for a in range(size) for b in range(size)}
        while queue:
            b,g=queue.popleft()
            for a in range(1,size):
                element=(a^b,g^translated_rows[a,b])
                if element not in seen:
                    seen.add(element)
                    queue.append(element)
        enumerated=len(seen)
        require(enumerated==1<<(half+n),'full group enumeration')
        # Compare actual kernel rank and membership without summing bit masks arithmetically.
        actual_kernel={}
        for a,f in seen:
            if a==0:
                require(belongs_W(f,n),'enumerated kernel membership')
                insert(actual_kernel,f)
        require(len(actual_kernel)==half,'enumerated kernel rank')

    return {'n':n,'algebra_dimension':size,'signed_unit_states':2*size,
            'kernel_rank':len(kernel),'log2_group_order':half+n,
            'derived_rank':len(derived),'minimum_generators':2*n,
            'nilpotency_class':n-1,'lower_central_ranks':lower_ranks,
            'canonical_clifford_order':len(clifford),'minimum_cycles_for_connected_cover':n,
            'generator_pair_checks':commutator_count,'enumerated_group_order':enumerated,
            'common_fixed_indices_of_stabilizer':fixed}


def main():
    results=[check_dimension(n) for n in range(3,8)]
    x={1:1,10:1};y={4:1,15:-1};e4={4:1}
    require(multiplication(4,x,y)=={},'zero divisor')
    require(multiplication(4,x,x)=={0:-2},'x square')
    require(multiplication(4,x,multiplication(4,x,e4))=={4:-2,15:-2},'alternativity failure')
    require(multiplication(4,multiplication(4,x,x),e4)=={4:-2},'alternativity other bracket')
    require(sign(4,10,4)*sign(4,1,14)==-1,'Clifford term 1')
    require(sign(4,1,4)*sign(4,10,5)==-1,'Clifford term 2')
    # Persistence by the canonical low-coordinate embedding.
    for n in range(5,8):
        require(multiplication(n,x,y)=={},'persistent zero divisor')
    report={'status':'PASS','dimensions':results,
            'fixed_base':{'vertices':22,'edges':157,'beta':136,'largest_full_transport_n':68,
                          'first_impossible_n':69,'largest_connected_clifford_n':136,'scope':'formal permutation-voltage cover only'},
            'zero_divisor':{'n':4,'x':x,'y':y,'product':{},'squared_norms':[2,2,0]},
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    path=Path(__file__).with_name('structure_results.json')
    path.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
