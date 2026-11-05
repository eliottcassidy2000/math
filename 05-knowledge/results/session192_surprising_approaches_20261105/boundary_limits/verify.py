#!/usr/bin/env python3
"""Exact controls for fixed-feedback boundary laws. Standard library only."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations, product
from math import factorial
from pathlib import Path
import hashlib, json, random


def require(ok, message):
    if not ok:
        raise ValueError(message)


def graph(core, word):
    r = len(core)
    n = r + len(word)
    rows = list(core) + [0] * len(word)
    for i, c in enumerate(word, r):
        for a in range(r):
            if c >> a & 1:
                rows[a] |= 1 << i
            else:
                rows[i] |= 1 << a
        for j in range(r, i):
            rows[j] |= 1 << i
    return rows


def hp(rows):
    n = len(rows)
    d = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        d[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            z = d[mask][v]
            if not z:
                continue
            avail = rows[v] & ~mask
            while avail:
                bit = avail & -avail
                avail -= bit
                d[mask | bit][bit.bit_length() - 1] += z
    return sum(d[-1])


def literal(core, word):
    r = len(core)
    rows = graph(core, word)
    total = 0
    for order in permutations(range(r)):
        for alloc in product(range(r + 1), repeat=len(word)):
            path = []
            for j in range(r + 1):
                path += [r + i for i, a in enumerate(alloc) if a == j]
                if j < r:
                    path.append(order[j])
            total += all(rows[a] >> b & 1 for a, b in zip(path, path[1:]))
    return total


def boundary_num(word, r, all_slots=True):
    # Integer number of successful allocations, denominator (r+1)^len(word).
    d = [0] * (1 << (r + 1))
    d[0] = 1
    for c in word:
        out = [0] * len(d)
        for seen, z in enumerate(d):
            if not z:
                continue
            for j in range(r + 1):
                if j == 0 or seen >> j & 1 or c >> (j - 1) & 1:
                    out[seen | 1 << j] += z
        d = out
    if all_slots:
        return d[-1]
    return d[-1] + d[-2]  # All constrained labels, free may be absent.


def reverse_complement(word, r):
    return tuple(((1 << r) - 1) ^ c for c in reversed(word))


def boundary_ifs(word):
    # U numerator in the 3-coordinate affine recursion; common denominator 3^m.
    a = b = u = 0
    den = 1
    for c in reversed(word):
        a, b, u = 2*a + (c & 1)*den, 2*b + ((c >> 1) & 1)*den, u + (c & 1)*b + ((c >> 1) & 1)*a
        den *= 3
    return a, b, u


# Independent exact path implementation: the inherited 36-state grammar.
TRANS = [[[] for _ in range(36)] for _ in range(4)]
ACCEPT = []
for o, x, y, z in product(range(2), range(3), range(3), range(2)):
    state = o*18 + x*6 + y*2 + z
    if x != 2 and (y == 1 or (y == 0 and o == 0)):
        ACCEPT.append(state)
    for c in range(4):
        p, q = (c >> o) & 1, (c >> (1-o)) & 1
        TRANS[c][state].append(o*18 + (2 if p else 1)*6 + y*2 + z)
        if y or p:
            TRANS[c][state].append(o*18 + x*6 + (2 if q else 1)*2 + z)
        if z or q:
            TRANS[c][state].append(o*18 + x*6 + y*2 + 1)
INIT = (1,) + (0,)*17 + (1,) + (0,)*17


def transfer(word):
    v = INIT
    for c in word:
        out = [0]*36
        for source, z in enumerate(v):
            if z:
                for dest in TRANS[c][source]:
                    out[dest] += z
        v = out
    return sum(v[j] for j in ACCEPT)


def all_cores(r):
    edges = [(i,j) for i in range(r) for j in range(i+1,r)]
    for mask in range(1 << len(edges)):
        rows = [0]*r
        for k,(i,j) in enumerate(edges):
            if mask >> k & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
        yield rows


def moment(r,q):
    @lru_cache(None)
    def solve(masks):
        if not any(masks):
            return F(1)
        total = F(0)
        stay = F(0)
        for choices in product(range(r+1), repeat=q):
            new = []
            requests = 0
            for mask,j in zip(masks,choices):
                active = (mask & (1 << (j-1))) if j else 0
                requests |= active
                new.append(mask ^ active)
            prob = F(1, (r+1)**q * 2**requests.bit_count())
            new = tuple(new)
            if new == masks:
                stay += prob
            else:
                total += prob*solve(new)
        require(stay < 1, 'nonabsorbing self-loop')
        return total/(1-stay)
    return solve(((1 << r)-1,)*q)


def finite_moment_direct_letters(r,q,steps):
    # Separate forward killed-chain iteration; explicitly reads every actual letter.
    states = list(product(range(1 << r), repeat=q))
    ix = {a:i for i,a in enumerate(states)}
    kernels = []
    for masks in states:
        counts = defaultdict(int)
        for choices in product(range(r+1), repeat=q):
            for letter in range(1 << r):
                new = []
                valid = True
                for mask,j in zip(masks,choices):
                    active = (mask & (1 << (j-1))) if j else 0
                    if active and not letter & active:
                        valid = False
                        break
                    new.append(mask ^ active)
                if valid:
                    counts[ix[tuple(new)]] += 1
        kernels.append(list(counts.items()))
    d = [0]*len(states)
    d[ix[((1 << r)-1,)*q]] = 1
    base = (r+1)**q * (1 << r)
    for _ in range(steps):
        new = [0]*len(states)
        for i,z in enumerate(d):
            if z:
                for j,c in kernels[i]:
                    new[j] += z*c
        d = new
    return F(d[ix[(0,)*q]],base**steps)


def tensor_path_moments(max_m):
    # Sum of H^2 over ALL words, by the tensor square, without moment recursion.
    kernel = []
    for i in range(36):
        for j in range(36):
            counts = defaultdict(int)
            for c in range(4):
                for a in TRANS[c][i]:
                    for b in TRANS[c][j]:
                        counts[a*36+b] += 1
            kernel.append(list(counts.items()))
    v = [a*b for a in INIT for b in INIT]
    accepted = [a*36+b for a in ACCEPT for b in ACCEPT]
    out = {0:sum(v[j] for j in accepted)}
    for m in range(1,max_m+1):
        nxt = [0]*1296
        for i,z in enumerate(v):
            if z:
                for j,c in kernel[i]:
                    nxt[j] += z*c
        v = nxt
        out[m] = sum(v[j] for j in accepted)
    return out


def main():
    receipt = {'scope':'finite exact controls, all-order proofs in PROOF.md'}
    physical = literal_checks = inequalities = affine = 0
    first_squares = {}
    # Complete two-apex words through m=5, independently counted on vertices.
    for m in range(6):
        sumh = sumh2 = 0
        for word in product(range(4), repeat=m):
            h = transfer(word)
            require(h == hp(graph([2,0], word)), 'two-apex subset DP')
            physical += 1
            if m <= 4:
                require(h == literal([2,0], word), 'literal two-apex paths')
                literal_checks += 1
            a,b,u = boundary_ifs(word)
            require(u == boundary_num(word,2,False), 'affine U recursion')
            affine += 1
            for p in range(m+1):
                q = m-p
                good = 2*boundary_num(word[:p],2)*boundary_num(reverse_complement(word[p:],2),2)
                error = h-good
                bound = 6*(2**p*3**q+3**p*2**q)
                require(0 <= error <= bound, 'rank-one all-split bound')
                inequalities += 1
            sumh += h
            sumh2 += h*h
        if m:
            require(F(sumh,4**m) == F(3**m+3*2**m+3,8), 'inherited exact mean')
        first_squares[m] = sumh2
    receipt['two_apex_full_words_by_vertex_DP'] = physical
    receipt['literal_two_apex_counts'] = literal_checks
    receipt['all_split_inequalities'] = inequalities
    receipt['affine_system_checks'] = affine
    # Different core orders/orientations, including zero-chain and cyclic controls.
    generic = 0
    for r,max_m in [(1,6),(2,3),(3,2)]:
        for core in all_cores(r):
            for m in range(max_m+1):
                for word in product(range(1 << r), repeat=m):
                    h = hp(graph(core,word))
                    require(h == literal(core,word), 'generic literal slots')
                    p = m//2
                    q = m-p
                    good = factorial(r)*boundary_num(word[:p],r)*boundary_num(reverse_complement(word[p:],r),r)
                    bound = factorial(r)*(r+1)*(r**p*(r+1)**q+(r+1)**p*r**q)
                    require(0 <= h-good <= bound, 'generic factorization')
                    generic += 1
    receipt['generic_all_core_physical_checks'] = generic
    # Nontrivial boundary blocks and longer physical controls.
    rng = random.Random(192031)
    longer = []
    for r in (2,3):
        for k in range(28):
            m = rng.randrange(2*(r+1),11)
            core = list(all_cores(r))[rng.randrange(1 << (r*(r-1)//2))]
            word = tuple(rng.randrange(1 << r) for _ in range(m))
            h = hp(graph(core,word))
            p = rng.randrange(r+1,m-r)
            q = m-p
            good = factorial(r)*boundary_num(word[:p],r)*boundary_num(reverse_complement(word[p:],r),r)
            bound = factorial(r)*(r+1)*(r**p*(r+1)**q+(r+1)**p*r**q)
            require(0 <= h-good <= bound, 'longer physical boundary bound')
            if r == 2 and core == [2,0]:
                require(h == transfer(word), 'larger transfer versus DP')
            # One named core-edge reversal, with all attachment incidences unchanged.
            flipped = list(core)
            flipped[0] ^= 1 << 1
            flipped[1] ^= 1 << 0
            h2 = hp(graph(flipped,word))
            require(abs(h-h2) <= factorial(r-1)*r**m, 'core-edge Lipschitz bound')
            longer.append({'r':r,'m':m,'word':list(word),'core':core,'H':h,'core_flip_H':h2,'good':good})
    receipt['larger_physical_controls'] = longer
    hostiles = []
    for m in range(1,11):
        word = (1,)*m
        a = hp(graph([2,0],word))
        b = hp(graph([0,1],word))
        require(a == 1 and b == 2**m+1, 'sharp lower-scale hostile')
        hostiles.append([m,a,b])
    receipt['sharp_core_flip_hostiles'] = hostiles
    # Exact rational infinite-boundary moments and independent finite-time method.
    moments = {}
    for r,max_q in [(1,4),(2,4),(3,3)]:
        moments[str(r)] = {}
        for q in range(1,max_q+1):
            z = moment(r,q)
            if r == 1:
                require(z == F(1,q+1), 'uniform one-apex control')
            if q <= 2:
                finite = finite_moment_direct_letters(r,q,60)
                require(0 <= z-finite <= q*r*F(r,r+1)**60, 'independent finite absorption moment')
            moments[str(r)][str(q)] = str(z)
    require(moments['2'] == {'1':'1/4','2':'29/320','3':'1177/29440','4':'237365069/11787776000'}, 'two-apex moment constants')
    require(moments['3'] == {'1':'1/8','2':'79/3360','3':'951449/167336960'}, 'three-apex moment constants')
    receipt['boundary_moments'] = moments
    mean = 2*moment(2,1)**2
    second = 4*moment(2,2)**2
    muplus = moment(2,1)*F(11,9)
    secondplus = moment(2,2)*F(35,27)
    strongmean = 2*muplus**2
    strongsecond = 4*secondplus**2
    require((mean,second,second-mean**2) == (F(1,8),F(841,25600),F(441,25600)), 'path limit constants')
    require((strongmean,strongsecond,strongsecond-strongmean**2) == (F(121,648),F(41209,746496),F(136625,6718464)), 'strong limit constants')
    receipt['two_apex_limits'] = {'mean':str(mean),'second_moment':str(second),'variance':str(second-mean**2),'strong_mean':str(strongmean),'strong_second_moment':str(strongsecond),'strong_variance':str(strongsecond-strongmean**2)}
    tensor = tensor_path_moments(100)
    for m,z in first_squares.items():
        require(z == tensor[m], 'tensor second moment versus complete words')
    tensor_rows = []
    for m in (20,40,60,80,100):
        val = F(tensor[m],36**m)
        p = m//2
        q = m-p
        # q_moment=2, maximum observable 2; multiply coupling bound by 4.
        bound = 48*(F(2,3)**p + F(2,3)**q)
        require(abs(val-second) <= bound, 'tensor-square limiting moment')
        tensor_rows.append({'m':m,'normalized_second_moment':str(val),'limit_error':str(val-second),'proved_error_bound':str(bound)})
    receipt['tensor_square_independent_checks'] = tensor_rows
    receipt['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    content = json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    Path(__file__).with_name('receipt.json').write_text(content)
    print(content,end='')

if __name__ == '__main__':
    main()
