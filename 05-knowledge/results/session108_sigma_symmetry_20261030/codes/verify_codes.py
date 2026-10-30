#!/usr/bin/env python3
"""Exact, dependency-free guards; checks remain active under python -O."""
from itertools import permutations, combinations
from collections import Counter
import json
from pathlib import Path

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def parity(x): return x.bit_count() % 2

def image(columns):
    return {sum(parity(x & z) << i for i, x in enumerate(columns)) for z in range(8)}

def syndrome(word, columns):
    out = 0
    for i, x in enumerate(columns):
        if word >> i & 1: out ^= x
    return out

def kernel(columns): return {w for w in range(1 << len(columns)) if syndrome(w, columns) == 0}
def permute(w, perm): return sum(((w >> i) & 1) << perm[i] for i in range(len(perm)))
def histogram(code): return dict(sorted(Counter(w.bit_count() for w in code).items()))
def sigma(n): return sum(d for d in range(1,n+1) if n % d == 0)
def squareclass(n):
    out = 1; p = 2
    while p*p <= n:
        e = 0
        while n % p == 0: n //= p; e ^= 1
        if e: out *= p
        p += 1
    return out*n

def main():
    bank = [2,5,7,19,29,59,89]
    columns = [2,3,1,4,7,6,5]
    classes = [squareclass(sigma(p)) for p in bank]
    require(classes == [3,6,2,5,30,15,10], 'inherited bank')
    for sc,x in zip(classes,columns):
        require(sc == (2 if x & 1 else 1)*(3 if x & 2 else 1)*(5 if x & 4 else 1), 'mask')
    C,H = image(columns),kernel(columns)
    require(histogram(C)=={0:1,4:7}, 'simplex enumerator')
    require(histogram(H)=={0:1,3:7,4:7,7:1}, 'Hamming enumerator')
    cols8 = [0]+columns
    C8,H8 = image(cols8),kernel(cols8)
    require(len(C8)==8 and len(H8)==32, 'zero extension dimensions')
    require(min(w.bit_count() for w in H8 if w)==1, 'not extended Hamming')
    require(all(not(w & 1) for w in C8), 'fixed zero response')
    # Enumerate all linear maps by images of the three standard basis vectors.
    maps=[]
    for a,b,c in permutations(range(1,8),3):
        T = tuple((a if x&1 else 0)^(b if x&2 else 0)^(c if x&4 else 0) for x in range(8))
        if len(set(T))==8: maps.append(T)
    require(len(maps)==168, 'GL order')
    auts = []
    for perm in permutations(range(7)):
        if {permute(w,perm) for w in H}==H: auts.append(perm)
    induced = {tuple(columns.index(T[x]) for x in columns) for T in maps}
    require(set(auts)==induced, 'all Hamming permutation automorphisms')
    pointed=[T for T in maps if T[4]==4]
    require(len(pointed)==24, 'point stabilizer')
    lines={frozenset((x,y,x^y)) for x,y in combinations(range(1,8),2)}
    require(len(lines)==7 and sum(len(L) for L in lines)==21, 'Fano 7/21 incidence count')
    flags={(x,L) for L in lines for x in L}
    x0,L0=next(iter(flags))
    require({(T[x0],frozenset(T[x] for x in L0)) for T in maps}==flags, 'flag transitivity')
    through=[L for L in lines if 4 in L]
    pairs=sorted(tuple(sorted(bank[columns.index(x)] for x in L if x != 4)) for L in through)
    require(pairs==[(2,59),(5,29),(7,89)], 'pointed19 pairs')
    other=[L for L in lines if 4 not in L]
    four_perms={tuple(other.index(frozenset(T[x] for x in L)) for L in other) for T in pointed}
    require(len(four_perms)==24, 'faithful S4 action')
    # RM(1,3) on all eight vector labels; NOT obtained from sigma(3)=4 alone.
    lin=image(list(range(8)))
    rm=lin|{w^255 for w in lin}
    require(histogram(rm)=={0:1,4:14,8:1}, 'RM enumerator')
    require(all(parity(u&v)==0 for u in rm for v in rm), 'RM self orthogonal')
    require(len(rm)==16, 'self dual dimension')
    affine={tuple(T[x]^t for x in range(8)) for T in maps for t in range(8)}
    require(len(affine)==1344, 'AGL order')
    full_auts=set()
    lin_auts=set()
    for perm in permutations(range(8)):
        if {permute(w,perm) for w in rm}==rm: full_auts.add(perm)
        if {permute(w,perm) for w in lin}==lin: lin_auts.add(perm)
    require(full_auts==affine, 'full RM group')
    require(lin_auts==set(maps), 'zero extended simplex group')
    require({permute(w,tuple(x^1 for x in range(8))) for w in lin}!=lin, 'translation hostile')
    require([sigma(p) for p in [3,5,19]]==[4,6,20], 'requested sigma')
    require([sigma(p)-p for p in [3,5,19]]==[1,1,1], 'aliquot collapse')
    # Basis (2,3,5,19), six rows direct then sigma.
    joint_columns=[2,4,8,0,3,4]
    joint={tuple(parity(z&v) for v in joint_columns) for z in range(16)}
    require(len(joint)==16, 'joint rank4')
    require(all(w[3]==0 and w[1]==w[5] for w in joint), 'joint complete constraints')
    require({tuple(w[i] for i in [0,1,2,4]) for w in joint}==set(__import__('itertools').product(range(2),repeat=4)), 'joint independent ports')
    require(squareclass(sigma(3*5))==squareclass(sigma(5)), 'coprime toggle')
    require(squareclass(sigma(3*3))!=squareclass(sigma(3)), 'noncoprime hostile')
    # Arithmetic parity completion on actual prime labels [3]+bank.
    even_kernel={w for w in H8 if parity(w)==0}
    rm_labelled=image(cols8)|{w^255 for w in image(cols8)}
    require(even_kernel==rm_labelled, 'arithmetic even kernel is self-dual RM')
    require(histogram(even_kernel)=={0:1,4:14,8:1}, 'arithmetic parity weights')
    fibers=Counter((parity(w),syndrome(w,cols8)) for w in range(256))
    require(len(fibers)==16 and set(fibers.values())=={16}, 'joint parity/sigma fibers')
    for w in range(256):
        require(parity(w^1)==(parity(w)^1) and syndrome(w^1,cols8)==syndrome(w,cols8), 'unique 3 toggle')
        primes=[p for i,p in enumerate([3]+bank) if w>>i&1]
        divisors=[1]
        for p in primes: divisors += [p*d for d in divisors]
        sig=sum(divisors)
        sc=squareclass(sig)
        mask=syndrome(w,cols8)
        expected=(2 if mask&1 else 1)*(3 if mask&2 else 1)*(5 if mask&4 else 1)
        require(sc==expected, 'all 256 arithmetic selections')
        require((w in even_kernel)==(sc==1 and len(primes)%2==0), 'even square-sigma arithmetic')
    for T in maps:
        for t in range(8):
            perm=tuple(cols8.index(T[x]^t) for x in cols8)
            for w in range(256):
                new=permute(w,perm)
                require(parity(new)==parity(w), 'affine preserves parity')
                require(syndrome(new,cols8)==(T[syndrome(w,cols8)]^(t if parity(w) else 0)), 'full affine syndrome equivariance')
    # Translation by mask1 swaps prime3 with prime7: odd square-sigma becomes nonsquare.
    trans=tuple(cols8.index(x^1) for x in cols8)
    require(permute(1,trans)==1<<cols8.index(1), '3 maps to 7 under translation')
    require(syndrome(1,cols8)==0 and syndrome(permute(1,trans),cols8)==1, 'old zero-fiber not preserved')
    receipt={'status':'PASS','bank':bank,'columns':columns,'Hamming_weights':histogram(H),'simplex_weights':histogram(C),'zero_extended_kernel_weights':histogram(H8),'zero_extended_image_weights':histogram(C8),'GL3_order':len(maps),'pointed19_stabilizer':len(pointed),'pointed19_pairs':pairs,'Fano_lines':len(lines),'Fano_flags':21,'RM_weights':histogram(rm),'AGL3_order':len(affine),'permutation_searches':[5040,40320],'checks_use_assert':False,'arithmetic_selections':256,'affine_syndrome_checks':1344*256,'parity_completion_weights':histogram(even_kernel),'joint_parity_sigma_fibers':len(fibers)}
    Path(__file__).with_name('verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__': main()
