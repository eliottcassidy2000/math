"""Full image witnesses for the new four-source trade, modulo 256.

Exact bit-plane dyadic elimination, adapted from Session44 rectangle-bank
solve_images.py. The checker validates output by literal edge multiplication.
"""
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p/'01-canon').is_dir())
N = 3528
n = N-1
bits_total = 8
modulus_total = 1 << bits_total
charges = [{1477: 1, 2837: 1, 1425: -1, 1921: -1}, {1965: 1, 2234: -1},
           {1477: 1, 1841: 1, 1394: -1, 2185: -1}]
larger = '--larger' in sys.argv
if larger:
    charges = [dict([(i,1) for i in (1480,1691,2071,2185,2581,3000)] +
                    [(i,-1) for i in (1638,1923,2176,2553,2835,2945)]),
               dict([(i,1) for i in (176,1477,1691,1974,2071,2087,2185,2905)] +
                    [(i,-1) for i in (1394,1425,1457,1638,1841,1921,2553,2581)])]
path = REPO/'05-knowledge/results/metagraph_order8_20261023/merged_edges.tsv.gz'
raw = path.read_bytes()
edges = [tuple(int(row[k]) for k in ('source','target','blue','black'))
         for row in csv.DictReader(gzip.decompress(raw).decode().splitlines(),delimiter='\t')]
rows = [[0]*bits_total for _ in range(n)]
degree = [0]*N
for i, j, b, k in edges:
    if i == j: continue
    w = b+k
    degree[i] += w; degree[j] += w
    value = -w % modulus_total
    for bit in range(bits_total):
        if value >> bit & 1 and i < n and j < n:
            rows[i][bit] |= 1 << j; rows[j][bit] |= 1 << i
for i in range(n):
    for bit in range(bits_total):
        if degree[i] >> bit & 1: rows[i][bit] |= 1 << i
for rhs, charge in enumerate(charges):
    for i, coefficient in charge.items():
        value = 2*coefficient % modulus_total
        for bit in range(bits_total):
            if value >> bit & 1: rows[i][bit] |= 1 << (n+rhs)
coefficient_mask = (1 << n)-1


def plus(x,y):
    carry = 0
    out = []
    for a,b in zip(x,y):
        z = a ^ b
        out.append(z ^ carry)
        carry = (a & b) | (z & carry)
    return tuple(out)


def times(x,c):
    k = len(x)
    out = (0,)*k
    for j in range(k):
        if c >> j & 1: out = plus(out, (0,)*j+tuple(x[:k-j]))
    return out


pivots = []
ranks = []
operations = 0
for bits in range(bits_total,0,-1):
    modulus = 1 << bits
    i = 0
    while i < len(rows):
        p = next((j for j in range(i,len(rows)) if rows[j][0] & coefficient_mask),None)
        if p is None: break
        rows[i],rows[p] = rows[p],rows[i]
        low = rows[i][0] & coefficient_mask
        col = (low & -low).bit_length()-1
        value = sum(((plane >> col)&1) << bit for bit,plane in enumerate(rows[i]))
        pivot = times(rows[i],pow(value,-1,modulus))
        rows[i] = pivot
        scaled = {}  # only scales actually requested by lower rows
        for j in range(i+1,len(rows)):
            value = sum(((plane >> col)&1) << bit for bit,plane in enumerate(rows[j]))
            if value:
                scale = -value % modulus
                if scale not in scaled: scaled[scale] = times(pivot,scale)
                rows[j] = plus(rows[j],scaled[scale]);operations += 1
        pivots.append((bits,col,pivot));i += 1
    ranks.append(i)
    print('bits',bits,'unit pivots',i,'remaining',len(rows)-i,flush=True)
    if any(row[0] for row in rows[i:]): raise RuntimeError('inconsistent residual RHS')
    rows = [tuple(row[1:]) for row in rows[i:]]
solutions = []
for rhs,charge in enumerate(charges):
    x = [0]*n
    for bits,col,row in reversed(pivots):
        modulus = 1 << bits
        target = sum(((plane >> (n+rhs))&1) << bit for bit,plane in enumerate(row))
        total = 0
        for bit,plane in enumerate(row):
            scan = plane & coefficient_mask & ~(1 << col)
            while scan:
                low = scan & -scan;scan -= low
                total += (1 << bit)*x[low.bit_length()-1]
        x[col] = (target-total) % modulus
    x.append(0)
    current = [0]*N
    for i,j,b,k in edges:
        value = (b+k)*(x[i]-x[j]);current[i] += value;current[j] -= value
    expected = [2*charge.get(i,0) for i in range(N)]
    if any((a-b) % modulus_total for a,b in zip(current,expected)):
        raise RuntimeError('literal full Laplacian certificate failed')
    solutions.append({'charge':[[i,v] for i,v in sorted(charge.items())], 'coefficient':2,
                      'potential_mod256':x,'integer_correction_div256_nonzero':
                      [[i,(a-b)//modulus_total] for i,(a,b) in enumerate(zip(current,expected)) if a != b]})
receipt = {'modulus':modulus_total,'input_edges_sha256':sha256(raw).hexdigest(),
           'unit_pivot_counts_by_modulus':dict(zip([1 << j for j in range(bits_total,0,-1)],ranks)),
           'row_operations':operations,'solutions':solutions}
(HERE/('larger_image.json' if larger else 'image_witnesses.json')).write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
print('PASS literal mod256 image witnesses',flush=True)
