#!/usr/bin/env python3
"""Independent exact audit of the guard-refined receipt. No verifier imports.

Source branch coefficients use a closed sum; residues use inverse rational
composition. All continuation certificates are rederived from integer orbits.
"""
from fractions import Fraction
from collections import Counter
from pathlib import Path
import hashlib
import json
import sys

DEFAULT_ROOT = Path(__file__).resolve().parent.parent
ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_ROOT
receipt = json.loads((ROOT / 'verification-receipt.json').read_text())
K = '110110110110110101011011101010'
F = '10'
H = '110101000'
V = K * 4 + F
banks = {'K': K, 'V': V, 'F': F, 'H': H}
checks = 0

def check(test, label):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(label)

def step(x):
    return x // 2 if not x % 2 else (3*x + 1) // 2

def independent_branch(word):
    """Closed affine constant, not recursive affine composition."""
    odd_indices = [i for i, bit in enumerate(word) if bit == '1']
    h = len(odd_indices)
    a = 3 ** h
    b = sum(2**i * 3**(h-k-1) for k, i in enumerate(odd_indices))
    return a, b, 2 ** len(word)

def independent_residue(word):
    """Pull endpoint zero back through rational inverse branches."""
    y = Fraction(0)
    for bit in word[::-1]:
        y = (2*y - 1)/3 if bit == '1' else 2*y
    M = 1 << len(word)
    r = y.numerator * pow(y.denominator, -1, M) % M
    t = r
    for bit in word:
        check(str(t % 2) == bit, 'source parity replay')
        t = step(t)
    a, b, M = independent_branch(word)
    check(t * M == a*r+b, 'closed branch against integer replay')
    return r

def continuation(r, d, limit):
    z = r + d
    a, b, M = 1, 0, 1
    bits = []
    for j in range(1, limit+1):
        parity = z & 1
        bits.append(str(parity))
        z = step(z)
        b = 3*b + M if parity else b
        a = a*3 if parity else a
        M *= 2
        if a < M:
            return j, a, b, M, ''.join(bits), z
    return None

def validate_certificate(source, r, d, data):
    j, a, b, M, word, endpoint = data
    check(independent_branch(word) == (a,b,M), 'continuation closed formula')
    check(endpoint*M == a*(r+d)+b, 'continuation integer replay')
    check(j <= len(source) and (1 << len(source)) % M == 0,
          'all-lift guard containment')
    check(a < M and r > 0 and (M-a)*r > a*d+b,
          'strict all-lift original-height inequality')
    return Fraction(a*d+b, M-a)

# The symbolic partition has seven K-started exits, six V-started exits,
# and exactly two surviving guards. Check all numerical facts it requires.
short = [(f'K^{n}{b}', K*n+banks[b]) for n in (1,2,3) for b in ('F','H')]
short += [('K^4H', K*4+H)]
short += [(f'VK^{n}{b}', V+K*n+banks[b]) for n in (0,1,2) for b in ('F','H')]
check(len(short) == 13, 'thirteen exits')
for (name, word), row in zip(short, receipt['short_descending_branches']):
    r = independent_residue(word)
    a,b,M = independent_branch(word)
    check(a < M and r > 0 and (M-a)*r > b, 'short all-lift descent')
    check(row == {'name':name,'steps':len(word),'residue':r,
                  'fixed_point':str(Fraction(b,M-a))}, 'short receipt')
for block in (K,V):
    for j in range(1,len(block)+1):
        a,b,M = independent_branch(block[:j])
        check(a > M, 'every nonempty expansion prefix expands')

C = Fraction(83829391903560839836419813590128166431,
             157489106280555890461927692448373459)
E = C+1
a,b,M = independent_branch(V)
check(C == Fraction(b,a-M), 'maximal expansion center')
check(13332 < 25*E < 13333, '25 carry budget')
check(533286 < 1000*E < 533287, '1000 carry budget')
sources = [('A',K*5),('B',V+K*3)]
sources = [(name,w,independent_residue(w)) for name,w in sources]
for name,w,r in sources:
    check(receipt['long_source_cylinders'][name] ==
          {'length':len(w),'residue':r,'modulus':1<<len(w)}, 'source receipt')

exceptions = []
max25 = Fraction(0)
maxbase = Fraction(0)
maxj25 = Counter()
maxj = 0
classes = {}
digest = hashlib.sha256()
for d in range(1,533287):
    failed = False
    for name, w, r in sources:
        result = continuation(r,d,len(w))
        if result is None:
            digest.update(f'{name},{d},NONE\n'.encode())
            failed = True
            continue
        j,a,b,M,bits,z = result
        check((M-a)*r > a*d+b, 'all offsets exact inequality')
        key = (name,bits)
        if key not in classes:
            validate_certificate(w,r,d,result)
            classes[key] = [d,d,1,M,a,b]
        else:
            classes[key][1] = d
            classes[key][2] += 1
        digest.update(f'{name},{d},{j},{a},{b}\n'.encode())
        maxj = max(maxj,j)
        threshold = Fraction(a*d+b,M-a)
        maxbase = max(maxbase,threshold)
        if d <= 13332:
            max25 = max(max25,threshold)
            maxj25[name] = max(maxj25[name],j)
    if failed:
        exceptions.append(d)
for (name,bits),(lo,hi,count,M,a,b) in classes.items():
    check((hi-lo)%M == 0 and count == (hi-lo)//M+1,
          'complete arithmetic progression per parity class')
check(exceptions == receipt['exceptional_offsets'], 'exact exceptions')
check(not any(d <= 13332 for d in exceptions), 'no main theorem exceptions')
check(dict(maxj25) == {'A':135,'B':183}, 'per-source main maximum depths')
check(maxj == 183, 'base maximum depth')
check(str(max25) == receipt['main_max_threshold'], 'main maximum threshold')
check(str(maxbase) == receipt['base_max_threshold'], 'base maximum threshold')
check(len(classes) == receipt['base_class_count'] == 33990, 'class count')
check(digest.hexdigest() == receipt['base_row_digest'], 'full offset digest')

# Audit the stored forests, not the published tree builder. Reconstruct every
# node, all parent links and children, first-crossing scope, and leaf claims.
node_count = 0
kind_counts = Counter()
tree_depth = 0
guard_depth = 0
for d in exceptions:
    tree = receipt['refinement_trees'][str(d)]
    check('K' in tree and 'V' in tree, 'both tree roots')
    for name,row in tree.items():
        node_count += 1
        check(name[0] in ('K','V'), 'valid root')
        check(all(b in banks for b in name), 'valid macroblock letters')
        if len(name)>1:
            check(name[:-1] in tree and tree[name[:-1]]['kind']=='branch',
                  'connected and parent nonterminal')
        word = ''.join(banks[b] for b in name)
        a,b,M = independent_branch(word)
        r = independent_residue(word)
        kind = row['kind']
        kind_counts[kind] += 1
        if kind == 'branch':
            check(a >= M, 'nonterminal is not an early crossing')
            check(all(name+s in tree for s in banks), 'four children present')
        else:
            check(not any(name+s in tree for s in banks), 'terminal has no children')
            if kind == 'source_descent':
                check(name[-1] in ('F','H'), 'first crossing at explicit reset')
                check(a<M and r>0 and (M-a)*r>b, 'source leaf all-lift descent')
                check(row['steps']==len(word), 'source leaf time')
            elif kind == 'offset_repair':
                check(a>=M, 'offset guard occurs strictly before first crossing')
                result = continuation(r,d,len(word))
                check(result is not None, 'offset leaf contraction exists')
                threshold = validate_certificate(word,r,d,result)
                check(row == {'kind':'offset_repair','steps':result[0],
                              'source_bits':len(word),'source_residue':r,
                              'threshold':str(threshold)}, 'offset leaf record')
                tree_depth = max(tree_depth,result[0])
                guard_depth = max(guard_depth,len(word))
            else:
                raise RuntimeError('unknown node type')
check(node_count == receipt['refinement_node_count'] == 838,'forest node count')
check(tree_depth == 273,'forest repair maximum')
check(guard_depth == receipt['refinement_max_source_bits'] == 394,'guard maximum')
check(node_count == 30 + 4*kind_counts['branch'], 'finite forest cardinality')

# Exact finite controls for the symbolic rational obstruction. Its unbounded
# residue argument and if-and-only-if are reviewed in the accompanying report.
c = Fraction(4912087879,88519643)
a,b,M = independent_branch(K)
check(-c == Fraction(-b,a-M),'K negative rational fixed point')
prior = 0
for L in range(1,17):
    word = K*L
    r = independent_residue(word)
    modulus = 1 << len(word)
    xr = (-c.numerator * pow(c.denominator,-1,modulus)) % modulus
    check(r == xr and r>=prior,'K residues are compatible rational reductions')
    prior = r
check(55<c<56,'positive rational displacement threshold')

print(json.dumps({'status':'PASS','independent_checks':checks,
 'per_source_main_depth':dict(maxj25),'max_base_depth':maxj,
 'maximum_main_threshold':str(max25),'maximum_base_threshold':str(maxbase),
 'class_count':len(classes),'exceptional_offsets':exceptions,
 'forest_node_count':node_count,'forest_kind_counts':dict(kind_counts),
 'max_refined_repair':tree_depth,'max_refined_guard':guard_depth,
 'base_row_digest':digest.hexdigest()},indent=2,sort_keys=True))
