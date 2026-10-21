"""Complete no-search C35 certificate consumer with exact compiled keys."""
from __future__ import annotations
import argparse
import base64
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import zlib

HERE = Path(__file__).resolve().parent
N = 35
FRAMES = 1 << 17
PAIRS = tuple((i, j) for i in range(N) for j in range(i+1, N))
WIDTH = len(PAIRS)
ALL = (1 << WIDTH)-1
GATES = 0

def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)

def integer(x, low, high):
    if type(x) is not int or not low <= x <= high:
        raise ValueError('literal integer range')

def word(text):
    if type(text) is not str or len(text) != N or set(text) != {chr(65+i) for i in range(N)}:
        raise ValueError('literal complete source permutation')
    return tuple(ord(c)-65 for c in text)

def source(mask):
    integer(mask, 0, FRAMES-1)
    connection = {35-d if mask >> (d-1) & 1 else d for d in range(1, 18)}
    return tuple(sum(1 << j for j in range(N) if (j-i) % N in connection) for i in range(N))

def key_literal(rows, p):
    return sum(((rows[p[i]] >> p[j]) & 1) << b for b, (i, j) in enumerate(PAIRS))

def key_compiler(p):
    base = 0
    columns = [0]*17
    for b, (i, j) in enumerate(PAIRS):
        d = (p[j]-p[i]) % N
        base |= int(d < 18) << b
        columns[min(d, N-d)-1] |= 1 << b
    return base, tuple(columns)

@lru_cache(None)
def _compile(text):
    p = word(text)
    base, columns = key_compiler(p)
    reverse_base, reverse_columns = key_compiler(p[::-1])
    constant = base | ((base ^ reverse_base ^ ALL) << WIDTH)
    cols = tuple(c | ((c ^ r) << WIDTH) for c, r in zip(columns, reverse_columns))
    demand_zero = demand_one = 0
    for a, b in zip(p, p[1:]):
        d = (b-a) % N
        if d < 18:
            demand_zero |= 1 << (d-1)
        else:
            demand_one |= 1 << (N-d-1)
    return constant, cols, demand_zero, demand_one

def compile_word(text):
    word(text)  # Public literal validation before any cache lookup.
    return _compile(text)

def evaluate(compiled, mask):
    value, columns, zero, one = compiled
    bits = mask
    while bits:
        bit = bits & -bits
        value ^= columns[bit.bit_length()-1]
        bits -= bit
    return value

def audit_words(mask, supplied):
    integer(mask, 0, FRAMES-1)
    if type(supplied) is not tuple or len(supplied) != 7:
        raise ValueError('seven supplied words')
    keys = set()
    for text in supplied:
        c = compile_word(text)
        _, _, zero, one = c
        if mask & zero or mask & one != one:
            raise ValueError('actual Hamilton arc rejected')
        value = evaluate(c, mask)
        if value >> WIDTH == 0:
            raise ValueError('blue word rejected')
        keys.add(value & ALL)
    if len(keys) != 7:
        raise ValueError('seven full Aut orbits required')
    return 7

def load():
    m = json.loads((HERE/'inputs35.json').read_text(encoding='utf-8'))
    need(type(m) is dict and set(m) == {'schema','frames','words','index_bytes','index_sha256','expanded_word_sha256','index_base85_zlib','discovery_only'}, 'input schema')
    integer(m['frames'], FRAMES, FRAMES)
    integer(m['index_bytes'], FRAMES*14, FRAMES*14)
    need(m['schema'] == 'C35-seven-black-v1', 'literal schema')
    need(type(m['words']) is list and len(set(m['words'])) == len(m['words']), 'distinct literal word bank')
    for w in m['words']:
        word(w)
    packed = base64.b85decode(m['index_base85_zlib'].encode('ascii'))
    dec = zlib.decompressobj()
    raw = dec.decompress(packed, FRAMES*14+1)
    need(len(raw) == FRAMES*14 and dec.eof and not dec.unused_data and not dec.unconsumed_tail, 'complete bounded index bank')
    need(hashlib.sha256(raw).hexdigest() == m['index_sha256'], 'source-ordered index hash')
    return m, raw

def refused(fn, *args):
    try:
        fn(*args)
    except (ValueError, TypeError):
        need(True, 'hostile rejected')
    else:
        need(False, 'hostile accepted')

def run():
    global GATES
    GATES = 0
    meta, raw = load()
    texts = tuple(meta['words'])
    compiled = tuple(compile_word(t) for t in texts)
    # Independently reconstruct all literal source matrices on an affine basis.
    # Equality on this basis verifies every coefficient of each affine key.
    basis = (0,) + tuple(1 << j for j in range(17))
    native_controls = 0
    for mask in basis:
        rows = source(mask)
        need(all(r.bit_count() == 17 for r in rows), 'source regularity')
        need(all(((rows[i] >> j) ^ (rows[j] >> i)) & 1 for i, j in PAIRS), 'source tournament')
        for text, c in zip(texts, compiled):
            p = word(text)
            key = key_literal(rows, p)
            reverse = key_literal(rows, p[::-1])
            expected = key | ((key ^ reverse ^ ALL) << WIDTH)
            need(evaluate(c, mask) == expected, 'all affine coefficients independently reconstructed')
            actual_path = all(rows[a] >> b & 1 for a, b in zip(p, p[1:]))
            need(actual_path == (not(mask & c[2]) and mask & c[3] == c[3]), 'literal path clause reconstruction')
            native_controls += actual_path
    expanded = hashlib.sha256()
    uses = [0]*len(texts)
    first = None
    for mask in range(FRAMES):
        ids = tuple(int.from_bytes(raw[14*mask+2*j:14*mask+2*j+2], 'little') for j in range(7))
        need(all(i < len(texts) for i in ids), 'literal word-bank indices')
        keys = set()
        for idx in ids:
            c = compiled[idx]
            need(not(mask & c[2]) and mask & c[3] == c[3], 'native Hamilton word on this exact mask')
            val = evaluate(c, mask)
            need(val >> WIDTH != 0, 'actual positional reversal is not anti')
            keys.add(val & ALL)
            expanded.update(texts[idx].encode('ascii'))
            uses[idx] += 1
        need(len(keys) == 7, 'seven distinct complete ordered keys')
        if first is None:
            first = tuple(texts[i] for i in ids)
        # Separate literal consumers on deterministic whole-source checkpoints.
        if mask % 8191 == 0 or mask == FRAMES-1:
            rows = source(mask)
            literal_keys = set()
            for i in ids:
                p = word(texts[i])
                need(all(rows[a] >> b & 1 for a, b in zip(p,p[1:])), 'direct matrix checkpoint path')
                k = key_literal(rows,p)
                need(k ^ key_literal(rows,p[::-1]) != ALL, 'direct matrix checkpoint blackness')
                literal_keys.add(k)
            need(len(literal_keys)==7, 'direct matrix checkpoint orbit keys')
    need(all(uses), 'no unused bank words')
    need(expanded.hexdigest() == meta['expanded_word_sha256'], 'lossless literal stream commitment')
    need(audit_words(0, first) == 7, 'public source-bound consumer')
    for bad in (False, 0.0, -1, FRAMES):
        refused(source,bad)
        refused(audit_words,bad,first)
    refused(audit_words,0,list(first))
    refused(audit_words,0,(first[0],)*7)
    refused(audit_words,FRAMES-1,first)
    refused(compile_word,first[0][:-1]+first[0][0])
    refused(compile_word,tuple(range(35)))
    refused(compile_word,False)
    # Full Aut rather than a chosen transitive subgroup is essential in the
    # lexicographic branch. Over F5, the augmentation code on seven blocks
    # has no nonzero constant vector and every shifted norm sum is zero.
    for i in range(6):
        v = [0]*7
        v[i], v[6] = 1,4
        total = [sum(v[(j-k)%7] for k in range(7))%5 for j in range(7)]
        need(total == [0]*7, 'augmentation shift has no order35 lift')
    need(all((7*a)%5 for a in range(1,5)), 'no nonzero constant augmentation vector')
    # Small permitted parameter theorem: neither cross prime divides p-1.
    need(all((p-1)%q for p,q in ((5,7),(7,5))), 'centralizing q-lift guard')
    after25 = [12,14,15,18,20,21,24,28,35]
    after35 = [n for n in after25 if n != 35]
    need(max(after35) == 28, 'own combined center-degree bound')
    return {'schema':'C35-complete-seven-black-v1','checks':GATES,
            'frames':FRAMES,'literal_word_positions':7*FRAMES,'distinct_literal_words':len(texts),
            'black_orbits_per_frame_at_least':7,'delta_even_hence_at_least':8,
            'affine_basis_sources':len(basis),'native_basis_controls':native_controls,
            'expanded_word_sha256':expanded.hexdigest(),'index_sha256':meta['index_sha256'],
            'input_lf_sha256':hashlib.sha256((HERE/'inputs35.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
            'discovery_only':meta['discovery_only'],
            'remaining_after_this_folder_alone':after35,
            'center_degree_at_most':28,'safe_recurrence_onset':41,
            'scope':'All VT35 sources by the proved reduction with CITED odd-order solvability. No isomorphism or full H/Aut counts.'}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--verify',action='store_true')
    args=p.parse_args()
    result=run()
    text=json.dumps(result,sort_keys=True,indent=2)+'\n'
    path=HERE/'receipt35.json'
    if args.verify:
        need(path.read_text(encoding='utf-8')==text,'saved receipt differs')
    else:
        path.write_text(text,encoding='utf-8',newline='\n')
    print(f"PASS {result['checks']} body checks; LF SHA256 {hashlib.sha256(text.encode()).hexdigest()}")

if __name__=='__main__':
    main()
