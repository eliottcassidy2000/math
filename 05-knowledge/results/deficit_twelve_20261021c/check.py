"""No-search consumer of the complete order-twelve black-score witness bank.

The independent C++ producer searches paths. This consumer only reconstructs
each literal source and replays the supplied words and complement targets.
"""
from __future__ import annotations
import argparse
import base64
from collections import Counter
import hashlib
import itertools as it
import json
from pathlib import Path
import zlib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
N = 12
FRAME_COUNT = 1 << 18
BALANCED_COUNT = 151040
BANK_SIZE = 8 + BALANCED_COUNT * 8 * 15
WORD_BANK_SIZE = 8 + BALANCED_COUNT * 8 * 6
GATES = 0
PAIRS = tuple(it.combinations(range(N), 2))
PAIR_BIT = [[0]*N for _ in range(N)]
for _i, (_a, _b) in enumerate(PAIRS):
    PAIR_BIT[_a][_b] = PAIR_BIT[_b][_a] = 1 << _i
ALL_EDGES = (1 << len(PAIRS))-1


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def integer(x, low=0, high=None):
    if type(x) is not int or x < low or high is not None and x > high:
        raise ValueError('exact integer range')
    return x


def sha(data):
    return hashlib.sha256(data).hexdigest()


def lf(path):
    return path.read_bytes().replace(b'\r\n', b'\n')


def source(t):
    if type(t) is not tuple or not t:
        raise ValueError('immutable nonempty source')
    n = len(t)
    for i, row in enumerate(t):
        integer(row, 0, (1 << n)-1)
        if row >> i & 1:
            raise ValueError('loop')
        for j in range(i):
            if (row >> j & 1)+(t[j] >> i & 1) != 1:
                raise ValueError('not a tournament')
    return t


def arc(t, a, b):
    return (t[a] >> b) & 1


def path(t, p):
    n = len(t)
    if type(p) is not tuple or len(p) != n or any(type(x) is not int for x in p) or sorted(p) != list(range(n)):
        raise ValueError('literal source permutation')
    if any(not arc(t, a, b) for a, b in zip(p, p[1:])):
        raise ValueError('non-native Hamilton word')
    return p


def graph_code(t):
    return sum(arc(t, a, b) << k for k, (a, b) in enumerate(it.combinations(range(len(t)), 2)))


def graph_decode(n, code):
    integer(n, 1)
    integer(code, 0, (1 << (n*(n-1)//2))-1)
    rows = [0]*n
    for k, (a, b) in enumerate(it.combinations(range(n), 2)):
        x, y = (a, b) if code >> k & 1 else (b, a)
        rows[x] |= 1 << y
    return tuple(rows)


def target_code(t, p):
    path(t, p)
    n = len(t)
    edge_index = {(a, b): k for k, (a, b) in enumerate(it.combinations(range(n), 2))}
    kept = sum(1 << edge_index[tuple(sorted((a, b)))] for a, b in zip(p, p[1:]))
    return graph_code(t) ^ (((1 << (n*(n-1)//2))-1) ^ kept)


def make_layout():
    """Independent per-edge canonical orbit reader, not producer's orbit sweep."""
    reps = {}
    descriptors = []
    for a, b in PAIRS:
        rep = min(tuple(sorted(((a+4*j)%N, (b+4*j)%N))) for j in range(3))
        c, d = rep
        forward = any(((c+4*j)%N, (d+4*j)%N) == (a, b) for j in range(3))
        reps.setdefault(rep, []).append((a, b, forward))
    kinds = {}
    free = 0
    for rep, edges in sorted(reps.items()):
        values = set()
        for a, b, f in edges:
            if b == a+1:
                values.add(int(f))
            if a == 0 and b == N-1:
                values.add(int(not f))
        if len(values) > 1:
            raise ValueError('inconsistent prescribed cycle')
        if values:
            kinds[rep] = ('fixed', values.pop())
        else:
            kinds[rep] = ('free', free)
            free += 1
    for a, b in PAIRS:
        rep = min(tuple(sorted(((a+4*j)%N, (b+4*j)%N))) for j in range(3))
        c, d = rep
        f = any(((c+4*j)%N, (d+4*j)%N) == (a, b) for j in range(3))
        descriptors.append((a, b, f, *kinds[rep]))
    return tuple(descriptors), tuple((rep, tuple(edges), kinds[rep]) for rep, edges in sorted(reps.items())), free


LAYOUT, ORBITS, FREE_BITS = make_layout()


def frame(mask):
    integer(mask, 0, FRAME_COUNT-1)
    rows = [0]*N
    for a, b, forward, kind, k in LAYOUT:
        bit = k if kind == 'fixed' else (mask >> k) & 1
        if bit == forward:
            rows[a] |= 1 << b
        else:
            rows[b] |= 1 << a
    return tuple(rows)


def balanced(t):
    hist = Counter(row.bit_count() for row in t)
    return all(hist[d] == hist[len(t)-1-d] for d in range(len(t)))


def score_ports(t, ports, required=8):
    """Authenticate a sufficient black-orbit lower bound, with full targets.

    This generic source API does not assume or certify self-converseness.
    It does not infer equality of path orbits from equality of score strings.
    """
    source(t)
    integer(required, 1)
    if type(ports) is not tuple or len(ports) != required:
        raise ValueError('exact supplied port count')
    n = len(t)
    degrees = tuple(row.bit_count() for row in t)
    seen = set()
    for item in ports:
        if type(item) is not tuple or len(item) != 2:
            raise ValueError('word and literal target code required')
        p, dest = item
        path(t, p)
        integer(dest, 0, (1 << (n*(n-1)//2))-1)
        sig = tuple(degrees[x] for x in p)
        if all(sig[i]+sig[n-1-i] == n-1 for i in range(n)):
            raise ValueError('score string does not certify blackness')
        if sig in seen:
            raise ValueError('score strings do not distinguish these ports')
        if target_code(t, p) != dest:
            raise ValueError('wrong complement target')
        seen.add(sig)
    return required


def audit_record(record):
    if type(record) is not dict or set(record) != {'mask', 'source', 'ports'}:
        raise ValueError('source record schema')
    t = source(record['source'])
    if t != frame(record['mask']):
        raise ValueError('source substitution')
    return score_ports(t, record['ports'])


def load_inputs():
    meta = json.loads((HERE/'inputs.json').read_text(encoding='utf-8'))
    if type(meta) is not dict or set(meta) != {'schema','frame_count','balanced_count','record_bytes','bank_bytes','bank_sha256','word_bank_bytes','word_bank_sha256','producer_sha256','pins','bank_base85_zlib'}:
        raise ValueError('input schema')
    need(meta['schema'] == 'D12SCR01', 'bank schema')
    for key, val in [('frame_count', FRAME_COUNT), ('balanced_count', BALANCED_COUNT), ('record_bytes', 15), ('bank_bytes', BANK_SIZE), ('word_bank_bytes', WORD_BANK_SIZE)]:
        integer(meta[key]); need(meta[key] == val, 'exact bank dimension')
    packed = base64.b85decode(meta['bank_base85_zlib'].encode('ascii'))
    dec = zlib.decompressobj()
    raw = dec.decompress(packed, WORD_BANK_SIZE+1)
    need(len(raw) == WORD_BANK_SIZE and dec.eof and not dec.unused_data and not dec.unconsumed_tail, 'complete bounded word bank')
    need(raw[:8] == b'D12SCR01' and sha(raw) == meta['word_bank_sha256'], 'raw source-ordered word bank')
    need(sha(lf(HERE/'census.cpp')) == meta['producer_sha256'], 'frozen producer pin')
    for p, expected in meta['pins'].items():
        need(sha(lf(ROOT/p)) == expected, 'inherited input pin '+p)
    return meta, raw


def unpack(raw, offset):
    code = int.from_bytes(raw[offset:offset+6], 'little')
    return tuple((code >> (4*j)) & 15 for j in range(N))


def refused(fn, *args):
    try:
        fn(*args)
    except (ValueError, TypeError):
        need(True, 'hostile rejected')
    else:
        need(False, 'hostile accepted')


def blue(t, p):
    path(t, p)
    rev = dict(zip(p, p[::-1]))
    return all(arc(t, rev[a], rev[b]) != arc(t, a, b) for a, b in it.combinations(range(len(t)), 2))


def literal_controls(first):
    # Positive boundary: the actual six-vertex deficit-six center.
    t = (36,17,10,19,37,14)
    hp = tuple(p for p in it.permutations(range(6)) if all(arc(t,a,b) for a,b in zip(p,p[1:])))
    au = tuple(p for p in it.permutations(range(6)) if all(arc(t,p[a],p[b]) == arc(t,a,b) for a,b in it.combinations(range(6),2)))
    b = sum(blue(t,p) for p in hp)
    witnesses = {}
    for p in hp:
        sig = tuple(t[x].bit_count() for x in p)
        if any(sig[i]+sig[5-i] != 5 for i in range(6)):
            witnesses.setdefault(sig, (p,target_code(t,p)))
    ports = tuple(witnesses.values())
    need((len(hp),len(au),b,len(witnesses)) == (45,3,27,6), 'independent six-center positive boundary')
    need(score_ports(t,ports,6) == 6, 'six-score certificate accepted')
    refused(score_ports,t,ports,8)
    refused(score_ports,t,ports+(ports[0],ports[1]),8)
    # Regular Paley7 has black paths, but its score observation is constant.
    p7 = tuple(sum(1 << ((i+d)%7) for d in (1,2,4)) for i in range(7))
    p = (0,1,2,4,5,6,3)
    need(not blue(p7,p), 'literal black regular hostile')
    refused(score_ports,p7,((p,target_code(p7,p)),),1)
    # Source/port/target typing is checked before arithmetic or hashing.
    need(audit_record(first)==8, 'public source receipt')
    bad = dict(first); bad['mask'] = float(first['mask']); refused(audit_record,bad)
    bad = dict(first); bad['source'] = tuple(float(x) for x in first['source']); refused(audit_record,bad)
    bad = dict(first); bad['source'] = frame(first['mask']^1); refused(audit_record,bad)
    ports = first['ports']; p,dest = ports[0]
    for malformed in (p[:-1]+(p[0],), p[:-1]+(-1,), tuple(float(x) for x in p), tuple(bool(x) for x in p), p[::-1]):
        refused(score_ports,first['source'],((malformed,dest),)+ports[1:])
    refused(score_ports,first['source'],((p,dest^1),)+ports[1:])
    refused(score_ports,first['source'],((p,float(dest)),)+ports[1:])
    refused(score_ports,first['source'],ports[:-1]+(ports[0],))
    refused(frame,False); refused(frame,1.0); refused(frame,-1); refused(frame,FRAME_COUNT)
    refused(audit_record,dict(first, score_claim='black'))
    return {'six_center':{'source':t,'H':45,'Aut':3,'blue_H':27,'delta':6,'score_port_count':6},
            'regular_black_hostile':{'source':p7,'word':p,'score_string':[3]*7},
            'first_full_record':first}


def verify():
    global GATES
    GATES = 0
    meta, raw = load_inputs()
    need(len(ORBITS)==22 and FREE_BITS==18 and sum(x[2][0]=='fixed' for x in ORBITS)==4, 'complete normalized edge universe')
    need(all(len(edges)==3 for _,edges,_ in ORBITS), 'odd edge orbits')
    # Direct mask injectivity: each free bit changes exactly its three pairs.
    zero = frame(0)
    for k in range(FREE_BITS):
        f = frame(1 << k)
        need(sum(arc(zero,a,b)!=arc(f,a,b) for a,b in PAIRS)==3, 'free orbit bit remains visible')
    # n=12, m odd, q=12/m <=Delta+2=8 leaves only (4,3).
    need([(12//m,m) for m in range(1,13,2) if 12%m==0 and 12//m<=8] == [(4,3)], 'complete cut-index reduction')
    asymmetric = balanced_n = 0
    offset = 8
    port_hash = hashlib.sha256(); source_hash = hashlib.sha256()
    expanded_hash = hashlib.sha256(b'D12SCR01')
    first = None
    degrees_hist = Counter()
    for mask in range(FRAME_COUNT):
        t = frame(mask)
        source_hash.update(b''.join(x.to_bytes(2,'little') for x in t))
        need(all(arc(t,i,(i+1)%N) for i in range(N)), 'literal Hamilton cycle and hence strongness')
        need(all(arc(t,a,b)==arc(t,(a+4)%N,(b+4)%N) for a,b in PAIRS), 'literal rotation automorphism')
        deg = tuple(x.bit_count() for x in t)
        hist = Counter(deg)
        if any(hist[d] != hist[N-1-d] for d in range(N)):
            asymmetric += 1
            continue
        balanced_n += 1
        degrees_hist[','.join(map(str,sorted(hist)))] += 1
        src_code = graph_code(t)
        sigs = set(); saved = []
        for j in range(8):
            p = unpack(raw,offset)
            encoded_word = raw[offset:offset+6]; offset += 6
            path(t,p)
            need(True, 'supplied native Hamilton permutation')
            sig = tuple(deg[x] for x in p)
            need(any(sig[i]+sig[N-1-i] != N-1 for i in range(N)), 'literal score blackness')
            need(sig not in sigs, 'distinct automorphism-invariant score words')
            sigs.add(sig)
            kept = sum(PAIR_BIT[a][b] for a,b in zip(p,p[1:]))
            dest = src_code ^ (ALL_EDGES ^ kept)
            expanded_hash.update(encoded_word); expanded_hash.update(dest.to_bytes(9,'little'))
            need(dest.bit_length() <= 66, 'exact target-code range')
            port_hash.update(bytes(p)); port_hash.update(dest.to_bytes(9,'little'))
            if first is None:
                saved.append((p,dest))
        if first is None:
            first = {'mask':mask,'source':t,'ports':tuple(saved)}
    need(offset == len(raw), 'no omitted or trailing source-ordered witnesses')
    need(expanded_hash.hexdigest()==meta['bank_sha256'], 'lossless target reconstruction matches independent producer bank')
    need((asymmetric,balanced_n)==(111104,151040), 'complete frame partition')
    controls = literal_controls(first)
    # Independent target construction on the displayed source (no bit-XOR).
    for p,dest in first['ports']:
        keep = {frozenset(e) for e in zip(p,p[1:])}
        rows = [0]*N
        for a,b in PAIRS:
            forward = arc(first['source'],a,b) if frozenset((a,b)) in keep else 1-arc(first['source'],a,b)
            x,y=(a,b) if forward else (b,a);rows[x]|=1<<y
        need(tuple(rows)==graph_decode(N,dest), 'independent literal target replay')
        need(not blue(first['source'],p), 'literal anti-map blackness independent of score test')
    return {'schema':'deficit-twelve-20261021c','checks':GATES,
            'status':'FINITE-EXACT complete necessary normalized universe; no order-twelve center',
            'frame_count':FRAME_COUNT,'nonselfconverse_score_rejections':asymmetric,
            'eight_black_score_certificates':balanced_n,'survivors':0,
            'literal_word_count':8*balanced_n,'raw_bank_sha256':meta['bank_sha256'],
            'source_rows_sha256':source_hash.hexdigest(),'word_target_sha256':port_hash.hexdigest(),
            'score_class_sets':dict(sorted(degrees_hist.items())),
            'remaining_orders':[14,15,18,20,21,24,25,28,35],
            'scope':'No Aut enumeration, isomorphism-class census, or classification of remaining orders.',
            'controls':controls,'pins':meta['pins'],'producer_sha256':meta['producer_sha256']}


def encoded(result):
    return (json.dumps(result,sort_keys=True,indent=2)+'\n').encode('utf-8')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--verify',action='store_true')
    ap.add_argument('--write',action='store_true')
    args=ap.parse_args()
    if args.verify and args.write:
        ap.error('choose immutable verification or write')
    data=encoded(verify())
    if args.write:
        (HERE/'receipt.json').write_bytes(data)
    if args.verify and data != lf(HERE/'receipt.json'):
        raise ValueError('saved receipt mismatch')
    print(f'{GATES:,} exact checks; LF SHA256 {sha(data)}')


if __name__=='__main__':
    main()
