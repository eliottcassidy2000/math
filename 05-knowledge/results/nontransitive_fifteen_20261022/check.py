"""No-search order-15 certificates and finite controls for the prime-shore lemma."""
from __future__ import annotations
import argparse
import base64
from collections import Counter
from functools import lru_cache
import hashlib
import itertools as it
import json
from pathlib import Path
import zlib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
N, Q, M, FREE = 15, 3, 5, 18
FRAMES = 1 << FREE
PAIRS = tuple(it.combinations(range(N), 2))
WIDTH = len(PAIRS)
ALL = (1 << WIDTH)-1
PAIR_INDEX = {p: i for i, p in enumerate(PAIRS)}
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def integer(x, lo=0, hi=None):
    if type(x) is not int or x < lo or hi is not None and x > hi:
        raise ValueError('literal integer range')
    return x


def word(text):
    if type(text) is not str or len(text) != N or set(text) != set('0123456789abcde'):
        raise ValueError('literal source permutation')
    return tuple(int(c, 16) for c in text)


def make_layout():
    # Canonical representative of each individual edge, rather than the
    # producer's visited-orbit sweep.
    groups = {}
    for a, b in PAIRS:
        rep = min(tuple(sorted(((a+Q*j)%N, (b+Q*j)%N))) for j in range(M))
        c, d = rep
        forward = any(((c+Q*j)%N, (d+Q*j)%N) == (a,b) for j in range(M))
        groups.setdefault(rep, []).append((a,b,forward))
    kinds, count = {}, 0
    for rep, edges in sorted(groups.items()):
        forced = {int(f) for a,b,f in edges if b == a+1}
        forced |= {int(not f) for a,b,f in edges if (a,b) == (0,N-1)}
        if len(forced) > 1:
            raise ValueError('inconsistent cycle gauge')
        if forced:
            kinds[rep] = ('fixed', forced.pop())
        else:
            kinds[rep] = ('free', count)
            count += 1
    layout = {}
    for rep, edges in groups.items():
        for a,b,f in edges:
            layout[a,b] = (f, *kinds[rep])
    if count != FREE or len(groups) != 21 or any(len(v) != M for v in groups.values()):
        raise ValueError('complete marked-cycle layout')
    return layout


LAYOUT = make_layout()


def frame(mask):
    integer(mask, 0, FRAMES-1)
    rows = [0]*N
    for a,b in PAIRS:
        f, kind, j = LAYOUT[a,b]
        bit = j if kind == 'fixed' else (mask >> j) & 1
        x,y = (a,b) if bit == f else (b,a)
        rows[x] |= 1 << y
    return tuple(rows)


def key(rows, p):
    return sum(((rows[p[i]] >> p[j]) & 1) << k for k,(i,j) in enumerate(it.combinations(range(len(rows)),2)))


def key_affine(p):
    base, cols = 0, [0]*FREE
    for pos,(i,j) in enumerate(PAIRS):
        u,v = p[i],p[j]
        a,b = sorted((u,v))
        f,kind,k = LAYOUT[a,b]
        value = int((k if kind == 'fixed' else 0) == f) ^ int(u > v)
        base |= value << pos
        if kind == 'free':
            cols[k] |= 1 << pos
    return base,tuple(cols)


def subset_table(cols):
    out = [0]*(1 << len(cols))
    for x in range(1,len(out)):
        b = x & -x
        out[x] = out[x-b] ^ cols[b.bit_length()-1]
    return tuple(out)


@lru_cache(None)
def _compile(text):
    p = word(text)
    base,cols = key_affine(p)
    rb,rc = key_affine(p[::-1])
    constant = base | ((base ^ rb ^ ALL) << WIDTH)
    packed = tuple(c | ((c ^ r) << WIDTH) for c,r in zip(cols,rc))
    zero = one = kept = 0
    impossible = False
    for u,v in zip(p,p[1:]):
        a,b = sorted((u,v))
        f,kind,k = LAYOUT[a,b]
        desired = int(f) ^ int(u > v)
        if kind == 'fixed':
            impossible |= k != desired
        elif desired:
            one |= 1 << k
        else:
            zero |= 1 << k
        kept |= 1 << PAIR_INDEX[a,b]
    impossible |= bool(zero & one)
    return constant,subset_table(packed[:9]),subset_table(packed[9:]),zero,one,impossible,kept


def compile_word(text):
    word(text)  # Validate before the private cache, including after warm calls.
    return _compile(text)


def value(c, mask):
    return c[0] ^ c[1][mask & 511] ^ c[2][mask >> 9]


def legal(c, mask):
    return not c[5] and not(mask & c[3]) and mask & c[4] == c[4]


def audit_ports(mask, ports):
    integer(mask,0,FRAMES-1)
    if type(ports) is not tuple or len(ports) != 7:
        raise ValueError('seven source-bound ports')
    source_code = key(frame(mask), tuple(range(N)))
    keys = set()
    for item in ports:
        if type(item) is not tuple or len(item) != 2:
            raise ValueError('word and full literal target')
        text,target = item
        integer(target,0,ALL)
        c = compile_word(text)
        if not legal(c,mask):
            raise ValueError('non-native Hamilton word')
        v = value(c,mask)
        if not(v >> WIDTH):
            raise ValueError('blue port')
        if target != source_code ^ ALL ^ c[6]:
            raise ValueError('wrong full complement target')
        keys.add(v & ALL)
    if len(keys) != 7:
        raise ValueError('seven distinct full Aut orbits required')
    return 7


def refused(fn,*args):
    try:
        fn(*args)
    except (ValueError,TypeError):
        need(True,'hostile rejected')
    else:
        need(False,'hostile accepted')


def load():
    meta = json.loads((HERE/'inputs_q3.json').read_text(encoding='utf-8'))
    need(meta['schema'] == 'D15-q3-full-key-v1','input schema')
    integer(meta['frames'],FRAMES,FRAMES)
    integer(meta['score_rejected'],91392,91392)
    integer(meta['seven_black'],170752,170752)
    integer(meta['index_bytes'],2652680,2652680)
    integer(meta['expanded_bytes'],26557960,26557960)
    need(type(meta['words']) is list and len(set(meta['words'])) == len(meta['words']),'distinct word bank')
    for text in meta['words']:
        word(text)
    dec = zlib.decompressobj()
    raw = dec.decompress(base64.b85decode(meta['index_base85_zlib'].encode('ascii')), meta['index_bytes']+1)
    need(len(raw) == meta['index_bytes'] and dec.eof and not dec.unused_data and not dec.unconsumed_tail,'complete bounded bank')
    need(raw[:8] == b'D15IDX01','index header')
    need(hashlib.sha256(raw).hexdigest() == meta['index_sha256'],'source-index commitment')
    return meta,raw


def small_controls():
    # Pure-blue and deficit-six controls use independent complete literal words.
    data = []
    for n,S,expected in ((3,{1},(3,1,1)),(5,{1,2},(15,3,3)),(7,{1,2,4},(189,9,3))):
        rows = tuple(sum(1 << j for j in range(n) if (j-i)%n in S) for i in range(n))
        paths = [p for p in it.permutations(range(n)) if all(rows[a] >> b & 1 for a,b in zip(p,p[1:]))]
        allkeys = {key(rows,p) for p in paths}
        bluekeys = {key(rows,p) for p in paths if key(rows,p) ^ key(rows,p[::-1]) == (1 << (n*(n-1)//2))-1}
        need((len(paths),len(allkeys),len(bluekeys)) == expected,'independent small Hamilton/full-key counts')
        data.append({'n':n,'H':len(paths),'h':len(allkeys),'blue':len(bluekeys)})
    # Nonconstant necklaces, with arbitrary fixed permutations of orbit labels.
    cases = 0
    for q in (3,5,7):
        for w in it.product(range(3),repeat=q):
            if len(set(w)) == 1:
                continue
            rotations = [w[s:]+w[:s] for s in range(q)]
            need(len(set(rotations)) == q,'prime nonconstant necklace has q rotations')
            for tau in it.permutations(range(3)):
                count = sum(r == tuple(tau[x] for x in r[::-1]) for r in rotations)
                need(count <= 1,'one fixed anti-label action allows at most one blue rotation')
            cases += 1
    # Composite q hostile: primitive necklace is essential in the general lemma.
    # The shorter-period hostile uses odd composite9 and a repeated3-word.
    w=(0,1,2)*3
    need(len({w[s:]+w[:s] for s in range(9)}) == 3,'composite repeated-necklace hostile')
    # A full singleton shore separates block patterns from cut patterns.
    for w in ((0,1,2,2,2),(0,1,2,3,4)):
        cuts = {w[s:]+w[:s] for s in range(5)}
        cuts = {r*3 for r in cuts}
        blocks = {tuple(c for x in w[s:]+w[:s] for c in (x,)*3) for s in range(5)}
        need(len(cuts) == len(blocks) == 5 and cuts.isdisjoint(blocks),'five source orbit-label block/cut patterns')
    # Literal two-family control: C3 substitution into a strong five-vertex
    # source with score classes1,3,2,2,2 along its marked cycle. This is an
    # actual source, not a relabelled abstract necklace.
    macro = [0]*5
    for a,b in ((0,1),(1,2),(1,3),(1,4),(2,3),(3,4),(4,2),(2,0),(3,0),(4,0)):
        macro[a] |= 1 << b
    rows=[]
    for x in range(15):
        i,j=x%5,x//5
        rows.append(sum(1 << y for y in range(15) if
                        (y%5 != i and macro[i] >> (y%5) & 1) or
                        (y%5 == i and (y//5-j)%3 == 1)))
    rows=tuple(rows)
    cutwords=tuple(tuple((s+i)%15 for i in range(15)) for s in range(5))
    blockwords=[]
    for s in range(5):
        v=s;p=[]
        for j in range(5):
            p.extend((v,(v+5)%15,(v+10)%15))
            v=(p[-1]+1)%15
        blockwords.append(tuple(p))
    family_keys=[];black_counts=[]
    for family in (cutwords,tuple(blockwords)):
        need(all(sorted(p)==list(range(15)) and all(rows[a] >> b & 1 for a,b in zip(p,p[1:])) for p in family),'literal native cut/block construction')
        keys={key(rows,p) for p in family}
        blacks={key(rows,p) for p in family if key(rows,p)^key(rows,p[::-1]) != ALL}
        need(len(keys)==5 and len(blacks)>=4,'actual five distinct orbits and four black ports')
        family_keys.append(keys);black_counts.append(len(blacks))
    need(family_keys[0].isdisjoint(family_keys[1]),'literal full-key family separation')
    return data,cases,{'rows':list(rows),'cut_words':[list(p) for p in cutwords],
                       'block_words':[list(p) for p in blockwords],'black_counts':black_counts}


def run():
    global GATES
    GATES = 0
    meta,raw = load()
    texts = tuple(meta['words'])
    compiled = tuple(compile_word(t) for t in texts)
    basis = (0,)+tuple(1 << j for j in range(FREE))
    for mask in basis:
        rows = frame(mask)
        need(all(rows[i] >> ((i+1)%N) & 1 for i in range(N)),'actual marked Hamilton cycle')
        for text,c in zip(texts,compiled):
            p = word(text)
            k,r = key(rows,p),key(rows,p[::-1])
            need(value(c,mask) == k | ((k ^ r ^ ALL) << WIDTH),'independent literal matrix verifies every affine coefficient')
            need(legal(c,mask) == all(rows[a] >> b & 1 for a,b in zip(p,p[1:])),'basis native-arc conditions')
    base_rows = frame(0)
    rowcols = tuple(tuple(x ^ y for x,y in zip(base_rows,frame(1 << j))) for j in range(FREE))
    source_base,source_cols = key_affine(tuple(range(N)))
    source_lo,source_hi = subset_table(source_cols[:9]),subset_table(source_cols[9:])
    current = list(base_rows)
    expanded = hashlib.sha256(b'D15Q302!')
    at,rejected,accepted = 8,0,0
    uses = [0]*len(texts)
    first = None
    for mask in range(FRAMES):
        if mask:
            changed = mask ^ (mask-1)
            while changed:
                bit=changed & -changed
                col=rowcols[bit.bit_length()-1]
                for i in range(3):
                    current[i] ^= col[i]
                changed -= bit
        # The three representative scores each occur five times by the
        # authenticated R^3 action; no full-Aut assumption is used here.
        scores = sorted(current[i].bit_count() for i in range(3))
        symmetric = scores[1] == 7 and scores[0]+scores[2] == 14
        tag=raw[at];at+=1
        need(tag in (0,1) and bool(tag) == symmetric,'complete source score filter, no skipped symmetric frame')
        expanded.update(bytes([tag]))
        if not tag:
            rejected += 1
            continue
        accepted += 1
        ids = tuple(int.from_bytes(raw[at+2*j:at+2*j+2],'little') for j in range(7))
        at += 14
        need(all(j < len(texts) for j in ids),'literal word-bank index range')
        source_code=source_base ^ source_lo[mask & 511] ^ source_hi[mask >> 9]
        keys=set(); ports=[]
        for idx in ids:
            c=compiled[idx]
            need(legal(c,mask),'actual Hamilton path on this literal source')
            val=value(c,mask)
            need(val >> WIDTH != 0,'full positional-reversal blackness')
            keys.add(val & ALL)
            target=source_code ^ ALL ^ c[6]
            encoded=sum(int(ch,16) << (4*i) for i,ch in enumerate(texts[idx]))
            expanded.update(encoded.to_bytes(8,'little'))
            expanded.update(target.to_bytes(14,'little'))
            ports.append((texts[idx],target));uses[idx]+=1
        need(len(keys)==7,'seven distinct full source-Aut orbit keys')
        if first is None:
            first=mask,tuple(ports)
        if mask % 4093 == 0 or mask == FRAMES-1:
            rows=frame(mask)
            need(scores == sorted(rows[i].bit_count() for i in range(3)),'independent direct-source checkpoint')
            direct=set()
            for text,target in ports:
                p=word(text);k=key(rows,p)
                need(all(rows[a] >> b & 1 for a,b in zip(p,p[1:])),'direct checkpoint Hamilton word')
                need(k ^ key(rows,p[::-1]) != ALL,'direct checkpoint black word')
                kept=sum(1 << PAIR_INDEX[tuple(sorted((a,b)))] for a,b in zip(p,p[1:]))
                need(target == key(rows,tuple(range(N))) ^ ALL ^ kept,'direct target reconstruction')
                direct.add(k)
            need(len(direct)==7,'direct checkpoint distinct full keys')
    need(at==len(raw) and rejected==91392 and accepted==170752,'exact complete frame universe')
    need(all(uses),'every bank word used')
    need(expanded.hexdigest()==meta['expanded_sha256'],'all full-source/word/target bytes match producer commitment')
    need(audit_ports(*first)==7,'public exact source-bound port consumer')
    mask,ports=first
    for bad in (False,0.0,-1,FRAMES):
        refused(frame,bad);refused(audit_ports,bad,ports)
    refused(audit_ports,mask,list(ports))
    refused(audit_ports,mask,(ports[0],)*7)
    refused(audit_ports,mask,((ports[0][0],ports[0][1]^1),)+ports[1:])
    refused(audit_ports,mask,((ports[0][0],float(ports[0][1])),)+ports[1:])
    refused(compile_word,ports[0][0][:-1]+ports[0][0][0])
    refused(compile_word,tuple(range(N)));refused(compile_word,False)
    small,necklaces,block_control=small_controls()
    pins={}
    for rel in ('center_transitivity_20261021c/VT15.md','center_transitivity_20261021c/receipt15.json','deficit_six_centers_20261021b/COMBINED_CUTOFF.md'):
        path=HERE.parent/rel
        pins[rel]=hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    return {'schema':'nontransitive-fifteen-20261022-v1','checks':GATES,
            'q3':{'frames':FRAMES,'score_rejected':rejected,'seven_black':accepted,
                  'word_positions':7*accepted,'distinct_words':len(texts),
                  'index_sha256':meta['index_sha256'],'full_target_stream_sha256':expanded.hexdigest()},
            'prime_necklace_controls':necklaces,'small_full_key_controls':small,
            'native_two_family_control':block_control,
            'inherited_pins':pins,'scope':'q3 complete finite exclusion; q5 non-VT analytic two-family theorem; inherited VT15 exclusion'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run()
    text=json.dumps(result,sort_keys=True,indent=2)+'\n'
    path=HERE/'receipt.json'
    if args.verify:
        need(path.read_bytes().replace(b'\r\n',b'\n')==text.encode(),'saved receipt equality')
    else:
        path.write_text(text,encoding='utf-8')
    print('PASS',GATES,'checks; LF receipt SHA256',hashlib.sha256(text.encode()).hexdigest())


if __name__=='__main__':
    main()
