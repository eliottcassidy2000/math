"""Complete order-fourteen normalized source and literal-port consumer."""
from __future__ import annotations
import argparse
import base64
from collections import Counter
import hashlib
import importlib.util
import itertools as it
import json
from pathlib import Path
import zlib

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SPEC=importlib.util.spec_from_file_location('deficit12',HERE.parent/'deficit_twelve_20261021c/check.py')
BASE=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
GATES=0
N=14
FRAMES=2048
SIZE=8+FRAMES*8*19


def need(ok,why):
    global GATES
    GATES+=1
    if not ok:
        raise ValueError(why)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def lf(path):
    return path.read_bytes().replace(b'\r\n',b'\n')


def frame(mask):
    """Shore/difference construction independent of producer edge-orbit sweep."""
    BASE.integer(mask,0,FRAMES-1)
    shore=[]
    for positions in ((0,2,4),(8,9,10)):
        shore.append({j if mask>>bit&1 else 7-j for j,bit in enumerate(positions,1)})
    cross=(1,*(int(mask>>bit&1) for bit in (1,3,5,6,7)),0)
    rows=[0]*N
    for a,b in it.combinations(range(N),2):
        if a%2==b%2:
            forward=(b//2-a//2)%7 in shore[a%2]
        elif a%2==0:
            forward=bool(cross[(b//2-a//2)%7])
        else:
            forward=not cross[(a//2-b//2)%7]
        x,y=(a,b) if forward else (b,a)
        rows[x]|=1<<y
    return tuple(rows)


def audit_record(record):
    if type(record) is not dict or set(record)!={'mask','source','ports'}:
        raise ValueError('exact source record')
    t=BASE.source(record['source'])
    if t!=frame(record['mask']):
        raise ValueError('source substitution')
    return BASE.score_ports(t,record['ports'],8)


def refused(fn,*args):
    try:
        fn(*args)
    except (ValueError,TypeError):
        need(True,'hostile rejected')
    else:
        need(False,'hostile accepted')


def load_bank():
    m=json.loads((HERE/'inputs.json').read_text(encoding='utf-8'))
    if type(m) is not dict or set(m)!={'schema','frames','bytes','sha256','producer_sha256','pins','bank_base85_zlib'}:
        raise ValueError('bank metadata schema')
    BASE.integer(m['frames']);BASE.integer(m['bytes'])
    need((m['schema'],m['frames'],m['bytes'])==('D14SCR01',FRAMES,SIZE),'complete declared universe')
    dec=zlib.decompressobj()
    raw=dec.decompress(base64.b85decode(m['bank_base85_zlib']),SIZE+1)
    need(len(raw)==SIZE and dec.eof and not dec.unused_data and not dec.unconsumed_tail,'complete bounded bank')
    need(raw[:8]==b'D14SCR01' and sha(raw)==m['sha256'],'raw bank commitment')
    need(sha(lf(HERE/'census.cpp'))==m['producer_sha256'],'frozen producer')
    for p,h in m['pins'].items():
        need(sha(lf(ROOT/p))==h,'pinned mathematical/API input '+p)
    return m,raw


def record(mask,raw):
    BASE.integer(mask,0,FRAMES-1)
    ports=[]
    offset=8+mask*8*19
    for j in range(8):
        at=offset+19*j
        word=int.from_bytes(raw[at:at+7],'little')
        p=tuple((word>>(4*k))&15 for k in range(N))
        ports.append((p,int.from_bytes(raw[at+7:at+19],'little')))
    return {'mask':mask,'source':frame(mask),'ports':tuple(ports)}


def verify():
    global GATES
    GATES=0
    meta,raw=load_bank()
    # Exact necessary cut-index parameters; no prime/VT assumption.
    need([(14//m,m) for m in range(1,15,2) if 14%m==0 and 14//m<=8]==[(2,7)],'order14 exact cycle-index reduction')
    edge_orbits=set()
    forced=set()
    for a,b in it.combinations(range(N),2):
        rep=min(tuple(sorted(((a+2*k)%N,(b+2*k)%N))) for k in range(7))
        edge_orbits.add(rep)
        if b==a+1 or (a,b)==(0,N-1):
            forced.add(rep)
    need(len(edge_orbits)==13 and len(forced)==2,'13 edge orbits,2 prescribed cycle orbits')
    rows_hash=hashlib.sha256()
    ports_hash=hashlib.sha256()
    sources=set();counts=Counter();samples=[]
    for mask in range(FRAMES):
        rec=record(mask,raw);t=rec['source'];sources.add(t)
        need(all(BASE.arc(t,i,(i+1)%N) for i in range(N)),'literal Hamilton cycle proves source strong')
        need(all(BASE.arc(t,a,b)==BASE.arc(t,(a+2)%N,(b+2)%N) for a,b in it.combinations(range(N),2)),'literal shift2 automorphism')
        degrees=tuple(x.bit_count() for x in t)
        need(len(set(degrees))==2 and all(degrees[i]==degrees[i%2] for i in range(N)) and degrees[0]+degrees[1]==13,'two distinct complementary score shores')
        counts[','.join(map(str,sorted(set(degrees))))]+=1
        need(audit_record(rec)==8,'eight source-bound black score orbits and exact targets')
        rows_hash.update(b''.join(x.to_bytes(2,'little') for x in t))
        for p,target in rec['ports']:
            need(not BASE.blue(t,p),'independent direct positional-reversal test')
            ports_hash.update(bytes(p));ports_hash.update(target.to_bytes(12,'little'))
        if mask in (0,FRAMES-1):samples.append(rec)
    need(len(sources)==FRAMES,'no duplicate labelled source masks')
    one=samples[0]
    bad=dict(one);bad['mask']=False;refused(audit_record,bad)
    bad=dict(one);bad['source']=frame(1);refused(audit_record,bad)
    bad=dict(one);bad['source']=tuple(float(x) for x in one['source']);refused(audit_record,bad)
    p,target=one['ports'][0]
    for wrong in (p[::-1],p[:-1]+(p[0],),p[:-1]+(-1,),tuple(float(x) for x in p)):
        bad=dict(one);bad['ports']=((wrong,target),)+one['ports'][1:];refused(audit_record,bad)
    bad=dict(one);bad['ports']=((p,target^1),)+one['ports'][1:];refused(audit_record,bad)
    bad=dict(one);bad['ports']=one['ports'][:-1]+(one['ports'][0],);refused(audit_record,bad)
    refused(frame,1.0);refused(frame,FRAMES);refused(frame,-1)
    # Reuse the exact independent n6 positive boundary, not a synthetic graph.
    t=(36,17,10,19,37,14);witnesses={}
    for p in it.permutations(range(6)):
        if not all(BASE.arc(t,a,b) for a,b in zip(p,p[1:])):continue
        s=tuple(t[x].bit_count() for x in p)
        if any(s[i]+s[5-i]!=5 for i in range(6)):
            witnesses.setdefault(s,(p,BASE.target_code(t,p)))
    ports=tuple(witnesses.values())
    need(len(ports)==6 and BASE.score_ports(t,ports,6)==6,'known deficit6 positive boundary')
    refused(BASE.score_ports,t,ports+(ports[0],ports[1]),8)
    return {'schema':'deficit-fourteen-20261021c','checks':GATES,
            'status':'FINITE-EXACT complete necessary normalized universe; no order14 center',
            'frames':FRAMES,'eight_black_orbit_certificates':FRAMES,'survivors':0,
            'literal_word_count':FRAMES*8,'score_shores':dict(sorted(counts.items())),
            'source_rows_sha256':rows_hash.hexdigest(),'word_target_sha256':ports_hash.hexdigest(),
            'raw_bank_sha256':meta['sha256'],'producer_sha256':meta['producer_sha256'],'pins':meta['pins'],
            'remaining_after_orders12_and14':[15,18,20,21,24,25,28,35],
            'scope':'Necessary normalized frames, not ordinary isomorphism classes or an exact deficit histogram.',
            'samples':samples,'six_center_score_certificate_count':6}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.verify and args.write:parser.error('choose immutable verification or write')
    data=(json.dumps(verify(),sort_keys=True,indent=2)+'\n').encode('utf-8')
    if args.write:(HERE/'receipt.json').write_bytes(data)
    if args.verify and data!=lf(HERE/'receipt.json'):raise ValueError('saved receipt mismatch')
    print(f'{GATES:,} exact checks; LF SHA256 {sha(data)}')


if __name__=='__main__':main()
