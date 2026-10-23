"""Independent direct no-search consumer for every vertex-transitive21 source."""
from __future__ import annotations
import argparse
import base64
import hashlib
import itertools as it
import json
from pathlib import Path
import zlib

HERE=Path(__file__).resolve().parent
GROUPS=('C21','C7:3')
N=21
PAIRS=tuple(it.combinations(range(N),2))
INDEX={p:i for i,p in enumerate(PAIRS)}
ALL=(1 << len(PAIRS))-1
GATES=0


def need(ok,why):
    global GATES
    GATES+=1
    if not ok:
        raise ValueError(why)


def integer(x,lo,hi):
    if type(x) is not int or not lo<=x<=hi:
        raise ValueError('literal integer range')
    return x


def group(name):
    if type(name) is not str or name not in GROUPS:
        raise ValueError('marked group')
    return name


def multiply(name,a,b):
    group(name);integer(a,0,20);integer(b,0,20)
    if name=='C21':
        return (a+b)%21
    x,y=a%7,a//7
    u,v=b%7,b//7
    return (x+pow(2,y,7)*u)%7+7*((y+v)%3)


def inverse(name,a):
    group(name);integer(a,0,20)
    if name=='C21':
        return (-a)%21
    x,y=a%7,a//7
    return (-pow(pow(2,y,7),-1,7)*x)%7+7*((-y)%3)


def source(name,mask):
    group(name);integer(mask,0,1023)
    reps=tuple(a for a in range(1,21) if a<inverse(name,a))
    S=tuple(inverse(name,a) if mask >> j & 1 else a for j,a in enumerate(reps))
    # Right multiplication by the selected connection elements, independent
    # of producer's inverse(x)*y membership implementation.
    return tuple(sum(1 << multiply(name,x,s) for s in S) for x in range(21))


def word(text):
    if type(text) is not str or len(text)!=21 or set(text)!=set('ABCDEFGHIJKLMNOPQRSTU'):
        raise ValueError('literal complete source word')
    return tuple(ord(c)-65 for c in text)


def key(rows,p):
    return sum(((rows[p[a]] >> p[b]) & 1) << i for i,(a,b) in enumerate(PAIRS))


def target(rows,p):
    kept=sum(1 << INDEX[tuple(sorted((a,b)))] for a,b in zip(p,p[1:]))
    return key(rows,tuple(range(21))) ^ ALL ^ kept


def audit(name,mask,ports):
    rows=source(name,mask)
    if type(ports) is not tuple or len(ports)!=7:
        raise ValueError('seven exact marked ports')
    keys=set()
    for item in ports:
        if type(item) is not tuple or len(item)!=2:
            raise ValueError('word/full-target pair')
        text,dest=item
        p=word(text);integer(dest,0,ALL)
        if not all(rows[a] >> b & 1 for a,b in zip(p,p[1:])):
            raise ValueError('non-native path')
        k=key(rows,p)
        if k ^ key(rows,p[::-1])==ALL:
            raise ValueError('blue port')
        if dest!=target(rows,p):
            raise ValueError('wrong full complement target')
        keys.add(k)
    if len(keys)!=7:
        raise ValueError('seven distinct full-Aut orbits')
    return 7


def refused(fn,*args):
    try:
        fn(*args)
    except (ValueError,TypeError):
        need(True,'hostile rejected')
    else:
        need(False,'hostile accepted')


def run():
    global GATES
    GATES=0
    meta=json.loads((HERE/'inputs21.json').read_text(encoding='utf-8'))
    need(meta['schema']=='VT21-two-groups-v1','schema')
    integer(meta['frames'],2048,2048);integer(meta['word_bytes'],301056,301056)
    dec=zlib.decompressobj()
    raw=dec.decompress(base64.b85decode(meta['word_base85_zlib']),301057)
    need(len(raw)==301056 and dec.eof and not dec.unused_data and not dec.unconsumed_tail,'complete bounded word stream')
    need(hashlib.sha256(raw).hexdigest()==meta['word_sha256'],'literal word stream hash')
    for name in GROUPS:
        for a in range(21):
            need(multiply(name,a,inverse(name,a))==multiply(name,inverse(name,a),a)==0,'two-sided inverse')
            for b in range(21):
                for c in range(21):
                    need(multiply(name,multiply(name,a,b),c)==multiply(name,a,multiply(name,b,c)),'complete21 group associativity')
    need(multiply('C7:3',7,1)==2+7 and multiply('C7:3',1,7)==1+7,'noncommutative group retained')
    for p,q in ((3,7),(7,3)):
        part=1;n=p-1
        while n%q==0:
            n//=q;part*=q
        need(part<=q,'exact conjugation q-part guard')
    full=hashlib.sha256(b'VT21KEY1')
    at=0;first=None;different_sources=0
    for name in GROUPS:
        for mask in range(1024):
            rows=source(name,mask)
            need(all(row.bit_count()==10 for row in rows),'source regularity')
            need(all(((rows[a] >> b)&1)+((rows[b] >> a)&1)==1 for a,b in PAIRS),'literal tournament')
            keys=set();ports=[]
            for j in range(7):
                text=raw[at:at+21].decode('ascii');at+=21;p=word(text)
                need(all(rows[a] >> b & 1 for a,b in zip(p,p[1:])),'every native Hamilton edge')
                k=key(rows,p)
                need(k ^ key(rows,p[::-1])!=ALL,'every full reversal test is black')
                keys.add(k)
                dest=target(rows,p)
                full.update(text.encode('ascii'));full.update(dest.to_bytes(27,'little'))
                ports.append((text,dest))
            need(len(keys)==7,'seven distinct full ordered keys')
            if first is None:
                first=(name,mask,tuple(ports))
            if name=='C7:3':
                different_sources+=rows!=source('C21',mask)
    need(at==len(raw),'every word position consumed')
    need(full.hexdigest()==meta['full_stream_sha256'],'every literal full target matches discovery commitment')
    need(different_sources>0,'group choice is not discarded')
    need(audit(*first)==7,'public source-bound exact port consumer')
    name,mask,ports=first
    for bad in (False,0.0,-1,1024):
        refused(source,name,bad);refused(audit,name,bad,ports)
    for bad in (False,0,'C3xC7','C7:7'):
        refused(source,bad,0)
    refused(audit,name,mask,list(ports));refused(audit,name,mask,(ports[0],)*7)
    refused(audit,name,mask,((ports[0][0],ports[0][1]^1),)+ports[1:])
    refused(audit,name,mask,((ports[0][0],False),)+ports[1:])
    refused(word,ports[0][0][:-1]+ports[0][0][0]);refused(word,tuple(range(21)))
    # Changing the source orientation invalidates the original native words.
    refused(audit,name,mask^1023,ports)
    return {'schema':'VT21-complete-two-groups-v1','checks':GATES,'groups':list(GROUPS),
            'frames_per_group':1024,'word_positions':14336,
            'full_target_stream_sha256':full.hexdigest(),'word_sha256':meta['word_sha256'],
            'same_mask_distinct_group_sources':different_sources,
            'scope':'Every VT21 tournament has a regular C21 or C7:3 subgroup under cited odd-order solvability; all2048 necessary frames excluded.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();text=json.dumps(result,sort_keys=True,indent=2)+'\n';path=HERE/'receipt21.json'
    if args.verify:
        need(path.read_bytes().replace(b'\r\n',b'\n')==text.encode(),'saved equality')
    else:
        path.write_text(text,encoding='utf-8')
    print('PASS',GATES,'checks; LF receipt SHA256',hashlib.sha256(text.encode()).hexdigest())


if __name__=='__main__':
    main()
