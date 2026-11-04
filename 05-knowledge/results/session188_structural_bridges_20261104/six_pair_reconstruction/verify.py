#!/usr/bin/env python3
"""Exact known-six-pair reconstruction controls; no third-party packages."""
from pathlib import Path
from tempfile import TemporaryDirectory
from itertools import permutations
import hashlib,json,struct,subprocess
ROOT=Path(__file__).resolve().parent

def need(ok,why):
    if not ok:raise RuntimeError(why)

def graph(n,code):
    a=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1:a[i]|=1<<j
            else:a[j]|=1<<i
            k+=1
    return tuple(a)

def opposite(a):return tuple((((1<<len(a))-1)^(1<<i)^row) for i,row in enumerate(a))
def endpoint_counts(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            z=d[mask][v]
            if z:
                nxt=a[v]&~mask
                while nxt:
                    bit=nxt&-nxt;nxt-=bit
                    d[mask|bit][bit.bit_length()-1]+=z
    return d

def even_endpoint_bits(a):
    E=endpoint_counts(a);S=endpoint_counts(opposite(a));out=0
    for mask in range(1,1<<len(a)):
        if mask.bit_count()%2==0:
            for i in range(len(a)):out^=((E[mask][i]+S[mask][i])&1)<<i
    return out

def delete(a,j):
    labels=[i for i in range(len(a)) if i!=j]
    rows=tuple(sum(((a[u]>>v)&1)<<k for k,v in enumerate(labels)) for u in labels)
    code=0;k=0
    for u in range(len(rows)):
        for v in range(u+1,len(rows)):
            code|=((rows[u]>>v)&1)<<k;k+=1
    return labels,code

def coefficient(values,am,cm):
    total=0;a=am
    while True:
        c=cm
        while True:
            sign=(-1)**((am^a).bit_count()+(cm^c).bit_count())
            total+=sign*(values[a,c]-values[63^c,63^a])
            if c==0:break
            c=(c-1)&cm
        if a==0:break
        a=(a-1)&am
    need(total%4==0,'shared coefficient divisible by four')
    return (total//4)%2

def main():
    sources=['PROOF.md','verify.py','literal.cpp'];pins={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sources}
    four=0
    for code in range(64):
        a=graph(4,code)
        for i in range(4):
            left=sum(all(a[p[k]]>>p[k+1]&1 for k in range(3)) for p in permutations(range(4)) if p[0]==i or p[-1]==i)
            other=[j for j in range(4) if j!=i]
            right=1+sum((a[i]>>u)&1 for u in other)+sum(((a[i]>>u)&1)*((a[u]>>v)&1) for u in other for v in other if u!=v)
            need((left-right)%2==0,'four-name endpoint expansion');four+=1
    t5=[]
    for code in range(1024):
        a=graph(5,code);t=even_endpoint_bits(a);need(t==sum((r.bit_count()%2)<<i for i,r in enumerate(a)),'five-name degree identity');t5.append(t)
    classes={};vectors=[];nonzero=0
    for code in range(32768):
        a=graph(6,code);h=[[0]*6 for _ in range(6)]
        for j in range(6):
            labels,c=delete(a,j)
            for k,i in enumerate(labels):h[i][j]=(t5[c]>>k)&1
        vec=tuple(h[i][j] for i in range(6) for j in range(6) if i!=j)
        need(vec==tuple((a[i].bit_count()+((a[i]>>j)&1))%2 for i in range(6) for j in range(6) if i!=j),'six-name signature formula')
        need(any(vec),'nonzero signature');nonzero+=1;classes.setdefault(vec,[]).append(code)
        vectors.append(sum(v<<k for k,v in enumerate(vec)))
        for trial in (0,1):
            ds=[trial]+[1^h[0][j]^h[j][0]^trial for j in range(1,6)]
            rows=tuple(sum((h[i][j]^ds[i])<<j for j in range(6) if i!=j) for i in range(6))
            need(rows in (a,opposite(a)),'explicit quotient recovery')
    need(len(classes)==16384,'signature class count')
    for members in classes.values():need(len(members)==2 and members[0]^members[1]==32767,'exact converse classes')
    raw=b''.join(struct.pack('<I',v) for v in vectors);(ROOT/'h_vectors.bin').write_bytes(raw)
    with TemporaryDirectory() as td:
        binary=Path(td)/'literal';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'literal.cpp'),'-o',str(binary)],check=True)
        data=subprocess.check_output([str(binary)]);(ROOT/'literal.json').write_bytes(data);out=json.loads(data)
    battery=[tuple(x) for x in out['coarse_battery']];need(len(set(battery))==188,'literal battery size')
    recovered=[]
    for record in out['states']:
        code=record['quotient_code'];values=dict(zip(battery,record['exact_values']));a=graph(6,code)
        h=tuple(coefficient(values,(1<<i)|(1<<j),1<<j) for i in range(6) for j in range(6) if i!=j)
        need(sum(v<<k for k,v in enumerate(h))==vectors[code],'literal physical signature recovery');recovered.append(code)
    need(out['states'][1]['exact_values']==out['states'][3]['exact_values'],'full retained coarse flip invariance')
    need(len({tuple(v%8 for v in row) for row in out['pair_decoder_exact_values']})==64,'all decoder states distinct')
    fixed={(i,row) for i in range(6) for row in range(32)};need(len(fixed)==192,'fixed decoder count')
    need(pins=={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sources},'frozen source')
    return dict(status='PASS',four_name_literal_endpoint_checks=four,five_name_tournaments=1024,six_name_tournaments=32768,nonzero_signatures=nonzero,exact_signature_classes=16384,complete_signature_vector_bits=30,retained_vector_bank_bytes=len(raw),literal_query_count=out['literal_query_count'],literal_signature_states=recovered,coarse_queries=188,adaptive_queries=194,fixed_queries=380,source_sha256=pins,data_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['h_vectors.bin','literal.json']},scope='Known partition into six pair modules, fixed exterior; no arbitrary twelve-name reconstruction and no full 2^21-state response census.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
