#!/usr/bin/env python3
"""Compose the incoming recursive B pairs with protected A and native gluing.

The universal argument uses the frozen transition contracts. This consumer
audits the whole cut interface, then independently stream-checks actual cycles.
"""
import hashlib
import json
from math import isqrt
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent
INCOMING=RESULTS/'session124_recursive_square_prefix_20261031'
OLD=RESULTS/'seeded_eventual_cycle_theorem_20261019'
NATIVE=RESULTS/'square_sum_spectrum_overnight_20261031a'
sys.path.insert(0,str(INCOMING))
from core import Core, prefix, threshold

P=(1,8,17,19,6,10,15,21,4,12,24)
S=[t for t in range(-1200,1201) if not(-24<=t<=24 or 74<=t<=122)]
checks=0


def need(x,m):
    global checks
    checks+=1
    if not x:raise RuntimeError(m)


def edge(a,b):return (a,b) if a<b else (b,a)


def table_contracts():
    paired_checks=0
    for parity in (0,1):
        for d in range(24,73):
            file=INCOMING/'paired'/f'B_p{parity}_d{d}.json'
            args=[sys.executable]+(['-O'] if sys.flags.optimize else [])
            receipt=json.loads(subprocess.check_output(args+[str(INCOMING/'verify_paired.py'),str(file)]))
            need(receipt['status']=='passed','incoming paired replay')
            paired_checks+=receipt['checks']
    cuts=set();states=0;retained_blocks=0
    for family in ('Aplus','Aminus'):
        for parity in (0,1):
            for d in range(24,73):
                a=json.loads((OLD/'certified-pair-transitions'/f'{family}_p{parity}_d{d}.json').read_text())
                need((a['family'],a['parity'],a['residue'])==(family,parity,d),'old A state')
                declared=set(a['split_prefix_blocks']);cuts|=declared
                need(declared<={0,-2},'outside uniform cut exclusion')
                for w in a['words']:
                    for c in range(-24,25):
                        if c in (0,-2):continue
                        occurrences=[z for z in w if z[0] in ('block','tail') and z[1]==c]
                        need(len(occurrences)==1 and occurrences[0][0]=='block','lost intact offset')
                        retained_blocks+=1
                states+=1
    need(cuts=={0,-2} and states==196,'full A interface census')
    need(len(S)==2303,'ray width')
    intervals={49*u+c for u in range(-24,25) if u not in (0,2) for c in range(-24,25)}
    need(intervals==set(S),'signed offset coverage')
    # No new final table: re-use the fully audited native table and all its
    # 2--7 cuts. Re-run its universal endpoint/ownership consumer.
    sys.path.insert(0,str(NATIVE))
    from verify_native_forest import check_table
    table=json.loads((NATIVE/'n24_forest_table.json').read_text())
    check_table(table,list(P))
    return table,{'incoming_pair_states':98,'incoming_paired_checks':paired_checks,
                  'old_A_states':states,'unaffected_block_occurrences':retained_blocks,
                  'excluded_offsets':sorted(cuts),'native_final_states':len(table)}


def old_A(pair,family,n):
    m=len(pair[0]);d=n-49*m
    record=json.loads((OLD/'certified-pair-transitions'/f'{family}_p{m%2}_d{d}.json').read_text())
    outputs=[]
    for word in record['words']:
        out=[]
        for a in word:
            if a[0]=='single':out.append(a[1])
            elif a[0]=='head':out.append(49-a[1])
            else:
                kind,c,i,o=a
                block=[49*v+(-1)**(j+1)*c for j,v in enumerate(pair[i])]
                if kind=='tail':block=block[1:]
                out.extend(block if o==1 else block[::-1])
        outputs.append(out)
    colors={x:j%2 for j,x in enumerate(outputs[0])}
    for i,out in enumerate(outputs):
        need(len(out)==n+i and sorted(out)==list(range(1,n+i+1)),'A coverage')
        need(out[:2]==[69,127],'A prefix')
        need(abs(out.index(2)-out.index(7))==1 and out.index(2)%2==0,'A cut/phase')
        for x,y in zip(out,out[1:]):need(isqrt(x+y)**2==x+y,'A square')
        if i:need(all(colors[x]==j%2 for j,x in enumerate(out) if x<=n),'A common parity')
    return outputs


def expand(N,l,table):
    need(N>=threshold(l+2),'joint threshold')
    n=(N-24)//49;d=N-49*n;m=(n-24)//49
    need(m>=threshold(l),'recursive B source threshold')
    builder=Core();pair=builder.build_B(m);small=prefix(l)
    need(all(p[:len(small)]==small for p in pair),'incoming unscaled prefix')
    j=int(d>=57);family='Aplus' if (n+j)%2==0 else 'Aminus'
    apair=old_A(pair,family,n)
    order=tuple(1 if p.index(2)<p.index(7) else -1 for p in apair)
    cuts=[min(p.index(2),p.index(7))+1 for p in apair]
    record=next(r for r in table if (r['parity'],r['residue'],r['order0'],r['order1'])==(n%2,d,*order))
    required={edge(a,b) for a,b in zip(P,P[1:])}|{edge(482,3362),edge(3362,6242)}
    forest_vertices=set()
    for t in S:
        ray=[2401*x+(-1)**j*t for j,x in enumerate(small)]
        need(not forest_vertices.intersection(ray),'rays are disjoint')
        forest_vertices.update(ray)
        required.update(edge(a,b) for a,b in zip(ray,ray[1:]))
    need(not forest_vertices.intersection(P),'native path disjoint from ray forest')
    missing=set(required);seen=bytearray(N+1);last=first=None;count=0;digest=hashlib.sha256()
    def vertices():
        for atom in record['word']:
            if atom[0]=='single':yield atom[1]
            else:
                _,c,i,part,o=atom;p=apair[i]
                lo,hi=(0,cuts[i]) if part==0 else (cuts[i],len(p))
                for k in (range(lo,hi) if o==1 else range(hi-1,lo-1,-1)):
                    yield 49*p[k]+(-1)**(k+1)*c
    for v in vertices():
        need(1<=v<=N and not seen[v],'stream vertex ownership')
        seen[v]=1;count+=1;digest.update((str(v)+'\n').encode())
        if first is None:first=v
        if last is not None:
            need(isqrt(last+v)**2==last+v,'stream square')
            missing.discard(edge(last,v))
        last=v
    need(count==N and sum(seen)==N,'stream full interval')
    need(isqrt(first+last)**2==first+last,'closing square')
    missing.discard(edge(first,last))
    need(not missing,'all prescribed paths retained')
    return {'N':N,'level':l,'threshold':threshold(l+2),'source_index':m,
            'ray_count':len(S),'vertices_per_ray':len(small),'prescribed_edges':len(required),
            'native_path':P,'seed':[482,3362,6242],
            'decimal_lines_sha256':digest.hexdigest(),'incoming_core_checks':builder.checks}


def main():
    table,interfaces=table_contracts()
    controls=[expand(128453,0,table),expand(128454,0,table),expand(6294221,1,table)]
    print(json.dumps({'status':'passed','interfaces':interfaces,'checks':checks,'controls':controls,
                      'scope':'All levels by the inherited recursive prefix and intact-segment proof; controls expand the three displayed cycles.'},indent=2,sort_keys=True))


if __name__=='__main__':main()
