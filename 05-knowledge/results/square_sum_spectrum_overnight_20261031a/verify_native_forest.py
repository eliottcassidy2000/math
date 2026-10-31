#!/usr/bin/env python3
"""Universal final-table consumer and concrete expansions, without a solver."""
import argparse
import gzip
import hashlib
import json
from math import isqrt
from pathlib import Path
import struct

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def need(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def square(n):
    return isqrt(n)**2 == n


def check_table(table, pattern):
    expected = {(p,d,a,b) for p in (0,1) for d in range(24,73)
                for a in (-1,1) for b in (-1,1)}
    seen = set()
    for record in table:
        p,d,o0,o1 = (record[k] for k in ('parity','residue','order0','order1'))
        state = p,d,o0,o1
        need(state in expected and state not in seen, 'state duplication or range')
        seen.add(state)
        chosen = {c:int(49+(-1)**(p+1)*c<=d) for c in range(-24,25)}
        endpoints = {}
        for c, i in chosen.items():
            first, last = ((98-c,343+c) if (o0,o1)[i] == 1 else (343+c,98-c))
            terminal = 10 if i == int(d>=57) else 8
            endpoints[c,0] = (3381-c,first)
            endpoints[c,1] = (last,49*terminal+(-1)**(p+i)*c)
        labels = list(range(1,25))+[v for e in endpoints.values() for v in e]
        need(len(labels)==len(set(labels))==220, 'port collision')
        used, singles, external, ends = set(), set(), set(), []
        for atom in record['word']:
            if atom[0] == 'single':
                need(len(atom)==2 and 1<=atom[1]<=24 and atom[1] not in singles,
                     'singleton coverage')
                singles.add(atom[1]); ends.append((atom[1],atom[1]))
            else:
                need(atom[0]=='segment' and len(atom)==5, 'segment schema')
                _,c,i,s,o = atom
                need((c,s) in endpoints and (c,s) not in used, 'segment coverage')
                need(i==chosen[c] and o in (-1,1), 'choice or direction')
                used.add((c,s)); pair = endpoints[c,s]
                ends.append(pair if o==1 else pair[::-1])
        need(used==set(endpoints) and singles==set(range(1,25)), 'missing atom')
        for first, second in zip(ends,ends[1:]+ends[:1]):
            a,b = first[1],second[0]
            need(square(a+b), 'nonsquare connector')
            external.add(frozenset((a,b)))
        for a,b in zip(pattern,pattern[1:]):
            need(frozenset((a,b)) in external, 'lost native path edge')
        need(frozenset((482,3362)) in external, 'lost inherited seed connector')
        need(3362+6242==98**2, 'inherited internal seed sum')
        # The prefix 69,127 is never cut at 2,7, so this port check retains
        # the internal 3362--6242 edge without knowing the rest of the input.
    need(seen==expected, 'incomplete universal table')


def base_pair(n):
    folder = ROOT/'square_sum_B_base_interval_20261019'
    manifest = json.loads((folder/'B_base_chunks.json').read_text())
    chunk = next(c for c in manifest if c['first']<=n<=c['last'])
    zipped = (folder/chunk['file']).read_bytes()
    need(hashlib.sha256(zipped).hexdigest()==chunk['sha256'], 'base integrity')
    raw = gzip.decompress(zipped)
    lo,hi = struct.unpack('<II',raw[8:16])
    need(raw[:8]==b'SSBAS1\0\0' and lo<=n<=hi, 'base range')
    offset = 16+(n-lo)*(n+lo-2)
    result = []
    for length in (n,n+1):
        path = [1]
        for r in raw[offset:offset+length-1]:
            path.append(r*r-path[-1])
        offset += length-1
        need(len(path)==length and sorted(path)==list(range(1,length+1)), 'base vertices')
        result.append(path)
    return result


def transition(pair, family, n):
    m = len(pair[0]); d = n-49*m
    folder = ROOT/'seeded_eventual_cycle_theorem_20261019/certified-pair-transitions'
    record = json.loads((folder/f'{family}_p{m%2}_d{d}.json').read_text())
    result=[]
    for word in record['words']:
        path=[]
        for atom in word:
            if atom[0]=='single': path.append(atom[1])
            elif atom[0]=='head': path.append(49-atom[1])
            else:
                kind,c,i,o=atom
                segment=[49*x+(-1)**(j+1)*c for j,x in enumerate(pair[i])]
                if kind=='tail': segment=segment[1:]
                path.extend(segment if o==1 else segment[::-1])
        result.append(path)
    return result


def B_pair(n):
    if n<=2620: return base_pair(n)
    m=(n-24)//49
    return transition(B_pair(m),'B',n)


def expand(N, table, pattern):
    n=(N-24)//49; d=N-49*n; j=int(d>=57)
    need(n>=2621,'eventual range')
    family='Aplus' if (n+j)%2==0 else 'Aminus'
    pair=transition(B_pair((n-24)//49),family,n)
    orders=tuple(1 if p.index(2)<p.index(7) else -1 for p in pair)
    record=next(r for r in table if (r['parity'],r['residue'],r['order0'],r['order1'])==(n%2,d,*orders))
    cuts=[min(p.index(2),p.index(7))+1 for p in pair]
    path=[]
    for atom in record['word']:
        if atom[0]=='single': path.append(atom[1])
        else:
            _,c,i,s,o=atom
            block=[49*x+(-1)**(k+1)*c for k,x in enumerate(pair[i])]
            segment=block[:cuts[i]] if s==0 else block[cuts[i]:]
            path.extend(segment if o==1 else segment[::-1])
    need(len(path)==N and sorted(path)==list(range(1,N+1)), 'full interval coverage')
    neighbors={x:set() for x in pattern+[482,3362,6242]}
    for a,b in zip(path,path[1:]+path[:1]):
        need(square(a+b), 'full cyclic square edge')
        if a in neighbors: neighbors[a].add(b)
        if b in neighbors: neighbors[b].add(a)
    for a,b in zip(pattern,pattern[1:]):need(b in neighbors[a], 'full native forest')
    need(neighbors[3362]=={482,6242}, 'full inherited seed')
    raw=('\n'.join(map(str,path))+'\n').encode()
    return {'N':N,'state':[n%2,d,*orders], 'cycle_sha256':hashlib.sha256(raw).hexdigest()}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--pattern',choices=('short','n24'),default='n24')
    args=parser.parse_args()
    pattern=[1,8,17] if args.pattern=='short' else [1,8,17,19,6,10,15,21,4,12,24]
    name='native_forest_table.json' if args.pattern=='short' else 'n24_forest_table.json'
    raw=(HERE/name).read_bytes(); table=json.loads(raw)
    check_table(table,pattern)
    expansions=[expand(N,table,pattern) for N in (128453,132063,132128)]
    print(json.dumps({'status':'passed','universal_states':len(table),'checks':checks,
                      'native_path':pattern,'protected_seed':[482,3362,6242],
                      'threshold':128453,'table_sha256':hashlib.sha256(raw).hexdigest(),
                      'expansions':expansions},indent=2,sort_keys=True))


if __name__=='__main__': main()
