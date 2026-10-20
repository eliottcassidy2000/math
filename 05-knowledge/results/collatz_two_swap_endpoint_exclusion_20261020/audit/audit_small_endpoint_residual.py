#!/usr/bin/env python3
"""Independent small residual using occupancy sets and a reverse carry scan."""
from hashlib import sha256
from math import gcd
from pathlib import Path
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def main():
    packet = Path(__file__).resolve().parent.parent
    producer = json.loads((packet/'small.receipt.json').read_text())
    need((packet/'small.receipt.json').read_bytes() == (packet/'small.receipt.optimized.json').read_bytes(), 'producer optimized receipt')
    need(sha256((packet/'verify_small_endpoint_residual.py').read_bytes()).hexdigest() == producer['script_sha256'], 'producer source hash')
    digest = sha256()
    rows = 0
    shapes = 0
    hits = []
    minimum_denominator = None
    for p in range(2,100):
        for r in ([1] if p==2 else range(p//2+1,p)):
            if gcd(p,r)!=1:
                continue
            q = p+r
            D = 2**q-3**p
            source = {i for i in range(q) if ((i+1)*p)//q > (i*p)//q}
            ordinary = sorted(source)
            label_positions = [ordinary[-1]]+ordinary[:-1]
            by_residue = {(-j*r)%p:location for j,location in enumerate(label_positions)}
            shapes += 1
            for family,direction,movers in [('plus',1,[p-1-r,p-1]),('minus',-1,[r,0])]:
                occupied = source.copy()
                edges = []
                for residue in movers:
                    old = by_residue[residue]
                    new = (old+direction)%q
                    need(old in occupied and new not in occupied, 'nonactual endpoint move')
                    occupied.remove(old)
                    occupied.add(new)
                    edges.append(old if direction==1 else new)
                word = ''.join('1' if i in occupied else '0' for i in range(q))
                B = 0
                suffix_weight = 1
                for position in range(q-1,-1,-1):
                    if position in occupied:
                        B += (1 << position)*suffix_weight
                        suffix_weight *= 3
                content = gcd(abs(D),B)
                formula = gcd(2*13**p-10**p,3*13**r-2*10**r) if family=='plus' else gcd(9**p-2*5**p,2*9**r-3*5**r)
                need(content==formula,'independent formula agreement')
                if p>=3:
                    need(content<abs(D),'independent nonintegrality')
                    denominator = abs(D)//content
                    minimum_denominator = denominator if minimum_denominator is None else min(denominator,minimum_denominator)
                else:
                    fixed = B//D
                    need(B % D == 0,'boundary integrality')
                    value = fixed
                    for bit in map(int,word):
                        need(value%2==bit,'boundary parity')
                        value = (3*value+1)//2 if bit else value//2
                    need(value==fixed,'boundary closure')
                    hits.append([family,p,r,word,fixed])
                row = {'family':family,'p':p,'r':r,'q':q,'D':D,'B':B,'content':content,'edges':edges,'word':word}
                digest.update(json.dumps(row,sort_keys=True,separators=(',',':')).encode()+b'\n')
                rows += 1
    need(shapes==1502 and rows==3004,'independent universe')
    need(digest.hexdigest()==producer['all_case_trace_sha256'],'every producer row independently matched')
    need(hits==[['plus',2,1,'101',-7],['minus',2,1,'110',-5]],'boundary classification')
    print(json.dumps({'status':'FINITE-EXACT independent complete small residual',
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'producer_script_sha256':producer['script_sha256'],
        'shapes':shapes,'paths':rows,'nonintegral_paths':rows-2,
        'minimum_reduced_denominator_for_p_ge_3':minimum_denominator,
        'boundary_hits':hits,'all_case_trace_sha256':digest.hexdigest(),
        'independence':'No producer implementation imported; occupancy-set moves and reverse carry scan.'},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
