"""Integral C3 coefficient split, its clock, and actual Paley7 inputs."""
from collections import Counter
from itertools import product
from math import gcd
from pathlib import Path
import json

GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def split(c):
    x, y, z = c
    return (x+y+z, x-z, z-y)


def norm(a, b):
    return a*a+a*b+b*b


def main():
    integers = 0
    for c in product(range(-12, 16), repeat=3):
        h, a, b = split(c)
        need((h-a+b) % 3 == 0, 'integral image guard')
        z = (h-a+b)//3
        need((z+a, z-b, z) == c, 'source-preserving inverse')
        need(split(c[1:]+c[:1]) == (h, -a-b, a), 'contragredient order3 clock')
        need(norm(a, b) == sum(x*x for x in c)-
             c[0]*c[1]-c[0]*c[2]-c[1]*c[2], 'Eisenstein norm of quotient coefficients')
        need(h*h+2*norm(a, b) == 3*sum(x*x for x in c), 'norm and trivial-character sum')
        integers += 1
    rows = []
    for m in (1,2,3,4,7,9,11,12,13,27):
        image = Counter(tuple(x % m for x in split(c)) for c in product(range(m), repeat=3))
        g = gcd(m,3)
        need(len(image) == m**3//g, 'complete modular image size')
        need(set(image.values()) == {g}, 'constant modular source fibers')
        for h, a, b in product(range(m), repeat=3):
            need(((h,a,b) in image) == ((h-a+b) % g == 0), 'complete image criterion')
        rows.append({'modulus':m,'source_states':m**3,'image_states':len(image),'fiber_size':g})
    need((1-0+0) % 3 != 0, 'rational split has an unrealizable integer state')
    p = Path(__file__).parent
    native = json.loads((p.parent/'overnight_rank_return_20261029/two_run.json').read_text())
    recovered = []
    for item in native['paley7_affine_convolution']:
        c = item['quotient_totals']
        h, a, b = split(c)
        need(h == item['trivial'] and norm(a,b) == item['cubic_norm'], 'actual native coefficient source')
        recovered.append({'sector':item['sector'],'coefficients':c,'source_sum':h,'phase':[a,b],'norm':norm(a,b)})
    need([x['coefficients'] for x in recovered] == [[15,12,12],[-11,2,2]], 'pinned native source triples')
    out = {'gates':GATES,'integer_sources':integers,'modular_censuses':rows,
           'actual_paley7_sources':recovered,'source_packet_commit':'4ff64c524c',
           'scope':'integer C3 coefficient lattice; not native attachment realizability'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
