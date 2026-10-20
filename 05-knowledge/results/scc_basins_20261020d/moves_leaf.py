"""Exact mixed-leaf and two uniform endpoint-hostile detour controls."""
from itertools import permutations
import json
from pathlib import Path
from moves import (need, CASES, R4, induced, join, is_path, blue,
                   complement, blue_tilings)

LEAF_CASES = [
    ((2,4,9,3), [(0,2,3,4,1,5),(4,1,5,0,2,3),(4,5,0,1,2,3)]),
    ((6,12,8,1), [(0,2,3,4,1,5),(4,1,5,0,2,3),(4,5,0,1,2,3)]),
    ((10,4,9,2), [(0,1,4,2,3,5),(1,4,2,3,5,0),(1,2,3,4,5,0)]),
]


def verify(source, paths):
    original = tuple(range(len(source)))
    need(blue(source,original), 'original presentation must be blue')
    state = source
    for p in paths:
        need(is_path(state,p), 'detour word is not a Hamiltonian path')
        need(not blue(state,p), 'detour word is blue')
        state = complement(state,p)
    need(state == complement(source,original), 'literal target mismatch')


def rotate(p, first):
    k = p.index(first)
    return p[k:]+p[:k]


def main():
    pure_four = tuple(sum(1 << v for v in range(u+1,4)) for u in range(4))
    blue_four = [b for b,_ in CASES if blue(b,tuple(range(4)))]
    need(blue_four == [b for b,_ in LEAF_CASES], 'three blue R4 presentations differ')
    empty_targets = []
    for b in blue_four:
        target = complement(b,tuple(range(4)))
        if any(induced(target,p) == R4 for p in permutations(range(4))):
            empty_targets.append('R4')
        elif any(induced(target,p) == pure_four for p in permutations(range(4))):
            empty_targets.append('transitive')
        else:
            raise RuntimeError('unexpected blue R4 neighbor')
    need(empty_targets.count('R4') == 2 and empty_targets.count('transitive') == 1,
         'empty-word blue loop/leaf boundary differs')
    for b,paths in LEAF_CASES:
        verify(join([(0,),b,(0,)]),paths)
    cyclic_records = []
    for m in range(2,8):
        count = 0
        for b in blue_tilings(m):
            c = list(join([(0,),b,(0,)]))
            c[0] ^= 1 << (m+1)
            c[m+1] ^= 1
            source = join([(0,),tuple(c),(0,)])
            a,v,w,z = 0,1,m+2,m+3
            block = tuple(range(2,m+2))
            ps = [(a,w,v)+block+(z,),
                  block+(z,a,w,v), block+(w,z,a,v)]
            verify(source,ps)
            count += 1
        cyclic_records.append(dict(middle_order=m,blue_presentations=count))
    for h in range(4,101):
        c = tuple(sum(1 << v for v in range(h) if v == u+1 or v < u-1)
                  for u in range(h))
        source = join([(0,),c,(0,)])
        p0 = tuple(range(h+2))
        p1 = (0,)+tuple(range(2,h+1))+(1,h+1)
        verify(source,[p1,rotate(p1,h),rotate(p0,h)])
    report = dict(singleton_wrapped_R4=LEAF_CASES,
                  empty_wrapper_blue_halfedge_targets=empty_targets,
                  cyclic_blowup_controls=cyclic_records,
                  reverse_transitive_path_orders=list(range(4,101)))
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: all three singleton-wrapped R4 blue presentations; empty wrapper;')
    print('all cyclic-blowup blue middles through order seven; reverse-path orders 4..100.')


if __name__ == '__main__':
    main()
