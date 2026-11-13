"""Exact neutral chirality/marked-response controls; no old broad census."""
from itertools import combinations, permutations, product
from pathlib import Path
from math import prod
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def tt(n):
    return tuple(((1 << n)-1) ^ ((1 << (i+1))-1) for i in range(n))


def cyclic(blocks):
    offsets = (0, len(blocks[0]), len(blocks[0])+len(blocks[1]))
    return tuple((row << offsets[i]) |
                 (((1 << len(blocks[(i+1)%3]))-1) << offsets[(i+1)%3])
                 for i, block in enumerate(blocks) for row in block)


def reverse(rows):
    return tuple(((1 << len(rows))-1) ^ (1 << i) ^ row
                 for i, row in enumerate(rows))


def inventory(rows):
    n = len(rows)
    ends = [[0]*n for _ in range(1 << n)]
    hs = [1]+[0]*((1 << n)-1)
    for mask in range(1, 1 << n):
        for v in range(n):
            if mask >> v & 1:
                rest = mask ^ (1 << v)
                ends[mask][v] = (1 if rest == 0 else
                    sum(ends[rest][u] for u in range(n) if rows[u] >> v & 1))
        hs[mask] = sum(ends[mask])
    return hs, ends


def profile(rows, root=None):
    hs, ends = inventory(rows)
    full = (1 << len(rows))-1
    _, starts = inventory(reverse(rows))
    cover = [1]+[0]*full
    cs = [0]
    ls = [0]
    for k in range(1, len(rows)+1):
        nxt = [0]*(full+1)
        for mask in range(1, full+1):
            sub = mask
            while sub:
                nxt[mask] += hs[sub]*cover[mask ^ sub]
                sub = (sub-1) & mask
        cs.append(nxt[full])
        if root is not None:
            ls.append(sum(starts[sub][root]*cover[full ^ sub]
                          for sub in range(1, full+1) if sub >> root & 1))
        cover = nxt
    state = None if root is None else (
        hs[full], starts[full][root], hs[full ^ (1 << root)],
        sum(ends[sub][root]*hs[full ^ sub] for sub in range(1, full+1)))
    return cs, ls, state, hs


def value(v, i):
    return v[i] if i < len(v) else 0


def triangles(rows):
    return tuple(t for t in combinations(range(len(rows)), 3)
                 if all(sum((rows[v] >> u)&1 for u in t if u != v) == 1
                        for v in t))


def flip_quotient(rows, sizes):
    owner = [i for i, size in enumerate(sizes) for _ in range(size)]
    return tuple(row ^ sum(1 << j for j in range(len(rows)) if owner[i] != owner[j])
                 for i, row in enumerate(rows))


def tournament(n, code):
    rows = [0]*n
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (i, j) if code >> bit & 1 else (j, i)
        rows[u] |= 1 << v
    return tuple(rows)


def physical_coefficients(positive, observer):
    """All maximal-run words on exterior r,a,b and an arbitrary M.

    A=TT2=(r,a). Positive orientation is A->M->b->A;
    negative is A->b->M->A. For F use two labelled paths,
    first nonempty and ending r, second possibly empty.
    """
    exterior = ('a','b') if observer == 'q' else ('r','a','b')
    slots = 2 if observer == 'F' else 1
    def edge(x, y):
        if x == y:
            return False  # consecutive M symbols would be one maximal run
        if x in ('r','a') and y in ('r','a'):
            return (x,y) == ('r','a')
        owner = {'r':0, 'a':0, 'M':1 if positive else 2,
                 'b':2 if positive else 1}
        return (owner[x]+1)%3 == owner[y]
    coefficients = []
    for k in range(1, len(exterior)+slots+1):
        total = 0
        for positions in combinations(range(k+len(exterior)), len(exterior)):
            for perm in permutations(exterior):
                word = ['M']*(k+len(exterior))
                for pos, label in zip(positions, perm):
                    word[pos] = label
                cuts = range(1,len(word)+1) if slots == 2 else (len(word),)
                for cut in cuts:
                    paths = (word[:cut],word[cut:]) if slots == 2 else (word,)
                    if observer == 'F' and paths[0][-1] != 'r':
                        continue
                    if observer == 's' and paths[0][0] != 'r':
                        continue
                    if all(all(edge(x,y) for x,y in zip(path,path[1:])) for path in paths):
                        total += 1
        coefficients.append(total)
    return coefficients


def main():
    bank = [tournament(n, code) for n in range(1,4)
            for code in range(1 << (n*(n-1)//2))]
    controls = []
    # All labelled B,C through3 and selected rooted A, with full subset DP.
    rooted = [(tt(1),0),(tt(2),0),(tt(2),1),((2,4,1),0)]
    for (A,r), B, C in product(rooted, bank, bank):
        cb,_,_,_ = profile(B)
        cc,_,_,_ = profile(C)
        _,la,_,_ = profile(A,r)
        plus = cyclic((A,B,C))
        minus = flip_quotient(plus, tuple(map(len,(A,B,C))))
        _,_,sp,hp = profile(plus,r)
        _,_,sm,hm = profile(minus,r)
        predicted = sum(value(la,k+1)*(value(cb,k+1)*value(cc,k)
                         -value(cb,k)*value(cc,k+1)) for k in range(len(A)))
        require(sp[1]-sm[1] == predicted, 'marked determinant mismatch')
        require(hp == hm, 'labelled induced-count neutrality mismatch')
        require(triangles(plus) == triangles(minus), 'cyclic triple changed')
        controls.append({'orders':list(map(len,(A,B,C))),'root':r,
                         'start_difference':predicted})
    coefficients = {str(sign):{key:physical_coefficients(sign,key)
                    for key in ('H','s','q','F')} for sign in (True,False)}
    require(coefficients['True'] == {'H':[5,3,0,0], 's':[2,1,0,0],
                                     'q':[3,1,0], 'F':[6,3,0,0,0]}, 'positive tensor')
    require(coefficients['False'] == {'H':[5,3,0,0], 's':[2,0,0,0],
                                      'q':[3,1,0], 'F':[6,7,2,0,0]}, 'negative tensor')
    special = []
    for M in bank + [tt(4),(2,4,9,3)]:
        cm,_,_,hs = profile(M)
        h = cm[1]
        p = 2*h+value(cm,2)
        w = 3*h+3*value(cm,2)+value(cm,3)
        plus = cyclic((tt(2),M,tt(1)))
        minus = flip_quotient(plus,(2,len(M),1))
        sp = profile(plus,0)[2]
        sm = profile(minus,0)[2]
        require(sp == (3*p-h,p,p+h,3*p), 'positive response family')
        require(sm == (3*p-h,2*h,p+h,2*w+p-2*h), 'negative response family')
        special.append({'M':M,'h_p_w':[h,p,w],'positive':sp,'negative':sm})
    # A smallest nontrivial old-port contrast; independent permutation path count.
    graph = cyclic((tt(2),tt(2),tt(1)))
    reflected = flip_quotient(graph,(2,2,1))
    perm_counts=[]
    for g in (graph,reflected):
        paths=[p for p in permutations(range(5)) if all(g[x]>>y&1 for x,y in zip(p,p[1:]))]
        perm_counts.append([len(paths),sum(p[0]==0 for p in paths)])
    require(perm_counts == [[11,4],[11,2]], 'permutation control')
    # Depth-two transitive-leaf hierarchy: degrees and starts cannot be one.
    inner = cyclic((tt(1),tt(1),tt(1)))
    tower = cyclic((inner,inner,inner))
    _,_,_,hs = profile(tower)
    starts = inventory(reverse(tower))[1][-1]
    require(min(starts)>=2 and min(8-row.bit_count() for row in tower)>=2,
            'deep hierarchy boundary')
    result={'status':'PASS','universe':'4 rooted A types times all 11 labelled B and C of orders1..3',
            'determinant_controls':len(controls),'physical_run_coefficients':coefficients,
            'special_families':special,'independent_permutation_controls':perm_counts,
            'depth_two_control':{'order':9,'H':hs[-1],'starts':starts},
            'scope':'new chirality formulas and controls; inherited full-Q neutrality is credited; no arbitrary-H selector',
            'reproduce':'python3 -B 05-knowledge/results/zeckendorf_target_selection_20261113/chirality.py'}
    Path(__file__).with_name('chirality.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(controls),'rooted determinant controls; complete run tensors and hierarchy boundary')


if __name__ == '__main__':
    main()
