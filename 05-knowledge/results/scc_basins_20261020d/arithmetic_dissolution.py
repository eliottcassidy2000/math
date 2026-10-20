"""Exact contextual R4 dissolution controls; no assertions disabled by -O.

Universe: every standard-Hamiltonian-path presentation of A of orders1..6.
No isomorphism pruning or sampling is used. The all-order proof is in
arithmetic.md; this finite check independently executes its literal paths.
"""
from itertools import combinations, permutations
from pathlib import Path
from hashlib import sha256
import json

HERE = Path(__file__).resolve().parent
B = (6, 4, 8, 3)
D = (10, 0, 3, 6)
FIRST = (3, 0, 1, 2)
SECOND = (3, 2, 0, 1)
# Paired-R4 theorem, presentation (6,12,8,1), moves.md section2.
DETOUR = [('A',1,2,3,0), (3,'A',1,2,0), (0,3,'A',1,2),
          (0,2,3,'A',1), (2,3,'A',0,1)]


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def arc(t, u, v):
    return bool(t[u] & (1 << v))


def paths(t):
    n = len(t)
    def visit(p, unused):
        if not unused:
            yield tuple(p)
        for v in sorted(unused):
            if not p or arc(t, p[-1], v):
                yield from visit(p + [v], unused - {v})
    yield from visit([], set(range(n)))


def is_path(t, p):
    return sorted(p) == list(range(len(t))) and all(
        arc(t,u,v) for u,v in zip(p,p[1:]))


def blue(t, p):
    n = len(p)
    return all(arc(t,p[i],p[j]) == arc(t,p[n-1-j],p[n-1-i])
               for i,j in combinations(range(n),2))


def k_move(t, p):
    need(is_path(t,p), 'K requires an actual Hamiltonian path')
    retained = {frozenset((u,v)) for u,v in zip(p,p[1:])}
    out = list(t)
    for u,v in combinations(range(len(t)),2):
        if frozenset((u,v)) not in retained:
            out[u] ^= 1 << v
            out[v] ^= 1 << u
    return tuple(out)


def join(a,b):
    m,n = len(a),len(b)
    right = ((1 << n)-1) << m
    return tuple(x | right for x in a) + tuple(x << m for x in b)


def presentations(n):
    free = [(u,v) for u,v in combinations(range(n),2) if v != u+1]
    for mask in range(1 << len(free)):
        t = [0] * n
        for u in range(n-1):
            t[u] |= 1 << (u+1)
        for i,(u,v) in enumerate(free):
            if mask >> i & 1:
                t[u] |= 1 << v
            else:
                t[v] |= 1 << u
        yield tuple(t)


def literal(a, p):
    m = len(a)
    source = join(a,B)
    p1 = p + tuple(m+i for i in FIRST)
    p2 = p + tuple(m+i for i in SECOND)
    middle = k_move(source,p1)
    target = k_move(middle,p2)
    need(target == join(a,D), 'two-K literal target mismatch')
    need(blue(middle,p2) == blue(target,p2), 'K blue invariance failed')
    return source,middle,target,p1,p2


def route(a):
    m = len(a)
    target = join(a,D)
    # This is independent of the macro: exhaust actual target HPs.
    target_words = list(paths(target))
    active = any(not blue(target,p) for p in target_words)
    a_words = list(paths(a))
    candidates = []
    for p in a_words:
        source,middle,dest,p1,p2 = literal(a,p)
        if not blue(middle,p2):
            candidates.append((p,source,middle,dest,p1,p2))
    need(bool(candidates) == active, 'fixed suffix lost black target HPs')
    if not active:
        return dict(active=False, A=list(a), target_paths=len(target_words))
    p,source,middle,dest,p1,p2 = candidates[0]
    first_blue = blue(source,p1)
    if first_blue:
        steps = []
        for word in DETOUR:
            normalized = tuple(v for token in word for v in
                               (range(m) if token == 'A' else [m+token]))
            steps.append(tuple(p1[i] for i in normalized))
        state = source
        for step in steps:
            need(not blue(state,step), 'paired-R4 repair is blue')
            state = k_move(state,step)
        need(state == middle, 'paired-R4 repair changed the intermediate')
    else:
        steps = [p1]
    steps.append(p2)
    state = source
    for step in steps:
        need(is_path(state,step), 'route has an absent path arc')
        need(not blue(state,step), 'route contains a blue step')
        state = k_move(state,step)
    need(state == target, 'route did not restore the complete prefix')
    standard = tuple(range(m))
    standard_second_blue = blue(literal(a,standard)[1],
                               literal(a,standard)[4])
    return dict(active=True, steps=len(steps), first_blue=first_blue,
                switched_path=(p != standard),
                standard_second_blue=standard_second_blue,
                A=list(a), prefix_path=list(p), route=[list(x) for x in steps])


def main():
    need(k_move(k_move(B,FIRST),SECOND) == D, 'empty-prefix identity')
    # Literal two-port identity is the formal ledger for arbitrary A.
    literal((2,0),(0,1))
    rows,witnesses = [],{}
    for n in range(1,7):
        stats = dict(order=n,presentations=0,active=0,pure_blue_target=0,
                     two_black=0,six_black=0,path_selection_needed=0)
        for a in presentations(n):
            r = route(a)
            stats['presentations'] += 1
            if not r['active']:
                stats['pure_blue_target'] += 1
                witnesses.setdefault('pure_blue_target',r)
                continue
            stats['active'] += 1
            stats['two_black' if r['steps']==2 else 'six_black'] += 1
            if r['standard_second_blue']:
                stats['path_selection_needed'] += 1
                witnesses.setdefault('path_selection_needed',r)
            if r['first_blue']:
                witnesses.setdefault('six_black_repair',r)
        need(stats['presentations'] == 2**((n-1)*(n-2)//2),
             'presentation universe incomplete')
        rows.append(stats)
    # A first active blue-standard target needs a mixed center beyond the
    # six-vertex prefix census. This eight-vertex prefix is explicit.
    padded = join(join((0,),(2,4,1)),(6,12,8,1))
    extra = route(padded)
    need(extra['active'] and extra['standard_second_blue'] and
         extra['switched_path'], 'chosen-path hostile did not exercise repair')
    witnesses['path_selection_needed'] = extra
    need(all(k in witnesses for k in ('pure_blue_target','path_selection_needed',
                                      'six_black_repair')),
         'a hostile or positive branch was never exercised')
    report = dict(scope='all standard-HP presentations of nonempty A through order6',
                  literal_B=list(B),literal_target_D=list(D),
                  first_suffix=list(FIRST),second_suffix=list(SECOND),
                  census=rows,extra_prefix_orders=[8],branch_witnesses=witnesses,
                  source_sha256_lf=sha256(Path(__file__).read_text().encode()).hexdigest())
    result = json.dumps(report,indent=2)+'\n'
    (HERE/'arithmetic_dissolution.json').write_text(result,encoding='utf-8')
    print(result,end='')


if __name__ == '__main__':
    main()
