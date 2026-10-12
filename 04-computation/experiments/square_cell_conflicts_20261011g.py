"""Exact half-open cell conflicts; optimal eight colours in the declared period.

The theorem concerns whole-cell colours with period 8a in both coordinates.
It does not cover diagonal cell splitting or a larger colour period.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product, combinations
from math import isqrt

CHECKS = 0


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def positive_rational(x):
    need(type(x) in (int, F) and x > 0, 'positive exact squared-distance threshold')
    return F(x)


@dataclass(frozen=True)
class Cell:
    x: int
    y: int
    width: int
    height: int


def audit_cell(cell):
    need(type(cell) is Cell and all(type(v) is int for v in
         (cell.x, cell.y, cell.width, cell.height)), 'exact integral cell coordinates')
    need(cell.width > 0 and cell.height > 0, 'positive cell sides')
    return cell


def translate(cell, dx, dy):
    audit_cell(cell)
    need(type(dx) is int and type(dy) is int, 'exact translation')
    return Cell(cell.x + dx, cell.y + dy, cell.width, cell.height)


def distance_window(left, right):
    """Open squared-distance range endpoints; zero may also be attained."""
    audit_cell(left)
    audit_cell(right)
    intervals = ((right.x-left.x-left.width, right.x+right.width-left.x),
                 (right.y-left.y-left.height, right.y+right.height-left.y))
    low = sum(max(lo, -hi, 0)**2 for lo, hi in intervals)
    high = sum(max(lo*lo, hi*hi) for lo, hi in intervals)
    return low, high


def conflicts(left, right, threshold):
    threshold = positive_rational(threshold)
    lo, hi = distance_window(left, right)
    return lo < threshold < hi


def motif():
    cells = []
    for i, j in product(range(2), repeat=2):
        width = 2 if (i+j) % 2 == 0 else 1
        for u, v in product(range(0, 4, width), repeat=2):
            cells.append(Cell(4*i+u, 4*j+v, width, width))
    return tuple(cells)


CELLS = motif()
SHIFTS = tuple(product(range(-1, 2), repeat=2))
WRAPPED = tuple(tuple(tuple(distance_window(a, translate(b, 8*p, 8*q))
                            for p, q in SHIFTS) for b in CELLS) for a in CELLS)
CLIQUES = ((0,1,2,3,4,6,12,14), (0,1,2,3,4,7,13,19),
           (0,1,2,3,5,17,24,27), (0,1,2,3,5,18,27,28),
           (0,1,2,3,4,6,13,15))
SCALE_COVER = ((0,8,10,True,False), (1,10,13,True,False),
               (2,13,17,True,False), (2,17,18,True,False),
               (3,18,25,True,True), (4,25,26,False,True),
               (0,26,36,False,True))


def compact_graph(threshold):
    """Complete infinite-period conflict graph on the compact scale window."""
    threshold = positive_rational(threshold)
    need(8 <= threshold <= 36, 'declared compact scale window')
    edges = set()
    loops = set()
    for i in range(40):
        for j in range(i, 40):
            if any(lo < threshold < hi for lo, hi in WRAPPED[i][j]):
                (loops if i == j else edges).add(i if i == j else (i,j))
    return frozenset(edges), frozenset(loops)


def clique_certificate(threshold):
    threshold = positive_rational(threshold)
    need(8 <= threshold <= 36, 'declared compact scale window')
    index = (0 if threshold < 10 or threshold > 26 else
             1 if threshold < 13 else 2 if threshold < 18 else
             3 if threshold <= 25 else 4)
    clique = CLIQUES[index]
    verify_clique(threshold, clique)
    return clique


def verify_clique(threshold, clique):
    threshold = positive_rational(threshold)
    need(type(clique) is tuple and len(clique) == 8 and
         all(type(v) is int and 0 <= v < 40 for v in clique) and len(set(clique)) == 8,
         'eight exact distinct motif-cell identifiers')
    edges, loops = compact_graph(threshold)
    need(not loops and all(tuple(sorted((i,j))) in edges for i,j in combinations(clique,2)),
         'all periodic pair conflicts retained')
    return True


def lifted_clique_edges(threshold, clique):
    """Each edge can use different actual period copies; retain every such lift."""
    verify_clique(threshold,clique)
    result=[]
    for i,j in combinations(clique,2):
        shift=next((p,q) for p,q in SHIFTS if
                   conflicts(CELLS[i],translate(CELLS[j],8*p,8*q),threshold))
        result.append((i,j,shift))
    return tuple(result)


def interval_clique_packets():
    """Fixed period-copy choice for every pair throughout each exact scale interval."""
    packets=[]
    for index,l,h,lc,hc in SCALE_COVER:
        pairs=[]
        for i,j in combinations(CLIQUES[index],2):
            hit=next((p,q,lo,hi) for (p,q),(lo,hi) in zip(SHIFTS,WRAPPED[i][j])
                     if (lo<l if lc else lo<=l) and (hi>h if hc else hi>=h))
            pairs.append((i,j,*hit))
        packets.append(((index,l,h,lc,hc),tuple(pairs)))
    return tuple(packets)


def repeated_cell_witness(threshold):
    threshold = positive_rational(threshold)
    need(threshold > 36, 'small-scale repeat regime')
    p = 1
    while 2*(8*p+2)**2 <= threshold:
        p += 1
    for q in range(p+1):
        if conflicts(CELLS[0], translate(CELLS[0],8*p,8*q),threshold):
            return p, q
    raise RuntimeError('proved overlapping-strip witness missing')


def obstruction(threshold):
    threshold = positive_rational(threshold)
    if threshold < 8:
        need(conflicts(CELLS[0], CELLS[0], threshold), 'interior self-conflict')
        return 'interior', 0
    if threshold > 36:
        return 'periodic self-conflict', repeated_cell_witness(threshold)
    return 'eight-clique', clique_certificate(threshold)


def eight_colouring():
    result = []
    for cell in CELLS:
        p, q = cell.x//2, cell.y//2
        result.append(4*((p//2+q//2) % 2)+2*(p % 2)+(q % 2))
    return tuple(result)


def full_shift_bound(threshold):
    threshold = positive_rational(threshold)
    floor_root = isqrt(threshold.numerator//threshold.denominator)
    ceil_root = floor_root + (floor_root*floor_root < threshold)
    return ceil_root//8 + 2


def verify_period_colouring(colours, thresholds):
    """One forbidden distance per colour, squared and divided by a^2.

    Complete source is this fixed40-cell motif with literal8a periods.
    """
    need(type(thresholds) is tuple and bool(thresholds), 'nonempty exact palette thresholds')
    thresholds = tuple(positive_rational(s) for s in thresholds)
    need(type(colours) is tuple and len(colours) == 40 and all(
         type(c) is int and 0 <= c < len(thresholds) for c in colours), 'complete exact cell colours')
    for i, left in enumerate(CELLS):
        for j in range(i,40):
            if colours[i] != colours[j]:
                continue
            s = thresholds[colours[i]]
            bound = full_shift_bound(s)
            for p,q in product(range(-bound,bound+1),repeat=2):
                need(not conflicts(left,translate(CELLS[j],8*p,8*q),s),
                     'same-colour forbidden distance in actual cells/period copies')
    return True


def palette_thresholds(a_squared, distances):
    a_squared = positive_rational(a_squared)
    need(type(distances) is tuple and bool(distances), 'exact palette distances')
    return tuple(positive_rational(d)**2/a_squared for d in distances)


def main():
    global CHECKS
    CHECKS = 0
    check(len(CELLS)==40 and sum(c.width*c.height for c in CELLS)==64,'complete40-cell motif')
    # Every open unit mesh square belongs to exactly one source cell.
    for x,y in product(range(8),repeat=2):
        check(sum(c.x <= x < c.x+c.width and c.y <= y < c.y+c.height for c in CELLS)==1,
              'literal partition, not just a cell count')
    strata = tuple(F(q,2) for q in range(16,73))
    for (index,l,h,lc,hc),pairs in interval_clique_packets():
        for i,j,p,q,lo,hi in pairs:
            check(distance_window(CELLS[i],translate(CELLS[j],8*p,8*q))==(lo,hi) and
                  (lo<l if lc else lo<=l) and (hi>h if hc else hi>=h),
                  'one fixed actual pair witness covers the entire real interval')
    for s in strata:
        clique = clique_certificate(s)
        for i,j in combinations(clique,2):
            check(any(lo<s<hi for lo,hi in WRAPPED[i][j]), 'exact clique pair')
        # Integer endpoints imply the same edge set throughout each open stratum.
        if s.denominator==2:
            left = F(s.numerator//2)
            check(compact_graph(left+F(1,3))==compact_graph(left+F(2,3)), 'open stratum invariance')
        check(obstruction(s)[0]=='eight-clique','complete compact scale branch')
        for i,j,(p,q) in lifted_clique_edges(s,clique):
            left=translate(CELLS[i],800,0)
            right=translate(CELLS[j],800+8*p,8*q)
            check(conflicts(left,right,s) and
                  min(left.x-left.y-left.height,right.x-right.y-right.height)>100,
                  'finite conflict witness moves entirely beyond a diagonal strip')
    for s in [F(q,3) for q in range(1,24)]:
        check(obstruction(s)==('interior',0),'all tested interior regime')
    for s in [F(q,3) for q in range(109,2000)]:
        p,q = repeated_cell_witness(s)
        check(conflicts(CELLS[0],translate(CELLS[0],8*p,8*q),s),'repeat witness actual metric')
    for p in range(1,101):
        check((8*p+6)**2 < 2*(8*p+2)**2,'successive strip overlap')
        for q in range(p):
            lo,_=distance_window(CELLS[0],translate(CELLS[0],8*p,8*(q+1)))
            _,hi=distance_window(CELLS[0],translate(CELLS[0],8*p,8*q))
            check(lo<hi,'within-strip strict overlap')
    colours=eight_colouring()
    check(verify_period_colouring(colours,(F(8),)*8),'constructive whole-plane8-colour certificate')
    check(set(colours)==set(range(8)),'all8colour classes used')
    check(not conflicts(CELLS[0],CELLS[0],8) and conflicts(CELLS[0],CELLS[0],7),
          'half-open diameter boundary is safe only at and above threshold')
    check(not conflicts(CELLS[0],translate(CELLS[0],8,0),36) and
          conflicts(CELLS[0],translate(CELLS[0],8,0),37),'unattained minimum boundary')
    check(not conflicts(CELLS[0],CELLS[6],8) and
          any(lo<8<hi for lo,hi in WRAPPED[0][6]),'period-copy edge invisible in one motif')
    check(conflicts(CELLS[0],CELLS[1],8) and F(4,8)!=1,'cell-centre-only hostile')
    typed=palette_thresholds(F(1,8),(1,1,1,1,F(1,2),F(1,2)))
    check(typed==(8,8,8,8,2,2),'colour-specific scale retained')
    for factor in (F(1,3),F(2),F(7,5)):
        check(palette_thresholds(F(1,8)*factor*factor,
              tuple(factor*d for d in (1,1,1,1,F(1,2),F(1,2))))==typed,
              'uniform dilation preserves all forbidden-distance ratios')
    hostile=0
    def reject(call):
        nonlocal hostile
        try: call()
        except ValueError:
            hostile+=1
            check(True,'malformed or false certificate rejected')
        else: raise RuntimeError('hostile accepted')
    for bad in (
        lambda: obstruction(True), lambda: obstruction(8.0), lambda: obstruction(0),
        lambda: compact_graph(37), lambda: verify_clique(8,(False,1,2,3,4,6,12,14)),
        lambda: verify_clique(8,(0,0,2,3,4,6,12,14)),
        lambda: distance_window(replace(CELLS[0],width=2.0),CELLS[1]),
        lambda: translate(CELLS[0],True,0),
        lambda: verify_period_colouring((0,)*40,(8,)),
        lambda: verify_period_colouring((True,)+colours[1:],(8,)*8),
        lambda: verify_period_colouring(colours,(F(9),)*8),
        lambda: verify_period_colouring(colours,(F(7),)*8),
        lambda: palette_thresholds(F(1,8),(1,1.0)),
    ): reject(bad)
    print('square_cell_conflicts_20261011g')
    print('model:40 half-open cells; constant per cell; colour period8a in both axes')
    print('all real scales: S<8 interior conflict;8<=S<=36 K8;S>36 periodic self-conflict')
    print('complete compact scale strata:',len(strata),'; retained clique templates:',len(CLIQUES),
          '; fixed-witness scale intervals:',len(SCALE_COVER))
    print('eight-colour positive control at a^2=1/8:',colours)
    print('critical K8:',CLIQUES[0])
    print('polychromatic thresholds:',typed)
    print('larger colour periods, split cells, and Stechkin historical tiling are not classified')
    print('hostiles:',hostile,'; exact checks:',CHECKS)


if __name__=='__main__':
    main()
