"""Exact basin, clock, and affine-specialization observers.

No Haar theorem, orbit search, or universal ROOT claim is used.  Cycles and
ROOT words are supplied finite receipts; every edge is authenticated exactly.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from math import isqrt

import collatz_square_norm_transfer_20261007f as norms

CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def integer(x, minimum=None):
    if type(x) is not int or (minimum is not None and x < minimum):
        raise ValueError('exact integer in declared domain required')
    return x


def prime(p):
    integer(p, 2)
    if any(p % d == 0 for d in range(2, isqrt(p)+1)):
        raise ValueError('prime base required')
    return p


def C(p, n):
    prime(p); integer(n, 1)
    i = n % p
    return n//p if i == 0 else ((p+1)*n+p-i)//p


def padic_residue(p, x):
    prime(p)
    if type(x) not in (int,F):
        raise ValueError('exact rational p-adic input required')
    x = F(x)
    if x.denominator % p == 0:
        raise ValueError('rational input is not p-adically integral')
    return x.numerator*pow(x.denominator,-1,p) % p


def rational_C(p, x):
    """Rational points of Z_p; distinct domain from positive ordinary C."""
    i = padic_residue(p,x)
    x = F(x)
    return x/p if i == 0 else ((p+1)*x+p-i)/p


def inverse_branches(p, target):
    padic_residue(p,target)
    target = F(target)
    return (p*target,)+tuple((p*target-p+i)/(p+1) for i in range(1,p))


@dataclass(frozen=True)
class Frame:
    p: int
    k: int
    A: int


def audit_frame(frame):
    if type(frame) is not Frame:
        raise ValueError('typed Frame required')
    prime(frame.p); integer(frame.k); integer(frame.A)
    return frame


def value(frame, reference):
    """Evaluate the marked formal affine relation at this exact source."""
    audit_frame(frame); integer(reference, 1)
    P = frame.p+1
    h = max(0, -frame.k)
    numerator = P**max(0, frame.k)*reference+frame.A
    denominator = P**h
    if numerator % denominator or numerator <= 0:
        raise ValueError('frame is not a positive integral source at this reference')
    return numerator//denominator


def equality_numerator(frame, reference):
    audit_frame(frame); integer(reference, 1)
    P = frame.p+1
    return (P**max(frame.k, 0)-P**max(-frame.k, 0))*reference+frame.A


def advance(frame, reference):
    """Rational branch derivation, independent of the inherited integer table."""
    source = value(frame, reference)
    p, P = frame.p, frame.p+1
    i, j = reference % p, source % p
    mi, mj = (P if i else 1), (P if j else 1)
    bi, bj = (p-i if i else 0), (p-j if j else 0)
    k1 = frame.k+int(j != 0)-int(i != 0)
    M1 = F(P)**k1
    e = F(frame.A, P**max(0, -frame.k))
    e1 = (mj*e+bj-M1*bi)/p
    A1 = e1*P**max(0, -k1)
    if A1.denominator != 1:
        raise ArithmeticError('integer frame lattice was not preserved')
    result = Frame(p, k1, A1.numerator)
    v1 = C(p, reference)
    if value(result, v1) != C(p, source):
        raise ArithmeticError('actual branch relation was not preserved')
    return result, v1


@dataclass(frozen=True)
class Path:
    p: int
    states: tuple


@dataclass(frozen=True)
class Cycle:
    p: int
    states: tuple


def typed_states(states):
    if type(states) is not tuple or not states:
        raise ValueError('nonempty tuple of exact positive states required')
    for n in states:
        integer(n, 1)


def audit_path(path):
    if type(path) is not Path:
        raise ValueError('typed Path required')
    prime(path.p); typed_states(path.states)
    if any(C(path.p, a) != b for a, b in zip(path.states, path.states[1:])):
        raise ValueError('unauthenticated path edge')
    return path


def audit_cycle(cycle):
    if type(cycle) is not Cycle:
        raise ValueError('typed Cycle required')
    prime(cycle.p); typed_states(cycle.states)
    if len(set(cycle.states)) != len(cycle.states):
        raise ValueError('primitive cycle receipt must have distinct states')
    if cycle.states[0] != min(cycle.states):
        raise ValueError('cycle must start at its least state')
    if any(C(cycle.p, a) != b for a, b in
           zip(cycle.states, cycle.states[1:]+cycle.states[:1])):
        raise ValueError('unauthenticated cycle edge')
    return cycle


def terminal_label(path, cycle):
    """Complete synchronous label on the supplied eventually periodic domain."""
    audit_path(path); audit_cycle(cycle)
    if path.p != cycle.p or path.states[-1] not in cycle.states:
        raise ValueError('path does not enter the supplied cycle')
    phase = (cycle.states.index(path.states[-1])-len(path.states)+1) % len(cycle.states)
    return cycle, phase


def reaches_one(path, cycle):
    terminal_label(path, cycle)
    return 1 in cycle.states


def synchronous(path_a, cycle_a, path_b, cycle_b):
    if path_a.p != path_b.p:
        raise ValueError('same dynamical map required')
    return terminal_label(path_a, cycle_a) == terminal_label(path_b, cycle_b)


def root_path(source, odd_word):
    """Expand a supplied first-hit odd word to its exact shortcut-T clock."""
    integer(source, 1)
    if source % 2 == 0 or type(odd_word) is not tuple:
        raise ValueError('positive odd source and tuple word required')
    for a in odd_word:
        integer(a, 1)
    states = [source]
    current = source
    for a in odd_word:
        if current == 1:
            raise ValueError('ROOT padding is not a first-hit word')
        numerator = 3*current+1
        if numerator % (1 << a) or (numerator >> a) % 2 == 0:
            raise ValueError('wrong exact odd valuation')
        for _ in range(a):
            current = C(2, current)
            states.append(current)
    if current != 1 or 1 in states[:-1]:
        raise ValueError('supplied word is not a strict ROOT receipt')
    return audit_path(Path(2, tuple(states)))


def finite_observer_collision(p, c, depth, multiplier=1):
    """Same frame and reference p-adic prefix; one equality, one inequality."""
    prime(p); integer(c, 1); integer(depth, 1); integer(multiplier, 1)
    frame = Frame(p, -1, p*c)
    reference = c+(p+1)*p**depth*multiplier
    return frame, c, reference


ROOT2 = Cycle(2, (1, 2))
BAD3 = Cycle(3, (7, 10, 14, 19, 26, 35, 47, 63, 21))
BAD11 = Cycle(11, (
    642,701,765,835,911,994,1085,1184,1292,1410,1539,1679,1832,
    1999,2181,2380,2597,2834,3092,3374,3681,4016,4382,4781,5216,
    5691,6209,6774,7390,8062,8795,9595,10468,11420,12459,13592,
    14828,1348,1471,1605,1751,1911,2085,2275,2482,2708,2955,
    3224,3518,3838,4187,4568,4984,5438,5933,6473,7062))


def rejected(thunk):
    try:
        thunk()
    except (ValueError, TypeError):
        return True
    return False


def main():
    global CHECKS
    CHECKS = 0
    for cycle in (ROOT2, BAD3, BAD11):
        audit_cycle(cycle)
        need(len(set(cycle.states)) == len(cycle.states), 'primitive supplied cycle')
    need(not reaches_one(Path(3, (30,10)), BAD3), 'non-ROOT merge p3')
    need(synchronous(Path(3, (30,10)), BAD3, Path(3, (7,10)), BAD3), 'p3 common cycle phase')
    need(not reaches_one(Path(11, (7711,701)), BAD11), 'non-ROOT merge p11')
    need(synchronous(Path(11, (7711,701)), BAD11, Path(11, (642,701)), BAD11), 'p11 common cycle phase')

    frame_controls = 0
    for p, u0, v0 in ((2,10,3), (3,30,7), (11,7711,642)):
        frame, v = Frame(p,0,u0-v0), v0
        frame, v = advance(frame, v)
        need(value(frame,v) == v, 'numeric merge after one step')
        need(frame.k == -1 and frame.A == p*v, 'nonidentity specialization')
        for _ in range(120):
            need(equality_numerator(frame,v) == 0, 'exact evaluation detects merge')
            need(frame.k == -1 and frame.A != 0, 'formal state does not absorb')
            frame, v = advance(frame,v)
            frame_controls += 1

    grid = 0
    for p in (2,3,11):
        for k in range(-3,4):
            for A in range(-12,13):
                frame = Frame(p,k,A)
                for v in range(1,31):
                    try:
                        u = value(frame,v)
                    except ValueError:
                        continue
                    need(equality_numerator(frame,v) == (p+1)**max(0,-k)*(u-v), 'exact numerator identity')
                    f1,v1 = advance(frame,v)
                    need(value(f1,v1) == C(p,u), 'independent direct edge check')
                    grid += 1

    observer_controls = 0
    for p in (2,3,11):
        for c in (1,5,7):
            for depth in range(1,13):
                frame, v0, v1 = finite_observer_collision(p,c,depth)
                need(value(frame,v0) == v0 and value(frame,v1) != v1, 'same formal frame, distinct exact event')
                need((v1-v0) % p**depth == 0, 'same reference finite residue')
                for _ in range(depth):
                    need(v0 % p == v1 % p, 'same branch itinerary prefix')
                    v0, v1 = C(p,v0), C(p,v1)
                observer_controls += 1

    inverse_edges = 0
    for p,max_depth in ((2,6),(3,5),(11,3)):
        layer = {F(1)}
        for depth in range(1,max_depth+1):
            next_layer = set()
            for y in layer:
                children = inverse_branches(p,y)
                need(tuple(padic_residue(p,x) for x in children) == tuple(range(p)), 'one inverse in each residue branch')
                for x in children:
                    need(rational_C(p,x) == y, 'inverse branch checked exactly')
                    next_layer.add(x)
                    inverse_edges += 1
            need(len(next_layer) == p**depth, 'exact root preimage layer size')
            layer = next_layer

    norm_controls = 0
    for bits in range(2,13):
        for depth in range(1,9):
            u = 5+3**depth*2**bits
            v = 5+3**(depth+1)*2**bits
            frame = Frame(2,-1,10)
            guard = norms.compile_guard(5 % (2**bits*3**depth), bits, depth)
            need(value(frame,v) == u and u != v, 'different integers in the same marked frame')
            need(norms.accepts(guard,5) and norms.accepts(guard,u) and norms.accepts(guard,v), 'both marked norm observers collide')
            need(norms.source_from_norm(norms.norm_coordinate(u)) == u, 'full positive norm remains lossless')
            norm_controls += 1

    a, b = root_path(151,(1,1,10)), root_path(75,(1,2,8))
    need(len(a.states)-1 == 12 and len(b.states)-1 == 11, 'factory shortcut clocks')
    need(reaches_one(a,ROOT2) and reaches_one(b,ROOT2), 'both supplied factory sources grounded')
    need(not synchronous(a,ROOT2,b,ROOT2), 'odd-clock coincidence does not fix shortcut phase')
    for p in (2,3,11):
        root = Cycle(p,tuple(range(1,p+1)))
        audit_cycle(root)
        for x in root.states:
            for y in root.states:
                need(synchronous(Path(p,(x,)),root,Path(p,(y,)),root) == (x == y), 'root-cycle phases stay distinct')
    # Terminal labels do not depend on where along the witnessed cycle we stop.
    for cycle in (ROOT2,BAD3,BAD11):
        for i in range(len(cycle.states)):
            path = Path(cycle.p,(cycle.states[i],))
            for _ in range(3):
                extended = Path(path.p,path.states+(C(path.p,path.states[-1]),))
                need(terminal_label(path,cycle) == terminal_label(extended,cycle), 'terminal label invariant under extension')
                path = extended

    hostile = (
        lambda: C(True,3), lambda: C(4,3), lambda: C(2,True),
        lambda: value(Frame(2,False,10),5), lambda: value(Frame(2,-1,10),6),
        lambda: audit_path(Path(2,(3,True))), lambda: audit_path(Path(2,(3,4))),
        lambda: audit_cycle(Cycle(2,(1,True))), lambda: audit_cycle(Cycle(2,(2,1))),
        lambda: audit_cycle(Cycle(2,(1,2,1))), lambda: root_path(1,(2,)),
        lambda: root_path(151,(1,1,9)), lambda: root_path(75,(1,2,8.0)),
        lambda: terminal_label(Path(2,(3,5)),ROOT2),
        lambda: synchronous(Path(2,(1,)),ROOT2,Path(3,(7,)),BAD3),
        lambda: finite_observer_collision(2,5,True),
        lambda: inverse_branches(3,True), lambda: inverse_branches(3,F(1,3)),
    )
    for test in hostile:
        need(rejected(test), 'malformed or false receipt rejected')
    print('Exact supplied cycles: p=3, minimum 7, length 9; p=11, minimum 642, length 57.')
    print('Merge 30/7 -> 10 is in the wrong p=3 basin; merge 10/3 -> 5 is in the p=2 ROOT basin.')
    print('Both actual merges retain formal debt k=-1; exact reference evaluation detects equality.')
    print(f'Frame controls: {grid} admissible grid states and {frame_controls} post-merge steps.')
    print(f'Finite observer collisions: {observer_controls}; mixed marked norm collisions: {norm_controls}.')
    print(f'Rational p-adic ROOT-preimage edges: {inverse_edges}; each depth-N layer has exactly p^N points.')
    print('Grounded factory 151/75: odd ROOT ranks 3/3, shortcut ROOT clocks 12/11, distinct synchronous phases.')
    print(f'Malformed/false receipt controls: {len(hostile)}.')
    print('Scope: supplied finite paths/cycles and observer theorems; no Haar proof audit or new ROOT coverage.')
    print(f'ALL {CHECKS} CHECKS PASSED')


if __name__ == '__main__':
    main()
