"""Exact Mersenne reset-two phase compiler; no orbit/ROOT lookup.

Imported debt portals are reused, not discovered or certified by paper headlines.
All guards use exact integers and all checks remain active under -O.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from collections import Counter

CHECKS = Counter()


def need(ok, text="invalid input"):
    if not ok:
        raise ValueError(text)


def check(ok, category):
    need(ok, category)
    CHECKS[category] += 1


def natural(x):
    need(type(x) is int and x >= 0, "exact natural")


def word(w):
    need(type(w) is tuple and all(type(a) is int and a > 0 for a in w),
         "exact positive valuation tuple")


def carrier(w):
    word(w)
    p = q = 1
    b = 0
    for a in w:
        p, q, b = 3*p, q << a, 3*b+q
    return p, q, b


def valuation(n):
    need(type(n) is int and n > 0, "positive exact valuation argument")
    return (n & -n).bit_length()-1


def replay(n, w):
    need(type(n) is int and n > 0 and n % 2, "positive odd source")
    word(w)
    for a in w:
        need(n != 1, "no ROOT padding")
        need(valuation(3*n+1) == a, "wrong actual valuation")
        n = (3*n+1) >> a
    return n


def odd_log3(target, precision):
    """K odd with 3**K=target mod2**precision; None if no such K.

    Returns its unique residue modulo2**(precision-2), precision>=3.
    """
    natural(target)
    need(type(precision) is int and precision >= 3, "precision>=3")
    target %= 1 << precision
    if target % 8 != 3:
        return None
    k = 1
    for j in range(1, precision-2):
        modulus = 1 << (j+3)
        if pow(3, k, modulus) != target % modulus:
            k += 1 << j
    need(pow(3, k, 1 << precision) == target, "lift invariant")
    return k


@dataclass(frozen=True)
class Phase:
    tail: tuple
    residue: int
    modulus: int
    least: int


def compile_phase(tail):
    """Exact first-hit tail at z_K=(3**K-1)/2, for odd K>=3."""
    word(tail)
    if not tail:
        return Phase((), 1, 2, 3)
    p, q, b = carrier(tail)
    cost = sum(tail)
    target = (p-2*b+2*q)*pow(p, -1, 4*q) % (4*q)
    k = odd_log3(target, cost+2)
    if k is None:
        return None
    modulus = 1 << cost
    # Proper-prefix values must exceed ROOT. All are positive affine in 3**K.
    lower = 3
    for i in range(1, len(tail)):
        pi, qi, bi = carrier(tail[:i])
        while pi*3**lower <= 2*qi+pi-2*bi:
            lower += 2
    least = k + modulus*max(0, (lower-k+modulus-1)//modulus)
    return Phase(tail, k, modulus, least)


def accepts(phase, k):
    phase_guard(phase)
    need(type(k) is int and k >= 3 and k % 2, "odd Mersenne exponent>=3")
    return k >= phase.least and k % phase.modulus == phase.residue


def phase_guard(phase):
    need(type(phase) is Phase and all(type(x) is int for x in
         (phase.residue, phase.modulus, phase.least))
         and phase == compile_phase(phase.tail), "canonical exact phase")


def endpoint_is_root(tail, k):
    p, q, b = carrier(tail)
    rhs = 2*q+p-2*b
    if rhs <= 0 or rhs % p:
        return False
    rhs //= p
    exponent = 0
    while rhs % 3 == 0:
        exponent += 1
        rhs //= 3
    return rhs == 1 and exponent == k


def last_valuation(tail, k):
    """Actual next exponent after the tail, using modular powers only."""
    phase = compile_phase(tail)
    need(phase is not None and accepts(phase, k), "tail guard")
    need(not endpoint_is_root(tail,k), "no next edge after ROOT")
    p, q, b = carrier(tail)
    c = 6*b-3*p+2*q
    precision = sum(tail)+3
    while True:
        residue = (3*p*pow(3, k, 1 << precision)+c) % (1 << precision)
        if residue:
            result = valuation(residue)-sum(tail)-1
            need(result >= 1, "positive next exponent")
            return result
        precision *= 2


@dataclass(frozen=True)
class Portal:
    s: int
    e: int
    source: tuple  # u without its variable final exponent
    child: tuple   # v without its shifted final exponent
    offset: int

    def tail(self):
        return (2,)*(self.e-1)+self.source


PORTALS = (
    Portal(1, 1, (1, 6), (1, 2, 1), 2),
    Portal(1, 2, (1, 6, 4), (2, 1, 1, 1, 3), 2),
    Portal(1, 3, (1, 2, 9), (1, 1, 1, 1, 1, 2), 4),
    Portal(1, 4, (1, 14), (2, 2, 1, 3, 3, 1), 2),
    Portal(2, 1, (6,), (1, 1), 2),
    Portal(2, 2, (10,), (3, 2, 1), 2),
    Portal(2, 3, (10,), (1, 1, 1, 3), 2),
    Portal(2, 4, (3, 1, 11), (1, 1, 2, 1, 1, 1, 2), 4),
)

# Parent's separately proved D=2 collision. This package only compiles its guard.
NEW_TAIL = (3, 2, 2, 2, 2, 1, 5, 1, 1, 2)
NEW_PARTNER = (2, 2, 2, 1, 2, 2, 2, 1, 1, 3, 1, 1, 1)


def residual_progression(phase):
    """Intersect with the inherited named cap-escape K=1 mod1458, K>1024."""
    phase_guard(phase)
    half = phase.modulus//2
    j = ((phase.residue-1)//2)*pow(729, -1, half) % half if half > 1 else 0
    period = 729*phase.modulus
    least = 1+1458*j
    lower = max(1025, phase.least)
    least += period*max(0, (lower-least+period-1)//period)
    return least, period


def cell_guard(cell):
    need(type(cell) is tuple and len(cell) == 2, "typed cell")
    r, m = cell
    need(type(r) is int and type(m) is int and m > 0 and m & (m-1) == 0
         and 0 <= r < m, "canonical dyadic cell")


def subtract_cell(cells, removed):
    """Exact finite complement refinement; no enumeration of a common modulus."""
    cell_guard(removed)
    result = []
    rr, rm = removed
    for cell in cells:
        cell_guard(cell)
        r, m = cell
        if (r-rr) % min(m, rm):
            result.append(cell)
        elif m >= rm:
            continue
        else:
            active = (r, m)
            while active[1] < rm:
                r, m = active
                children = ((r, 2*m), (r+m, 2*m))
                active = children[rr % (2*m) != r]
                result.append(children[rr % (2*m) == r])
    return tuple(sorted(result, key=lambda c: (c[1], c[0])))


def ledger(phases):
    """Dyadic cylinder union and asymptotic mass; finite least-cuts stay separate.

    All named rows here start at their positive residue, so on odd K>=3 this
    is also their pointwise complement. Generic callers must apply accepts.
    """
    cells = ((0, 1),)
    for phase in phases:
        phase_guard(phase)
        cells = subtract_cell(cells, ((phase.residue-1)//2, phase.modulus//2))
    remainder = sum((F(1, m) for _, m in cells), F())
    return cells, 1-remainder


def main():
    # Independent exhaustive inversion of the finite exponential maps.
    for precision in range(3, 13):
        direct = {pow(3, k, 1 << precision): k for k in range(1, 1 << (precision-2), 2)}
        for c in range(1 << precision):
            check(odd_log3(c, precision) == direct.get(c), "odd_power_phase_iff")
    for depth in range(1, 11):
        m = 1 << depth
        images = {(3*(pow(9, j, 8*m)-1)//8) % m for j in range(m)}
        check(len(images) == m, "normalized_chart_bijection")
    for length in range(4):
        for tail in product(range(1, 5), repeat=length):
            phase = compile_phase(tail)
            for k in range(3, 130, 2):
                try:
                    endpoint = replay((3**k-1)//2, tail)
                    actual = True
                except ValueError:
                    actual = False
                check(actual == (phase is not None and accepts(phase, k)), "literal_tail_iff")
                if actual and endpoint != 1:
                    check(last_valuation(tail, k) == valuation(3*endpoint+1), "modular_next_valuation")
    expected = (None, (5589,8192), (32945,65536), (1035585,2097152),
                (31,64), (1289,4096), (1889,16384), (1522049,2097152))
    phases = []
    for portal, want in zip(PORTALS, expected):
        phase = compile_phase(portal.tail())
        check((None if phase is None else (phase.residue,phase.modulus)) == want, "portal_phases")
        for a in range(1, 8):
            # Identity at one initial1; appending more initial1s commutes with half-child.
            source = (1,)+(2,)*portal.e+portal.source+(a,)
            child = (2*portal.e+portal.s,)+portal.child+(a+portal.offset,)
            p,q,b = carrier(source)
            pp,qq,bb = carrier(child)
            check(p == pp and q == 2*qq and 2*bb-b == p, "inherited_portal_carry")
        if phase is None:
            continue
        phases.append(phase)
        start, period = residual_progression(phase)
        for j in (0,1,17,10**30):
            k = start+period*j
            check(accepts(phase,k) and k%1458 == 1 and k>1024, "residual_CRT")
            check(last_valuation(phase.tail,k)>=1, "large_symbolic_receipt")
    # Actual full route replay for modest exponents only; no huge-source expansion.
    for index in (1,4,5,6):
        portal=PORTALS[index]
        phase=compile_phase(portal.tail())
        k=phase.least
        a=last_valuation(portal.tail(),k)
        sw=(1,)*(k-1)+(2,)*portal.e+portal.source+(a,)
        cw=(1,)*(k-2)+(2*portal.e+portal.s,)+portal.child+(a+portal.offset,)
        check(replay((1<<k)-1,sw)==replay((1<<(k-1))-1,cw), "actual_whole_join")
    new=compile_phase(NEW_TAIL)
    check((new.residue,new.modulus)==(1459,2097152), "parent_D2_phase")
    # The D2 collision's one-run base uses source4x+3 and partner x.
    for a in range(1,8):
        p,q,b=carrier((1,1,2)+NEW_TAIL+(a,))
        pp,qq,bb=carrier(NEW_PARTNER+(a+2,))
        check(4*p*qq==pp*q and (3*p+b)*qq==bb*q, "parent_D2_affine_identity")
    k=1459;a=last_valuation(NEW_TAIL,k)
    check(replay((1<<k)-1,(1,)*(k-1)+(2,)+NEW_TAIL+(a,)) ==
          replay((1<<(k-2))-1,(1,)*(k-3)+NEW_PARTNER+(a+2,)), "actual_parent_D2_join")
    cells,mass=ledger(tuple(phases))
    check(mass==F(16849,524288), "disjoint_portal_mass")
    full,fullmass=ledger(tuple(phases)+(new,))
    check(fullmass==F(33699,1048576), "full_union_mass")
    check(ledger(tuple(phases)+(new,new))[1]==fullmass, "duplicate_not_new_coverage")
    for i,c in enumerate(full):
        check(all((c[0]-d[0])%min(c[1],d[1]) for d in full[i+1:]), "complement_disjoint")
    for k in range(3,8194,2):
        covered=any(k%p.modulus==p.residue for p in phases+[new])
        t=(k-1)//2
        check(covered != any(t%m==r for r,m in full), "ledger_partition")
    # Arbitrarily small residual measure may retain a designated positive integer.
    for depth in range(1,25):
        m=1<<depth
        resolved=subtract_cell(((0,1),),(1459%m,m))
        check(not any(1459%q==r for r,q in resolved)
              and sum((F(1,q) for r,q in resolved),F())==1-F(1,m),
              "singleton_survival_hostile")
        # Exact repetition of one event retains mass1/2, not(1/2)**depth.
        repeated=F(sum(all(x==1 for _ in range(depth)) for x in (0,1)),2)
        check(repeated==F(1,2) and (depth==1 or repeated>F(1,2)**depth),
              "repeated_event_not_independent")
    for operation in (lambda:compile_phase([2]),lambda:compile_phase((True,)),
                      lambda:odd_log3(3,True),lambda:odd_log3(3.0,4),
                      lambda:accepts(new,True),lambda:accepts(new,1459.0),
                      lambda:accepts(Phase((2,),True,4,5),5),
                      lambda:subtract_cell(((0,1),),(True,2)),
                      lambda:replay(1,(2,)),lambda:last_valuation(NEW_TAIL,31),
                      lambda:last_valuation((3,2),3)):
        try: operation()
        except ValueError: check(True,"malformed_or_wrong_guard")
        else: check(False,"malformed_or_wrong_guard")
    print("Exact reset-two Mersenne phase compiler; supplied four paper headlines accepted")
    print("K odd>=3; z=(3^K-1)/2; normalized exponent chart is a 2-adic isometry")
    for portal in PORTALS:
        p=compile_phase(portal.tail())
        print('portal', (portal.s,portal.e), 'tail',portal.tail(), 'phase',
              None if p is None else (p.residue,p.modulus), 'named residual',
              None if p is None else residual_progression(p))
    print('parent D2 tail',NEW_TAIL,'phase',(new.residue,new.modulus),
          'named residual',residual_progression(new))
    print('seven portals odd-K mass',mass,'six beyond old31/64',mass-F(1,32))
    print('with separately proved parent D2 row',fullmass,'uncovered',1-fullmass)
    print('complement dyadic leaves',len(full),'largest modulus',max(m for _,m in full))
    print('Whole original sources preserved; smaller dependencies are not ROOT certificates')
    for name,count in sorted(CHECKS.items()): print(name,count)
    print('CHECKS',sum(CHECKS.values()))


if __name__=='__main__':
    main()
