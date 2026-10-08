"""Complete bounded ordinary inverse-word bank at the original Mersenne source.

All 252 labelled entries are retained; antichains are counting views only.
Production does not materialize the Mersenne source or discover ROOT words.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache

import collatz_complement_guard_fusion_20261007e as prior

routes = prior.routes
BASE, STRIDE = prior.BASE, prior.STRIDE
MAX_LENGTH = 8
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def natural(n, least=0):
    if type(n) is not int or n < least:
        raise ValueError('exact integer in declared domain required')
    return n


def compositions(total, length):
    natural(total, 1); natural(length, 1)
    if total < length:
        return
    if length == 1:
        yield (total,)
        return
    for first in range(1, total-length+2):
        for tail in compositions(total-first, length-1):
            yield (first,)+tail


def word_universe():
    """Every positive word of length 1..8 with 2**cost < 3**length."""
    answer = []
    for length in range(1, MAX_LENGTH+1):
        p = 3**length
        for cost in range(length, p.bit_length()):
            if 1 << cost < p:
                answer.extend(compositions(cost, length))
    return tuple(answer)


@dataclass(frozen=True)
class Entry:
    word: tuple
    P: int
    Q: int
    B: int
    least_source: int
    ordinary_period: int
    exponent_residue: int
    exponent_period: int
    parameter_residue: int
    parameter_period: int
    minimum_parameter: int


@lru_cache(maxsize=None)
def _compile(word):
    p, q, b = routes.carrier(word)
    if q >= p or word[-1] % 2:
        raise ValueError('contracting word with even final valuation required')
    least = b*pow(q, -1, p) % p
    if least % 2 == 0:
        least += p
    if not (least > 1 and 0 < (q*least-b)//p < least):
        raise ArithmeticError('general positive native inverse boundary failed')
    log, period = prior.forms.two_log_three_power((q+b) % p, len(word))
    exponent = (log-sum(word)) % period
    if exponent % 2 != 1:
        raise ArithmeticError('odd Mersenne exponent phase lost')
    modulus = period//2
    residue = ((exponent-BASE)//2)*pow(STRIDE//2, -1, modulus) % modulus
    # The least exponent in this parameter cell is already far above every
    # finite ordinary-source threshold; compare bit lengths, never 2**BASE.
    if BASE+STRIDE*residue <= least.bit_length():
        raise ArithmeticError('fixed chart requires an explicit height cut')
    return Entry(word, p, q, b, least, 2*p, exponent, period,
                 residue, modulus, 0)


def compile_word(word):
    routes.letters(word)
    if not 1 <= len(word) <= MAX_LENGTH:
        raise ValueError('bounded bank length must be 1..8')
    return _compile(word)


def audit(entry):
    if type(entry) is not Entry:
        raise ValueError('exact Entry required')
    routes.letters(entry.word)
    for name in ('P', 'Q', 'B', 'least_source', 'ordinary_period',
                 'exponent_residue', 'exponent_period', 'parameter_residue',
                 'parameter_period', 'minimum_parameter'):
        natural(getattr(entry, name))
    if compile_word(entry.word) != entry:
        raise ValueError('noncanonical word/guard/height packet')
    return entry


def entries():
    return tuple(compile_word(w) for w in word_universe() if w[-1] % 2 == 0)


def parameter_contains(entry, parameter):
    audit(entry); natural(parameter)
    return (parameter >= entry.minimum_parameter and
            (parameter-entry.parameter_residue) % entry.parameter_period == 0)


def ordinary_source(entry, parameter):
    audit(entry); natural(parameter)
    return entry.least_source+entry.ordinary_period*parameter


def receipt(entry, source):
    """Validate a supplied ordinary odd source; export a paid dependency."""
    audit(entry); routes.odd(source)
    if source < entry.least_source or (source-entry.least_source) % entry.ordinary_period:
        raise ValueError('supplied source outside native inverse guard')
    child = (entry.Q*source-entry.B)//entry.P
    return routes.audit(routes.Receipt(source, child, (), entry.word, source))


def discharge(entry, source, supplied_child_root_word):
    return routes.discharge(receipt(entry, source), supplied_child_root_word)


def symbolic_child_mod(entry, parameter, modulus):
    """Read the same child's residue without materializing its source."""
    audit(entry); natural(parameter); natural(modulus, 1)
    if not parameter_contains(entry, parameter):
        raise ValueError('supplied parameter outside the exact native phase')
    precision = entry.P*modulus
    source = (pow(2, BASE+STRIDE*parameter, precision)-1) % precision
    numerator = entry.Q*source-entry.B
    if numerator % entry.P:
        raise ArithmeticError('lost ternary denominator precision')
    return (numerator//entry.P) % modulus


def _cell(cell):
    if type(cell) is not tuple or len(cell) != 2:
        raise ValueError('canonical residue/period pair required')
    r, p = cell
    natural(r); natural(p, 1)
    q = p
    while q % 3 == 0:
        q //= 3
    if q != 1 or r >= p:
        raise ValueError('canonical ternary-power cell required')
    return cell


def cell_union(cells):
    """Counting-only antichain; this never deletes labelled entries."""
    if type(cells) is not tuple:
        raise ValueError('tuple of cells required')
    for cell in cells:
        _cell(cell)
    result = []
    for r, p in sorted(set(cells), key=lambda cell: (cell[1], cell[0])):
        if not any((r-s) % q == 0 for s, q in result):
            result.append((r, p))
    return tuple(result)


def coverage_cells(bank=None):
    if bank is None:
        bank = entries()
    if type(bank) is not tuple:
        raise ValueError('tuple of labelled entries required')
    for entry in bank:
        audit(entry)
    return cell_union(tuple((e.parameter_residue, e.parameter_period) for e in bank))


def cell_mass(cells):
    return sum((F(1, p) for _, p in cell_union(cells)), F(0))


def combined_cells():
    old = tuple((e.parameter_residue, e.parameter_period) for e in prior.entries(12))
    return cell_union(coverage_cells()+old)


def main():
    universe = word_universe()
    need(len(universe) == 953 and len(set(universe)) == 953, 'complete word count')
    # Independent composition generator: choose cuts between cost many ones.
    alternate = set()
    for cost in range(1, 13):
        for mask in range(1 << (cost-1)):
            w, block = [], 1
            for j in range(cost-1):
                if mask >> j & 1:
                    w.append(block); block = 1
                else:
                    block += 1
            w.append(block)
            if len(w) <= 8 and 1 << cost < 3**len(w):
                alternate.add(tuple(w))
    need(set(universe) == alternate, 'independent complete composition universe')
    bank = entries()
    need(len(bank) == 252, 'all even-terminal labels retained')
    need(sum(w[-1] % 2 for w in universe) == 701, 'all missing phases explained')
    for w in universe:
        p, q, b = routes.carrier(w)
        need(((q+b) % 3 != 0) == (w[-1] % 2 == 0), 'terminal parity iff unit target')
        family = routes.compile_family(w, 0)
        residue = b*pow(q, -1, p) % p
        if residue % 2 == 0:
            residue += p
        need(family.least == residue, 'no lost positive native head, all 953 words')
    literal = 0
    for e in bank:
        for j in range(3):
            n = ordinary_source(e, j)
            proof = receipt(e, n)
            need(proof.child < n and proof.endpoint == n, 'immutable source payment')
            # Independent backwards rational arithmetic, not another forward replay.
            z = F(n)
            for a in reversed(e.word):
                z = (2**a*z-1)/3
                need(z.denominator == 1 and z.numerator > 0 and z.numerator % 2,
                     'every inverse stage positive odd')
            need(z == proof.child, 'independent inverse child')
            literal += 1
        for shift in (0, 1, 2):
            t = e.parameter_residue+shift*e.parameter_period
            need(parameter_contains(e, t), 'phase lift')
            exponent = BASE+STRIDE*t
            need((e.Q*(pow(2, exponent, e.P)-1)-e.B) % e.P == 0,
                 'direct modular native guard')
            for modulus in (1, 3, 19, 64):
                child = symbolic_child_mod(e, t, modulus)
                source = (pow(2, exponent, e.P*modulus)-1) % (e.P*modulus)
                need((e.P*child+e.B-e.Q*source) % (e.P*modulus) == 0,
                     'source-owned child modular readout')
    core = coverage_cells(bank)
    need(len(core) == 12 and cell_mass(core) == F(284, 729), 'exact bounded union')
    census = 0
    for t in range(3**7):
        labelled = any((t-e.parameter_residue) % e.parameter_period == 0 for e in bank)
        counted = any((t-r) % p == 0 for r, p in core)
        need(labelled == counted, 'counting view preserves guard union')
        census += labelled
    need(census == 852, 'independent full modulus census')
    oldmass = prior.coverage_mass(prior.entries(12))
    unionmass = cell_mass(combined_cells())
    need(unionmass == F(150929272, 387420489), 'exact old/new union')
    need(unionmass-oldmass == F(332, 6561), 'exact newly counted ternary mass')
    need(not any(parameter_contains(e, 23) for e in bank), 't23 remains missing')
    need(not any((23-r) % p == 0 for r, p in combined_cells()), 't23 misses combined bank')
    need(all((23-r) % min(27, p) != 0 for r, p in combined_cells()),
         'entire t23 mod27 cell remains outside the combined ternary bank')
    need(not any(parameter_contains(e, 0) for e in bank), 'least missing parameter is zero')
    e12 = compile_word((1, 2))
    need(receipt(e12, 31).child == 27, 'paid dependency is not a ROOT certificate')
    need(discharge(e12, 13, (1, 2, 3, 4)) == (3, 4),
         'supplied child ROOT proof is cut at the same source endpoint')
    need(routes.replay(3, (1, 4)) == (3, 5, 1) and routes.replay(5, (4,)) == (5, 1),
         'all positive odd possible smaller children of7 exclude7')
    need(routes.replay(7, (1, 1, 2, 3, 4))[-1] == 1,
         'forward ROOT certificate exists outside the pure inverse payment scheme')
    e_old = compile_word((1, 1, 2, 2))
    need((e_old.P, e_old.Q, e_old.B, e_old.least_source) == (81, 64, 73, 91),
         'inherited F91 recognized, not claimed new')
    malformed = (
        lambda: compile_word((1, True)), lambda: compile_word([1, 2]),
        lambda: compile_word(()), lambda: compile_word((2,)),
        lambda: compile_word((1, 1)), lambda: receipt(e12, 1),
        lambda: receipt(e12, 3), lambda: receipt(e12, True),
        lambda: discharge(e12, 13, (4,)),
        lambda: parameter_contains(e12, 1.0),
        lambda: symbolic_child_mod(e12, 23, 19),
        lambda: audit(replace(e12, minimum_parameter=False)),
        lambda: cell_union(((True, 3),)), lambda: cell_union(((1, 6),)),
    )
    for call in malformed:
        try:
            call()
        except (ValueError, TypeError):
            pass
        else:
            raise ValueError('malformed or non-native input accepted')
        need(True, 'typed hostile rejected')
    print('Complete contracting ordinary-word universe: 953; lengths 1..8; costs <=12')
    print('No Mersenne phase (odd final valuation): 701; retained labels: 252')
    print('Independent ordinary-source receipts:', literal)
    print('Maximum least ordinary source:', max(e.least_source for e in bank))
    print('Counting cells (residue, period, one shortest word):')
    for r, p in core:
        e = min((e for e in bank if (e.parameter_residue, e.parameter_period) == (r, p)),
                key=lambda e: (len(e.word), sum(e.word), e.word))
        print(' ', r, p, e.word)
    print('Standalone ternary mass:', cell_mass(core), '; full census 852/2187')
    print('Old fusion mass:', oldmass)
    print('Combined ternary mass:', unionmass)
    print('Exact increment over old fusion:', unionmass-oldmass)
    print('Entire t=23 mod27 cell: no matching row; t=0: no matching row')
    print('Ordinary source7: no smaller inverse-word child at any depth; supplied forward proof checks')
    print('No ROOT search; children remain source-owned external obligations')
    print('Exact checks:', CHECKS)


if __name__ == '__main__':
    main()
