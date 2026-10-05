"""Finite seed obligations, checked receipts, and infinite controller families.

No search for a seed's orbit is used by any production API.  No import-time
experiment or output-file mutation.  The companion note proves the schema.
"""
from dataclasses import dataclass, replace
from functools import lru_cache
from itertools import product
import translation_phase_decoder_20261005 as controller
import collatz_recursive_dependency_kernel_20261004 as kernel

codec = kernel.codec
AWARDS = {'H': 2, 'L': 4, 'A': 1, 'B': 1, 'G': -1}
SEED_WORDS = {1: (), 3: (1, 4), 5: (4,), 7: (1, 1, 2, 3, 4),
              11: (1, 2, 3, 4), 13: (3, 4)}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def integer(x, minimum=0):
    need(type(x) is int and x >= minimum, 'exact integer in range')


def odd(x):
    integer(x, 1)
    need(x % 2 == 1, 'positive odd integer')


def literal(n, word, root=False):
    """Strict first-hit replay; the sole literal verifier used by the checker."""
    return kernel.literal(n, word, require_root=root)


@dataclass(frozen=True)
class Node:
    source: int
    kind: str
    child: object = None
    source_word: tuple = ()
    child_word: tuple = ()


def obligations(nodes, requested, seeds):
    """Validate a finite dependency DAG and expose every assumed constant.

    A join proves Root(source) from Root(child), by an independently checked
    common future.  A smaller join additionally checks the natural-number rank.
    Empty assumptions and cycles never create a ROOT proof.
    """
    need(type(nodes) is dict and type(requested) is tuple and requested,
         'finite proof graph and nonempty request tuple')
    need(type(seeds) is frozenset, 'explicit finite assumption set')
    for seed in seeds:
        odd(seed)
        need(seed != 1, 'ROOT is checked, not an assumption')
    need(all(type(key) is str and key for key in nodes), 'named proof nodes')
    active, memo = set(), {}

    def visit(key):
        need(type(key) is str and key in nodes, 'existing named dependency')
        need(key not in active, 'circular dependencies are not proofs')
        if key in memo:
            return memo[key]
        active.add(key)
        node = nodes[key]
        need(type(node) is Node, 'typed proof node')
        odd(node.source)
        kernel.word_guard(node.source_word)
        kernel.word_guard(node.child_word)
        if node.kind in ('root', 'seed'):
            need(node.child is None and not node.source_word and not node.child_word,
                 'leaf contains no hidden dependency or receipt')
            if node.kind == 'root':
                need(node.source == 1, 'the checked ROOT leaf is one')
                result = frozenset()
            else:
                need(node.source in seeds, 'assumption explicitly declared')
                result = frozenset((node.source,))
        else:
            need(node.kind in ('actual', 'join', 'smaller'), 'known proof rule')
            result = visit(node.child)
            child = nodes[node.child]
            need(node.source_word or node.child_word, 'nontrivial checked receipt')
            if node.kind == 'actual':
                need(not node.child_word, 'actual dependency has empty child word')
            if node.kind == 'smaller':
                need(child.source < node.source, 'strict well-founded dependency')
            need(literal(node.source, node.source_word) ==
                 literal(child.source, node.child_word), 'same exact common future')
        active.remove(key)
        memo[key] = result
        return result

    return frozenset().union(*(visit(key) for key in requested))


def discharge(nodes, requested, seeds, supplied):
    """Replace assumed leaves by independently checked supplied ROOT words."""
    needed = obligations(nodes, requested, seeds)
    need(type(supplied) is dict and needed <= supplied.keys(), 'every seed obligation supplied')
    for seed in supplied:
        odd(seed)
    for seed in needed:
        literal(seed, supplied[seed], root=True)
    memo = {}

    def compile_node(key):
        if key in memo:
            return memo[key]
        node = nodes[key]
        if node.kind == 'root':
            word = ()
        elif node.kind == 'seed':
            word = supplied[node.source]
        else:
            child_route = compile_node(node.child)
            cut = len(node.child_word)
            need(child_route[:cut] == node.child_word, 'actual child-prefix agreement')
            word = node.source_word + child_route[cut:]
        literal(node.source, word, root=True)
        memo[key] = word
        return word

    result = {key: compile_node(key) for key in requested}
    for key, word in result.items():
        cert = kernel.codec_from_word(word)
        need(codec.expand(cert, codec.bit_bounds(cert)[1]) == nodes[key].source,
             'independent inverse certificate retains source')
    return result


@dataclass(frozen=True)
class Family:
    word: str
    seed: int


def family_guard(family):
    need(type(family) is Family, 'typed family descriptor')
    controller.symbols(family.word)
    need(family.word and not family.word.endswith('L') and 'LB' not in family.word,
         'nonempty native word with a non-L terminal')
    odd(family.seed)
    need(family.seed % 3 != 0, 'seed has an odd predecessor')
    credit = 0
    for letter in family.word:
        credit += AWARDS[letter]
        need(credit >= 0, 'conservative funding at every prefix')


@lru_cache(None)
def _parameters(word, seed):
    P, Q, B = controller.carrier(word)
    R = sum(controller.LETTERS[c][0] for c in word)
    A = Q.bit_length()-1
    modulus = 3*P
    epsilon = 0 if seed % 3 == 1 else 1
    target = ((1+3*B*pow(Q, -1, modulus)) *
              pow((2**epsilon)*seed, -1, modulus)) % modulus
    phase, phase_period = kernel.principal_log4(target, R+1)
    need(phase_period == P, 'full principal-unit exponent period')
    a0, period = epsilon+2*phase, 2*P
    # A sufficient positive-source cutoff, avoiding expansion of a huge 2^a.
    cutoff = max(4, max(1, Q+3*B).bit_length())
    if a0 < cutoff:
        a0 += ((cutoff-a0+period-1)//period)*period
    need((Q*seed*pow(2, a0, modulus)-Q-3*B) % modulus == 0,
         'unique ternary exponent phase')
    return P, Q, B, R, A, a0, period


def parameters(family):
    family_guard(family)
    return _parameters(family.word, family.seed)


def family_obligations(family):
    parameters(family)
    return frozenset() if family.seed == 1 else frozenset((family.seed,))


def exponent(family, parameter):
    integer(parameter)
    *_, a0, period = parameters(family)
    return a0+period*parameter


def source_residue(family, parameter, modulus):
    integer(modulus, 1)
    P, Q, B, *_ = parameters(family)
    a = exponent(family, parameter)
    denominator = 3*P
    numerator = (Q*family.seed*pow(2, a, denominator*modulus)-Q-3*B) % (denominator*modulus)
    need(numerator % denominator == 0, 'retained division precision')
    return numerator//denominator


def source(family, parameter, bit_cap=20000):
    integer(bit_cap, 1)
    P, Q, B, _, A, *_ = parameters(family)
    a = exponent(family, parameter)
    need(a+A+family.seed.bit_length()+2 <= bit_cap, 'declared expansion bit cap')
    numerator = Q*family.seed*(1 << a)-Q-3*B
    need(numerator % (3*P) == 0, 'expanded exact source division')
    n = numerator//(3*P)
    odd(n)
    return n


def source_word(family, parameter):
    """Symbolic route to the still-assumed seed; never observes its orbit."""
    P, Q, B, R, A, *_ = parameters(family)
    a = exponent(family, parameter)
    route = (a,)
    for letter in reversed(family.word):
        if letter in controller.ACTUAL:
            route = controller.ACTUAL[letter]+route
        elif letter == 'H':
            need(route, 'H has its terminal edge')
            route = kernel.V+(route[0]+2,)+route[1:]
        else:
            need(len(route) >= 3 and route[:2] == (1, 2), 'L child prefix retained')
            route = (1, 2, 1, 1, route[2]+2)+route[3:]
    need(len(route) == R+1 and sum(route) == A+a, 'relative odd and halving costs')
    # The only potentially enormous exponent stays last.  Its affine carry
    # is checked without ever constructing the huge full denominator.
    carry = cost = 0
    for i, valuation in enumerate(route):
        need(cost <= A+2, 'finite prefix carry before the large final exponent')
        carry = 3*carry+(1 << cost)
        if i+1 < len(route):
            cost += valuation
    need(carry == Q+3*B, 'exact source-to-seed affine identity')
    return route


def complete_family(family, parameter, supplied_seed_word):
    """Discharge one constant obligation; return a compact inherited ROOT AST."""
    parameters(family)
    literal(family.seed, supplied_seed_word, root=True)
    route = source_word(family, parameter)+supplied_seed_word
    cert = kernel.codec_from_word(route)
    R = parameters(family)[3]
    A = parameters(family)[4]
    a = exponent(family, parameter)
    need(codec.ranks(cert) == (R+1+len(supplied_seed_word),
                              R+1+len(supplied_seed_word)+A+a+sum(supplied_seed_word)),
         'completed exact first-hit ranks')
    return cert


def rejected(call):
    try:
        call()
    except (ValueError, TypeError, KeyError):
        return 1
    raise ValueError('hostile accepted')


def experiment():
    for seed, word in SEED_WORDS.items():
        literal(seed, word, root=True)
    nodes = {'s3': Node(3, 'seed'), 's5': Node(5, 'seed'),
             'n7': Node(7, 'smaller', 's3', (1, 1, 2, 3), (1,)),
             'n9': Node(9, 'actual', 'n7', (2,)),
             'n11': Node(11, 'actual', 's5', (1, 2, 3))}
    seeds, requests = frozenset((3, 5)), ('n7', 'n9', 'n11')
    need(obligations(nodes, requests, seeds) == seeds, 'two exact obligations, reused DAG')
    routes = discharge(nodes, requests, seeds, {3: SEED_WORDS[3], 5: SEED_WORDS[5]})
    need(routes['n7'] == SEED_WORDS[7] and routes['n11'] == SEED_WORDS[11],
         'independent checked seed discharge')
    cycle = {'x': Node(7, 'join', 'y', (1, 1, 2, 3), (1,)),
             'y': Node(3, 'join', 'x', (1,), (1, 1, 2, 3))}
    hostiles = rejected(lambda: obligations(cycle, ('x',), frozenset()))
    hostiles += rejected(lambda: discharge(nodes, requests, seeds, {3: SEED_WORDS[3]}))
    hostiles += rejected(lambda: discharge(nodes, requests, seeds, {3: (4,), 5: (4,)}))
    hostiles += rejected(lambda: discharge(nodes, requests, seeds, {3.0: (1, 4), 5: (4,)}))
    hostiles += rejected(lambda: obligations(nodes, requests, frozenset((5,))))
    forged = dict(nodes, n7=replace(nodes['n7'], source_word=(1, 2, 3)))
    hostiles += rejected(lambda: obligations(forged, requests, seeds))
    hostiles += rejected(lambda: obligations({'x': Node(True, 'root')}, ('x',), frozenset()))
    hostiles += rejected(lambda: discharge({'x': Node(5, 'seed')}, ('x',),
                                           frozenset((5,)), {5: (4, 2)}))

    programs = []
    for length in range(1, 4):
        for letters in product(AWARDS, repeat=length):
            word = ''.join(letters)
            try:
                family_guard(Family(word, 1))
            except ValueError:
                continue
            programs.append(word)
    plans = literal_count = phase_checks = codecs = 0
    max_a = max_depth = 0
    for word in programs:
        for seed in (1, 5, 7, 11, 13):
            family = Family(word, seed)
            need(family_obligations(family) == (frozenset() if seed == 1 else frozenset((seed,))),
                 'one finite obligation independently of program length or parameter')
            P, Q, B, R, A, a0, period = parameters(family)
            for t in (0, 1):
                a = exponent(family, t)
                need(a >= 4 and period == 2*3**R, 'all-height phase parameters')
                need((Q*seed*pow(2, a, 3*P)-Q-3*B) % (3*P) == 0,
                     'independent modular integrality')
                route = source_word(family, t)
                cert = complete_family(family, t, SEED_WORDS[seed])
                for depth in (1, 2, 5, 11, 31):
                    need(codec.mod2(cert, depth) == source_residue(family, t, 2**depth),
                         'independent binary source reader')
                    need(codec.mod3(cert, depth) == source_residue(family, t, 3**depth),
                         'independent ternary source reader')
                    phase_checks += 2
                guard = controller.native_guard(controller.encode(word))
                need(source_residue(family, t, guard[1]) == guard[0], 'entire native cylinder')
                if a+A+seed.bit_length()+2 <= 20000:
                    n = source(family, t)
                    need(literal(n, route) == seed, 'literal source-to-seed receipt')
                    literal(n, route+SEED_WORDS[seed], root=True)
                    _, states, paid = controller.recover_source(n, controller.encode(word), True)
                    y = (seed*(1 << a)-1)//3
                    need(paid and states[-1] == y < n, 'immutable-source payment')
                    need(codec.expand(cert, codec.bit_bounds(cert)[1]) == n,
                         'literal and inverse AST identity')
                    literal_count += 1
                plans += 1
                codecs += 1
                max_a = max(max_a, a)
                max_depth = max(max_depth, R+1)

    for bad in (Family('L', 5), Family('LB', 5), Family('G', 5),
                Family('LG', 3), Family('H', True), Family('H', 1.0), Family('', 1)):
        hostiles += rejected(lambda bad=bad: parameters(bad))
    good = Family('LG', 5)
    hostiles += rejected(lambda: source_word(good, -1))
    hostiles += rejected(lambda: complete_family(good, 0, SEED_WORDS[7]))
    hostiles += rejected(lambda: source(good, 100, bit_cap=10))
    hostiles += rejected(lambda: source_residue(good, 0, False))
    # A finite assumption is not secretly discharged by a cycle or by itself.
    need(obligations({'x': Node(7, 'seed')}, ('x',), frozenset((7,))) == frozenset((7,)),
         'self-assumption remains an explicit obligation')

    print('FINITE-EXACT: six explicitly supplied seed words checked; no orbit search used.')
    print('Conditional DAG requests7,9,11 expose exactly{3,5}; both independent seed receipts discharge all three.')
    print('Conservative funded native programs length1..3 ending non-L:', len(programs))
    print('Family universe: seeds1,5,7,11,13; parameters0,1:', plans, 'plans;', codecs, 'completed inherited ASTs.')
    print('Independent source residue comparisons:', phase_checks, '; expanded literal controls under20000bits:', literal_count)
    print('Maximum symbolic terminal exponent:', max_a, '; maximum source-to-seed odd depth:', max_depth)
    for seed in (1, 5, 7):
        fam = Family('LG', seed)
        *_, a0, period = parameters(fam)
        print('LG seed', seed, ': a=', a0, '+', period, '*t; n=(128*seed*2^a-287)/243; route1211(a+2).')
    print('Typed, malformed, circular, missing/wrong seed, root-padding and budget hostiles rejected:', hostiles)
    print('PROVED: each fixed admissible program/3-unit seed has infinitely many constructed sources, with one constant ROOT obligation.')
    print('No arbitrary-source coverage claim. Finite seed constants and universal parameterized schemas are distinct obligations.')


if __name__ == '__main__':
    experiment()
