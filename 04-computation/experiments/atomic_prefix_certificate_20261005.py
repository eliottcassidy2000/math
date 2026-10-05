"""Prefix-free input scheduling and exact atomic coverage of checked ROOTs.

Standard library only.  Normal/-O checks are explicit; no import-time run or
default file mutation.  Coverage weights belong to sources, not proof codes.
"""
from dataclasses import dataclass, field
from collections import deque
from fractions import Fraction as F
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def integer(n, lower=0):
    need(type(n) is int and n >= lower, 'exact integer in range')


def odd(n):
    integer(n, 1)
    need(n % 2 == 1, 'positive odd source')


def gamma(n):
    integer(n, 1)
    bits = bin(n)[2:]
    return '0'*(len(bits)-1)+bits


def read_gamma(bits, offset=0):
    need(type(bits) is str and all(c in '01' for c in bits), 'binary string')
    integer(offset)
    need(offset < len(bits), 'gamma start present')
    first = offset
    while first < len(bits) and bits[first] == '0':
        first += 1
    need(first < len(bits), 'gamma leading one present')
    width = first-offset+1
    end = first+width
    need(end <= len(bits), 'complete gamma payload')
    n = int(bits[first:end], 2)
    need(gamma(n) == bits[offset:end], 'canonical gamma field')
    return n, end


def atom(source):
    odd(source)
    index = (source+1)//2
    return F(1, 1 << (2*(index.bit_length()-1)+1))


def prefix_tail(count):
    """Exact mass outside the first count odd inputs; count zero is allowed."""
    integer(count)
    if count == 0:
        return F(1)
    k = count.bit_length()-1
    # Complete blocks0..k-1, then count-2^k+1 entries of block k.
    return F(1, 1 << k)-F(count-(1 << k)+1, 1 << (2*k+1))


def step(n):
    odd(n)
    need(n != 1, 'ROOT has no certificate outgoing edge')
    value = 3*n+1
    exponent = (value & -value).bit_length()-1
    return value >> exponent, exponent


def replay(source, word):
    odd(source)
    need(type(word) is tuple and all(type(a) is int and a > 0 for a in word),
         'tuple of positive exact valuations')
    current = source
    for expected in word:
        current, actual = step(current)
        need(actual == expected, 'actual valuation')
    need(current == 1, 'completed first-hit ROOT certificate')
    return True


def encode_receipt(source, word):
    replay(source, word)
    return gamma((source+1)//2)+gamma(len(word)+1)+''.join(gamma(a) for a in word)


def decode_receipt(bits):
    index, offset = read_gamma(bits)
    length_code, offset = read_gamma(bits, offset)
    word = []
    for _ in range(length_code-1):
        valuation, offset = read_gamma(bits, offset)
        word.append(valuation)
    need(offset == len(bits), 'no trailing suffix or second receipt')
    source, word = 2*index-1, tuple(word)
    replay(source, word)
    return source, word


class Ledger:
    def __init__(self):
        self.words = {}
        self.mass = F(0)
        self.verification_edges = 0
        self.submit(1, ())

    def submit(self, source, word):
        replay(source, word)
        self.verification_edges += len(word)
        if source in self.words:
            need(self.words[source] == word, 'canonical first-hit identity')
            return False
        self.words[source] = word
        self.mass += atom(source)
        need(self.mass < 1, 'finite full-support ledger has a positive tail')
        return True

    def submit_code(self, bits):
        source, word = decode_receipt(bits)
        return self.submit(source, word)

    def residual(self):
        return 1-self.mass

    def certifies_prefix_by_mass(self, count):
        integer(count, 1)
        # This ledger is finite: any missing requested atom has additional
        # unlisted positive atoms beside it. Equality therefore suffices.
        return self.residual() <= atom(2*count-1)

    def missing_prefix(self, count):
        integer(count)
        return tuple(2*m-1 for m in range(1, count+1) if 2*m-1 not in self.words)

    def prefix_residual(self, count):
        integer(count)
        return 1-prefix_tail(count)-sum((atom(n) for n in self.words if n <= 2*count-1), F(0))


@dataclass
class Task:
    source: int
    states: list = field(default_factory=list)
    word: list = field(default_factory=list)
    finished: bool = False

    def __post_init__(self):
        odd(self.source)
        self.states = [self.source]


class Scheduler:
    """Finite resumable stages; each task is a literal Collatz proof search.

    A quantum observes at most one odd edge.  It is not asserted to have
    constant bit complexity.  Prefix allocation never presumes termination.
    """
    def __init__(self, suffix_memory):
        need(type(suffix_memory) is bool, 'explicit suffix-memory switch')
        self.share = suffix_memory
        self.ledger = Ledger()
        self.tasks = {}
        self.stage = self.queries = self.quanta = self.reserved = 0
        self.observations = {}

    def quantum(self, task):
        self.quanta += 1
        known = self.ledger.words if self.share else {1: ()}
        current = task.states[-1]
        if current not in known:
            target, valuation = step(current)
            self.queries += 1
            if current in self.observations:
                need(self.observations[current] == (target, valuation), 'retained edge identity')
            self.observations[current] = (target, valuation)
            task.word.append(valuation)
            task.states.append(target)
            current = target
        if current in known:
            tail = known[current]
            route = tuple(task.word)+tail
            if self.share:
                for i, source in enumerate(task.states):
                    self.ledger.submit(source, tuple(task.word[i:])+tail)
            else:
                self.ledger.submit(task.source, route)
            task.finished = True

    def advance(self):
        self.stage += 1
        T = self.stage
        largest_index = (1 << ((T+1)//2))-1
        stage_reserved = 0
        for index in range(1, largest_index+1):
            source = 2*index-1
            task = self.tasks.setdefault(source, Task(source))
            length = len(gamma(index))
            need(length <= T, 'activated prefix-code length')
            budget = 1 << (T-length)
            stage_reserved += budget
            for _ in range(budget):
                if task.finished:
                    break
                self.quantum(task)
        self.reserved += stage_reserved
        need(stage_reserved < 1 << T, 'Kraft reservation bound')
        return self.snapshot()

    def snapshot(self):
        count = (1 << ((self.stage+1)//2))-1
        residual = self.ledger.residual()
        certified_by_mass = 0
        # This is the guarantee from the scalar residual alone, not direct lookup.
        while residual <= atom(2*(certified_by_mass+1)-1):
            certified_by_mass += 1
        return {'stage': self.stage, 'active_inputs': count,
                'active_missing': list(self.ledger.missing_prefix(count)),
                'verified_sources': len(self.ledger.words),
                'covered_mass': str(self.ledger.mass), 'residual': str(residual),
                'prefix_guaranteed_by_global_mass': certified_by_mass,
                'odd_edge_queries': self.queries, 'distinct_observed_edges': len(self.observations),
                'used_quanta': self.quanta, 'reserved_quanta': self.reserved,
                'certificate_verification_edges': self.ledger.verification_edges}


def grounded_observation_ledger(observations):
    """Inherited least ROOT closure, applied only to the retained exact edges."""
    need(type(observations) is dict, 'explicit finite observed-edge table')
    parents = {}
    for source, edge in observations.items():
        need(type(edge) is tuple and len(edge) == 2, 'typed observed-edge pair')
        target, valuation = edge
        odd(target)
        integer(valuation, 1)
        need(step(source) == edge, 'independent observed-edge verification')
        parents.setdefault(target, []).append((source, valuation))
    ledger, queue = Ledger(), deque((1,))
    while queue:
        target = queue.popleft()
        for source, valuation in parents.get(target, ()):
            if source not in ledger.words:
                ledger.submit(source, (valuation,)+ledger.words[target])
                queue.append(source)
    return ledger


def reject(call):
    try:
        call()
    except (ValueError, TypeError):
        return 1
    raise ValueError('hostile accepted')


def experiment():
    gamma_codes = [gamma(m) for m in range(1, 1024)]
    for m, code in enumerate(gamma_codes, 1):
        need(read_gamma(code) == (m, len(code)), 'exact input-code round trip')
        need(len(code) == 2*(m.bit_length()-1)+1, 'exact gamma length')
    for i, code in enumerate(gamma_codes):
        need(not any(other.startswith(code) for other in gamma_codes[i+1:]), 'prefix-free finite control')
    need(sum((atom(2*m-1) for m in range(1, 1024)), F(0)) == 1-F(1, 1024),
         'dyadic-block partial Kraft equality')
    for count in range(1024):
        need(prefix_tail(count) == 1-sum((atom(2*m-1) for m in range(1, count+1)), F(0)),
             'independent arbitrary-cut tail formula')
    need(sum((F(1, 4**(k+1))*2**k for k in range(10)), F(0)) == F(1, 2)*(1-F(1, 1024)),
         'uncorrected quarter-power weights have total one-half')

    results, final_searches = {}, {}
    receipt_count = 0
    for share in (False, True):
        search, snapshots = Scheduler(share), []
        for T in range(1, 15):
            row = search.advance()
            if T in (6, 8, 10, 12, 14):
                snapshots.append(row)
        for source, word in search.ledger.words.items():
            bits = encode_receipt(source, word)
            need(decode_receipt(bits) == (source, word), 'independent prefix receipt verification')
            before = search.ledger.mass
            need(not search.ledger.submit_code(bits) and search.ledger.mass == before,
                 'repeated proofs do not duplicate input mass')
            receipt_count += 1
        need(all(not row['active_missing'] or row['residual'] != '0' for row in snapshots),
             'finite missing tasks retain positive mass')
        for count in (1, 7, 15, 31, 63, 127):
            if search.ledger.certifies_prefix_by_mass(count):
                need(not search.ledger.missing_prefix(count), 'global mass finite-coverage certificate')
            need((search.ledger.prefix_residual(count) == 0) ==
                 (not search.ledger.missing_prefix(count)), 'localized exact finite criterion')
        results['suffix_memory' if share else 'independent_tasks'] = snapshots
        final_searches[share] = search

    need(final_searches[False].observations == final_searches[True].observations,
         'same final distinct observed-edge evidence')
    closure = grounded_observation_ledger(final_searches[True].observations)
    need(set(final_searches[True].ledger.words) <= set(closure.words), 'closure retains every proved suffix')
    closure_result = {'same_evidence_distinct_edges': len(final_searches[True].observations),
                      'rooted_sources': len(closure.words), 'covered_mass': str(closure.mass),
                      'residual': str(closure.residual()), 'active_missing': list(closure.missing_prefix(127)),
                      'postprocessing_certificate_verification_edges': closure.verification_edges}

    # A fixed unresolved integer is an atom; a density-one complement does not erase it.
    need(atom(27) == F(1, 128), '27 retains its exact positive atom')
    root_only = Ledger()
    need(root_only.residual() == atom(1) and root_only.certifies_prefix_by_mass(1),
         'finite-ledger equality boundary certifies its requested prefix')
    need(1-(1-atom(27)) == atom(27),
         'an infinite complement-of-one set has equality despite its missing atom')
    need(len(encode_receipt(1, ())) == 2 and atom(1) == F(1, 2),
         'ROOT proof-code weight1/4 differs from source weight1/2')
    for J in range(1, 10):
        threshold = atom(2*((1 << J)-1)-1)
        need(prefix_tail((1 << (2*J-1))-1) == threshold,
             'finite prefix scalar criterion succeeds at equality')
        need(prefix_tail((1 << (2*J-2))-1) > threshold,
             'one fewer dyadic depth does not meet the scalar sufficient criterion')
    hostiles = 0
    for call in (lambda: gamma(True), lambda: read_gamma('001'), lambda: read_gamma('000'),
                 lambda: decode_receipt('111'), lambda: encode_receipt(3, ()),
                 lambda: encode_receipt(1, (2,)), lambda: encode_receipt(3, (1, 3)),
                 lambda: atom(2), lambda: Ledger().submit(True, ()),
                 lambda: Scheduler(1), lambda: grounded_observation_ledger({3: (5.0, True)})):
        hostiles += reject(call)
    print(json.dumps({'status': 'PASS', 'gamma_round_trips': 1023,
                      'tail_formula_controls': 1024, 'verified_receipt_round_trips': receipt_count,
                      'hostiles_rejected': hostiles, 'experiments': results,
                      'same_evidence_grounded_closure': closure_result,
                      'scope': 'Exact atomic proof coverage and fair conditional discovery. No universal Collatz termination, no optimal universal-machine claim.'},
                     sort_keys=True, indent=2))


if __name__ == '__main__':
    experiment()
