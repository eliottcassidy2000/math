"""Five-letter translation receipts, exact native guards, and phase packets.

No import-time experiment. Run normally and with -O; no output-file mutation.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
import adaptive_credit_potential_20261004 as prior
import paid_guard_budget_20261005 as inherited


def need(ok, message):
    if not ok:
        raise ValueError(message)


def integer(x, minimum=0):
    need(type(x) is int and x >= minimum, 'exact integer in range')


# P=3^r, Q=2^a; entries are (r,a,B,native residue,native modulus).
LETTERS = {'H': (6, 10, 669, 155, 2048),
           'G': (2, 3, 5, 11, 16),
           'A': (4, 7, 85, 187, 256),
           'B': (4, 7, 73, 7, 256),
           'L': (2, 4, -3, 219, 256)}
RESIDUES = {3: 'H', 4: 'G', 2: 'A', 5: 'B', 6: 'L'}
ACTUAL = {'G': (1, 2), 'A': (1, 2, 1, 3), 'B': (1, 1, 2, 3)}


def symbols(word):
    need(type(word) is str and all(x in LETTERS for x in word), 'five-letter word')


@dataclass(frozen=True)
class Translation:
    numerator: int
    cost: int


def carrier(word):
    symbols(word)
    P = Q = 1
    B = 0
    for c in word:
        r, a, b, _, _ = LETTERS[c]
        P, Q, B = 3**r*P, 2**a*Q, 3**r*B+b*Q
    return P, Q, B


def encode(word):
    _, Q, B = carrier(word)
    return Translation(B, Q.bit_length()-1)


def validate(code):
    need(type(code) is Translation, 'translation type')
    need(type(code.numerator) is int, 'exact signed numerator')
    integer(code.cost)
    need((code.cost == 0 and code.numerator == 0) or
         (code.cost > 0 and code.numerator % 2 == 1), 'reduced dyadic code')


def decode(code):
    validate(code)
    B, a_total = code.numerator, code.cost
    reverse = []
    while a_total:
        residue = B*pow(2**a_total, -1, 9) % 9
        need(residue in RESIDUES, 'unknown last-letter residue')
        c = RESIDUES[residue]
        r, a, b, _, _ = LETTERS[c]
        need(a <= a_total, 'letter exceeds remaining denominator cost')
        numerator = B-b*2**(a_total-a)
        need(numerator % 3**r == 0, 'exact peeling divisibility')
        B, a_total = numerator//3**r, a_total-a
        need((a_total == 0 and B == 0) or (a_total > 0 and B % 2 == 1),
             'peeled prefix is reduced or empty')
        reverse.append(c)
    word = ''.join(reversed(reverse))
    need(encode(word) == code, 'exact reconstruction')
    return word


def phase_packet(code):
    decode(code)
    m = 2*code.cost//3
    modulus = 3**m
    return code.cost, code.numerator*pow(2**code.cost, -1, modulus) % modulus


def decode_phase(cost, residue):
    integer(cost)
    integer(residue)
    m = 2*cost//3
    need(residue < 3**m, 'canonical phase representative')
    original = cost, residue
    reverse = []
    while cost:
        need(m >= 2 and residue % 9 in RESIDUES, 'phase resolves next letter')
        c = RESIDUES[residue % 9]
        r, a, b, _, _ = LETTERS[c]
        need(a <= cost and r <= m, 'phase and binary budgets')
        z = (2**a*residue-b) % 3**m
        need(z % 3**r == 0, 'phase peeling divisibility')
        cost, m, residue = cost-a, m-r, z//3**r
        reverse.append(c)
    need(residue == 0, 'empty prefix phase')
    word = ''.join(reversed(reverse))
    need(phase_packet(encode(word)) == original, 'phase packet reconstruction')
    return word


def real_phase(q, word):
    """Real dyadic translation modulo one, not reduction modulo an odd prime."""
    integer(q, 1)
    need(q % 2 == 1, 'positive odd multiplier')
    need(type(word) is tuple and len(word) > 0 and
         all(type(a) is int and a >= 1 for a in word), 'nonempty positive valuation word')
    Q, B = 1, 0
    for a in word:
        Q, B = Q*2**a, q*B+Q
    return Translation(B % Q, Q.bit_length()-1)


def decode_real_phase(q, length, code):
    """Total fixed-length membership test for the real phase Phi_q modulo1."""
    integer(q, 1)
    need(q % 2 == 1, 'positive odd multiplier')
    integer(length, 1)
    validate(code)
    N, A, m = code.numerator, code.cost, length
    need(A >= m and 0 < N < 2**A, 'canonical nonzero real dyadic phase')
    word = []
    while m > 1:
        Z = (N-pow(q, m-1, 2**A)) % 2**A
        need(Z != 0, 'nonzero odd-tail residue')
        a = (Z & -Z).bit_length()-1
        need(a >= 1 and A-a >= m-1, 'positive first valuation and enough remaining cost')
        word.append(a)
        N, A, m = Z >> a, A-a, m-1
    need(N == 1 and A >= 1, 'one-letter final numerator')
    word = tuple(word)+(A,)
    need(real_phase(q, word) == code, 'exact real-phase reconstruction')
    return word


def native_guard(code):
    word = decode(code)
    if 'LB' in word:
        return None
    P, Q, B = carrier(word)
    t, target = (16, 11) if word.endswith('L') else (2, 1)
    modulus = t*Q
    residue = ((target*Q-B)*pow(P, -1, modulus)) % modulus
    return residue, modulus


def operate(n, c):
    prior.positive_odd(n)
    r, a, b, residue, modulus = LETTERS[c]
    need(n % modulus == residue, 'native source guard')
    y = (3**r*n+b)//2**a
    if c in ACTUAL:
        need(prior.actual_word(n, ACTUAL[c]) == y, 'actual macro replay')
    elif c == 'H':
        need(prior.action(n, 'H') == y, 'H common-future replay')
    else:
        need(inherited.local_l_receipt(n)[0] == y, 'L common-future replay')
    return y


def recover_source(n, code, require_paid=False):
    """Recover history and guard from its translation; source is never replaced."""
    prior.positive_odd(n)
    word = decode(code)
    guard = native_guard(code)
    need(guard is not None and n % guard[1] == guard[0], 'whole native domain')
    states = [n]
    credit = funding = 0
    paid = True
    for c in word:
        x = states[-1]
        award = -1 if c == 'G' else 2 if c == 'H' else 4 if c == 'L' else 1 if x == 7 else 2
        credit += award
        funding += c != 'G'
        paid &= credit >= 0
        y = operate(x, c)
        if paid:
            need(6 <= prior.potential(y, credit) <= (n+5)*prior.ADAPTIVE_ETA**funding,
                 'immutable-source potential')
            need(y < n and prior.root_rank(y) < prior.root_rank(n), 'immutable-source payment')
        states.append(y)
    paid = bool(word) and paid
    if require_paid:
        need(paid, 'nonempty paid program required')
    P, Q, B = carrier(word)
    need((P*n+B)//Q == states[-1], 'decoded endpoint')
    return word, tuple(states), paid


def compile_join(n, code):
    word, states, _ = recover_source(n, code)
    terminal = states[-1]
    endpoint = terminal
    route = tail = ()
    for c in reversed(word):
        required = 1 if c == 'H' else 3 if c == 'L' else 0
        while len(route) < required:
            need(endpoint != 1, 'no root padding')
            endpoint, a = prior.step(endpoint)
            route += (a,)
            tail += (a,)
        if c in ACTUAL:
            route = ACTUAL[c]+route
        elif c == 'H':
            route = prior.HEAD+(route[0]+2,)+route[1:]
        else:
            need(route[:2] == (1, 2), 'L required child prefix')
            route = (1, 2, 1, 1, route[2]+2)+route[3:]
    need(prior.actual_word(n, route) == endpoint, 'source first-hit-safe receipt')
    need(prior.actual_word(terminal, tail) == endpoint, 'terminal first-hit-safe receipt')
    need(sum(route)-sum(tail) == code.cost, 'denominator is relative halving cost')
    need(len(route)-len(tail) == sum(LETTERS[c][0] for c in word), 'numerator is relative odd depth')
    return route, tail


def rejected(call):
    try:
        call()
    except (ValueError, TypeError):
        return 1
    raise ValueError('hostile accepted')


def real_phase_controls():
    count = moments = 0
    for q in (1, 3, 5, 7, 9, 11):
        for m in range(1, 6):
            phases = {}
            for word in product(range(1, 5), repeat=m):
                code = real_phase(q, word)
                need(decode_real_phase(q, m, code) == word, 'fixed-length real-phase decoder')
                phase = F(code.numerator, 2**code.cost)
                need(phase not in phases, 'fixed-length real-phase injectivity')
                phases[phase] = F(1, 2**sum(word))
                count += 1
            need(all((phase+F(1, 2)) % 1 not in phases for phase in phases),
                 'no fixed-length antipodal pair')
            energy = sum((weight**2 for weight in phases.values()), F(0))
            need(energy == F(85, 256)**m, 'finite truncated second moment')
            moments += 1
    collision = Translation(1, 2)
    need(real_phase(3, (2,)) == real_phase(3, (1, 1)) == collision,
         'length is necessary even when dyadic cost is retained')
    need(decode_real_phase(3, 1, collision) == (2,) and
         decode_real_phase(3, 2, collision) == (1, 1), 'cross-length decoder distinction')
    phase_pairs = []
    for w in ((1,), (3,), (1, 2), (1, 8)):
        code = real_phase(3, w)
        real = F(code.numerator, 2**code.cost)
        P = 3**len(w)
        # Re-encode the full actual carry; its integer part was discarded by real_phase.
        Q, B = 1, 0
        for a in w:
            Q, B = Q*2**a, 3*B+Q
        modular = F(B*pow(Q, -1, P) % P, P)
        phase_pairs.append((w, str(real), str(modular)))
    need(phase_pairs == [((1,), '1/2', '2/3'), ((3,), '1/8', '2/3'),
                         ((1, 2), '5/8', '4/9'), ((1, 8), '5/512', '4/9')],
         'real phase and modular target phase have different meanings')
    bad = sum((rejected(lambda: decode_real_phase(2, 1, Translation(1, 1))),
               rejected(lambda: decode_real_phase(3, True, Translation(1, 1))),
               rejected(lambda: decode_real_phase(3, 2, Translation(3, 2))),
               rejected(lambda: decode_real_phase(3, 3, Translation(1, 2))),
               rejected(lambda: decode_real_phase(3, 1, Translation(3, 2))),
               rejected(lambda: decode_real_phase(3, 1, Translation(5, 2)))))
    return count, moments, phase_pairs, bad


def experiment():
    for c, (r, a, b, _, _) in LETTERS.items():
        need(3**r % 9 == 0 and b*pow(2**a, -1, 9) % 9 in RESIDUES,
             'terminal residue mechanism')
        need(RESIDUES[b*pow(2**a, -1, 9) % 9] == c, 'terminal residue identity')
    seen = {}
    formal = legal = empty = phase_checks = 0
    inherited_ops = {c: inherited.Guarded(3**r, 2**a, b, s, m)
                     for c, (r, a, b, s, m) in LETTERS.items()}
    for length in range(7):
        for letters in product(LETTERS, repeat=length):
            w = ''.join(letters)
            code = encode(w)
            need(decode(code) == w, 'translation-only round trip')
            need(code not in seen, 'finite no-collision control')
            seen[code] = w
            formal += 1
            packet = phase_packet(code)
            need(decode_phase(*packet) == w, 'bounded phase packet round trip')
            phase_checks += 1
            op = inherited.I
            try:
                for c in w:
                    op = inherited.compose(op, inherited_ops[c])
            except ValueError:
                need('LB' in w and native_guard(code) is None, 'only forbidden adjacency')
                empty += 1
                continue
            need(native_guard(code) == (op.residue, op.modulus), 'independent native composition')
            legal += 1
    replayed = receipts = funded = 0
    for length in range(1, 5):
        for letters in product(LETTERS, repeat=length):
            w = ''.join(letters)
            code = encode(w)
            guard = native_guard(code)
            if guard is None:
                continue
            for lift in (0, 1, 7):
                n = guard[0]+lift*guard[1]
                _, _, paid = recover_source(n, code)
                compile_join(n, code)
                replayed += 1
                receipts += 1
                funded += paid
    # Canonical actual-word source/target pairs and exact Fourier residues.
    pairs = []
    for w in ((1,), (3,), (1, 2), (1, 8)):
        P, Q, B = prior.IDENTITY
        for a in w:
            P, Q, B = 3*P, Q*2**a, 3*B+Q
        c = ((Q-B)*pow(P, -1, 2*Q)) % (2*Q)
        d = prior.actual_word(c, w)
        need(0 < c < 2*Q and 0 < d < 2*P, 'canonical positive lifts')
        need(d % P == B*pow(Q, -1, P) % P, 'target Fourier summand')
        pairs.append((w, c, d))
    need(pairs == [((1,),3,5),((3,),13,5),((1,2),11,13),((1,8),739,13)], 'same-target payment hostile')
    need(F(1, 2)/3 != F(pow(2, -1, 3), 3),
         'rational division is not modular reduction before Fourier exponentiation')
    # Fixed precision cannot retain arbitrarily long history, even at equal cost.
    aliases = 0
    for m in range(1, 13):
        s = (m+1)//2
        t = (s+1)//2
        u, v = 'H'*t+'HG'+'G'*s, 'H'*t+'GH'+'G'*s
        a, b = encode(u), encode(v)
        modulus = 3**m
        need(a.cost == b.cost and a != b, 'distinct equal-clock words')
        need((a.numerator-b.numerator) % modulus == 0, 'same finite phase')
        for w in (u, v):
            c, _ = native_guard(encode(w))
            recover_source(c, encode(w), require_paid=True)
        aliases += 1
    need(encode('LG') == Translation(53, 7), 'K is LG translation collision')
    need(inherited.compose(inherited.L, inherited.G) == inherited.K, 'K is LG guarded collision')
    need(decode(Translation(53, 7)) == 'LG', 'five-letter convention remains unique')
    hostiles = sum((rejected(lambda: decode(Translation(True, 3))),
                    rejected(lambda: decode(Translation(5, 3.0))),
                    rejected(lambda: decode(Translation(10, 4))),
                    rejected(lambda: decode(Translation(1, 1))),
                    rejected(lambda: decode(Translation(7, 3))),
                    rejected(lambda: recover_source(27, encode('L'))),
                    rejected(lambda: recover_source(219, encode('LB'))),
                    rejected(lambda: decode_phase(3, 3**2)),
                    rejected(lambda: recover_source(1, encode(''), require_paid=True))))
    real_count, moment_count, real_pairs, real_bad = real_phase_controls()
    print('Alphabet H/G/A/B/L translations modulo9:3/4/2/5/6; exact translation recovers formal word.')
    print('All words length0..6:', formal, 'translation round trips;', phase_checks, 'bounded-phase round trips.')
    print('Native domains:', legal, '; empty exactly when LB occurs:', empty)
    print('Native guard modulus2Q unless endingL, then16Q; full lost-three-bit sidecar reconstructed.')
    print('Live words length1..4,lifts0,1,7:', replayed, 'source replays;', receipts, 'actual join receipts;', funded, 'paid instances.')
    print('Phase packet(cost A, translation modulo3^floor(2A/3)) reconstructs the same exact word.')
    print('Canonical-target hostiles:', pairs)
    print('Fixed-phase equal-clock aliases m1..12:', aliases, '; native funded positive controls retained.')
    print('Adding sixth letter K destroys word injectivity:K=LG as complete guarded maps.')
    print('Typed/reduced/precision/source/empty-domain hostiles rejected:', hostiles)
    print('Incoming real dyadic phase: q1,3,5,7,9,11; lengths1..5; valuations1..4:', real_count, 'fixed-length round trips.')
    print('No antipodal pairs; exact truncated second moments:', moment_count, '; cross-length collision(2)/(1,1) at q3 retained.')
    print('Real phase versus modular target phase:', real_pairs)
    print('Real-phase type/length/nonzero-tail/base/normalization hostiles rejected:', real_bad)
    print('Compression recovers history and native guard; coverage and arbitrary-child ROOT remain OPEN.')


if __name__ == '__main__':
    experiment()
