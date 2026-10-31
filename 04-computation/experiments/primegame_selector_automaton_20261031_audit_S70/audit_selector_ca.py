"""Independent audit (S70 audit) of THM-4852 (3): FRACTRAN with numerators/denominators dividing a base B as a
digit-local ("selector") update on base-B digits.  Written from scratch; the true map is always computed with
Python big integers and the digit formulas are tested against it.

  T1  digit formulas: y = x/q has y_i = floor(x_i/q) + (x_{i+1} mod q) B/q; z = p y has z_i = (p y_i mod B) +
      floor(p y_{i-1}/B), with every z_i < B ("no further carry"); random B (primorial and non-primorial),
      random p, q | B, random multi-digit x with q | x.
  T2  locality of the TRUE map x -> p x / q: perturbing digits at distance >= 2 from position i never changes
      digit i of p x / q (radius 1), and radius 0 is not enough in general.
  T3  selector: the first applicable fraction is determined by the units digit x_0 alone.
  T4  random squarefree FRACTRAN programs (reduced fractions over random prime sets): the digit update with the rule
      chosen from x_0 reproduces the big-integer orbit; PRIMEGAME (both lists) from large random 29-smooth starts.
  T5  the composite map is NOT a cellular automaton in the strict sense: it does not commute with the shift
      x -> B x, and digit i of F(x) depends on x_0 at unbounded distance i.
  T6  caveats: (a) fractions must be in lowest terms (an unreduced squarefree 6/10 behaves differently);
      (b) squarefreeness is not needed: any p, q | B work (e.g. B = (30)^2, p/q = 4/9).
"""
import random
from math import gcd, prod

rng = random.Random(20261031)


def digits(x, B):
    d = []
    while x:
        d.append(x % B)
        x //= B
    return d or [0]


def undigits(d, B):
    return sum(v * B ** i for i, v in enumerate(d))


def local_step(d, p, q, B):
    """digit formulas of the note: divide by q (right neighbour), multiply by p (left neighbour)."""
    x = d + [0, 0]
    n = len(d) + 2
    y = [x[i] // q + (x[i + 1] % q) * (B // q) if i + 1 < n else x[i] // q for i in range(n)]
    z = []
    for i in range(n):
        v = (p * y[i]) % B + ((p * y[i - 1]) // B if i > 0 else 0)
        assert 0 <= v < B, "further carry"
        z.append(v)
    assert all(v < B for v in y)
    while len(z) > 1 and z[-1] == 0:
        z.pop()
    return z


def rand_divisor(primes):
    return prod(p for p in primes if rng.random() < 0.5)


def say(*a):
    print(*a, flush=True)


PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]

# ---------------------------------------------------------------- T1
say('T1  digit formulas (random B, p, q | B, multi-digit x)')
n_tests = 0
for trial in range(4000):
    S = sorted(rng.sample(PRIMES, rng.randint(2, 10)))
    B = prod(S)
    p, q = rand_divisor(S), rand_divisor(S)
    g = gcd(p, q)
    p, q = p // g, q // g
    x = rng.randrange(1, B ** rng.randint(1, 25))
    x = x * q  # make q | x
    d = digits(x, B)
    z = local_step(d, p, q, B)
    assert undigits(z, B) == p * x // q, (B, p, q, x)
    # carry bound used in the proof
    y = digits(x // q, B)
    for yi in y:
        assert (p * yi) % B <= B - p and (p * yi) // B <= p - 1
    n_tests += 1
say(f'    {n_tests} random cases: formulas reproduce p*x/q exactly; every z_i < B; p*y_i mod B <= B-p, carry <= p-1')

# ---------------------------------------------------------------- T2
say('T2  locality of the true map x -> p x / q')
viol = 0
radius0_fail = 0
checks = 0
for trial in range(3000):
    S = sorted(rng.sample(PRIMES, rng.randint(2, 8)))
    B = prod(S)
    p, q = rand_divisor(S), rand_divisor(S)
    g = gcd(p, q)
    p, q = p // g, q // g
    L = rng.randint(6, 14)
    d = [rng.randrange(B) for _ in range(L)]
    d[0] -= d[0] % q  # q | x_0 <=> q | x
    x = undigits(d, B)
    z = digits(p * x // q, B) + [0] * 3
    i = rng.randint(2, L - 3)
    # perturb every digit at distance >= 2 from i (keep x_0 mod q = 0 by perturbing only multiples of q at 0)
    d2 = list(d)
    for j in range(L):
        if abs(j - i) >= 2:
            d2[j] = rng.randrange(B) if j > 0 else rng.randrange(0, B, q)
    z2 = digits(p * undigits(d2, B) // q, B) + [0] * (L + 3)
    checks += 1
    if z2[i] != z[i]:
        viol += 1
    # radius 0: perturb only the neighbours i-1, i+1
    d3 = list(d)
    d3[i - 1] = rng.randrange(B)
    d3[i + 1] = rng.randrange(B)
    z3 = digits(p * undigits(d3, B) // q, B) + [0] * (L + 3)
    if z3[i] != z[i]:
        radius0_fail += 1
say(f'    {checks} perturbation tests at distance >= 2: {viol} changes of the digit (radius 1 holds)')
say(f'    perturbing only the two neighbours changed the digit in {radius0_fail} of {checks} cases (radius 0 fails)')
assert viol == 0 and radius0_fail > 0

# ---------------------------------------------------------------- T3, T4
OWNER = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11),
         (15, 2), (1, 7), (55, 1)]
GUY = OWNER[:11] + [(15, 14), (15, 2), (55, 1)]


def run_both(prog, B, x, steps):
    d = digits(x, B)
    for t in range(steps):
        ks = [k for k, (p, q) in enumerate(prog) if x % q == 0]
        kd = [k for k, (p, q) in enumerate(prog) if d[0] % q == 0]
        if not ks:
            assert not kd
            return t, True
        assert ks[0] == kd[0]
        p, q = prog[ks[0]]
        x = x * p // q
        d = local_step(d, p, q, B)
        assert undigits(d, B) == x
    return steps, False


say('T3/T4  random squarefree FRACTRAN programs (reduced fractions), rule chosen by the units digit')
tot_steps = 0
halts = 0
for trial in range(400):
    S = sorted(rng.sample(PRIMES[:12], rng.randint(3, 7)))
    B = prod(S)
    prog = []
    for _ in range(rng.randint(3, 12)):
        p, q = rand_divisor(S), rand_divisor(S)
        g = gcd(p, q)
        prog.append((p // g, q // g))
    x = prod(rng.choice(S) ** rng.randint(0, 6) for _ in range(6)) * rng.choice([1, B, B * B])
    t, h = run_both(prog, B, x, 400)
    tot_steps += t
    halts += h
say(f'    400 random programs, {tot_steps} steps (halting runs: {halts}): digit update == big-integer orbit, selector == x_0')

B29 = prod(PRIMES[:10])
assert B29 == 6469693230
for name, prog in (('owner/Conway 1987', OWNER), ('Guy 1983', GUY)):
    tot = 0
    for trial in range(20):
        x = prod(rng.choice(PRIMES[:10]) ** rng.randint(0, 40) for _ in range(12))
        t, h = run_both(prog, B29, x, 3000)
        tot += t
    t, h = run_both(prog, B29, 2, 20000)
    say(f'    PRIMEGAME {name}: 20 large random 29-smooth starts x 3000 steps ({tot} steps) + 20,000 steps from 2: match')

# ---------------------------------------------------------------- T5
say('T5  the composite map is not a CA in the strict sense')


def F(x, prog):
    for p, q in prog:
        if x % q == 0:
            return x * p // q
    return None


x = 2
lhs, rhs = F(B29 * x, OWNER), B29 * F(x, OWNER)
say(f'    shift test: F(B*2) = {lhs} (rule 17/91) but B*F(2) = {rhs} (rule 15/2): equal? {lhs == rhs}')
assert lhs != rhs
# x and x' differ only in the units digit; digit 6 of the image differs (rule 15/2 versus rule 55/1)
hi = 123456789 * B29 ** 6
xa, xb = hi + 2, hi + 1
za, zb = digits(F(xa, OWNER), B29), digits(F(xb, OWNER), B29)
say(f'    x = c*B^6 + 2 vs c*B^6 + 1 (digits differ only at position 0): image digit 6 = {za[6]} vs {zb[6]}')
assert za[6] != zb[6]

# ---------------------------------------------------------------- T6
say('T6  caveats')
# (a) unreduced 6/10: FRACTRAN multiplies by the rational 6/10 = 3/5 whenever the product is integral (5 | x)
x = 5 * 7
fr_val = x * 6 // 10 if (x * 6) % 10 == 0 else None
ca_sel = (digits(x, 30)[0] % 10 == 0)
say(f'    6/10 on x = 35: FRACTRAN product integral -> {fr_val}; unreduced selector "10 | x_0" -> {ca_sel}  (needs lowest terms)')
assert fr_val == 21 and not ca_sel
# (b) non-squarefree p, q | B
cnt = 0
DIVS = {B: [d for d in range(1, B + 1) if B % d == 0] for B in (900, 30 ** 3, 2 ** 3 * 3 ** 2 * 5 * 7, (2 * 3 * 5 * 7) ** 2)}
for trial in range(2000):
    B = rng.choice(sorted(DIVS))
    divs = DIVS[B]
    p, q = rng.choice(divs), rng.choice(divs)
    g = gcd(p, q)
    p, q = p // g, q // g
    x = rng.randrange(1, B ** rng.randint(1, 12)) * q
    z = local_step(digits(x, B), p, q, B)
    assert undigits(z, B) == p * x // q
    cnt += 1
say(f'    {cnt} cases with non-squarefree p, q | B (B = 900, 27000, 2520, 44100): same formulas exact')
say('ALL CA CHECKS PASSED')
