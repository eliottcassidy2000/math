#!/usr/bin/env python3
"""collatz_sqdbl_20260925_foundry.py -- the squares/doubles transport foundry.

Session collatz-squares-doubles-20260925 (opus). Owner's seed: "multiplication
relates to the squares the way addition relates to the doubles" -- the diagonal
of a binary operation (x+x, x*x) determines the operation (polarization), and
the exponential carries (add, doubles, halving, 3n+1) to (multiply, squares,
square root, q*X^3).

This generator enumerates TRANSPORT x INGREDIENT cells, types each cell by a
rule table (what the transport preserves, whether it can see the 3n-1 SHEET,
the 5n+1 DRIFT, planted DEFECTs, or is a density-zero object), and prints the
curated cards (statement, status, cheapest decisive test) for the non-empty
cells. Statuses are those proved or computed in the companion note
05-knowledge/results/collatz_sqdbl_20260925_squares_doubles_foundry.md and in
collatz_sqdbl_20260925_probes.py. The typing is a modelling judgement, not a
theorem; PROVED/FINITE-EXACT labels refer to the companion note.
"""
import sys
from collections import Counter

# ------------------------------------------------------------------ grammar
TRANSPORTS = {
    # key: (name, kind, ambient-structure-gained, sees SHEET?, sees DRIFT?, sees DEFECT?)
    'EXP':    ('exponential conjugation n -> q^n (doubles->squares, halving->sqrt, 3n+1 -> q X^3)', 'lossless', None, 'inherits', 'inherits', 'inherits'),
    'DLOG':   ('discrete-log realisation on a cyclic 2-group (Fermat primes p = 2^k + 1)', 'lossless-residue', 'field addition of F_p', 'no (residue)', 'no (residue)', 'no (residue)'),
    'MONOID': ('extension of EXP to the multiplicative monoid N (exponent vectors, gate = all even)', 'embedding', 'other prime coordinates', 'inherits', 'inherits', 'inherits'),
    'POLAR':  ('polarization: binary operation n o m = (nm+1)/2^v with slice m=3 and diagonal (n^2+1)/2', 'lossless', 'the whole slice family qn+1', 'inherits', 'yes (slice index)', 'inherits'),
    'QSQ':    ('quarter-square identity 3n = ((n+3)/2)^2 - ((n-3)/2)^2 (odd n), Q(x) = floor(x^2/4)', 'lossless', None, 'inherits', 'inherits', 'inherits'),
    'DIFFSQ': ('difference of squares of the clock: 2^(2a) - 3^(2b) = (2^a - 3^b)(2^a + 3^b)', 'lossy (clock parity)', 'divisor 2^a + 3^b', 'yes (sign of 2^K-3^L)', 'yes (base 3 vs 4)', 'no'),
    'QCHAR':  ('quadratic characters mod 8 (chi_{-4}, chi_8) as strategies and as valuation bits', 'residue', None, 'no (word function)', 'no', 'no'),
    'JACOBI': ('Jacobi symbols (3/m), (m/3) along Syracuse orbits, quadratic reciprocity', 'residue', 'reciprocity sign', 'no (word function)', 'yes (3 vs 5 mod 4 changes the law shape)', 'no'),
    'LOCSQ':  ('local squares: 2-adic (m = 1 mod 8) and 3-adic (m = 1 mod 3) square units', 'residue', None, 'no (word function)', 'no', 'no'),
    'GLOBSQ': ('global squares s^2 in the orbit; figurate numbers; Pell x^2 - 3y^2 = 1', 'sparse', 'positivity', 'yes (squares are positive)', 'no', 'no (density zero)'),
    'HILB':   ('Hilbert symbols / product formula: sign(2^K - 3^L) as a product of local symbols', 'sign transfer', 'finite places', 'yes (moves the real place)', 'no', 'no'),
    'REALVAL': ('real value R(d) of the Bernstein series vs its 2-adic value; orbit product prod(1 +- 1/(3m_l))', 'two-place', 'archimedean size', 'yes (the +- in the product)', 'yes (3 -> q in the product)', 'yes (orbitwise)'),
    'THETA3': ('cubic theta polarization (k+1)^3 - k^3 = 3k^2 + 3k + 1: two-variable unipotent q-system', 'functional equation', 'second variable', 'n/a', 'n/a', 'n/a'),
    'GMFAIR': ('geometric-mean-fair pairings (multiplicative analogue of AM-fair, THM-4470)', 'family', None, 'no', 'no', 'no'),
    'HEIGHT': ('elliptic transport [n]G: canonical height quadratic, parallelogram law = polarization', 'lossless', 'Mordell-Weil group', 'inherits', 'inherits', 'inherits'),
    'SUMSQ':  ('sums of two/four squares as positivity witnesses (Fermat, Lagrange), norms from Z[i]', 'sparse/sign', 'positivity', 'yes (odd sums of two squares are 1 mod 4)', 'no', 'no (density zero)'),
}

INGREDIENTS = {
    'GATE':  'parity gate (is n a double?)',
    'HALF':  'halving (undo the double)',
    'AFF':   'the affine step 3n+1',
    'VAL':   'the Syracuse valuation v_2(3n+1)',
    'CYC':   'the cycle equation m(2^K - 3^L) = B and the clock (K, L)',
    'WORD':  'the parity word / Bernstein series (divergence half)',
    'SIGN':  'the sign law (3n+1 vs 3n-1 sheet)',
    'DRIFT': 'the drift log 3 < 2 log 2',
}

# ------------------------------------------------------------------ curated cards
# (transport, ingredient) -> (statement, status, verdict class, cheapest decisive test)
CARDS = {
    ('EXP', 'GATE'): ('n even <=> X = q^n is a square in <q>', 'PROVED (trivial)', 'COSMETIC', 'none needed'),
    ('EXP', 'HALF'): ('n -> n/2 becomes X -> sqrt(X); the two lifts mod 2^k are the two square roots', 'PROVED', 'COSMETIC', 'P2'),
    ('EXP', 'AFF'): ('3n+1 becomes X -> q X^3: cube then multiply by the generator; 3n-1 is X -> q^-1 X^3', 'PROVED', 'COSMETIC', 'P2'),
    ('EXP', 'CYC'): ('a cycle is X^(3^L/2^K - 1) = q^(-c): X = q^(c 2^K/(2^K - 3^L)); the clock is the exponent (3/2^v) word', 'PROVED (restatement)', 'COSMETIC', 'none'),
    ('EXP', 'SIGN'): ('the sheet is q vs q^-1 in the cube step; positivity of n is X > 1 for q > 1', 'PROVED (restatement)', 'COSMETIC', 'none'),
    ('EXP', 'DRIFT'): ('the drift is the exponent 3/4 = 3/2^2 of "cube then two square roots": 3 < 2^2', 'PROVED (restatement)', 'COSMETIC', 'none'),
    ('DLOG', 'GATE'): ('at level 2^k = p - 1, node s is even iff y = g^s is a square in F_p^*', 'PROVED', 'EXACT MODEL', 'P2'),
    ('DLOG', 'HALF'): ('the parity graph G_+ mod 2^k (THM-4474) is the +-sqrt graph y -> +-sqrt(y) (squares), y -> +-sqrt(g y^3) (non-squares); verified p = 17, 257, 65537', 'PROVED + FINITE-EXACT', 'EXACT MODEL', 'P2 (done)'),
    ('DLOG', 'CYC'): ('parity-graph cycles = square-root cycles in F_p^*; THM-4474 Theorem A reads: provable iff every sqrt-cycle has non-square density < log_3 2', 'PROVED (restatement)', 'EXACT MODEL', 'Karp on the sqrt graph = existing engine'),
    ('DLOG', 'SIGN'): ('3n-1 is the same graph with g^-1; negation x -> -x is y -> y^-1 (THM-4474 E)', 'PROVED (restatement)', 'EXACT MODEL', 'none'),
    ('DLOG', 'DRIFT'): ('the field addition of F_p is new structure invisible in Z/2^k; no Collatz predicate uses it', 'OPEN (no map)', 'AMBIENT ONLY', 'define a predicate mixing y+1 with the sqrt graph; none found'),
    ('MONOID', 'GATE'): ('"X is a square" = all exponents even: the gate becomes a conjunction over primes', 'PROVED', 'NEW OBJECT', 'P1'),
    ('MONOID', 'HALF'): ('sqrt halves every exponent at once (synchronised halving)', 'PROVED', 'NEW OBJECT', 'P1'),
    ('MONOID', 'AFF'): ('rad(X) X^3 sends every nonzero exponent to 3e+1 and keeps zero exponents zero', 'PROVED', 'NEW OBJECT', 'P1'),
    ('MONOID', 'WORD'): ('Theorem 2: the orbit of X reaches rad(X) iff X = m^e (m squarefree) and e reaches 1; unequal exponents desynchronise forever and diverge (Terras injectivity)', 'PROVED + FINITE-EXACT (X <= 10^5, both sheets)', 'NEW OBJECT', 'P1 (done)'),
    ('MONOID', 'SIGN'): ('3n-1 version: rad(X)^-1 X^3; convergent set = {m^e: e reaches 1}, cycles at e in {5,7,10} and the 17-cycle', 'PROVED + FINITE-EXACT', 'NEW OBJECT', 'P1 (done)'),
    ('MONOID', 'DRIFT'): ('equivalent to Collatz on the diagonal, provably divergent off it: the square gate is a synchronisation constraint, not a new mechanism', 'PROVED', 'NEW OBJECT / NO MECHANISM', 'none'),
    ('POLAR', 'AFF'): ('Syracuse = slice m = 3 of the commutative operation n o m = (nm+1)/2^(v_2(nm+1)); qn+1 maps are the other slices', 'PROVED (definition)', 'FAMILY VIEW', 'none'),
    ('POLAR', 'VAL'): ('the diagonal ("square") is n o n = (n^2+1)/2 with v_2(n^2+1) = 1 always; it is monotone, values 5 mod 8 after one step', 'PROVED', 'DEGENERATE (no balance)', 'none'),
    ('POLAR', 'DRIFT'): ('polarization nm = ((n+m)^2 - (n-m)^2)/4 needs addition and halving, the two Collatz operations; the slice index q is where the drift lives', 'OBSERVATION', 'FAMILY VIEW', 'none'),
    ('QSQ', 'AFF'): ('(3n+1)/2 = (Q(n+3) - Q(n-3) + 1)/2 for odd n; Q(x) = floor(x^2/4) = Turan number t(x,3)', 'PROVED (identity)', 'COSMETIC', 'none'),
    ('QSQ', 'HALF'): ('n/2 = sqrt(Q(n)) for even n', 'PROVED (identity)', 'COSMETIC', 'none'),
    ('DIFFSQ', 'CYC'): ('at a doubled clock (2a, 2b) the denominator factors; a positive cycle needs 2^a + 3^b | B (the digit sum) and 2^a > 3^b', 'PROVED (elementary)', 'SIDE CONSTRAINT', 'search for cycles at doubled clocks: none exist below the Hercher bound anyway'),
    ('DIFFSQ', 'SIGN'): ('for L even, 2^K - 3^L = 7 mod 8 when positive: never a sum of three squares (Legendre); for L odd it is 5 mod 8', 'PROVED (mod 8)', 'RESIDUE LAW', 'none'),
    ('QCHAR', 'VAL'): ('Lemma 3: v_2(3n+s) >= 2 iff s = chi_{-4}(n); >= 3 iff also chi_8(n) = -1. Higher bits are not quadratic (only four quadratic characters mod 2^k)', 'PROVED + FINITE-EXACT (n < 2^20)', 'RESIDUE LAW', 'P3 (done)'),
    ('QCHAR', 'SIGN'): ('the unique class-(i) level-2 strategy of THM-4474 is chi_{-4}; the unique class-(ii) is -chi_{-4}; Collatz and 3n-1 are the two constant characters', 'PROVED + FINITE-EXACT', 'RESIDUE LAW', 'P3 (done)'),
    ('QCHAR', 'GATE'): ('the down class of T_+ is {(-1/n) = +1}, of T_- is {(-1/n) = -1}: the sheet sign is the value of chi_{-4} on the down class', 'PROVED', 'RESIDUE LAW (word function)', 'none'),
    ('JACOBI', 'VAL'): ('Lemma 4: (3/m_i) = (-1)^(v_i + [v_(i+1) = 1]) on both sheets; for 5n+-1, (5/m_i) = (-1)^(v_i). The reciprocity sign (3 = 3 mod 4) couples adjacent valuations', 'PROVED + FINITE-EXACT (n < 2*10^5, 12 steps, 4 maps)', 'RESIDUE LAW', 'P4 (done)'),
    ('JACOBI', 'CYC'): ('over a cycle, prod_i (3/m_i) = (-1)^(K + L_1) with L_1 the number of v = 1 steps; equals (3/prod B_i) (3/D)^L, all residue data', 'PROVED', 'RESIDUE LAW', 'none'),
    ('JACOBI', 'SIGN'): ('for positive m the law is a function of the word on both sheets: sheet-blind, as the word-function theorem predicts', 'PROVED', 'RESIDUE LAW', 'none'),
    ('LOCSQ', 'VAL'): ('Lemma 5b: an odd iterate is a 2-adic square iff v_next = 2, a 3-adic square unit iff v_prev is even; a global square forces both', 'PROVED + FINITE-EXACT', 'RESIDUE LAW', 'P5 (done)'),
    ('GLOBSQ', 'AFF'): ('Lemma 5a: (2t+1)^2 -> 3t(t+1) + 1 = 6 T_t + 1; every odd square is a down point on the plus sheet and an up point on the minus sheet', 'PROVED', 'SPARSE', 'P5 (done)'),
    ('GLOBSQ', 'CYC'): ('the successor of an odd square is again a square iff 2t+1 = y_j, j odd, of x^2 - 3y^2 = 1 (s = 1, 15, 209, 2911, 40545, 564719, ...); no chain of length 3 except at 1', 'PROVED + FINITE-EXACT (s <= 2*10^6)', 'SPARSE', 'P5 (done)'),
    ('GLOBSQ', 'WORD'): ('the trunk (4^i - 1)/3 = the n with 3n+1 a square of a power of two; the squares in an orbit are a density-zero set with no closed dynamics', 'PROVED (trunk) / OBSERVATION', 'SPARSE', 'none'),
    ('GLOBSQ', 'SIGN'): ('"the orbit contains a square" is sign-aware (squares are positive) but has no measure content', 'OBSERVATION', 'SPARSE', 'none'),
    ('HILB', 'SIGN'): ('prod_v (2^K - 3^L, -1)_v = 1: 2^K > 3^L iff the odd-multiplicity primes 3 mod 4 of |2^K - 3^L| have even count when L is even and odd count when L is odd', 'PROVED (product formula)', 'SIGN TRANSFER', 'needs the factorisation of 2^K - 3^L: no gain over Baker'),
    ('HILB', 'CYC'): ('every prime of the denominator divides the digit sum B; their quadratic characters are constrained by the sign: a residue condition on B per prime, no global leverage found', 'OBSERVATION', 'SIGN TRANSFER', 'test on the four known cycles: consistent, vacuous'),
    ('REALVAL', 'WORD'): ('Proposition 6: R_L(d) = n(prod_(l<L)(1 + 1/(3m_l)) - 1) on the plus sheet and n(1 - prod(1 - 1/(3m_l))) on the minus sheet; R(d) < infinity iff sum 1/m_l < infinity; c = lim m_L 2^(d_L)/3^L = n +- R(d)', 'PROVED + FINITE-EXACT (n <= 2000, both sheets, exact rationals)', 'TWO-PLACE', 'P6 (done)'),
    ('REALVAL', 'SIGN'): ('the sheet enters the divergence half exactly once: R(d) <= n on the minus sheet for EVERY positive orbit (equality iff sum 1/m_l = infinity), no constraint on the plus sheet; the 2-adic identities are n = -R_2(d) (plus) and n = R_2(d) (minus)', 'PROVED', 'TWO-PLACE', 'P6 (done)'),
    ('REALVAL', 'DRIFT'): ('same identities with 3 -> q; on 5n+-1 the presumed divergent orbits have c > 0 (full rate): 37.374..., 47.839..., 8.528..., 10.660...', 'FINITE-EXACT (4000 odd steps)', 'TWO-PLACE', 'P6 (done)'),
    ('REALVAL', 'CYC'): ('for eventually periodic words R(d) = R_2(d) = n (the Eliahou identity 2^K/3^L = prod(1 + b/(3n_i)) is the cycle case; S4 of the counterexample portrait)', 'PROVED (known) + FINITE-EXACT (1000 orbits)', 'TWO-PLACE', 'P6 (done)'),
    ('REALVAL', 'VAL'): ('HYP-9139: a divergent 3n-1 orbit (equivalently a negative divergent 3n+1 orbit) with sum 1/m_l = infinity would make the real and 2-adic values of a non-periodic word coincide at an integer; conjecturally impossible; open beyond bounded critical discrepancy', 'OPEN (new hypothesis)', 'TWO-PLACE', 'none decisive (needs a divergent orbit); partial: prove it for discrepancy O(log l)'),
    ('THETA3', 'WORD'): ('Theta_3(x,y) = sum q^(k^3) x^(k^2) y^k satisfies Theta_3(x,y) = 1 + qxy Theta_3(q^3 x, q^3 x^2 y): the shift acts on cubes through squares (polarization); HYP-9127 already records this system as a unipotent Mahler system, inadmissible', 'PROVED (identity); route BLOCKED (HYP-9127 update)', 'BLOCKED', 'none'),
    ('GMFAIR', 'AFF'): ('no GM-fair pairing exists for T: (a/2)(3b+1)/2 = ab forces b = 1; AM-fairness (THM-4470) transports to GM-fairness only inside the EXP image', 'PROVED (one line)', 'DEGENERATE', 'none'),
    ('HEIGHT', 'AFF'): ('on a rank-one curve, hat-h(3P + G) = (9n^2 + 6n + 1) hat-h(G): the parallelogram law is the polarization of the seed; nothing beyond n', 'PROVED (restatement)', 'COSMETIC', 'none'),
    ('SUMSQ', 'GATE'): ('odd sums of two squares are 1 mod 4: all are down points of T_+ and none are down points of T_-; the norm form of Z[i] is positive definite, so this is an order statement encoded algebraically', 'PROVED', 'SPARSE / SIGN', 'none'),
    ('SUMSQ', 'SIGN'): ('n > 0 iff n is a sum of four squares (Lagrange): positivity is algebraic globally but every 2-adic integer is a sum of four squares, so nothing local sees it', 'PROVED (classical)', 'SPARSE / SIGN', 'none'),
}


def cell_class(t, i):
    if (t, i) in CARDS:
        return CARDS[(t, i)][2]
    kind = TRANSPORTS[t][1]
    if kind in ('lossless', 'lossless-residue'):
        return 'EMPTY (lossless transport adds nothing to this ingredient)'
    if kind == 'residue':
        return 'EMPTY (residue transport does not touch this ingredient)'
    if kind in ('sparse', 'sparse/sign'):
        return 'EMPTY (density-zero object)'
    return 'EMPTY'


def main():
    verbose = '--cards' in sys.argv or True
    print("# squares/doubles transport foundry -- session collatz-squares-doubles-20260925")
    print("transports: %d, ingredients: %d, cells: %d, curated cards: %d"
          % (len(TRANSPORTS), len(INGREDIENTS), len(TRANSPORTS) * len(INGREDIENTS), len(CARDS)))
    print()
    print("## transport typing")
    print("| key | transport | kind | ambient structure gained | sees SHEET | sees DRIFT | sees DEFECT |")
    print("|---|---|---|---|---|---|---|")
    for k, (name, kind, amb, sh, dr, de) in TRANSPORTS.items():
        print("| %s | %s | %s | %s | %s | %s | %s |" % (k, name, kind, amb or '-', sh, dr, de))
    print()
    print("## grid (verdict class per cell)")
    print("| transport \\ ingredient | " + " | ".join(INGREDIENTS) + " |")
    print("|---|" + "---|" * len(INGREDIENTS))
    classes = Counter()
    for t in TRANSPORTS:
        row = []
        for i in INGREDIENTS:
            c = cell_class(t, i)
            classes[c.split(' (')[0]] += 1
            row.append(c.split(' (')[0])
        print("| %s | %s |" % (t, " | ".join(row)))
    print()
    print("## class counts over the %d cells" % (len(TRANSPORTS) * len(INGREDIENTS)))
    for c, n in classes.most_common():
        print("- %s: %d" % (c, n))
    print()
    print("## curated cards")
    print("| # | transport | ingredient | statement | status | class | cheapest decisive test |")
    print("|---|---|---|---|---|---|---|")
    for idx, ((t, i), (stmt, status, cls, test)) in enumerate(CARDS.items(), start=1):
        print("| %d | %s | %s | %s | %s | %s | %s |" % (idx, t, i, stmt, status, cls, test))
    print()
    print("## ranking of the non-cosmetic classes (modelling judgement)")
    print("1. TWO-PLACE (REALVAL): the only class that sees the sheet AND the drift AND is orbitwise (defect-aware); it names one new open statement (HYP-9139) and no new mechanism.")
    print("2. NEW OBJECT (MONOID): an exact equivalent of Collatz on the diagonal of the exponent lattice with a proved divergence theorem off it; the square gate is synchronisation, not a mechanism.")
    print("3. EXACT MODEL (DLOG): the parity graph is the +-sqrt graph on F_p^* for Fermat primes; useful as a picture of THM-4474, adds no predicate.")
    print("4. RESIDUE LAW (QCHAR, JACOBI, LOCSQ): exact word-function lemmas; sheet-blind by construction.")
    print("5. SPARSE, SIGN TRANSFER, SIDE CONSTRAINT, FAMILY VIEW: exact but without measure or leverage.")
    print("6. COSMETIC, DEGENERATE, BLOCKED: recorded so that they are not regenerated.")


if __name__ == '__main__':
    main()
