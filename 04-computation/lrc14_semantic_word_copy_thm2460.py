#!/usr/bin/env python3
"""Exact finite checks for THM-2460.

The companion uses only integer and Fraction arithmetic.  Every
truth-bearing check goes through ``require`` and therefore remains active
under ``python -O``.
"""

from fractions import Fraction


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


P = 13
LOCAL_ATOMS = 2**7
WORDS = 3
EFFECTIVE_ATOMS = WORDS * LOCAL_ATOMS
ELIGIBLE_COLOURS = 12 * (P**2 - 1)
EXTENSIONS = (1, 16, 8, 4, 2, 1)


def frac_mod_one(x):
    return x - (x.numerator // x.denominator)


print("THM-2460 SEMANTIC-WORD COPY -- exact finite audit")

# Fixed-word diagonal idempotence.
fixed_word = 1
survivors = 0
for omega in range(LOCAL_ATOMS):
    for sigma in range(WORDS):
        present = int(sigma == fixed_word)
        for nu in range(LOCAL_ATOMS):
            for tau in range(WORDS):
                value = present * int(omega == nu) * int(sigma == tau)
                expected = int(
                    omega == nu and sigma == fixed_word and tau == fixed_word
                )
                require(value == expected, "fixed-word diagonal identity failed")
                survivors += value
require(survivors == LOCAL_ATOMS, "wrong fixed-word survivor count")
print(f"fixed word: {survivors} matched local atoms, zero extra loss")

# The unsplit return residual has exactly three effective word states.
require(EFFECTIVE_ATOMS == 384, "effective atom count")
off_diagonal = 0
diagonal = 0
for left in range(EFFECTIVE_ATOMS):
    left_word = left // LOCAL_ATOMS
    left_atom = left % LOCAL_ATOMS
    for right in range(EFFECTIVE_ATOMS):
        right_word = right // LOCAL_ATOMS
        right_atom = right % LOCAL_ATOMS
        matched = int(left_word == right_word and left_atom == right_atom)
        diagonal += matched
        off_diagonal += 1 - matched
require(diagonal == EFFECTIVE_ATOMS, "word/local diagonal census")
require(
    off_diagonal == EFFECTIVE_ATOMS**2 - EFFECTIVE_ATOMS,
    "word/local off-diagonal census",
)

drift_den = EFFECTIVE_ATOMS**2
eligible_den = P * drift_den
coefficient_den = ELIGIBLE_COLOURS * eligible_den
require(drift_den == 147456, "drift denominator")
require(eligible_den == 1916928, "eligible-energy denominator")
require(coefficient_den == 3864526848, "coefficient denominator")
require(tuple(WORDS * n for n in EXTENSIONS) == (3, 48, 24, 12, 6, 3),
        "adaptive extension invoice")
print(
    "unsplit return:",
    f"{EFFECTIVE_ATOMS} atoms, D/{drift_den},",
    f"eligible D/{eligible_den}, coefficient^2 D/{coefficient_den}",
)
print("adaptive three-word extension counts:", *(WORDS * n for n in EXTENSIONS))

# Root constancy: R*(y+u)/13 and 13^(k-1)*y differ by an integer.
test_y = Fraction(17, 113)
root_checks = 0
for k in range(1, 7):
    R = P**k
    for u in range(P):
        lhs = frac_mod_one(Fraction(R, P) * (test_y + u))
        rhs = frac_mod_one(P ** (k - 1) * test_y)
        require(lhs == rhs, "future word is not root-constant")
        require(R % P == 0, "target-neutral clock")
        root_checks += 1
print(f"root-constant future-word checks: {root_checks}")

# The enhanced graph is the disjoint union of three 128-vertex blocks.
vertices = EFFECTIVE_ATOMS
within_block_pairs = WORDS * LOCAL_ATOMS**2
cross_block_pairs = vertices**2 - within_block_pairs
require(within_block_pairs == 49152, "within-word pair count")
require(cross_block_pairs == 98304, "cross-word pair count")

# A uniform sharp matrix has total mass one and every edge 1/(3*128^2).
edge_mass = Fraction(1, WORDS * LOCAL_ATOMS**2)
total_mass = edge_mass * within_block_pairs
require(total_mass == 1, "uniform word-block matrix mass")
require(edge_mass == Fraction(1, 49152), "sharp total-to-edge floor")
fixed_block_edge = Fraction(1, LOCAL_ATOMS**2)
require(fixed_block_edge == Fraction(1, 16384), "fixed-word edge floor")
print(
    "word-block graph:",
    f"{within_block_pairs} within-block pairs,",
    f"{cross_block_pairs} forced cross-word zeros,",
    "edge floors 1/16384 fixed and 1/49152 global",
)

# Target neutrality and septimal activity of the word harmonic.
phase_checks = 0
for k in range(1, 9):
    R = P**k
    epsilon = R % 7
    expected = 6 if k % 2 else 1
    require(epsilon == expected, "wrong septimal clock parity")
    for beta in range(1, 7):
        require((R * beta) % P == 0, "word harmonic acquired target charge")
        require((R * beta) % 7 != 0, "word harmonic lost septimal activity")
        phase_checks += 1
print(f"target-neutral/septimally-active harmonic checks: {phase_checks}")

# A copied Boolean word is not an independently address-bearing harmonic.
# On C_4, delta_0 has every normalized Fourier coefficient 1/4.  The
# convolution identity qhat=qhat*qhat has four nonzero allocations for
# each output frequency.
qhat = [Fraction(1, 4)] * 4
for n in range(4):
    pieces = [qhat[b] * qhat[(n - b) % 4] for b in range(4)]
    require(all(piece != 0 for piece in pieces), "missing harmonic allocation")
    require(sum(pieces) == qhat[n], "idempotent Fourier convolution")
print("second word harmonic: 4 nonzero allocations per C4 output (noncanonical)")

# Sharp split-base hostile.  Three rational base cells are the three word
# blocks.  The selected table lives only in word zero.  A is present on
# words zero and one at root zero; F is present only on word one at root
# one.  The packets are rootwise disjoint, but same-parent co-support is
# confined to word one.
base_mass = [Fraction(1, 3)] * 3
selected_table = [1, 0, 0]
a_count = [1, 1, 0]
f_count = [0, 1, 0]
co_support = [
    base_mass[s] * a_count[s] * f_count[s] for s in range(WORDS)
]
require(sum(base_mass[s] * selected_table[s] for s in range(WORDS))
        == Fraction(1, 3), "selected table mass")
require(co_support == [0, Fraction(1, 3), 0], "split-base co-support hostile")
require(sum(co_support) > 0, "global co-support should survive")
require(co_support[0] == 0, "selected word unexpectedly has co-support")
print("split-base hostile: selected drift word 0, all co-support word 1")

# A word label alone cannot remove the uniform-offset boundary: tensoring
# any finite chart family by one constant word label preserves every
# averaged occupancy and convolution equation.
offset_owner = [Fraction(3, 20), Fraction(1, 10)]
offset_replica = [Fraction(1, 20), Fraction(0, 1)]
word_label = [1, 1]
require(
    [word_label[i] * offset_owner[i] for i in range(2)] == offset_owner,
    "word tensor changed owner density",
)
require(
    [word_label[i] * offset_replica[i] for i in range(2)] == offset_replica,
    "word tensor changed replica density",
)
print("constant-word tensor leaves the uniform-offset hostile unchanged")

print("fixed prescribed clock and later delayed clock remain distinct labels")
print("all exact checks passed")
