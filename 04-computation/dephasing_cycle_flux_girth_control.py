#!/usr/bin/env python3
"""Exact controls for the first Wilson-flux Schur word on unit cycles."""

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def superoperator_blocks(h):
    """Return B=P K Q, C=Q K P, D=Q K Q for K=-i[H,-]."""
    n = h.rows
    basis = [(a, b) for a in range(n) for b in range(n)]
    index = {pair: k for k, pair in enumerate(basis)}
    k_super = sp.zeros(n * n)
    for a, b in basis:
        column = index[a, b]
        for vertex in range(n):
            k_super[index[vertex, b], column] += -sp.I * h[vertex, a]
            k_super[index[a, vertex], column] += sp.I * h[b, vertex]
    p_indices = [index[j, j] for j in range(n)]
    q_indices = [index[a, b] for a, b in basis if a != b]
    return (
        k_super.extract(p_indices, q_indices),
        k_super.extract(q_indices, p_indices),
        k_super.extract(q_indices, q_indices),
    )


def raw_schur_word(blocks, length):
    """The length-m word A_m=B D^(m-2) C, for m>=2."""
    b_block, c_block, d_block = blocks
    return sp.simplify(b_block * d_block ** (length - 2) * c_block)


def unit_cycle(n, z):
    h = sp.zeros(n)
    for vertex in range(n - 1):
        h[vertex, vertex + 1] = 1
        h[vertex + 1, vertex] = 1
    h[n - 1, 0] = z
    h[0, n - 1] = 1 / z
    return h


def forward_shift(n):
    """S acts on columns by (S p)_j=p_(j+1), indices modulo n."""
    shift = sp.zeros(n)
    for vertex in range(n):
        shift[vertex, (vertex + 1) % n] = 1
    return shift


z = sp.symbols("z", nonzero=True)

print("UNIT-CYCLE STRONG-DEPHASING GIRTH CONTROL")
print("convention=K=-i[H,-], H_(g-1,0)=z, H_(0,g-1)=z^-1")
print("raw_word=A_m=B D^(m-2) C")

for g in range(3, 7):
    h_cycle = unit_cycle(g, z)
    blocks = superoperator_blocks(h_cycle)

    for length in range(2, g):
        lower_delta = sp.simplify(
            raw_schur_word(blocks, length)
            - raw_schur_word(blocks, length).subs(z, 1)
        )
        require(
            lower_delta == sp.zeros(g),
            f"C{g} acquired flux before commutator length {g}",
        )

    first_delta = sp.simplify(
        raw_schur_word(blocks, g) - raw_schur_word(blocks, g).subs(z, 1)
    )
    shift = forward_shift(g)
    binomial_circulant = sp.expand((sp.eye(g) - shift) ** g)
    if g % 2:
        wilson_factor = sp.I * (-1) ** ((g - 1) // 2) * (z - 1 / z)
        parity = "odd:sine/oriented/skew"
    else:
        wilson_factor = (-1) ** (g // 2) * (z + 1 / z - 2)
        parity = "even:cosine/reciprocal/symmetric"
    expected = sp.simplify(wilson_factor * binomial_circulant)

    require(first_delta == expected, f"C{g} Wilson formula failed")
    require(first_delta != sp.zeros(g), f"C{g} first Wilson coefficient vanished")
    require(
        binomial_circulant.T
        == ((-1) ** g) * binomial_circulant,
        f"C{g} transpose parity failed",
    )
    require(
        binomial_circulant.rank() == g - 1,
        f"C{g} binomial circulant has unexpected rank",
    )
    require(
        [sp.simplify(sum(first_delta.row(i))) for i in range(g)] == [0] * g,
        f"C{g} population conservation failed",
    )
    print(
        f"C{g}: lower_words_phase_blind=True; parity={parity}; "
        f"rank_base={g - 1}"
    )
    print(f"C{g}: base_matrix={(binomial_circulant).tolist()}")
    print(f"C{g}: wilson_factor={wilson_factor}")

# Weighted C4 / K2,2 control.  A spanning-tree gauge makes three hopping
# amplitudes positive; z carries the remaining Wilson phase.
r0, r1, r2, r3 = sp.symbols("r0 r1 r2 r3", positive=True)
h_weighted_c4 = sp.zeros(4)
for vertex, magnitude in enumerate((r0, r1, r2)):
    h_weighted_c4[vertex, vertex + 1] = magnitude
    h_weighted_c4[vertex + 1, vertex] = magnitude
h_weighted_c4[3, 0] = r3 * z
h_weighted_c4[0, 3] = r3 / z
weighted_blocks = superoperator_blocks(h_weighted_c4)
weighted_delta = sp.simplify(
    raw_schur_word(weighted_blocks, 4)
    - raw_schur_word(weighted_blocks, 4).subs(z, 1)
)
c4_base = (sp.eye(4) - forward_shift(4)) ** 4
cycle_magnitude = r0 * r1 * r2 * r3
weighted_expected = sp.simplify(
    cycle_magnitude * (z + 1 / z - 2) * c4_base
)
require(
    sp.simplify(weighted_delta - weighted_expected) == sp.zeros(4),
    "weighted C4 scale failed",
)
require(
    sp.simplify(weighted_delta[0, 0] - 2 * cycle_magnitude * (z + 1 / z - 2))
    == 0,
    "weighted C4 diagonal readout failed",
)
print("weighted_C4_delta=R*(z+z^-1-2)*(I-S)^4; R=r0*r1*r2*r3")
print("K2,2_readout=Delta_A4_00=4*(Re(ad*conj(bc))-|abcd|)")
print("det_reconstruction=|ad|^2+|bc|^2-2|abcd|-Delta_A4_00/2")

# Tree hostile control: a formal phase on one edge of P6 is gauge removable.
n_tree = 6
h_path = sp.zeros(n_tree)
for vertex in range(n_tree - 1):
    value = z if vertex == 2 else sp.Integer(1)
    h_path[vertex, vertex + 1] = value
    h_path[vertex + 1, vertex] = 1 / value
path_blocks = superoperator_blocks(h_path)
for length in range(2, 9):
    path_delta = sp.simplify(
        raw_schur_word(path_blocks, length)
        - raw_schur_word(path_blocks, length).subs(z, 1)
    )
    require(path_delta == sp.zeros(n_tree), "tree phase survived vertex gauge")
print("P6_tree_control=phase_blind_for_A2_through_A8")

# Chord hostile control: on C5 plus chord 0-2, the marked edge belongs to a
# four-cycle, so its phase is already visible at length four, not length five.
h_chord = unit_cycle(5, z)
h_chord[0, 2] = 1
h_chord[2, 0] = 1
chord_blocks = superoperator_blocks(h_chord)
chord_a3_delta = sp.simplify(
    raw_schur_word(chord_blocks, 3) - raw_schur_word(chord_blocks, 3).subs(z, 1)
)
chord_a4_delta = sp.simplify(
    raw_schur_word(chord_blocks, 4) - raw_schur_word(chord_blocks, 4).subs(z, 1)
)
require(chord_a3_delta == sp.zeros(5), "marked edge entered through wrong triangle")
require(chord_a4_delta != sp.zeros(5), "four-cycle containing marked edge was missed")
print("C5_plus_chord_0-2=marked_phase_blind_at_A3_nonzero_at_A4")

print("VERDICT=formal Wilson data first occurs at raw length g on chordless C_g")
print("SLOW_ORDER=Gamma^(-(g-2)) after leading-rate time normalization")
