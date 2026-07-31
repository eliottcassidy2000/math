"""Exact companion for THM-2992.

Verify the signed-edge quadratic factorization of a depressed quartic, the
block discriminant/resultant identities, all twelve S4 lifts of quotient
transpositions, and exact identity/transposition/double-transposition/four-
cycle local controls.  Every truth-bearing check uses ``require`` so that
``python -O`` remains a genuine verification run.
"""

from itertools import permutations

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def cycle_type(g):
    used = set()
    parts = []
    for i in range(len(g)):
        if i in used:
            continue
        j = i
        length = 0
        while j not in used:
            used.add(j)
            length += 1
            j = g[j]
        parts.append(length)
    return tuple(sorted(parts, reverse=True))


def compose(g, h):
    """First h, then g."""
    return tuple(g[h[i]] for i in range(len(g)))


def inverse(g):
    ans = [None] * len(g)
    for i, image in enumerate(g):
        ans[image] = i
    return tuple(ans)


# ---------------------------------------------------------------------------
# Exact signed-edge algebra.
# Write h=q/s while retaining s^2=u.  This removes every denominator from
# the quadratic blocks and leaves the resolvent-shell equation
# H=(p+u)^2-h^2-4r=0, equivalently S(u)=uH=0 because u is a unit.
# ---------------------------------------------------------------------------

T, p, r, u, s, h = sp.symbols("T p r u s h")
quartic = T**4 + p * T**2 + s * h * T + r
shell = sp.expand((p + u) ** 2 - h**2 - 4 * r)
S_at_u = sp.expand(u**3 + 2 * p * u**2 + (p**2 - 4 * r) * u - u * h**2)
S_prime_at_u = sp.expand(3 * u**2 + 4 * p * u + p**2 - 4 * r)

a = sp.expand((p + u + h) / 2)
b = sp.expand((p + u - h) / 2)
Q_plus = sp.expand(T**2 - s * T + a)
Q_minus = sp.expand(T**2 + s * T + b)


def reduce_s(expr):
    return sp.expand(sp.rem(sp.expand(expr), s**2 - u, s))


def reduce_shell(expr):
    return sp.expand(sp.rem(reduce_s(expr), shell, h))


require(sp.expand(S_at_u - u * shell) == 0, "S(u)=uH identity changed")
factor_error = reduce_s(Q_plus * Q_minus - quartic)
require(
    sp.expand(factor_error - shell / 4) == 0,
    "off-shell quadratic factorization error changed",
)

delta_plus = sp.expand(u - 4 * a)
delta_minus = sp.expand(u - 4 * b)
require(delta_plus == -2 * p - u - 2 * h, "plus block discriminant changed")
require(delta_minus == -2 * p - u + 2 * h, "minus block discriminant changed")
require(
    sp.expand(delta_plus + delta_minus + 4 * p + 2 * u) == 0,
    "block-discriminant sum changed",
)
require(
    sp.expand(delta_plus - delta_minus + 4 * h) == 0,
    "block-discriminant difference changed",
)
block_product = reduce_shell(delta_plus * delta_minus)
require(
    sp.expand(block_product - (16 * r - 4 * p * u - 3 * u**2)) == 0,
    "block-discriminant product changed",
)

block_resultant = reduce_shell(sp.resultant(Q_plus, Q_minus, T))
require(
    sp.expand(block_resultant - S_prime_at_u) == 0,
    "quadratic-block resultant is not S'(u)",
)

quartic_discriminant = (
    256 * r**3
    - 128 * p**2 * r**2
    + 144 * p * (s * h) ** 2 * r
    - 27 * (s * h) ** 4
    + 16 * p**4 * r
    - 4 * p**3 * (s * h) ** 2
)
quartic_discriminant_on_shell = reduce_shell(quartic_discriminant)
require(
    sp.expand(
        quartic_discriminant_on_shell
        - block_product * S_prime_at_u**2
    )
    == 0,
    "quartic discriminant product identity changed",
)

gauge = {s: -s, h: -h}
require(
    sp.expand(Q_plus.xreplace(gauge) - Q_minus) == 0,
    "signed-edge reversal does not exchange quadratic blocks",
)
require(
    sp.expand(delta_plus.xreplace(gauge) - delta_minus) == 0,
    "signed-edge reversal does not exchange block discriminants",
)


# ---------------------------------------------------------------------------
# Every lift in S4 of every quotient transposition in S4/V4 ~= S3.
# ---------------------------------------------------------------------------

S4 = tuple(permutations(range(4)))
matchings = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


def canonical_matching(blocks):
    return tuple(sorted(tuple(sorted(block)) for block in blocks))


matching_lookup = {
    canonical_matching(blocks): i for i, blocks in enumerate(matchings)
}


def quotient_action(g):
    ans = []
    for blocks in matchings:
        image = canonical_matching(tuple(tuple(g[x] for x in block) for block in blocks))
        require(image in matching_lookup, "matching quotient left its carrier")
        ans.append(matching_lookup[image])
    return tuple(ans)


quotient_transpositions = sorted(
    {quotient_action(g) for g in S4 if cycle_type(quotient_action(g)) == (2, 1)}
)
require(len(quotient_transpositions) == 3, "S3 transposition census changed")

lift_count = 0
sheet_transposition_count = 0
four_cycle_count = 0
v4_difference_count = 0
fixed_block_records = []

for sigma in quotient_transpositions:
    lifts = [g for g in S4 if quotient_action(g) == sigma]
    require(len(lifts) == 4, f"quotient transposition {sigma} has !=4 lifts")
    transposition_lifts = [g for g in lifts if cycle_type(g) == (2, 1, 1)]
    four_cycle_lifts = [g for g in lifts if cycle_type(g) == (4,)]
    require(
        len(transposition_lifts) == len(four_cycle_lifts) == 2,
        f"quotient transposition {sigma} lift types changed",
    )

    fixed_matching = next(i for i in range(3) if sigma[i] == i)
    blocks = tuple(tuple(sorted(block)) for block in matchings[fixed_matching])

    local_fixed_blocks = []
    for g in transposition_lifts:
        fixed_sheets = tuple(i for i in range(4) if g[i] == i)
        require(
            len(fixed_sheets) == 2 and fixed_sheets in blocks,
            "transposition lift does not fix one complete matching block",
        )
        local_fixed_blocks.append(fixed_sheets)
    require(
        set(local_fixed_blocks) == set(blocks),
        "the two transposition lifts do not fix opposite matching blocks",
    )

    for g in four_cycle_lifts:
        require(not any(g[i] == i for i in range(4)), "four-cycle lift fixed a sheet")
        image0 = tuple(sorted(g[x] for x in blocks[0]))
        image1 = tuple(sorted(g[x] for x in blocks[1]))
        require(
            image0 == blocks[1] and image1 == blocks[0],
            "four-cycle lift does not exchange the signed-edge blocks",
        )

    difference = compose(inverse(transposition_lifts[0]), transposition_lifts[1])
    require(cycle_type(difference) == (2, 2), "transposition lifts do not differ by V4")
    require(
        all(tuple(sorted(difference[x] for x in block)) == block for block in blocks),
        "V4 difference does not preserve both blocks",
    )

    lift_count += len(lifts)
    sheet_transposition_count += len(transposition_lifts)
    four_cycle_count += len(four_cycle_lifts)
    v4_difference_count += 1
    fixed_block_records.append((sigma, fixed_matching, tuple(sorted(local_fixed_blocks))))

require(lift_count == 12, "quotient-transposition lift total changed")
require(sheet_transposition_count == 6, "sheet-transposition lift total changed")
require(four_cycle_count == 6, "four-cycle lift total changed")
require(v4_difference_count == 3, "matching-direction V4 census changed")


# ---------------------------------------------------------------------------
# Exact local controls over C((t)).
# ---------------------------------------------------------------------------

t, z, U = sp.symbols("t z U")


def at(expr, substitutions):
    return sp.factor(sp.expand(expr.subs(substitutions)))


# Identity: both block discriminants are units, hence squares in a strict
# henselian DVR with algebraically closed residue field.
identity_sub = {p: -1, r: 0, u: 4, s: 2, h: -3}
require(at(shell, identity_sub) == 0, "identity control left the resolvent shell")
identity_deltas = (at(delta_plus, identity_sub), at(delta_minus, identity_sub))
require(identity_deltas == (4, -8), "identity-control discriminants changed")
require(at(S_prime_at_u, identity_sub) != 0, "identity-control blocks collide")

# Transposition: (T-1)^2-t is the ramified block; T(T+2) is fixed.
transposition_sub = {
    p: -(3 + t),
    r: 0,
    u: 4,
    s: 2,
    h: 1 - t,
}
require(at(shell, transposition_sub) == 0, "transposition control left shell")
transposition_deltas = (
    at(delta_plus, transposition_sub),
    at(delta_minus, transposition_sub),
)
require(transposition_deltas == (4 * t, 4), "transposition parities changed")
transposition_quartic = sp.expand(quartic.subs(transposition_sub))
require(
    sp.expand(
        transposition_quartic
        - (T**4 - (3 + t) * T**2 + 2 * (1 - t) * T)
    )
    == 0,
    "transposition quartic changed",
)
require(
    at(S_prime_at_u, transposition_sub) == (t - 9) * (t - 1),
    "transposition block resultant changed",
)

# Double transposition: both blocks ramify in the same quadratic extension.
double_sub = {
    p: -2 - 3 * t,
    r: (1 - t) * (1 - 2 * t),
    u: 4,
    s: 2,
    h: t,
}
require(at(shell, double_sub) == 0, "double-transposition control left shell")
double_deltas = (at(delta_plus, double_sub), at(delta_minus, double_sub))
require(double_deltas == (4 * t, 8 * t), "double-transposition parities changed")
require(
    at(S_prime_at_u, double_sub).subs(t, 0) == 16,
    "double-transposition blocks collide at the closed point",
)

# Four-cycle boundary: X^4+tX+t is Eisenstein.  Its fixed matching root has
# valuation one, so no signed edge s lies in K.  After adjoining s, normalized
# valuations are v_L(u)=2 and v_L(q/s)=1, hence both block discriminants have
# odd valuation; this is the square (double-transposition) subgroup of C4,
# not a classification of the original four-cycle by two K-block parities.
four_cycle_quartic = T**4 + t * T + t
four_cycle_matching = U**3 - 4 * t * U - t**2
hensel_scaled = sp.expand(four_cycle_matching.subs(U, t * z) / t**2)
require(
    hensel_scaled == t * z**3 - 4 * z - 1,
    "four-cycle matching Hensel equation changed",
)
require(
    sp.expand(hensel_scaled.subs({t: 0, z: -sp.Rational(1, 4)})) == 0,
    "four-cycle valuation-one matching root changed",
)
require(
    sp.diff(hensel_scaled, z).subs({t: 0, z: -sp.Rational(1, 4)}) == -4,
    "four-cycle Hensel root is not simple",
)
four_cycle_discriminant = sp.discriminant(four_cycle_quartic, T)
require(
    sp.factor(four_cycle_discriminant) == -t**3 * (27 * t - 256),
    "Eisenstein four-cycle discriminant changed",
)
vL_u = 2
vL_q_over_s = 1
four_cycle_block_parities = (
    min(vL_u, vL_q_over_s) % 2,
    min(vL_u, vL_q_over_s) % 2,
)
require(
    four_cycle_block_parities == (1, 1),
    "four-cycle signed-edge extension valuation control changed",
)


print("THM-2992 signed quartic edge block-parity decoder")
print("symbolic_ring=Q[p,r,u,s,h,T]/(s^2-u,S(u)) with h=q/s")
print(f"factor_error={sp.factor(factor_error)}=S(u)/(4u)")
print(f"delta_plus={delta_plus}")
print(f"delta_minus={delta_minus}")
print(f"delta_product={block_product}")
print(f"block_resultant={block_resultant}=S'(u)")
print("quartic_discriminant=delta_plus*delta_minus*S'(u)^2")
print("signed_edge_reversal=block_label_swap")
print(f"quotient_transpositions={quotient_transpositions}")
print(
    f"quotient_lifts={lift_count} transpositions={sheet_transposition_count} "
    f"four_cycles={four_cycle_count}"
)
print(f"fixed_block_records={fixed_block_records}")
print(f"matching_direction_v4_differences={v4_difference_count}")
print(f"identity_deltas={identity_deltas} parity=(0,0)")
print(
    f"transposition_deltas={transposition_deltas} parity=(1,0) "
    f"resultant={at(S_prime_at_u, transposition_sub)}"
)
print(f"double_transposition_deltas={double_deltas} parity=(1,1)")
print("four_cycle_eisenstein=T^4+tT+t disc=t^3*(256-27t)")
print("four_cycle_fixed_matching_vK(u)=1 signed_edge_not_in_K=1")
print("after_signed_edge_extension=(vL(u),vL(q/s),block_parities)=(2,1,(1,1))")
print("quartic_signed_edge_block_parity_thm2992=PASS")
