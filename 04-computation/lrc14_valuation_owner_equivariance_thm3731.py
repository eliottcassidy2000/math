"""Exact owner-equivariance and owner-loop-drift interface probe.

Universe 1 is the labelled lift of the 165 THM-2349 valuation profiles.
The S_3 action moves blocker labels.  Universe 2 is one explicit distinct-
coefficient lift of every profile.  The final control replays the THM-2478
nonnegative diagonal-zero tensor hostile at the literal THM-2365 projection.
"""

from fractions import Fraction as F
from itertools import permutations, product


P = 13
LABELS = tuple(range(3))
PERMS = tuple(permutations(LABELS))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def act_tuple(values, perm):
    """Move old label i to new label perm[i]."""
    moved = [None] * len(values)
    for i, value in enumerate(values):
        moved[perm[i]] = value
    return tuple(moved)


def act_set(indices, perm):
    return frozenset(perm[i] for i in indices)


def valuation_13(value):
    require(isinstance(value, int) and value > 0, "valuation input")
    exponent = 0
    while value % P == 0:
        value //= P
        exponent += 1
    return exponent


def admissible_owners(profile):
    minimum = min(profile)
    return frozenset(i for i, value in enumerate(profile) if value == minimum)


def coefficient_owner(coefficients):
    return min(LABELS, key=lambda i: (valuation_13(coefficients[i]), coefficients[i]))


def profile_equivariance_probe():
    quotient = tuple(
        (1, middle, deepest)
        for deepest in range(5, 20)
        for middle in range(1, deepest)
    )
    strict = tuple(profile for profile in quotient if profile[1] > 1)
    repeated = tuple(profile for profile in quotient if profile[1] == 1)
    labelled = frozenset(act_tuple(profile, perm) for profile in quotient for perm in PERMS)

    require(len(quotient) == 165, "165-row quotient size")
    require(len(strict) == 150 and len(repeated) == 15, "strict/repeated split")
    require(len(labelled) == 945, "labelled profile lift size")
    require(
        frozenset(tuple(sorted(profile)) for profile in labelled) == frozenset(quotient),
        "the labelled lift has exactly 165 S3 orbits",
    )

    # The set-valued minimum-owner field is an equivariant positive control.
    for profile in labelled:
        for perm in PERMS:
            require(
                admissible_owners(act_tuple(profile, perm))
                == act_set(admissible_owners(profile), perm),
                "set-valued owner equivariance",
            )

    # A deterministic equivariant selector on an orbit exists precisely when
    # one admissible owner is fixed by the representative's stabilizer.
    orbit_pass = []
    orbit_fail = []
    failure_certificates = []
    for profile in quotient:
        stabilizer = tuple(perm for perm in PERMS if act_tuple(profile, perm) == profile)
        admissible = admissible_owners(profile)
        fixed = tuple(
            owner
            for owner in admissible
            if all(perm[owner] == owner for perm in stabilizer)
        )
        if fixed:
            orbit_pass.append(profile)
        else:
            orbit_fail.append(profile)
            swaps = tuple(
                perm
                for perm in stabilizer
                if any(perm[owner] != owner for owner in admissible)
            )
            require(swaps, "failed orbit has a moving stabilizer element")
            failure_certificates.append((profile, admissible, swaps[0]))

    require(tuple(orbit_pass) == strict, "all and only strict profiles admit a selector")
    require(tuple(orbit_fail) == repeated, "all and only repeated profiles obstruct")
    require(
        all(len(admissible_owners(profile)) == 1 for profile in strict),
        "strict unique-min positive control",
    )

    # The minimal repair is the marked lift (profile, chosen admissible owner).
    marked = tuple(
        (profile, owner)
        for profile in labelled
        for owner in admissible_owners(profile)
    )
    require(len(marked) == 990, "marked owner lift size")
    for profile, owner in marked:
        for perm in PERMS:
            moved_profile = act_tuple(profile, perm)
            moved_owner = perm[owner]
            require(moved_owner in admissible_owners(moved_profile), "marked lift closure")
            require(moved_owner == perm[owner], "marked projection equivariance")

    return quotient, strict, repeated, labelled, marked, tuple(failure_certificates)


def coefficient_sidecar_probe(quotient):
    # Literal coefficient representatives from the full-owner hostile in
    # THM-2367 (20).  That theorem explicitly says they are not scalar-cover
    # rows; here they are only a positive control for the coefficient key.
    q_two = 2940
    c_one = P * 7 * q_two
    samples = []
    for one, middle, deepest in quotient:
        require(one == 1, "normalized first depth")
        c_two = c_one * 7 * P ** (middle - 1)
        c_three = c_two * 7 * P ** (deepest - middle)
        base = (c_one, c_two, c_three)
        require(len(set(base)) == 3, "distinct coefficient realization")
        require(
            tuple(valuation_13(value) for value in base) == (1, middle, deepest),
            "THM-2367 representative valuation profile",
        )
        for perm in PERMS:
            samples.append(act_tuple(base, perm))

    require(len(samples) == 990, "coefficient sample count")
    require(len(set(samples)) == 990, "coefficient samples are distinct")

    dilations = (1, 2, 5, P, 2 * P, P**2)
    for coefficients in samples:
        profile = tuple(valuation_13(value) for value in coefficients)
        owner = coefficient_owner(coefficients)
        require(owner in admissible_owners(profile), "coefficient key chooses a minimum-depth owner")
        for perm in PERMS:
            moved = act_tuple(coefficients, perm)
            require(
                coefficient_owner(moved) == perm[owner],
                "coefficient key is S3-equivariant",
            )
        for scale in dilations:
            dilated = tuple(scale * value for value in coefficients)
            require(
                coefficient_owner(dilated) == owner,
                "coefficient key is invariant under common positive dilation",
            )

    # On every labelled repeated-depth profile, the coefficient sidecar has
    # two realizations selecting opposite depth-one labels.  Thus this repair
    # cannot factor through the valuation profile alone.
    owners_by_profile = {}
    for coefficients in samples:
        profile = tuple(valuation_13(value) for value in coefficients)
        owners_by_profile.setdefault(profile, set()).add(coefficient_owner(coefficients))
    repeated_labelled = tuple(
        profile
        for profile, owners in owners_by_profile.items()
        if len(admissible_owners(profile)) == 2
    )
    require(len(repeated_labelled) == 45, "labelled repeated-profile count")
    require(
        all(owners_by_profile[profile] == set(admissible_owners(profile)) for profile in repeated_labelled),
        "coefficient sidecar resolves both stabilizer-related owners",
    )
    return tuple(samples), repeated_labelled


def complete_atom_probe():
    atoms = tuple(product((0, 1), repeat=7))
    source_atoms = tuple(
        atom
        for atom in atoms
        if atom[:5] == (0, 0, 0, 0, 0) and atom[5:] == (1, 0)
    )
    require(len(atoms) == 128, "complete atom bank size")
    require(source_atoms == ((0, 0, 0, 0, 0, 1, 0),), "unique semantic source atom")
    return atoms, source_atoms[0]


def projection(table):
    projected = {}
    for r, s, t in table:
        projected[r, s, t] = sum(
            table[(r + b) % P, (s + a) % P, (t + b) % P]
            for a in range(P)
            for b in range(P)
        ) / (P * P)
    return projected


def drift(table):
    projected = projection(table)
    return sum((table[key] - projected[key]) ** 2 for key in table) / P**3


def owner_loop_hostile_probe():
    keys = tuple(product(range(P), repeat=3))
    owner = {
        (r, s, t): F(1, 4) if (r - t) % P == 1 else F(0)
        for r, s, t in keys
    }
    nonowner = {
        (r, s, t): F(1, 4) if (r, s, t) == (1, 1, 0) else F(0)
        for r, s, t in keys
    }
    aggregate = {key: owner[key] + nonowner[key] for key in keys}

    owner_projection = projection(owner)
    require(owner_projection == owner, "owner component is circulant")
    require(
        all(aggregate[t, s, t] == 0 for s in range(P) for t in range(P)),
        "aggregate diagonal zero",
    )
    require(sum(nonowner[r, 0, 0] for r in range(P)) == 0, "nonowner untwisted mass zero")

    owner_drift = drift(owner)
    nonowner_drift = drift(nonowner)
    aggregate_drift = drift(aggregate)
    require(owner_drift == 0, "owner-loop drift vanishes")
    require(nonowner_drift == aggregate_drift, "all aggregate drift lies off the owner")
    require(aggregate_drift == F(21, 742586), "THM-2478 exact drift value")
    return owner_drift, nonowner_drift, aggregate_drift


def main():
    quotient, strict, repeated, labelled, marked, certificates = profile_equivariance_probe()
    samples, repeated_labelled = coefficient_sidecar_probe(quotient)
    atoms, source_atom = complete_atom_probe()
    owner_drift, nonowner_drift, aggregate_drift = owner_loop_hostile_probe()

    profile, admissible, swap = certificates[0]
    print("LRC14 semantic-arrival owner equivariance probe")
    print("status=PASS")
    print(f"valuation_quotient_rows={len(quotient)}")
    print(f"strict_rows={len(strict)} repeated_rows={len(repeated)}")
    print(f"labelled_profiles={len(labelled)} S3_orbits={len(quotient)}")
    print(f"deterministic_owner_orbits_pass={len(strict)} fail={len(repeated)}")
    print(f"first_failure_profile={profile} admissible={tuple(sorted(admissible))} stabilizer_swap={swap}")
    print(f"set_valued_owner_equivariance=PASS marked_lift_size={len(marked)}")
    print(f"coefficient_samples={len(samples)} coefficient_key_equivariance=PASS")
    print("coefficient_source=THM-2367-equation-20-noncover-hostile-representatives")
    print(f"repeated_labelled_profiles_resolved_by_sidecar={len(repeated_labelled)}")
    print("coefficient_key=argmin_i(nu13(c_i),c_i)")
    print(f"complete_atoms={len(atoms)} unique_source_atom={source_atom}")
    print(f"owner_drift={owner_drift}")
    print(f"nonowner_drift={nonowner_drift}")
    print(f"aggregate_drift={aggregate_drift}")
    print("aggregate_positive_owner_zero_interface_hostile=PASS")
    print("scope=no scalar-cover row excluded; no semantic arrival proved")


if __name__ == "__main__":
    main()
