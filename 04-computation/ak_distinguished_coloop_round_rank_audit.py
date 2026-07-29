#!/usr/bin/env python3
"""
Exact round-by-round rank audit for arithmetic-Kakeya forcing certificates.

For a live vertex set of size u, use (s,d)=(a+b,a-b) and write the generator
matrix as A=[B | D_rho B].  If F is the simultaneous firing set, every d_v
with v in F is a coloop.  Deleting those d-columns loses |F| rank; deleting
the corresponding s-columns loses a further p in [0,|F|].  Thus

    r_next = r - f - p,       H_next = H - p,

where r=rank(A), f=|F|, and H=r-u.  In a successful honest graph certificate,
each live component is pinned, so rank(B)=u and f<=H.

This script independently replays four frozen records:

  * strict/verifiable 13/7;
  * per-suffix, merge-free 7/4;
  * legal-slot quotient 12/7;
  * legal-slot quotient 9/5.

The mode-three witnesses are used only as quotient labelled multigraphs.
Passing this audit proves the rank identities for the resulting quotient
rows; it does not prove the still-owed same-H/stagewise constructibility.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "ak_projective_matroid_audit.py":
        "f5159f97717a09dcd49fccd2a28e39a1c535c46309f23e0dbf483eed93fb7f05",
    "ak_mode3_v2.py":
        "5a85b0b9de423028d929d3ce1d86a1a0a90cab5286f6434300d1104cc44cef55",
}
for name, expected in PINNED.items():
    payload = (HERE / name).read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(payload).hexdigest() == expected,
        f"pinned dependency changed: {name}",
    )


from ak_projective_matroid_audit import rref, strict_13_over_7  # noqa: E402
from ak_mode3_v2 import MERGE, Mode3Instance  # noqa: E402


def rank(rows, ncols):
    if ncols == 0 or not rows:
        return 0
    return len(rref(rows, ncols))


def bridge_count(brows, u):
    """Rows whose deletion lowers the full incidence-block rank."""
    if u == 0:
        return 0
    return sum(
        rank(brows[:i] + brows[i + 1:], u) < u
        for i in range(len(brows))
    )


def two_forest_partition(brows, u):
    """Brute exact partition into two B-independent row sets, if one exists."""
    e = len(brows)
    if e == 0:
        return ((), ())
    if e > 2 * u:
        return None
    indices = tuple(range(e))
    # Swap symmetry lets us require row 0 to lie in the first part.
    for left_size in range(max(1, e - u), min(u, e) + 1):
        for tail in combinations(indices[1:], left_size - 1):
            left = (0,) + tail
            left_set = set(left)
            right = tuple(i for i in indices if i not in left_set)
            if rank([brows[i] for i in left], u) != len(left):
                continue
            if rank([brows[i] for i in right], u) != len(right):
                continue
            return left, right
    return None


def topology_defect_witness(brows, u, target):
    """Smallest S with |S|-2 rank(B[S]) equal to the requested defect."""
    indices = tuple(range(len(brows)))
    for size in range(len(brows) + 1):
        for subset in combinations(indices, size):
            rr = rank([brows[i] for i in subset], u)
            if size - 2 * rr == target:
                return subset, rr
    return None


class QuotientCertificate:
    """The only data the rank audit needs: paid rows, live classes, and cost."""

    def __init__(self, name, n, t0, generators, cost):
        self.name = name
        self.n = int(n)
        self.t0 = frozenset(t0)
        self.generators = tuple(generators)
        self.cost = int(cost)

    @classmethod
    def from_strict(cls, name, instance):
        generators = tuple(
            sparse for _, sparse in instance.labelled_generators("strict")
        )
        t0 = {instance.index[v] for v in instance.t0}
        return cls(name, instance.n, t0, generators, instance.m + instance.r)

    @classmethod
    def from_mode3(cls, name, instance):
        return cls(
            name,
            instance.n,
            instance.T0,
            instance.generators(),
            instance.m() + instance.r(),
        )

    def matrices(self, live):
        """Return the active B and interleaved [B|D B] rows over Q."""
        order = tuple(sorted(live))
        pos = {v: j for j, v in enumerate(order)}
        brows = []
        arows = []
        for sparse in self.generators:
            brow = [Fraction(0)] * len(order)
            arow = [Fraction(0)] * (2 * len(order))
            for v, (a, b) in sparse.items():
                if v not in pos:
                    continue
                j = pos[v]
                brow[j] = Fraction(a + b)
                arow[2 * j] = Fraction(a + b)
                arow[2 * j + 1] = Fraction(a - b)
            if any(arow):
                brows.append(brow)
                arows.append(arow)
        return order, brows, arows


def fired_vertices(order, arows):
    if not arows:
        return ()
    reduced = rref(arows, 2 * len(order))
    fired = []
    for pivot, row in reduced:
        if pivot % 2 != 1:
            continue
        if sum(x != 0 for x in row) == 1:
            fired.append(order[pivot // 2])
    return tuple(sorted(fired))


def audit_profile(cert):
    live = set(range(cert.n)) - set(cert.t0)
    u0 = len(live)
    rounds = []
    while live:
        order, brows, arows = cert.matrices(live)
        u = len(live)
        rb = rank(brows, u)
        r = rank(arows, 2 * u)
        e = len(arows)
        c = e - u
        sigma = e - r
        z = bridge_count(brows, u)
        forest_partition = two_forest_partition(brows, u)
        defect_witness = topology_defect_witness(brows, u, sigma)
        fired = fired_vertices(order, arows)
        require(fired, f"{cert.name}: forcing stopped with {u} live vertices")
        f = len(fired)
        next_live = live - set(fired)
        _, next_brows, next_arows = cert.matrices(next_live)
        r_next = rank(next_arows, 2 * len(next_live))
        rb_next = rank(next_brows, len(next_live))
        p = r - f - r_next
        H = r - u
        H_next = r_next - len(next_live)
        nu = 2 * u - r
        nu_next = 2 * len(next_live) - r_next

        require(rb == u, f"{cert.name}: live incidence rank {rb} != {u}")
        require(c >= 0, f"{cert.name}: negative cycle rank")
        require(sigma >= max(0, 2 * c - (e - z)),
                f"{cert.name}: cycle/bicirculation lower bound failed")
        require(H == c - sigma, f"{cert.name}: H cycle identity failed")
        require(H <= min(c, u - z),
                f"{cert.name}: cycle/bridge headroom bound failed")
        require(sigma != 0 or forest_partition is not None,
                f"{cert.name}: full row rank without two-forest topology")
        require(defect_witness is not None,
                f"{cert.name}: slope-specific cancellation detected")
        require(rb_next == len(next_live),
                f"{cert.name}: next incidence rank is not full")
        require(0 <= p <= f, f"{cert.name}: p={p} not in [0,{f}]")
        require(r_next == r - f - p, f"{cert.name}: rank transition failed")
        require(H_next == H - p, f"{cert.name}: H transition failed")
        require(nu_next == nu - f + p,
                f"{cert.name}: nullity transition failed")
        require(r >= u + f, f"{cert.name}: coloop/full-B inequality failed")
        require(f <= H, f"{cert.name}: round size exceeds rank headroom")

        rounds.append({
            "u": u,
            "e": e,
            "rB": rb,
            "r": r,
            "H": H,
            "nu": nu,
            "c": c,
            "sigma": sigma,
            "bridges": z,
            "two_forest": forest_partition is not None,
            "dense_core_size": len(defect_witness[0]),
            "dense_core_rank": defect_witness[1],
            "f": f,
            "p": p,
            "r_next": r_next,
            "H_next": H_next,
            "fired": fired,
        })
        live = next_live

    r0 = rounds[0]["r"]
    H0 = rounds[0]["H"]
    sizes = tuple(rec["f"] for rec in rounds)
    ps = tuple(rec["p"] for rec in rounds)
    require(sum(sizes) == u0, f"{cert.name}: firing sizes do not sum to u")
    require(H0 == sum(ps), f"{cert.name}: H0 != sum p_j")
    require(r0 == u0 + sum(ps), f"{cert.name}: total rank identity failed")
    require(max(sizes) <= H0, f"{cert.name}: max-round bound failed")
    require(cert.cost >= r0, f"{cert.name}: paid rows below row rank")
    require(cert.cost >= u0 + max(sizes),
            f"{cert.name}: global round-profile bound failed")

    q = len(rounds)
    require(cert.cost * q >= u0 * (q + 1),
            f"{cert.name}: q-round score bound failed")

    if q == 2:
        f1, f2 = sizes
        p1, p2 = ps
        require(p2 == f2, f"{cert.name}: terminal p2=f2 failed")
        require(r0 == u0 + f2 + p1,
                f"{cert.name}: two-round rank identity failed")
        require(2 * cert.cost >= 3 * u0,
                f"{cert.name}: two-round 3/2 lower bound failed")
        if 3 * cert.cost <= 5 * u0:
            lo = 2 * u0 - cert.cost
            hi = cert.cost - u0
            require(lo <= f1 <= hi and lo <= f2 <= hi,
                    f"{cert.name}: <=5/3 balance window failed")

    return {
        "u": u0,
        "g": cert.cost,
        "score": Fraction(cert.cost, u0),
        "r0": r0,
        "slack": cert.cost - r0,
        "H0": H0,
        "sizes": sizes,
        "ps": ps,
        "rounds": rounds,
    }


def per_suffix_7_over_4():
    slots = [
        {
            ((1,), (1, 1)): (1, 2),
            ((1,), (1, 2)): (2, 1),
            ((1,), (2, 2)): (1, 1),
            ((1,), (2, 1)): (1, 2),
        },
        {
            ((1, 1), (2,)): (2, 1),
            ((2, 1), (2,)): (0, 1),
        },
        {
            ((2, 2, 1), ()): (2, 1),
            ((2, 1, 1), ()): (1, 2),
            ((1, 1, 1), ()): (1, 1),
            ((1, 2, 1), ()): (0, 1),
        },
    ]
    seeds = [
        ((1, 1, 2), (1, 0)),
        ((1, 1, 1), (0, 1)),
        ((1, 2, 1), (1, 1)),
        ((2, 1, 1), (1, 0)),
    ]
    return Mode3Instance([2, 2, 2], slots, [], seeds)


def quotient_12_over_7():
    base = per_suffix_7_over_4()
    slots = [dict(layer) for layer in base.slots]
    slots[0][((1,), (1, 1))] = MERGE
    seeds = list(base.Rbase[:-1])
    return Mode3Instance([2, 2, 2], slots, [], seeds)


def quotient_9_over_5():
    slots = [
        {
            ((1,), (1,)): (2, 1),
            ((2,), (1,)): (1, 2),
            ((2,), (2,)): MERGE,
            ((1,), (2,)): MERGE,
        },
        {
            ((1, 1), ()): (0, 1),
            ((3, 2), ()): (0, 1),
            ((3, 1), ()): MERGE,
            ((2, 2), ()): (1, 1),
            ((1, 2), ()): MERGE,
        },
    ]
    seeds = [
        ((1, 1), (1, 0)),
        ((1, 2), (1, 1)),
        ((2, 3), (1, 0)),
        ((3, 3), (1, 0)),
    ]
    return Mode3Instance([3, 3], slots, [], seeds)


def print_mode3_quotient(label, mi):
    classes = {}
    for v in mi.baseverts:
        classes.setdefault(mi.cid[mi.rep[v]], []).append(v)
    print(f"{label} quotient_classes={tuple(tuple(classes[c]) for c in sorted(classes))}")


def main():
    strict = strict_13_over_7()
    seven_four = per_suffix_7_over_4()
    twelve_seven = quotient_12_over_7()
    nine_five = quotient_9_over_5()

    require(seven_four.score() == Fraction(7, 4), "7/4 witness misencoded")
    require(twelve_seven.score() == Fraction(12, 7),
            "12/7 witness misencoded")
    require(nine_five.score() == Fraction(9, 5), "9/5 witness misencoded")

    print_mode3_quotient("7/4", seven_four)
    print_mode3_quotient("12/7", twelve_seven)
    print_mode3_quotient("9/5", nine_five)

    certs = (
        QuotientCertificate.from_strict("strict-13/7", strict),
        QuotientCertificate.from_mode3("per-suffix-7/4", seven_four),
        QuotientCertificate.from_mode3("quotient-12/7", twelve_seven),
        QuotientCertificate.from_mode3("quotient-9/5", nine_five),
    )
    for cert in certs:
        profile = audit_profile(cert)
        print(
            f"{cert.name}: score={profile['score']} u={profile['u']} "
            f"g={profile['g']} r0={profile['r0']} "
            f"row_slack={profile['slack']} H0={profile['H0']} "
            f"f={profile['sizes']} p={profile['ps']}"
        )
        for j, rec in enumerate(profile["rounds"], 1):
            print(
                f"  round{j}: u={rec['u']} rB={rec['rB']} r={rec['r']} "
                f"e={rec['e']} c={rec['c']} sigma={rec['sigma']} "
                f"bridges={rec['bridges']} H={rec['H']} nu={rec['nu']} "
                f"twoforest={'yes' if rec['two_forest'] else 'no'} "
                f"densecore={rec['dense_core_size']}/{rec['dense_core_rank']} "
                f"f={rec['f']} p={rec['p']} "
                f"r_next={rec['r_next']} fired={rec['fired']}"
            )
    print("ROUND_RANK_PROFILE_AUDIT_OK")


if __name__ == "__main__":
    main()
