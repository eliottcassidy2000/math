#!/usr/bin/env python3
"""Exact finite referee for the THM-2681 A4 resolvent corollary.

This companion distinguishes three objects which all have a three-letter
shadow:

* the A4 action on the three 2+2 matchings of four sheets;
* the S4 action on those matchings; and
* the actual THM-1310 cubic root action, whose Galois closure is S3 and whose
  root stabilizer is a transposition.

It also checks the mod-two divisor vectors of the tempting local Kummer
classes r1/r3 and r2/r3.  The geometric facts that the divisors (ri=0) are
prime and that the global ordered-root normalization is a UFD are proved in
THM-2681, not inferred from this finite calculation.
"""

from itertools import permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compose(p, q):
    """p after q."""
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    ans = [0] * len(p)
    for i, j in enumerate(p):
        ans[j] = i
    return tuple(ans)


def order(p):
    identity = tuple(range(len(p)))
    x = identity
    for n in range(1, 100):
        x = compose(p, x)
        if x == identity:
            return n
    raise RuntimeError("permutation order exceeded finite bound")


def parity(p):
    inversions = sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))
    return inversions % 2


def cycle_lengths(p):
    seen = set()
    lengths = []
    for i in range(len(p)):
        if i in seen:
            continue
        j = i
        length = 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = p[j]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def normalize_matching(blocks):
    return tuple(sorted(tuple(sorted(block)) for block in blocks))


def act_matching(p, matching):
    return normalize_matching(tuple(tuple(p[i] for i in block) for block in matching))


def action_image(group, objects, action):
    index = {obj: i for i, obj in enumerate(objects)}
    return {
        tuple(index[action(g, obj)] for obj in objects)
        for g in group
    }


def stabilizer(group, obj, action):
    return {g for g in group if action(g, obj) == obj}


def is_normal(subgroup, group):
    return all(
        compose(compose(g, h), inverse(g)) in subgroup
        for g in group
        for h in subgroup
    )


def core(subgroup, group):
    return {
        h
        for h in subgroup
        if all(compose(compose(g, h), inverse(g)) in subgroup for g in group)
    }


S4 = set(permutations(range(4)))
A4 = {p for p in S4 if parity(p) == 0}
identity4 = tuple(range(4))
V4 = {p for p in S4 if p == identity4 or cycle_lengths(p) == (2, 2)}

MATCHINGS = tuple(
    normalize_matching(m)
    for m in (
        ((0, 1), (2, 3)),
        ((0, 2), (1, 3)),
        ((0, 3), (1, 2)),
    )
)

image_s4 = action_image(S4, MATCHINGS, act_matching)
image_a4 = action_image(A4, MATCHINGS, act_matching)
kernel_s4 = {g for g in S4 if all(act_matching(g, m) == m for m in MATCHINGS)}
kernel_a4 = kernel_s4 & A4

require(len(S4) == 24 and len(A4) == 12 and len(V4) == 4, "group orders")
require(kernel_s4 == V4 and kernel_a4 == V4, "matching kernel is not V4")
require(len(image_s4) == 6, "S4/V4 matching image is not S3")
require(len(image_a4) == 3, "A4/V4 matching image is not C3")
require({order(p) for p in image_a4} == {1, 3}, "A4 quotient is not cyclic order three")

stab_a4 = stabilizer(A4, MATCHINGS[0], act_matching)
stab_s4 = stabilizer(S4, MATCHINGS[0], act_matching)
require(stab_a4 == V4 and is_normal(stab_a4, A4), "A4 resolvent root stabilizer")
require(len(stab_s4) == 8 and not is_normal(stab_s4, S4), "S4 resolvent root stabilizer")
require(core(stab_s4, S4) == V4, "S4 matching stabilizer core")

# The actual THM-1310 cubic action proved in THM-2681 is the natural S3
# action with root stabilizer a transposition.  Recheck its field-type
# fingerprint independently of the polynomial algebra.
S3 = set(permutations(range(3)))
root_stab = {p for p in S3 if p[0] == 0}
require(len(root_stab) == 2, "cubic root stabilizer order")
require(not is_normal(root_stab, S3), "THM-1310 cubic root field unexpectedly normal")
require(core(root_stab, S3) == {tuple(range(3))}, "THM-1310 root stabilizer core")

# Mod-two valuations along D1,D2,D3 of r1/r3 and r2/r3.  Their span is the
# even-weight standard plane and is stable under every coordinate permutation.
v13 = (1, 0, 1)
v23 = (0, 1, 1)
local_plane = {
    (0, 0, 0),
    v13,
    v23,
    tuple(a ^ b for a, b in zip(v13, v23)),
}
even_plane = {v for v in ((a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)) if sum(v) % 2 == 0}
require(local_plane == even_plane and len(local_plane) == 4, "local Kummer plane rank")
for p in S3:
    require({tuple(v[p[i]] for i in range(3)) for v in local_plane} == local_plane,
            "local Kummer plane is not S3-stable")
require(any(v != (0, 0, 0) for v in local_plane), "local Kummer hostile missing")

print("THM-2681 A4 RESOLVENT FIELD-TYPE REFEREE")
print(f"orders S4/A4/V4={len(S4)}/{len(A4)}/{len(V4)}")
print(f"matching images S4/A4={len(image_s4)}/{len(image_a4)}")
print(f"matching kernels S4/A4={len(kernel_s4)}/{len(kernel_a4)}")
print(f"A4 matching stabilizer={len(stab_a4)}, normal={int(is_normal(stab_a4, A4))}")
print(f"S4 matching stabilizer={len(stab_s4)}, core={len(core(stab_s4, S4))}")
print(f"THM1310 root stabilizer={len(root_stab)}, normal={int(is_normal(root_stab, S3))}, core={len(core(root_stab, S3))}")
print(f"local divisor Kummer plane size={len(local_plane)}")
print("VERDICT: A4 gives a cyclic cubic; the actual THM-1310 cubic is nonnormal S3")
print("ALL CHECKS PASSED")
