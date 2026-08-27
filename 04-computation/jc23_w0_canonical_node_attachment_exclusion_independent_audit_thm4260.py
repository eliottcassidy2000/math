#!/usr/bin/env python3
"""Hostile clean-room referee for the JC 1512 canonical-node packet.

The class universe is rebuilt from the standard-library THM-4249 audit, not
from the packet's primary projector module.  Reciprocal denominators are then
computed with a small native polynomial/rational-function implementation,
not SymPy's rational-function field used by the packet.
"""

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import subprocess
import sys
import types

REV = "d526261faf8fb33a3b210b19c7d0039b8ef18ca8"
REPO = None
INDEPENDENT_PATH = (
    "04-computation/"
    "jc23_w0_cyclic_projector_squeeze_thm4249_independent_audit.py"
)


def require(ok, msg):
    if not ok:
        raise AssertionError(msg)


def load_independent():
    src = subprocess.run(
        ["git", "-C", REPO, "show", f"{REV}:{INDEPENDENT_PATH}"],
        check=True, stdout=subprocess.PIPE,
    ).stdout
    require(
        sha256(src).hexdigest()
        == "1c1ae0d47f5218af5978cb840c0f6f9c564a6df338a7b650700cbca774e5e3c4",
        "clean-room dependency hash changed",
    )
    mod = types.ModuleType("jc1512_ref_cleanroom_lattice")
    exec(compile(src, f"{REV}:{INDEPENDENT_PATH}", "exec"), mod.__dict__)
    return mod


# Polynomials are low-to-high tuples over F_q.
Q = 0


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] % Q == 0:
        a.pop()
    return tuple(x % Q for x in a)


def padd(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                 for i in range(n)])


def pneg(a):
    return trim([-x for x in a])


def psub(a, b):
    return padd(a, pneg(b))


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % Q
    return trim(out)


def pscale(a, c):
    return trim([c * x for x in a])


def pdivmod(a, b):
    a, b = list(trim(a)), trim(b)
    require(b != (0,), "polynomial division by zero")
    if len(a) < len(b):
        return (0,), tuple(a)
    out = [0] * (len(a) - len(b) + 1)
    inv = pow(b[-1], -1, Q)
    while len(a) >= len(b) and trim(a) != (0,):
        k = len(a) - len(b)
        c = a[-1] * inv % Q
        out[k] = c
        for j, x in enumerate(b):
            a[j + k] = (a[j + k] - c * x) % Q
        a = list(trim(a))
    return trim(out), trim(a)


def pgcd(a, b):
    a, b = trim(a), trim(b)
    while b != (0,):
        _, r = pdivmod(a, b)
        a, b = b, r
    return pscale(a, pow(a[-1], -1, Q)) if a != (0,) else (0,)


class RF:
    __slots__ = ("n", "d")

    def __init__(self, n=0, d=1):
        n = n.n if isinstance(n, RF) else ((n,) if isinstance(n, int) else n)
        d = d.d if isinstance(d, RF) else ((d,) if isinstance(d, int) else d)
        n, d = trim(n), trim(d)
        require(d != (0,), "rational denominator zero")
        if n == (0,):
            self.n, self.d = (0,), (1,)
            return
        g = pgcd(n, d)
        n, rn = pdivmod(n, g)
        d, rd = pdivmod(d, g)
        require(rn == (0,) and rd == (0,), "gcd division failed")
        inv = pow(d[-1], -1, Q)
        self.n, self.d = pscale(n, inv), pscale(d, inv)

    @staticmethod
    def coerce(x):
        return x if isinstance(x, RF) else RF(x)

    def __add__(self, other):
        other = RF.coerce(other)
        return RF(padd(pmul(self.n, other.d), pmul(other.n, self.d)),
                  pmul(self.d, other.d))

    __radd__ = __add__

    def __neg__(self):
        return RF(pneg(self.n), self.d)

    def __sub__(self, other):
        return self + (-RF.coerce(other))

    def __rsub__(self, other):
        return RF.coerce(other) - self

    def __mul__(self, other):
        other = RF.coerce(other)
        return RF(pmul(self.n, other.n), pmul(self.d, other.d))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = RF.coerce(other)
        require(other.n != (0,), "rational division by zero")
        return RF(pmul(self.n, other.d), pmul(self.d, other.n))

    def __rtruediv__(self, other):
        return RF.coerce(other) / self

    def __pow__(self, exponent):
        require(exponent >= 0, "negative RF power")
        answer, base, e = RF(1), self, exponent
        while e:
            if e & 1:
                answer = answer * base
            base = base * base
            e //= 2
        return answer

    def __eq__(self, other):
        other = RF.coerce(other)
        return self.n == other.n and self.d == other.d


def point_neg(point):
    return None if point is None else (point[0], -point[1])


def point_double(point, t, u):
    if point is None:
        return None
    a, b = point
    if b == 0:
        return None
    a2 = 9 * t * a**4 / (2 * (u - 1) * b**2) - 2 * a
    b2 = 3 * t * a**2 * (a - a2) / ((u - 1) * b) - b
    return a2, b2


def point_add(left, right, t, u):
    if left is None:
        return right
    if right is None:
        return left
    a1, b1 = left
    a2, b2 = right
    if a1 == a2:
        if b1 == -b2:
            return None
        require(b1 == b2, "equal X with incompatible Y")
        return point_double(left, t, u)
    slope = (b2 - b1) / (a2 - a1)
    c = (u - 1) / (2 * t)
    a3 = c * slope**2 - a1 - a2
    b3 = slope * (a1 - a3) - b1
    return a3, b3


def integer_multiple(k, point, t, u):
    if k < 0:
        return point_neg(integer_multiple(-k, point, t, u))
    answer, summand = None, point
    while k:
        if k & 1:
            answer = point_add(answer, summand, t, u)
        summand = point_double(summand, t, u)
        k //= 2
    return answer


def eisenstein_multiple(alpha, point, omega, t, u):
    m, n = alpha
    omega_point = (omega * point[0], point[1])
    return point_add(integer_multiple(m, point, t, u),
                     integer_multiple(n, omega_point, t, u), t, u)


def reciprocal(poly):
    n = len(poly) - 1
    out = [0] * (n + 1)
    for k, c in enumerate(poly):
        out[n - k] = c * (-1 if k % 2 else 1)
    return trim(out)


def build_representatives(ind):
    shells = ind.enumerate_shells()
    final_sets = {}
    for degree, shell in shells.items():
        vectors = {
            v for v in shell
            if v[0] == ind.ZERO
            and (v[3] == ind.ZERO or ind.e_norm(v[3]) >= 3)
            and ind.split_coordinates(v)[2] not in (1, 2)
            and not (degree == 42 and v[3] == ind.ZERO)
        }
        if degree == 42:
            vectors = {
                v for v in vectors
                if (ind.split_coordinates(v)[1], ind.split_coordinates(v)[2])
                != (3, 13)
            }
        final_sets[degree] = vectors
    reps = {degree: ind.orbit_partition(vectors)[0]
            for degree, vectors in final_sets.items()}
    require({d: len(v) for d, v in reps.items()} == {34: 176, 42: 104},
            "independent refined orbit universe changed")

    # Classwise second-node firewall, checked on the clean-room representatives.
    for degree in (34, 42):
        for v in reps[degree]:
            _, b, c, d = v
            rho = (ind.ZERO, ind.e_add(b, ind.e_mul(ind.OMEGA2, d)),
                   ind.e_add(c, d), ind.e_neg(d))
            t2 = ind.tau(ind.tau(v))
            rhs = ind.vector_scale(ind.e_neg(ind.OMEGA2), t2)
            require(rho == rhs, "rho=(-omega^2)T^2 failed")
    return reps


def audit_prime(ind, reps, embedding):
    global Q
    Q, z, p, scale = embedding
    require(pow(z, 4, Q) - pow(z, 2, Q) + 1 == 0, "bad zeta embedding")
    sqrt3 = (2 * z - z**3) % Q
    require((p*p - (1 + sqrt3)*p + 1) % Q == 0, "bad quartic embedding")
    sd = (2*p**3 + 3*p*p - 1) % Q
    require(sd and (pow(scale, 6, Q)*sd - 4) % Q == 0, "bad scale embedding")

    t = RF((0, 1)); u = t*t
    inv2 = pow(2, -1, Q)
    ascale = scale*scale*inv2 % Q
    bscale = scale**3*inv2 % Q
    omega = pow(z, 4, Q)
    f = (ascale*(u-p*p)/t, bscale*(u+p**3)/t)
    g = (ascale*(z*z % Q)*(p*p*u-1)/t,
         (-bscale*pow(z, 3, Q))*(1+p**3*u)/t)

    grouped = defaultdict(list)
    all_gcds = Counter()
    for degree in (34, 42):
        for vector in reps[degree]:
            _, _, k, a, b = ind.split_coordinates(vector)
            point = point_add(eisenstein_multiple(a, f, omega, t, u),
                              eisenstein_multiple(b, g, omega, t, u), t, u)
            require(point is not None, "hidden projection became zero")
            d = point[0].d
            require(len(d)-1 == 4*k-1, "denominator missed sharp bound")
            require(d[0] == 0 and d[1] != 0
                    and all(d[j] == 0 for j in range(0, len(d), 2)),
                    "denominator is not tD(t^2)")
            gcd = pgcd(d, reciprocal(d))
            all_gcds[gcd] += 1
            require(gcd == ((Q-1) % Q, 0, 1), "non-wall reciprocal gcd")
            # Packet digests list monic coefficients high-to-low.
            grouped[(degree, k)].append(",".join(map(str, reversed(d))))

    require(all_gcds == Counter({((Q-1) % Q, 0, 1): 280}),
            "gcd census changed")
    digests = {key: sha256(";".join(rows).encode("ascii")).hexdigest()
               for key, rows in grouped.items()}
    return digests


def main():
    global REPO
    parser = argparse.ArgumentParser(
        description="Independent clean-room audit for THM-4260.",
    )
    parser.add_argument(
        "--repo", type=Path, default=Path("."),
        help="repository root (default: current directory)",
    )
    args = parser.parse_args()
    REPO = args.repo.resolve()

    ind = load_independent()
    reps = build_representatives(ind)
    profiles = Counter()
    for degree in (34, 42):
        for v in reps[degree]:
            nd, k = ind.split_coordinates(v)[1:3]
            profiles[(degree, nd, k)] += 1
    print(f"classes=d34:{len(reps[34])} d42:{len(reps[42])} total:{sum(map(len,reps.values()))}")
    print(f"profiles={dict(sorted(profiles.items()))}")
    for embedding in ((397,157,161,27),(577,57,224,25)):
        digests = audit_prime(ind, reps, embedding)
        print(f"q={embedding[0]} groups={len(digests)} gcd=t^2-1 classes=280 PASS")
        for key in sorted(digests):
            print(f"digest q={embedding[0]} degree={key[0]} K={key[1]} {digests[key]}")
    print("CLEANROOM_REFEREE_PASS")


if __name__ == "__main__":
    main()
