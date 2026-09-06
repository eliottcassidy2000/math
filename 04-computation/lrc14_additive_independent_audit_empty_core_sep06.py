#!/usr/bin/env python3
"""Independent exact referee for the additive norm-three 6/55 family.

Universe: positive primitive a<b<c=a+b, all ternary units. The complete
min-projection head is c<=59; the complete physical-only head is c<=33.
No import of the concurrent additive proof producer is used. Native sheet
contacts come from the hash-pinned original Sep05 one-ray interval engine.
Run: python3 -B 04-computation/lrc14_additive_independent_audit_empty_core_sep06.py
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import gcd
import subprocess
import types

REVISION = "a80f2e90aa37de41ba93d206ae58125dfa1ee9c5"
FROZEN_PATH = "04-computation/lrc14_one_ray_overnight_hexagon_sep05.py"
FROZEN_HASH = "6b41a879700632aa934650f27dafe9d076c051ddcee3262fabc07556a7aaf875"
CHECKS = 0


def need(value, label):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(label)


def frozen_native():
    source = subprocess.run(["git", "show", REVISION + ":" + FROZEN_PATH],
                            check=True, capture_output=True).stdout
    need(sha256(source).hexdigest() == FROZEN_HASH, "frozen native source hash")
    module = types.ModuleType("frozen_sep05_native_for_additive_referee")
    exec(compile(source, "git:" + FROZEN_PATH, "exec"), module.__dict__)
    return module


def curves(t):
    r = Q(3, 14)
    def A(x):
        return min(2*r, (r*(2-t)-x)/(1-t))
    def B(x):
        return min(2*r, (r*(1+t)-x)/t)
    def C(x):
        return min(2*r, (r-x)/(t*(1-t)))
    def physical(x):
        return min(A(x), C(x))
    return A, B, C, physical


def integrate(curve, knots):
    return sum(((right-left)*(curve(left)+curve(right))/2
                for left, right in zip(knots, knots[1:])), Q(0))


def main():
    native = frozen_native()
    print("INDEPENDENT ADDITIVE REFEREE: exact native/ray/continuum audit")
    print("frozen_native_revision", REVISION)
    print("frozen_native_sha256", FROZEN_HASH)
    count = physical_head = 0
    min_equalities = []
    physical_equalities = []
    best_min = (Q(0), None)
    best_physical = (Q(0), None)
    semantic = sha256()
    target = Q(6, 55)
    r = Q(3, 14)
    for c in range(3, 60):
        for a in range(1, (c+1)//2):
            b = c-a
            if a >= b or gcd(a,b) != 1 or any(v % 3 == 0 for v in (a,b,c)):
                continue
            w = a,b,c
            count += 1
            physical_head += c <= 33
            ks = [k for k in range(1,c) if 14*k < 3*c and k % 3]
            predicted = {(sign*k,sign*k,-sign*k) for k in ks for sign in (-1,1)}
            raw = native.carriers(w)
            need(raw == predicted, "complete full carrier reduction " + str(w))
            t = Q(a,c)
            A,B,C,physical = curves(t)
            E = tuple(2*sum((curve(Q(k,c)) for k in ks),Q(0))/c for curve in (A,B,C))
            mass = 2*sum((physical(Q(k,c)) for k in ks),Q(0))/c
            literal = native.literal_projection_data(w)
            need(literal == (E,mass), "native pair-then-third interval equality " + str(w))
            need(native.projection_data(w,raw) == literal, "raw/native independent equation agreement")
            need(E[0] <= E[1] and mass <= min(E), "projection/physical type order")
            need(min(E) <= target and mass <= target, "exact finite-head target")
            if min(E) == target:
                min_equalities.append(w)
            if mass == target:
                physical_equalities.append(w)
            if min(E) > best_min[0]:
                best_min = min(E),w
            if mass > best_physical[0]:
                best_physical = mass,w

            # Every curve is linear between these exact knots. Include both
            # the capped-C kink and the A/C physical crossover separately.
            knots = sorted({Q(0),r*t,r*(1-t),r*(1-2*t+2*t*t),r})
            baselines = tuple(Q(4,3)*integrate(curve,knots) for curve in (A,B,C,physical))
            expected = ((9+3*t)/98,(12-3*t)/98,(12-12*t+12*t*t)/98,Q(9,98))
            need(baselines == expected, "independent exact piecewise integrals")
            for value, baseline in zip((*E,mass),baselines):
                need(value < baseline + Q(4,7*c), "strict residue quadrature on head")
            need(min(baselines[0],baselines[2]) <= Q(39,392), "continuum min bound")
            semantic.update(repr((w,tuple(sorted(raw)),E,mass,baselines)).encode())

    need(count == 136 and physical_head == 42, "complete additive head counts")
    need(min_equalities == [(1,10,11)] and physical_equalities == [(1,10,11)],
         "complete equality loci")
    need(Q(6,55)-(Q(39,392)+Q(4,7*60)) == Q(1,12936), "strict min tail at c60")
    need(Q(6,55)-(Q(9,98)+Q(4,7*34)) == Q(41,91630), "strict physical tail at c34")
    print("COMPLETE_MIN_HEAD c<=59 rows",count,"max",best_min,"equalities",min_equalities)
    print("COMPLETE_PHYSICAL_HEAD c<=33 rows",physical_head,"max",best_physical,
          "equalities",physical_equalities)
    print("CONTINUUM min_E_max=39/392 at a/c=1/4; physical=9/98 independently of a/c")
    print("TAIL_MARGINS min_E at c60=1/12936; physical at c34=41/91630")
    E,mass = native.literal_projection_data((4,7,11))
    need(E == (Q(6,49),Q(41,308),Q(223,2156)) and mass == Q(215,2156),
         "strict physical/projection distinction")
    need(mass < min(E), "physical is not minimum projection")
    print("TYPE_HOSTILE (4,7,11) E",tuple(map(str,E)),"physical",mass,
          "physical_strictly_below_minE",mass < min(E))
    print("SEMANTIC_SHA256",semantic.hexdigest())
    print("CHECKS",CHECKS,"frozen_native_checks",native.CHECKS)
    print("PASS: supports, exact heads, piecewise integrals, strict tail margins, unique equality")
    print("SCOPE: additive primitive positive ternary-unit triples only; general mixed parity remains open")


if __name__ == "__main__":
    main()
