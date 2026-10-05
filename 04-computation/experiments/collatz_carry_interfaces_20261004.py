#!/usr/bin/env python3
"""Exact audits of carry interfaces; standard library, no convergence oracle.

Run with python3 -B and python3 -O -B. Both write the same JSON and text.
The accompanying note states the all-parameter proofs and finite universes.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import factorial, gcd
from pathlib import Path
import json

CHECKS = Counter()


def check(test, group):
    CHECKS[group] += 1
    if not test:
        raise RuntimeError(f"failed {group}, check {CHECKS[group]}")


@dataclass(frozen=True)
class Affine:
    p: int = 1
    q: int = 1
    b: int = 0

    def then(self, other):
        return Affine(other.p * self.p, other.q * self.q,
                      other.p * self.b + other.b * self.q)

    def at(self, x):
        return F(self.p * x + self.b, self.q)

    def ports(self):
        c = ((self.q - self.b) * pow(self.p, -1, 2 * self.q)) % (2 * self.q)
        d = self.at(c)
        check(d.denominator == 1, "port_integrality")
        return c, int(d)


LETTERS = {
    "H": Affine(729, 1024, 669),
    "G": Affine(9, 8, 5),
    "A": Affine(81, 128, 85),
    "B": Affine(81, 128, 73),
    "L": Affine(9, 16, -3),
}
VALUATIONS = {"G": (1, 2), "A": (1, 2, 1, 3), "B": (1, 1, 2, 3)}
AWARD = {"H": 2, "G": -1, "A": 1, "B": 1, "L": 4}
TAGS = {3: "H", 4: "G", 2: "A", 5: "B", 6: "L"}
H_PREFIX = (1, 2, 1, 1, 1, 2)


def encode(word):
    out = Affine()
    for letter in word:
        out = out.then(LETTERS[letter])
    return out


def decode_translation(t):
    """Exact finite recognition as well as decoding, using dyadic cost."""
    if not isinstance(t, F):
        raise ValueError("exact rational required")
    reverse = []
    while t:
        b, q = t.numerator, t.denominator
        if q & (q - 1):
            raise ValueError("not dyadic")
        tag = b * pow(q, -1, 9) % 9
        if tag not in TAGS:
            raise ValueError("no last-letter tag")
        letter = TAGS[tag]
        op = LETTERS[letter]
        if q % op.q:
            raise ValueError("insufficient denominator cost")
        old_q = q // op.q
        old_b, rem = divmod(b - op.b * old_q, op.p)
        if rem:
            raise ValueError("invalid carry predecessor")
        reverse.append(letter)
        t = F(old_b, old_q)
    return "".join(reversed(reverse))


def valuation(z):
    if z == 0:
        raise ValueError("zero has no finite valuation")
    return (z & -z).bit_length() - 1


def replay(x, word, first_hit=False):
    for a in word:
        if first_hit and x == 1:
            raise ValueError("root padding")
        z = 3 * x + 1
        if x % 2 != 1 or valuation(z) != a:
            raise ValueError("wrong exact valuation")
        x = z // (1 << a)
    return x


def macro_replay(n, word):
    x = n
    for letter in word:
        op = LETTERS[letter]
        c, _ = op.ports()
        if x % (2 * op.q) != c:
            raise ValueError("wrong macro guard")
        y = op.at(x)
        if y.denominator != 1 or int(y) % 2 != 1:
            raise ValueError("wrong macro endpoint")
        y = int(y)
        if letter == "H":
            z = replay(x, H_PREFIX)
            check(z == 4 * y + 1, "H_dependency")
            check(replay(z, (valuation(3 * z + 1),)) ==
                  replay(y, (valuation(3 * y + 1),)), "H_dependency")
        else:
            check(replay(x, VALUATIONS[letter]) == y, "native_word")
        x = y
    return x


def compiled_common_future(word, y):
    tail = (valuation(3 * y + 1),)
    for letter in reversed(word):
        if letter == "H":
            tail = H_PREFIX + (tail[0] + 2,) + tail[1:]
        else:
            tail = VALUATIONS[letter] + tail
    return tail


def compose_ports(u, v):
    c, d = u.ports()
    e, f = v.ports()
    t = ((e - d) // 2 * pow(u.p, -1, v.q)) % v.q
    s = (d + 2 * u.p * t - e) // (2 * v.q)
    return c + 2 * u.q * t, f + 2 * v.p * s


def words(alphabet, depth):
    yield ""
    for size in range(1, depth + 1):
        for w in product(alphabet, repeat=size):
            yield "".join(w)


def audit_macros():
    translations = {}
    funded, cases = 0, 0
    for w in words("HGAB", 6):
        op = encode(w)
        t = F(op.b, op.q)
        check(decode_translation(t) == w, "translation_decoder")
        check(t not in translations, "translation_injective")
        translations[t] = w
        c, d = op.ports()
        check(0 < c < 2 * op.q and c % 2 == 1 and
              0 < d < 2 * op.p and d % 2 == 1, "canonical_ports")
        check(op.b == op.q * d - op.p * c, "port_carry")
        for lift in (0, 1):
            cases += 1
            n = c + 2 * op.q * lift
            y = macro_replay(n, w)
            check(y == d + 2 * op.p * lift, "port_progression")
            check(macro_replay(n - 2 * op.q * (lift + 1), w) == d - 2 * op.p,
                  "negative_port_control")
            if n == 1 and not w:
                check(y == 1, "root_identity")
            else:
                receipt = compiled_common_future(w, y)
                check(replay(n, receipt, first_hit=True) ==
                      replay(y, (valuation(3 * y + 1),), first_hit=True),
                      "compiled_common_future")
            credit, valid, x = 0, True, n
            for letter in w:
                old_e = (x + 5) * F(9, 8) ** credit
                credit += AWARD[letter]
                valid = valid and credit >= 0
                x = int(LETTERS[letter].at(x))
                new_e = (x + 5) * F(9, 8) ** credit
                bound = F(1) if letter == "G" else F(15, 16)
                check(new_e <= bound * old_e, "inherited_credit_interface")
                if valid:
                    check(x < n and new_e < n + 5, "paid_prefix")
            if w and valid:
                funded += 1
    short = [encode(w) for w in words("HGAB", 3)]
    for u in short:
        for v in short:
            check(compose_ports(u, v) == u.then(v).ports(), "port_composition")
    invalid = [F(-1), F(1, 3), F(1), F(1, 2), F(3, 8), F(1, 1024)]
    for t in invalid:
        try:
            decode_translation(t)
        except ValueError:
            CHECKS["decoder_hostile"] += 1
        else:
            raise RuntimeError(f"accepted hostile {t}")
    check(85 * pow(128, -1, 3) % 3 == 73 * pow(128, -1, 3) % 3,
          "mod3_tag_collision")
    counts = Counter()
    for w in words("HGAB", 8):
        balance, minimum = 0, 0
        for l in w:
            balance += AWARD[l]
            minimum = min(minimum, balance)
        if balance == 0 and minimum == 0:
            counts[w.count("H"), w.count("A") + w.count("B")] += 1
    for (h, p), count in counts.items():
        expected = (2 ** p * factorial(3 * h + 2 * p) //
                    (factorial(h) * factorial(p) * factorial(2 * h + p + 1)))
        check(count == expected, "mixed_tree_count")
    return {"words": len(translations), "positive_sources": cases,
            "funded_nonempty_sources": funded,
            "port_composition_pairs": len(short) ** 2,
            "tags_mod9": TAGS,
            "clock_collision": {w: encode(w).__dict__ for w in ("H", "GA", "GB", "AG", "BG")}}


def ordinary(word, q=3):
    out = Affine()
    for a in word:
        out = out.then(Affine(q, 1 << a, 1))
    return out


@dataclass(frozen=True)
class GuardedAP:
    c: int = 1
    m: int = 2
    d: int = 1
    ell: int = 2

    def then(self, other):
        g = gcd(self.ell, other.m)
        if (other.c - self.d) % g:
            return None
        period = other.m // g
        t = ((other.c - self.d) // g * pow(self.ell // g, -1, period)) % period
        s = (self.d + self.ell * t - other.c) // other.m
        return GuardedAP(self.c + self.m * t, self.m * period,
                         other.d + other.ell * s, self.ell * other.ell // g)

    def at(self, n):
        if (n - self.c) % self.m:
            raise ValueError("outside native guard")
        return self.d + self.ell * ((n - self.c) // self.m)

    def affine(self):
        slope = F(self.ell, self.m)
        shift = F(self.d) - slope * self.c
        return slope, shift


NATIVE = {"H": GuardedAP(155, 2048, 111, 1458),
          "G": GuardedAP(11, 16, 13, 18),
          "A": GuardedAP(187, 256, 119, 162),
          "B": GuardedAP(7, 256, 5, 162),
          "L": GuardedAP(219, 256, 123, 144)}


def guarded_encode(word):
    ap = GuardedAP()
    for letter in word:
        if ap is None:
            return None
        ap = ap.then(NATIVE[letter])
    return ap


def extended_native_step(x, letter):
    y = NATIVE[letter].at(x)
    check(LETTERS[letter].at(x) == y, "extended_native_affine")
    if letter != "L":
        check(macro_replay(x, letter) == y, "extended_native_replay")
    else:
        k = replay(y, (1, 2), first_hit=True)
        a = valuation(3 * k + 1)
        check(replay(x, (1, 2, 1, 1, a + 2), first_hit=True) ==
              replay(y, (1, 2, a), first_hit=True), "L_local_common_future")
    return y


def audit_extended_library():
    translations = set()
    legal, empty, funded = 0, 0, 0
    formal_balanced, legal_balanced = Counter(), Counter()
    for w in words("HGABL", 6):
        op = encode(w)
        shift = F(op.b, op.q)
        check(decode_translation(shift) == w, "five_letter_decoder")
        check(shift not in translations, "five_letter_injective")
        translations.add(shift)
        ap = guarded_encode(w)
        balance, minimum = 0, 0
        for letter in w:
            balance += AWARD[letter]
            minimum = min(minimum, balance)
        if balance == 0 and minimum == 0:
            index = (w.count("H"), w.count("A") + w.count("B"), w.count("L"))
            formal_balanced[index] += 1
            if ap is not None:
                legal_balanced[index] += 1
        if ap is None:
            empty += 1
            continue
        legal += 1
        check(ap.affine() == (F(op.p, op.q), shift), "native_port_affine")
        check(0 < ap.c < ap.m and 0 < ap.d < ap.ell, "native_port_sign_boundary")
        for lift in (0, 1):
            n = ap.c + ap.m * lift
            x, credit, valid = n, 0, True
            for letter in w:
                old_e = (x + 5) * F(9, 8) ** credit
                x = extended_native_step(x, letter)
                credit += AWARD[letter]
                valid = valid and credit >= 0
                new_e = (x + 5) * F(9, 8) ** credit
                factor = F(1) if letter == "G" else F(15, 16)
                check(new_e <= factor * old_e, "five_letter_credit_interface")
                if valid:
                    check(x < n, "five_letter_paid_prefix")
            check(x == ap.at(n), "native_port_endpoint")
            if w and valid:
                funded += 1
    for (h, p, l), count in formal_balanced.items():
        expected = (2 ** p * factorial(3 * h + 2 * p + 5 * l) //
                    (factorial(h) * factorial(p) * factorial(l) * factorial(2 * h + p + 4 * l + 1)))
        check(count == expected, "five_letter_tree_count")
    check(guarded_encode("LB") is None and guarded_encode("LBGGGGG") is None,
          "formal_tree_native_empty")
    k = guarded_encode("LG")
    check(k == GuardedAP(219, 256, 139, 162), "LG_equals_K_guarded")
    check(encode("LG") == Affine(81, 128, 53), "LG_equals_K_affine")
    check(LETTERS["L"].at(27) == 15 and (27 - NATIVE["L"].c) % NATIVE["L"].m != 0,
          "cancelled_guard_hostile")
    small = [guarded_encode(w) for w in words("HGABL", 2)]
    def comp(a, b):
        return None if a is None or b is None else a.then(b)
    for a in small:
        for b in small:
            for c in small:
                check(comp(comp(a, b), c) == comp(a, comp(b, c)), "native_port_associativity")
    return {"formal_words": len(translations), "nonempty_native_guards_including_identity": legal,
            "empty_native_guards": empty, "funded_source_replays": funded,
            "native_associativity_triples": len(small) ** 3,
            "new_letter": {"name": "L", "affine": LETTERS["L"].__dict__,
                           "native_guard": NATIVE["L"].__dict__, "tag_mod9": 6},
            "explicit_empty_program": "LBGGGGG", "relation_if_K_is_added": "K=LG"}


def compositions(total, length):
    for cuts in combinations(range(1, total), length - 1):
        points = (0,) + cuts + (total,)
        yield tuple(b - a for a, b in zip(points, points[1:]))


def audit_fourier():
    tested = 0
    for q in (3, 5, 7, 9):
        for m in range(1, 6):
            for cost in range(m, 15):
                addresses = set()
                for w in compositions(cost, m):
                    op = ordinary(w, q)
                    c, d = op.ports()
                    r = (-op.b * pow(op.p, -1, op.q)) % op.q
                    check(r not in addresses, "fixed_cost_address_injective")
                    addresses.add(r)
                    check(0 < d < 2 * op.p and gcd(d, q) == 1,
                          "general_q_target_port")
                    for u in (1, 2, 5):
                        phase = F((u * op.b * pow(op.q, -1, op.p)) % op.p, op.p)
                        check(phase == F(u * d, op.p) % 1, "target_phase")
                        check(phase == (F(u * r, op.q) + F(u * op.b, op.p * op.q)) % 1,
                              "source_phase_bridge")
                    # Expand the incoming Fourier recursion from its top level.
                    rw = tuple(reversed(w))
                    s = 0
                    freq = F(0)
                    for j, a in enumerate(rw, 1):
                        s += a
                        freq += F(1, q ** (m - j + 1) * 2 ** s)
                    check(freq == F(op.b, op.p * op.q), "reversed_phase_path")
                    bound = (F(1, 2 ** m) - F(1, q ** m)) / (q - 2)
                    check(0 < freq <= bound, "uniform_phase_correction")
                    check((freq == bound) == all(a == 1 for a in w),
                          "phase_correction_equality")
                    tested += 1
    # Independent exact probability recursion, finite valuation cutoff a<=3.
    for q in (3, 5, 7):
        law = {0: F(1)}
        for m in range(1, 5):
            modulus = q ** m
            new = defaultdict(F)
            for y, mass in law.items():
                for a in range(1, 4):
                    new[((q * y + 1) * pow(2 ** a, -1, modulus)) % modulus] += mass / 2 ** a
            brute = defaultdict(F)
            for w in product(range(1, 4), repeat=m):
                op = ordinary(w, q)
                _, d = op.ports()
                brute[d % modulus] += F(1, op.q)
            check(dict(new) == dict(brute), "Fourier_law_independent_DP")
            check(sum(new.values()) == F(7, 8) ** m, "truncated_law_mass")
            law = new
    collision = []
    for w in ((1, 7), (7, 1)):
        op = ordinary(w)
        c, d = op.ports()
        collision.append({"word": w, "carrier": op.__dict__, "source": c, "target": d})
    check(collision[0]["target"] == collision[1]["target"] and
          collision[0]["source"] != collision[1]["source"], "phase_loses_source")
    w = (4, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3)
    op = ordinary(w)
    c, d = op.ports()
    check((c, d) == (165, 167) and op.q > op.p, "height_hostile")
    check(replay(c, w, True) == d and d > c, "height_hostile")
    check(replay(c + 2 * op.q, w, True) == d + 2 * op.p < c + 2 * op.q,
          "height_hostile")
    return {"word_cases": tested, "universe": "q=3,5,7,9; lengths 1..5; total cost <=14",
            "proved_phase_error_bound": "2*pi*abs(u)*(2^(-m)-q^(-m))/(q-2)",
            "phase_collision": collision,
            "height_hostile": {"word": w, "carrier": op.__dict__,
                               "unpaid": [c, d], "paid": [c + 2 * op.q, d + 2 * op.p]}}


def divide_bits(x, divisor, width, state=0):
    out = 0
    for d in range(width):
        b = (x >> d) & 1
        bp = (b - state) % 2
        state = (state + divisor * bp - b) // 2
        check(0 <= state < divisor, "division_state_range")
        out |= bp << d
    return out, state


def audit_phase_grid():
    for q in (3, 5, 7, 9, 11):
        for seed in (-17, -5, -1, 0, 1, 7, 19):
            for n in range(20):
                for d in range(1, 33):
                    M = 1 << d
                    r = seed * pow(q, -n, M) % M
                    s = seed * pow(q, -n - 1, M) % M
                    c = (q * s - r) // M
                    rr = seed * pow(q, -n, 2 * M) % (2 * M)
                    ss = seed * pow(q, -n - 1, 2 * M) % (2 * M)
                    cc = (q * ss - rr) // (2 * M)
                    b, bp = (rr - r) // M, (ss - s) // M
                    check(c + q * bp == b + 2 * cc, "carry_square")
                    check(c == -r * pow(M, -1, q) % q, "root_branch_digit")
                r = seed * pow(q, -n, 1 << 32) % (1 << 32)
                s, c = divide_bits(r, q, 32)
                check(s == seed * pow(q, -n - 1, 1 << 32) % (1 << 32),
                      "phase_transducer")
    counts = {}
    for q in (3, 5, 7):
        for depth in range(1, 5):
            M = q ** depth
            width = (2 * M).bit_length()
            modulus = 1 << width
            states = set()
            for x in range(1, modulus, 2):
                y = x * pow(M, -1, modulus) % modulus
                states.add((M * y - x) // modulus)
            check(states == set(range(M)), "all_composite_states_reachable")
            suffix_width = M.bit_length()
            suffix_modulus = 1 << suffix_width
            outputs = {(-c * pow(M, -1, suffix_modulus)) % suffix_modulus for c in states}
            check(len(outputs) == M, "composite_states_distinguishable")
            counts[f"{q}^{depth}"] = len(states)
    table = []
    for c in range(3):
        for b in (0, 1):
            bp = (b - c) % 2
            table.append([c, b, bp, (c + 3 * bp - b) // 2])
    return {"q3_table_columns": ["carry", "input_bit", "output_bit", "next_carry"],
            "q3_table": table, "minimal_state_counts_controls": counts,
            "boundary": "finite-state division at one level does not give bounded memory for all levels"}


CENTRE = (F(1, 3),) * 3


def fan(c, x):
    a, b = (c + 1) % 3, (c + 2) % 3
    out = [x[2] / 3] * 3
    out[a] += x[0]
    out[b] += x[1]
    return tuple(out)


def flag(sigma, x):
    a, b, c = sigma
    out = [F(0)] * 3
    out[a] = x[0] + x[1] / 2 + x[2] / 3
    out[b] = x[1] / 2 + x[2] / 3
    out[c] = x[2] / 3
    return tuple(out)


def point_word(w, family):
    x = CENTRE
    for letter in reversed(w):
        x = fan(letter, x) if family == "fan" else flag(letter, x)
    return x


def point_step(x, family):
    if family == "fan":
        if x.count(min(x)) != 1:
            raise ValueError("no unique fan minimum")
        c = x.index(min(x))
        a, b = (c + 1) % 3, (c + 2) % 3
        return c, (x[a] - x[c], x[b] - x[c], 3 * x[c])
    sigma = tuple(sorted(range(3), key=lambda i: x[i], reverse=True))
    a, b, c = sigma
    if not x[a] > x[b] > x[c]:
        raise ValueError("no strict flag order")
    return sigma, (x[a] - x[b], 2 * (x[b] - x[c]), 3 * x[c])


def decode_point(x, family, budget):
    w = []
    while x != CENTRE:
        if len(w) >= budget:
            raise ValueError("depth budget exhausted")
        letter, x = point_step(x, family)
        if min(x) <= 0 or sum(x) != 1:
            raise ValueError("outside the strict interior")
        w.append(letter)
    return tuple(w)


def audit_geometric_codes():
    counts = {}
    for family, alphabet, depth in (("fan", tuple(range(3)), 7),
                                    ("flag", tuple(permutations(range(3))), 5)):
        seen = {}
        for size in range(depth + 1):
            for w in product(alphabet, repeat=size):
                point = point_word(w, family)
                check(point not in seen, "marked_centroid_injective")
                seen[point] = w
                check(decode_point(point, family, size) == w, "marked_centroid_decoder")
        counts[family] = len(seen)
    cycle = [(F(1, 9), F(2, 9), F(6, 9)),
             (F(1, 9), F(5, 9), F(3, 9)),
             (F(4, 9), F(2, 9), F(3, 9))]
    for i, x in enumerate(cycle):
        check(point_step(x, "fan")[1] == cycle[(i + 1) % 3], "triadic_cycle_hostile")
    try:
        decode_point(cycle[0], "fan", 20)
    except ValueError:
        CHECKS["triadic_cycle_rejected"] += 1
    else:
        raise RuntimeError("periodic point accepted as finite address")
    # The exact mixed-history matrix collision also collapses marked centroids.
    def bisect(x):
        return (x[0] + x[1] / 2, x[1] / 2, x[2])
    def mixed(w):
        x = CENTRE
        for l in reversed(w):
            x = bisect(x) if l == "E" else fan(int(l), x)
        return x
    check(mixed("101EE") == mixed("EE011"), "mixed_history_collision")
    return {"centroid_word_controls": counts,
            "triadic_period_three": [[str(t) for t in row] for row in cycle],
            "mixed_collision_centroid": [str(t) for t in mixed("101EE")]}


def main():
    result = {"scope": "PROVED identities in the note; finite exact controls here; Collatz and H1 remain open",
              "macros": audit_macros(), "extended_library": audit_extended_library(),
              "Fourier": audit_fourier(),
              "phase_grid": audit_phase_grid(), "geometry": audit_geometric_codes()}
    result["checks"] = dict(sorted(CHECKS.items()))
    result["total_checks"] = sum(CHECKS.values())
    lines = ["Exact carry interface audit",
             f"Four-letter translation codes: {result['macros']['words']}",
             f"Positive native source replays: {result['macros']['positive_sources']}",
             f"Funded nonempty source replays: {result['macros']['funded_nonempty_sources']}",
             f"Independent port composition pairs: {result['macros']['port_composition_pairs']}",
             f"Five-letter codes: {result['extended_library']['formal_words']}; empty native guards: {result['extended_library']['empty_native_guards']}",
             f"Ordinary-word phase and source controls: {result['Fourier']['word_cases']}",
             f"Marked-centroid codes: {result['geometry']['centroid_word_controls']}",
             f"Total explicit checks: {result['total_checks']}",
             "All checks passed. No universal coverage, deterministic cancellation, or root completion inferred."]
    root = Path(__file__).resolve().parents[2]
    dest = root / "05-knowledge" / "results"
    stem = "collatz_carry_interfaces_20261004"
    (dest / (stem + ".json")).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    output = "\n".join(lines) + "\n"
    (dest / (stem + ".out")).write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
