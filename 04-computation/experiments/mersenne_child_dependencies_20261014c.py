"""Source-owned Mersenne dependency transport and fixed inverse-word obstruction.

Default execution discovers no Collatz orbit. Two existing compressed ROOT
words are imported and strictly replayed; derived words are exact splices.
The inverse-word census decides all exponents by a one-candidate theorem.
"""
from dataclasses import dataclass, replace
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json

import collatz_child_normal_forms_20261007 as inverse
import collatz_two_twos_grounding_20261007c as frozen


CHECKS = 0
BIT_CAP = 100000


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def natural(n, least=0):
    need(type(n) is int and n >= least, "exact integer in declared range")
    return n


def letters(word):
    need(type(word) is tuple and all(type(a) is int and a > 0 for a in word),
         "exact tuple of positive valuations")
    return word


def mersenne(E, bit_cap=BIT_CAP):
    natural(E, 1)
    natural(bit_cap, 1)
    need(E <= bit_cap, "explicit materialization cap")
    return (1 << E)-1


def v2(n):
    natural(n, 1)
    return (n & -n).bit_length()-1


def carrier(word):
    letters(word)
    P = Q = 1
    B = 0
    for a in word:
        P, Q, B = 3*P, Q*(1 << a), 3*B+Q
    return P, Q, B


@dataclass(frozen=True)
class Root:
    exponent: int
    word: tuple


def replay(E, word, root=False):
    n = mersenne(E)
    letters(word)
    need(type(root) is bool, "exact ROOT flag")
    peak = E
    for a in word:
        need(n != 1, "strict first-hit boundary")
        z = 3*n+1
        need((z & -z).bit_length()-1 == a, "actual source valuation")
        n = z >> a
        peak = max(peak, n.bit_length())
    need(not root or n == 1, "supplied word reaches ROOT")
    return n, peak


def audit_root(proof):
    need(type(proof) is Root, "exact ROOT record")
    replay(proof.exponent, proof.word, True)
    return proof


@dataclass(frozen=True)
class Join:
    source_exponent: int
    child_exponent: int
    left: tuple
    right: tuple
    endpoint: int


def audit_join(join):
    need(type(join) is Join, "exact join record")
    natural(join.source_exponent, 2)
    natural(join.child_exponent, 1)
    natural(join.endpoint, 1)
    need(join.source_exponent > join.child_exponent, "strict original-source payment")
    need(replay(join.source_exponent, join.left)[0] == join.endpoint,
         "same source-owned endpoint")
    need(replay(join.child_exponent, join.right)[0] == join.endpoint,
         "same child-owned endpoint")
    return join


def even_join(E):
    """Inherited reset-switch specialization; E=2 is a direct ROOT boundary."""
    natural(E, 4)
    need(E % 2 == 0, "even exponent required")
    mersenne(E)
    s = v2(E)
    left = (1,)*(E-1)+(3+s,)
    right = (1,)*(E-2)+(2, 1+s)
    endpoint = (3**E-1) >> (2+s)
    return audit_join(Join(E, E-1, left, right, endpoint))


def transport(join, proof):
    """Both ROOT directions are valid; payment is still source -> smaller child."""
    audit_join(join)
    audit_root(proof)
    if proof.exponent == join.source_exponent:
        old, new, E = join.left, join.right, join.child_exponent
    elif proof.exponent == join.child_exponent:
        old, new, E = join.right, join.left, join.source_exponent
    else:
        raise ValueError("supplied ROOT source is neither marked port")
    need(proof.word[:len(old)] == old, "same actual port prefix")
    return audit_root(Root(E, new+proof.word[len(old):]))


def from_legacy(receipt):
    need(type(receipt) is frozen.Receipt, "authenticated inherited receipt type")
    exponents = []
    for n in (receipt.source, receipt.child):
        natural(n, 1)
        need((n & (n+1)) == 0, "actual Mersenne port")
        exponents.append(n.bit_length())
    return audit_join(Join(*exponents, receipt.source_word,
                           receipt.child_word, receipt.endpoint))


def strict_profile(proof):
    audit_root(proof)
    return len(proof.word), sum(proof.word)-proof.exponent


def strict_affine_join_possible(source, child):
    """Complete criterion for the fixed bit-deletion affine ansatz, not all joins."""
    audit_root(source)
    audit_root(child)
    need(source.exponent > child.exponent, "ordered exponents")
    return (len(source.word) == len(child.word)
            and sum(source.word)-sum(child.word) == source.exponent-child.exponent)


def two_bit_heads():
    seed, u, v, period = frozen.load_exit()
    need(seed == 18269 and v[0] == 5, "retained even-child reset")
    partner = (2, 3)+v[1:]
    P, Q, B = carrier(u)
    p, q, b = carrier(partner)
    need(p == 9*P and Q == 4*q and b-p == B-P+q,
         "complete drop-two sibling identity, including carry")
    return seed, period, u, partner


def two_bit_exit(E):
    """Same inherited native phase, composed odd-to-odd smaller obligation."""
    natural(E, 2)
    mersenne(E)
    seed, period, u, v = two_bit_heads()
    need((E-seed) % period == 0, "unchanged exact exponent phase")
    x = 2*3**(E-1)-1
    P, Q, B = carrier(u)
    need((P*x+B) % Q == 0, "same original source head")
    Z = (P*x+B)//Q
    c = v2(3*Z+1)
    return audit_join(Join(E, E-2, (1,)*(E-1)+u+(c,),
                           (1,)*(E-3)+v+(c+2,), (3*Z+1) >> c))


def inverse_mersenne_pair(word):
    """Complete all-exponent test for one fixed contracting inverse word.

    Return (parent exponent K, child exponent F), or None. The returned
    forward word is word itself: M_F --word--> M_K, so child < parent.
    No search or arbitrary exponent cutoff occurs in this API.
    """
    letters(word)
    need(bool(word), "nonempty inverse word")
    need(sum(word) < 2*len(word), "necessary cost bound for inverse contraction")
    P, Q, B = carrier(word)
    need(Q < P, "contracting inverse slope")
    R = Q+B-P
    if R == 0:
        return None
    F = (abs(R) & -abs(R)).bit_length()-1
    z = P*(1 << F)+R
    if F < 2 or z <= 0 or z & (z-1):
        return None
    K = z.bit_length()-1-sum(word)
    if K <= F:
        return None
    # This finite word forces small F,K relative to its own exact carrier;
    # replay independently rejects padding and authenticates all valuations.
    n = (1 << F)-1
    for a in word:
        need(n != 1, "candidate cannot leave ROOT")
        z = 3*n+1
        need((z & -z).bit_length()-1 == a, "candidate's actual valuation")
        n = z >> a
    need(n == (1 << K)-1, "candidate's original-source identity")
    return K, F


def compositions(total, length, prefix=()):
    if length == 1:
        yield prefix+(total,)
    else:
        for a in range(1, total-length+2):
            yield from compositions(total-a, length-1, prefix+(a,))


def inverse_words(max_length):
    natural(max_length, 1)
    for length in range(1, max_length+1):
        cost = length
        while (1 << cost) < 3**length:
            yield from compositions(cost, length)
            cost += 1


NAMES = (6125, 6126, 6129, 6130, 18267, 18268, 18269, 18270, 18273, 18274)


def named_root(E):
    natural(E, 1)
    need(E in NAMES, "declared source-specific recipe")
    return _named_root(E)


@lru_cache(maxsize=None)
def _named_root(E):
    if E in (6125, 18269):
        return audit_root(Root(E, frozen.root_word(E)))
    if E in (6129, 18273):
        u, v = ((1, 10), (1, 1, 2, 2, 3)) if E == 6129 else ((10,), (3, 1, 1, 3))
        join = from_legacy(frozen.parent_receipt(E, u, v, 1))
        return transport(join, named_root(E-4))
    if E == 18268:
        return transport(from_legacy(frozen.exit_receipt(18269)), named_root(18269))
    if E == 18267:
        return transport(even_join(18268), named_root(18268))
    return transport(even_join(E), named_root(E-1))


def summary(proof):
    endpoint, peak = replay(proof.exponent, proof.word, True)
    encoded = json.dumps(proof.word, separators=(",", ":")).encode()
    return {"exponent": proof.exponent, "rank": len(proof.word),
            "cost": sum(proof.word), "cost_minus_exponent": sum(proof.word)-proof.exponent,
            "peak_bits": peak, "endpoint": endpoint,
            "word_json_sha256": sha256(encoded).hexdigest()}


def reject(fn, *args):
    try:
        fn(*args)
    except (ValueError, TypeError):
        need(True, "malformed/false receipt rejected")
    else:
        need(False, "invalid receipt accepted")


def main():
    # Exact reset-pair controls, including the separate E=2 ROOT boundary.
    for E in range(4, 130, 2):
        j = even_join(E)
        P, Q, B = carrier(j.left)
        p, q, b = carrier(j.right)
        need(P == p and Q == 2*q and B == 2*b-p,
             "full half-child affine identity")
        need(len(j.left) == len(j.right) == E and sum(j.left)-sum(j.right) == 1,
             "retained rank/cost labels")
    audit_root(Root(1, ()))
    audit_root(Root(2, (1, 4)))
    reject(even_join, 2)
    reject(even_join, 3)
    reject(even_join, True)
    reject(replay, 1, (2,), True)

    # Complete fixed-word universes, with no exponent cutoff.
    ordinary = tuple(inverse_words(8))
    need(len(ordinary) == 953, "inherited complete contracting-word universe")
    need(all(inverse_mersenne_pair(w) is None for w in ordinary),
         "all-exponent ordinary-bank Mersenne closure obstruction")
    generator_count = 0
    for length in range(1, 9):
        for labels in product((1, 5, 17), repeat=length):
            word = inverse.forward_word(labels)
            record = inverse.encode(labels)
            P, Q, B = carrier(word)
            need((P, Q, B) == (3**record.B, 1 << record.A, record.C),
                 "independent reverse-word carrier")
            need(inverse_mersenne_pair(word) is None,
                 "all-exponent three-generator-bank obstruction")
            generator_count += 1
    need(generator_count == 9840, "complete length-eight generator universe")
    # Positive native inverse receipt leaves the Mersenne category.
    need(inverse.apply(inverse.encode((5,)), 31) == 27
         and inverse.replay(27, (1, 2)) == 31, "positive off-category dependency")
    for w in ordinary:
        P, Q, B = carrier(w)
        need(Q+B-P >= 0 and ((Q+B-P == 0) == all(a == 1 for a in w)),
             "genuine words have nonnegative R, with exact all-ones boundary")
        for K in range(2, 33):
            numerator = Q*((1 << K)-1)-B
            if numerator <= 0 or numerator % P:
                continue
            h = numerator//P
            need(h < (1 << K)-1 and h & (h+1) != 0,
                 "independent bounded direct substitution")

    # Strict profiles obstruct the prescribed affine ansatz, not ROOT itself.
    three = Root(2, (1, 4))
    seven = Root(3, (1, 1, 2, 3, 4))
    need(not strict_affine_join_possible(seven, three), "7/3 all-length strict hostile")
    imported = {E: frozen.frozen_record(E) for E in (6125, 18269, 46941, 98477)}
    need(len({r["rank"] for r in imported.values()}) == 4,
         "retained four seed profiles are distinct (metadata, not fresh replay)")
    proofs = {E: named_root(E) for E in NAMES}
    rows = [summary(proofs[E]) for E in NAMES]
    direct = transport(two_bit_exit(18269), proofs[18269])
    need(direct == proofs[18267], "direct head equals the two-stage source-owned splice")
    seed, period, _, _ = two_bit_heads()
    for t in (0, 1, 9, 10**30):
        E = seed+period*t
        need(frozen.exit_guard(E) and E % 8 == 5 and (E-1) % 8 == 4,
             "symbolic all-height phase and fixed even-reset valuation")
    for E in (6126, 6130, 18270, 18274):
        need(strict_profile(proofs[E]) == strict_profile(proofs[E-1]),
             "reset pair has identical marked profile")
    for E in (18267, 18268, 18270, 18273, 18274):
        need(strict_profile(proofs[E]) == strict_profile(proofs[18269]),
             "first frozen seed component remains source-marked")
    need(not strict_affine_join_possible(proofs[18269], proofs[6125]),
         "different grounded components cannot be silently merged")
    # The derived word is independently reconstructed backward, not discovered.
    frozen.reverse_verify(mersenne(18267), proofs[18267].word)
    j = even_join(6126)
    reject(transport, j, proofs[18269])
    reject(audit_join, replace(j, child_exponent=6123))
    reject(audit_join, replace(j, right=j.right[:-1]+(j.right[-1]+1,)))
    reject(audit_root, Root(6125, proofs[6125].word+(2,)))
    reject(named_root, 6127)
    reject(named_root, True)
    reject(named_root, 6125.0)
    reject(two_bit_exit, 18271)
    reject(inverse_mersenne_pair, (2,))
    reject(inverse_mersenne_pair, ())
    reject(inverse_mersenne_pair, (True,))
    reject(mersenne, 100001)
    print(json.dumps({
        "status": "PROVED scoped source-preserving transport and all-exponent fixed-word test",
        "checks": CHECKS,
        "imported_seed_records": {E: {k: imported[E][k] for k in ("rank", "cost", "peak_bits")}
                                  for E in imported},
        "seeds_reauthenticated_this_run": [6125, 18269],
        "derived_with_no_discovery": [6126, 6130, 18267, 18270, 18274],
        "source_owned_root_words": rows,
        "inverse_universes": {"ordinary_length_le_8": 953,
                              "generator_length_le_8": 9840,
                              "mersenne_to_smaller_mersenne_pairs": 0},
        "open": "Universal guard coverage and ROOT for arbitrary phase members; no orbit search used"
    }, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
