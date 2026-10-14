#!/usr/bin/env python3
"""Length-21 univalence from a tiny exact shifted-numerator certificate.

No full frontier is enumerated. The 55 pruned subtrees plus 11 terminal
candidates form an explicit prefix-free cover of the only possible high branch.
"""
from functools import lru_cache
from collections import Counter
import hashlib
import json

checks = 0

def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise ArithmeticError(message)

L = 21
P3 = [3 ** i for i in range(L + 1)]
budget = [q.bit_length() - 1 for q in P3]
q = P3[L]
A_FINAL = budget[L]
shift = 4 * q
minimum = q - 2 ** L
maximum = sum(P3[L - 1 - i] * 2 ** budget[i] for i in range(L))
high_threshold = minimum + shift
low_ceiling = maximum - shift
require(shift <= maximum - minimum < 2 * shift, "only possible nonzero shift is 4*3^21")


def affine(word, native=True):
    A = B = 0
    for i, a in enumerate(word, 1):
        require(a >= 1, "positive valuation")
        B, A = 3 * B + 2 ** A, A + a
        if native:
            require(A <= budget[i], "prefix slope bound")
    return A, B


def upper_at_cost(A):
    return sum(P3[L - 1 - i] * 2 ** min(budget[i], A - L + i) for i in range(L))

possible_costs = [A for A in range(L, A_FINAL + 1)
                  if upper_at_cost(A) - minimum >= shift]
require(possible_costs == [32, 33], "only costs 32 and 33 meet the range")
# A low word beginning (1,2) has its remaining valuations minimized by ones.
low_branch_minimum = P3[L - 2] * 5 + 2 ** 3 * (P3[L - 2] - 2 ** (L - 2))
require(low_branch_minimum > low_ceiling, "the high word, not low word, begins (1,2)")

# Tail bound: at every later index j, cumulative cost is at most budget[j].
upper_tail = [0] * (L + 1)
for j in range(L - 1, -1, -1):
    upper_tail[j] = upper_tail[j + 1] + P3[L - 1 - j] * 2 ** budget[j]


def completion_upper(i, A, B):
    return P3[L - i] * B + P3[L - 1 - i] * 2 ** A + upper_tail[i + 1]


@lru_cache(None)
def descendants(i, A):
    """Independent cost-only count of all length-20 continuations."""
    if i == L - 1:
        return 1
    return sum(descendants(i + 1, A + a)
               for a in range(1, budget[i + 1] - A + 1))


def reject_low(original):
    """The unique possible affine decoder; return its first exact obstruction."""
    B, used, word = original, 0, []
    for remaining in range(L, 1, -1):
        base = P3[remaining - 1]
        remainder = B - base
        if remainder <= 0:
            return {"reason": "nonpositive_remainder", "decoded_prefix": word,
                    "next_position": len(word) + 1, "remaining_length": remaining,
                    "residual_B": B, "required_base": base, "remainder": remainder}
        a = (remainder & -remainder).bit_length() - 1
        require(a >= 1, "candidate numerator parity during decoding")
        word.append(a)
        used += a
        B = remainder >> a
        if used > budget[len(word)]:
            return {"reason": "prefix_budget", "decoded_prefix": word,
                    "position": len(word), "cumulative_cost": used,
                    "allowed_cost": budget[len(word)], "residual_B": B}
    require(B == 1, "decoder terminal B")
    last = A_FINAL - used
    require(last >= 1, "decoder final cost available")
    return {"reason": "collision", "word": word + [last]}


visited = 0
pruned = []
candidates = []
terminal_prefixes = []


def visit(word, A, B):
    global visited
    visited += 1
    i = len(word)
    bound = completion_upper(i, A, B)
    if bound < high_threshold:
        pruned.append({"prefix": word, "cost": A, "B_prefix": B,
                       "completion_upper": bound, "covered_prefixes": descendants(i, A)})
        terminal_prefixes.append(tuple(word))
        return
    if i == L - 1:
        high = 3 * B + 2 ** A
        require(high == bound and high >= high_threshold, "terminal candidate range")
        low = high - shift
        full_word = word + [A_FINAL - A]
        require(affine(full_word) == (A_FINAL, high), "candidate exact affine data")
        require(affine(word + [32 - A]) == (32, high), "same candidate at cost32")
        rejection = reject_low(low)
        require(rejection["reason"] != "collision", "no length21 collision")
        carry = high * pow(2 ** A_FINAL, -1, q) % q
        require(carry == low * pow(2 ** A_FINAL, -1, q) % q, "shifted candidate equal carry")
        candidates.append({"high_word_at_cost33": full_word, "prefix_cost": A,
                           "B_high": high, "B_low_candidate": low,
                           "carry_at_cost33": carry, "low_rejection": rejection})
        terminal_prefixes.append(tuple(word))
        return
    for a in range(1, budget[i + 1] - A + 1):
        visit(word + [a], A + a, 3 * B + 2 ** A)


visit([1, 2], 3, 5)
require(visited == 168 and len(pruned) == 55 and len(candidates) == 11,
        "small certificate shape")
# Verify every bound and that every pruning bound is attained by a legal completion.
for record in pruned:
    word = record["prefix"]
    A, B = affine(word)
    require((A, B) == (record["cost"], record["B_prefix"]), "pruning prefix data")
    require(completion_upper(len(word), A, B) == record["completion_upper"], "pruning bound formula")
    maximum_word = list(word)
    for j in range(len(word) + 1, L):
        maximum_word.append(budget[j] - A)
        A = budget[j]
    maximum_word.append(A_FINAL - A)
    require(affine(maximum_word) == (A_FINAL, record["completion_upper"]), "pruning maximum attained")
    require(record["completion_upper"] < high_threshold, "all pruned branches impossible")

# Prefix-free cover plus exact descendant counts verifies exhaustion independently.
for i, word in enumerate(terminal_prefixes):
    for other in terminal_prefixes[i + 1:]:
        common = min(len(word), len(other))
        require(word[:common] != other[:common], "certificate terminal branches prefix-free")
covered = sum(record["covered_prefixes"] for record in pruned) + len(candidates)
require(covered == descendants(2, 3), "certificate covers every possible high prefix")
require(len({record["B_high"] for record in candidates}) == 11, "candidate numerators distinct")
failures = dict(sorted(Counter(record["low_rejection"]["reason"] for record in candidates).items()))
require(failures == {"nonpositive_remainder": 10, "prefix_budget": 1}, "exact rejection summary")

result = {
    "status": "EXACT length21 no-collision certificate; all-length hypothesis remains open",
    "length": L, "canonical_cost": A_FINAL, "modulus": q, "shift": shift,
    "B_minimum": minimum, "B_maximum": maximum,
    "high_threshold": high_threshold, "low_ceiling": low_ceiling,
    "possible_costs_before_normalization": possible_costs,
    "low_branch_12_minimum": low_branch_minimum,
    "visited_prefix_states": visited, "pruned_subtrees": len(pruned),
    "covered_length20_prefixes_starting_12": covered,
    "candidate_count": len(candidates), "rejection_counts": failures,
    "candidates": candidates, "pruning_certificate": pruned,
    "checks": checks,
    "limitations": ["This certifies length21, not every length.",
                    "Maximal-cost normalization is an all-length equivalence, not a proof of injectivity.",
                    "No claim about convergence of a fixed positive integer source."]}
payload = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
result["payload_sha256"] = hashlib.sha256(payload).hexdigest()
print(json.dumps(result, sort_keys=True, indent=2))
