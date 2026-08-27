"""Exact F4 recurrence audit for the THM-4264 visible-gated observer."""

from collections import Counter
from itertools import product


# F_4=F_2[w]/(w^2+w+1), encoded as a+2b for a+bw.
ZERO, ONE, W, W2 = 0, 1, 2, 3


def add(left, right):
    return left ^ right


def mul(left, right):
    a0, a1 = left & 1, (left >> 1) & 1
    b0, b1 = right & 1, (right >> 1) & 1
    return ((a0 * b0) ^ (a1 * b1)) | (
        ((a0 * b1) ^ (a1 * b0) ^ (a1 * b1)) << 1
    )


def power(value, exponent):
    answer = ONE
    for _ in range(exponent):
        answer = mul(answer, value)
    return answer


def cubic_word(initial):
    # P_0(S)=S^3+w^2*S^2+w*S+1 in characteristic two.
    answer = list(initial)
    for index in range(9):
        answer.append(add(
            add(mul(W2, answer[index + 2]), mul(W, answer[index + 1])),
            answer[index],
        ))
    return tuple(answer)


def shift(word, index):
    return word[index % 12]


def p0_zero(word):
    return all(add(
        add(shift(word, j + 3), mul(W2, shift(word, j + 2))),
        add(mul(W, shift(word, j + 1)), shift(word, j)),
    ) == ZERO for j in range(12))


def visible_core(word):
    # Omit the unit in P_v=(S+w)(S^2+w).
    quadratic = tuple(add(shift(word, j + 2), mul(W, shift(word, j)))
                      for j in range(12))
    return tuple(add(shift(quadratic, j + 1),
                     mul(W, shift(quadratic, j)))
                 for j in range(12))


def order_two_recurrence(word):
    # (S+w^2)^2=S^2+w.
    return all(add(shift(word, j + 2), mul(W, shift(word, j))) == ZERO
               for j in range(12))


def telescope(word):
    answer = ZERO
    for value in word:
        answer = add(answer, value)
    return answer == ZERO


def period(word):
    for candidate in (1, 2, 3, 4, 6, 12):
        if all(word[j] == word[(j + candidate) % 12] for j in range(12)):
            return candidate
    raise RuntimeError("period divisor of twelve not found")


def name(value):
    return ("0", "1", "w", "w2")[value]


def show(word):
    return ",".join(name(value) for value in word)


def main():
    ambient = [cubic_word(initial) for initial in product(range(4), repeat=3)]
    assert len(set(ambient)) == 64
    assert all(p0_zero(word) and telescope(word) for word in ambient)
    ambient_periods = Counter(period(word) for word in ambient)
    assert ambient_periods == Counter({1: 1, 3: 3, 6: 12, 12: 48})

    projected = [word for word in ambient
                 if not any(visible_core(word))]
    assert len(projected) == 16
    assert all(order_two_recurrence(word) and telescope(word)
               for word in projected)
    projected_periods = Counter(period(word) for word in projected)
    assert projected_periods == Counter({1: 1, 3: 3, 6: 12})
    assert all(not any(word) for word in projected
               if word[0] == word[1] == ZERO)

    one_edge_hostile = (
        ZERO, ONE, ZERO, W, ZERO, W2,
        ZERO, ONE, ZERO, W, ZERO, W2,
    )
    assert one_edge_hostile in projected and one_edge_hostile[0] == ZERO

    ambient_hostile = tuple(
        mul(power(W2, index),
            ONE if (index * (index - 1) // 2) % 2 else ZERO)
        for index in range(12)
    )
    expected_ambient = (
        ZERO, ZERO, W, ONE, ZERO, ZERO,
        ONE, W2, ZERO, ZERO, W2, W,
    )
    expected_visible = (
        W, ONE, W2, W, ONE, W2,
        W, ONE, W2, W, ONE, W2,
    )
    assert ambient_hostile == expected_ambient
    assert ambient_hostile in ambient and ambient_hostile not in projected
    assert ambient_hostile[0] == ambient_hostile[1] == ZERO
    assert ambient_hostile[2] != ZERO
    assert visible_core(ambient_hostile) == expected_visible

    print("VISIBLE-GATED TWO-EDGE F4 RECURRENCE AUDIT")
    print(f"ambient_words={len(ambient)} "
          f"period_profile={dict(sorted(ambient_periods.items()))}")
    print("ambient_two_prefix_hostile=" + show(ambient_hostile))
    print("ambient_hostile_visible_projector=" + show(expected_visible))
    print(f"projected_words={len(projected)} "
          f"period_profile={dict(sorted(projected_periods.items()))}")
    print("projected_order_two_recurrence=PASS")
    print("projected_first_two_zero_implies_zero=PASS")
    print("one_edge_hostile=" + show(one_edge_hostile))
    print("scope=cyclic_module_only; geometric_1512_incidence_values_not_evaluated")
    print("VISIBLE_TWO_EDGE_RECURRENCE_AUDIT: PASS")


if __name__ == "__main__":
    main()
