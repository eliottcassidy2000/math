#!/usr/bin/env python3
"""Exact F4 certificate for THM-4230's genus-two j=0 gluing condition.

F4 is represented as F2[t]/(t^2+t+1), with elements a+b*t encoded by
a+2*b.
"""


def add(first, second):
    return first ^ second


def multiply(first, second):
    a0, a1 = first & 1, (first >> 1) & 1
    b0, b1 = second & 1, (second >> 1) & 1
    # t^2=t+1 in characteristic two.
    constant = (a0 * b0) ^ (a1 * b1)
    linear = (a0 * b1) ^ (a1 * b0) ^ (a1 * b1)
    return constant | (linear << 1)


def square(value):
    return multiply(value, value)


def main():
    assert [square(value) for value in range(4)] == [0, 1, 3, 2]
    assert all(square(square(value)) == value for value in range(4))
    admissible = []
    for a in range(4):
        for b in range(4):
            if all(add(multiply(a, value), multiply(b, square(value))) == 0
                   for value in range(4)):
                admissible.append((a, b))
    assert admissible == [(0, 0)]
    print("F4_FROBENIUS", [square(value) for value in range(4)])
    print("DESCENT_RESIDUES", admissible)
    print("CONCLUSION a_mod_2=b_mod_2=0")
    print("PASS")


if __name__ == "__main__":
    main()
