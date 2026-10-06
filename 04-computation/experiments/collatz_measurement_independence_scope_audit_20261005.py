"""Small exact regression audit of the incoming measurement-independence note.

Abstract bank completions below preserve bank atoms and normalization only;
they are not alternative canonical Collatz measures. No broad census.
"""
from fractions import Fraction as F
from math import factorial
import json
from collatz_backward_measurement_compiler_20261005 import deadline_floor

CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(label)


def odd_step(n):
    value = 3*n+1
    a = (value & -value).bit_length()-1
    return value >> a, a


def root_word(n):
    word = []
    while n != 1:
        n, a = odd_step(n)
        word.append(a)
        if len(word) > 512:
            raise ValueError("declared finite audit cap")
    return tuple(word)


def exact_weight(word):
    if not word:
        return F(1)
    l, k = len(word)-1, sum((a-1)//2 for a in word)
    return F(2*factorial(k)*factorial(l+1), factorial(l+k+2))


def shortcut_prefix(source, length):
    n, p, q, carry = source, 1, 1, 0
    for _ in range(length):
        if n % 2:
            n = (3*n+1)//2
            p, carry = 3*p, 3*carry+q
        else:
            n //= 2
        q *= 2
    check(q*n == p*source+carry, "exact shortcut affine identity")
    return n, p, q, carry


def h(m, j):
    gap = abs(m-j)
    t = 1 << gap
    return F(4*t, (1+t)**2)


def q_kernel(m, degree, j):
    x = h(m, j)
    return (9*x-8)*x**degree


def readout(law, m, degree):
    return sum((mass*q_kernel(m, degree, j) for j, mass in law.items()), F(0))


def main():
    report = {}
    ray = []
    for j in range(1, 13):
        n = (4**(j+1)-1)//3
        word = root_word(n)
        actual = exact_weight(word)
        check(word == (2*j+2,), "all ray sources have odd ROOT time one")
        check(actual == F(2, (j+1)*(j+2)), "exact vanishing weight on the ray")
        bound = deadline_floor(n, 1)["floor"]
        check(0 < bound <= actual, "source-sensitive reverse bound survives")
        if j >= 2:
            check(actual < F(1, 3), "time-only claimed floor fails")
        ray.append({"j":j, "source":n, "tau":1, "weight":str(actual),
                    "source_sensitive_lower":str(bound)})
    report["one_step_ray"] = ray

    endpoint, p, denominator, carry = shortcut_prefix(1, 2)
    check((endpoint,p,denominator,carry) == (1,3,4,1),
          "ROOT leaves coefficient cone but does not strictly descend")
    check(denominator > p and (denominator-p)*1 == carry,
          "exact equality boundary requires the carry")
    forty_one = root_word(41)
    endpoint,p,denominator,carry = shortcut_prefix(41, 2)
    check((endpoint,p,denominator,carry) == (31,3,4,1),
          "41 leaves first coefficient cone at shortcut depth two")
    check(len(forty_one) == 40 and len(forty_one) > 2,
          "cone-exit depth is not a ROOT deadline")

    current, p, denominator, carry, word = 165, 1, 1, 0, []
    for _ in range(17):
        current, a = odd_step(current)
        word.append(a)
        p, carry, denominator = 3*p, 3*carry+denominator, denominator*(1 << a)
    check(current == 167 and p < denominator, "contracting actual prefix may grow")
    check((denominator-p)*165 < carry, "carry explains this exact growth")
    report["clock_and_carry"] = {
        "root_boundary":{"source":1,"shortcut_depth":2,"endpoint":1,"slope":"3/4"},
        "first_exit_not_ROOT_time":{"source":41,"shortcut_exit":2,"odd_ROOT_time":40},
        "later_contracting_growth":{"source":165,"endpoint":167,"word":word,
                                   "P":p,"Q":denominator,"carry":carry},
        "scope":"The 165 witness is a later prefix, not a counterexample to first-exit equality above ROOT."}

    bank = {0:F(1,2)}
    positive = {0:F(1,2),1:F(1,2)}
    zero = {0:F(1,2),2:F(1,2)}
    for law in (positive, zero):
        check(sum(law.values()) == 1 and law[0] == bank[0],
              "bank atoms and normalization agree")
    for d in range(9):
        separate_lower = sum((mass*q_kernel(1,d,j) for j,mass in bank.items()), F(0))-4
        check(separate_lower == -4, "the separate interval lower bound is negative")
        check(readout(positive,1,d) == F(1,2), "actual readout can be positive")
        check(readout(zero,1,d) == 0, "same bank admits a zero-target completion")
    minimizing = {0:F(1,2),4:F(1,2)}
    check(readout(minimizing,1,1) == -F(640,729), "coupled residual witness")
    check(-F(640,729) > -4, "separate moment extrema are not jointly sharp")
    check(h(0,1) == F(8,9) and h(0,3) == F(32,81),
          "a gapped bank's tail maximum need not occur after its largest index")
    report["finite_bank"] = {
        "bank":{"0":"1/2"},"target":1,"separate_lower":"-4",
        "positive_completion_readout":"1/2","zero_completion_readout":"0",
        "degree_one_coupled_witness":"-640/729",
        "scope":"These completions preserve bank data and normalization, not every canonical Collatz identity."}

    check(root_word(5) == (4,) and root_word(21) == (6,),
          "two distinct sources reach ROOT")
    check(odd_step(21)[0] == 1 and odd_step(5)[0] == 1,
          "21 reaches ROOT while avoiding the singleton bank {5}")
    report["bank_forward_closure"] = {
        "bank":[5], "member":5, "its_successor":1,
        "successor_is_in_reach_bank":False,
        "rooted_source_avoiding_bank":21,
        "repair":"Close a certified bank under its known suffixes before using forward closure."}
    report["checks"] = CHECKS
    report["scope"] = ("Quantitative deadline converse repaired with source height; "
                       "coefficient/actual/ROOT clocks separated; bank-only information scoped.")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
