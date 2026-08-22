"""Exact pure four-way ANOVA hostile for the THM-3679 chart marginal."""

from collections import defaultdict


P = 13
W = (1, 2, 3, 4, 5, 6, 0, 0, 0)
SOURCE = 6
TARGET_A = 7
TARGET_B = 8
CHARTS = ((0, 1), (2, 3), (4, 5), (1, 0))
ZERO = (0, 0)

VECTORS = (
    (1, 8, 8, 1, 8, 1, 12, 8, 1),
    (1, 1, 7, 1, 1, 1, 12, 1, 1),
    (1, 1, 1, 1, 2, 1, 12, 1, 1),
    (8, 7, 8, 7, 8, 7, 12, 8, 7),
)
LOCAL_VALUES = ((7, 6), (7, 0), (12, 0), (1, 12))

gates = 0


def require(condition, message):
    global gates
    if not condition:
        raise RuntimeError(message)
    gates += 1


def dot_w(address):
    return sum(x * y for x, y in zip(address, W)) % P


def chart(address, index):
    k, ell = CHARTS[index]
    return (
        (address[TARGET_A] - address[k]) % P,
        (address[TARGET_B] - address[ell]) % P,
    )


def joint(address):
    return tuple(chart(address, index) for index in range(4))


def normalized_subset(mask):
    address = [
        sum(VECTORS[index][coordinate] for index in range(4) if mask >> index & 1)
        % P
        for coordinate in range(9)
    ]
    # The source speed is zero, and no chart reads the source coordinate.
    address[SOURCE] = 1
    return tuple(address)


for index, vector in enumerate(VECTORS):
    expected = tuple(
        LOCAL_VALUES[index] if other == index else ZERO for other in range(4)
    )
    require(dot_w(vector) == 0, ("vector outside relation kernel", index))
    require(joint(vector) == expected, ("vector is not chart-local", index))

measure = {}
for mask in range(1, 1 << 4):
    address = normalized_subset(mask)
    coefficient = 1 if mask.bit_count() % 2 else -1
    expected = tuple(
        LOCAL_VALUES[index] if mask >> index & 1 else ZERO for index in range(4)
    )
    require(dot_w(address) == 0, ("normalized address outside kernel", mask))
    require(all(coordinate != 0 for coordinate in address), ("nonunit", mask))
    require(joint(address) == expected, ("wrong joint label", mask))
    require(address not in measure, ("duplicate address", mask))
    measure[address] = coefficient

require(len(measure) == 15, "wrong support size")
require(sum(measure.values()) == 1, "wrong total mass")
require(all(joint(address) != (ZERO,) * 4 for address in measure), "joint-zero atom")

for index in range(4):
    pushed = defaultdict(int)
    for address, coefficient in measure.items():
        pushed[chart(address, index)] += coefficient
    pushed = {label: mass for label, mass in pushed.items() if mass}
    require(pushed == {ZERO: 1}, ("wrong chart pushforward", index, pushed))
    print(f"chart_{index}_pushforward={pushed}")

print("support_size=15")
print("total_mass=1")
print("joint_zero_mass=0")
print(f"RESULT=PASS;gates={gates}")
