"""THM-4950 (5): the strong base-3 spectrum below 28,981.

Strong H_3 values of orders <= 9 come from S86's census (alpha_n.txt: code H strong SC a1 a2 a3 a4 ...).  Order-10
values come from witnesses_order10.txt (one witness per value, re-checked by verify_witness10.py).  Since
H <= H_3 for every tournament and every strong tournament of order n >= 11 has H >= f(11) = 225 (Busch; THM-1370 (3):
the strong floor is nondecreasing), the strong values below 225 are decided exactly by orders <= 10.
usage: python strong_spectrum10.py CENSUS_DIR
"""
import sys
import os

here = os.path.dirname(os.path.abspath(__file__))


def main():
    d = sys.argv[1]
    strong = set()
    for n in range(1, 10):
        for line in open(os.path.join(d, f"alpha_{n}.txt")):
            f = line.split()
            if f[2] == "1":
                a = [int(x) for x in f[4:8]]
                strong.add(1 + 3 * a[0] + 9 * a[1] + 27 * a[2] + 81 * a[3])
    v10 = set(int(line.split()[0]) for line in open(os.path.join(here, "witnesses_order10.txt")))
    print(f"strong H_3 values: orders <= 9: {len(strong)}, order 10: {len(v10)}")
    strong |= v10
    exact = [v for v in range(1, 225) if v % 3 == 1 and v not in strong]
    print(f"EXACT (Busch floor f(11) = 225): values v = 1 (mod 3), v < 225, with no strong realiser: {exact}")
    miss = [v for v in range(225, 32000) if v % 3 == 1 and v not in strong]
    print(f"first value >= 225 with no strong realiser of order <= 10: {miss[0]}; "
          f"every v = 1 (mod 3) in [225, {miss[0] - 3}] has one")


if __name__ == "__main__":
    main()
