"""Independent replay for the finite mod-2^43 signed-Collatz certificate.

No construction/search algorithm is invoked.  Run next to cap4_certificate.txt.
"""
from fractions import Fraction
from pathlib import Path

HOLES = {7, 21}
SEEDS = {1, 5}
SMALL = {
    3: "-", 9: "-+", 11: "+++", 13: "+", 15: "-", 17: "+",
    19: "+-+++--+--", 23: "-", 25: "+", 27: "++-", 29: "-+++--+-",
    31: "-", 33: "+", 35: "-", 37: "--+", 39: "-", 41: "+",
    43: "+++", 45: "+", 47: "-", 49: "+", 51: "-", 53: "---",
    55: "-", 57: "+", 59: "-", 61: "+", 63: "-",
}


def require(cond):
    if not cond:
        raise ValueError("cap-4 certificate failed")


def v2(x):
    return (x & -x).bit_length()-1


def verify():
    entries = []
    for line in Path(__file__).with_name("cap4_certificate.txt").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        r, k, word = line.split()
        entries.append((int(k), int(r, 16), word))

    disjoint = set()
    mass = Fraction(0)
    max_contraction_threshold = max_endpoint_threshold = Fraction(0)
    max_len = max_k = 0
    min_val = max_val = 1
    peak = Fraction(0)
    for k, r, word in sorted(entries):
        require(1 <= k <= 43 and 0 < r < 2**k and r % 2 == 1)
        require(1 <= len(word) <= 26 and set(word) <= set("01234567"))
        require(all((d, r % 2**d) not in disjoint for d in range(1, k+1)))
        disjoint.add((k, r))
        mass += Fraction(1, 2**k)

        # For n ≡ r (mod 2^k), the current value is (p*n+b)/2^A.
        p, b, A = 1, 0, 0
        for char in word:
            digit = int(char)
            s, a = (1 if digit >= 4 else -1), 1+digit % 4
            max_val = max(max_val, a)
            peak = max(peak, Fraction(3*(65*p+max(0,b)),65*2**A))
            p, b = 3*p, 3*b+s*2**A
            A += a
            require(A < k)
            # Exact valuation of p*n+b; congruence of the source forces it.
            require((p*r+b) % 2**(A+1) == 2**A)
            # Every endpoint above 12 avoids all holes/ladders.
            max_endpoint_threshold = max(max_endpoint_threshold,
                                         Fraction(12*2**A-b,p))
        require(p < 2**A)
        max_contraction_threshold = max(max_contraction_threshold,Fraction(b,2**A-p))
        max_len, max_k = max(max_len,len(word)), max(max_k,k)
    require(mass == Fraction(1,2))
    require((len(entries),max_len,max_k,max_val) == (1068,26,43,4))
    require(max_contraction_threshold == Fraction(144967391,5077565) < 29)
    require(max_endpoint_threshold == Fraction(102660143,1594323) < 65)

    # Low odd sources use separate exact words, checked without sampling.
    require(set(SMALL) == set(range(1,65,2)) - SEEDS - HOLES)
    for start, word in SMALL.items():
        x = start
        for char in word:
            s = 1 if char == "+" else -1
            q = 3*x+s
            a = v2(q)
            m = q >> a
            require(1 <= a <= 4 and m > 0 and m not in HOLES)
            require(all((m<<i)-s not in HOLES for i in range(1,a+1)))
            x = m
        require(x < start)
    print("PASS: 1068 prefix-free cylinders, odd mass 1/2; 26 steps, v2 at most 4")
    print("thresholds:", max_contraction_threshold, max_endpoint_threshold,
          "coarse expanded peak/n bound:", peak)


if __name__ == "__main__":
    verify()
