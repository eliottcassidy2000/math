"""Exact positive digit-current certificate; no optimizer or external package.

Universe: six stated shortcut edges, powers/Mersenne edges through 256 bits,
all six padding templates at M=5..80, and every ordinary edge n=1..4095
lifted above its own digit height at three offsets. Checks survive python -O.
"""
from fractions import Fraction


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def step(n):
    return (3*n+1)//2 if n & 1 else n//2


def bits(n):
    return {s for s in range(n.bit_length()) if (n >> s) & 1}


def current(n):
    out = {}
    for s in bits(step(n)) | bits(n):
        v = int(s in bits(step(n))) - int(s in bits(n))
        if v:
            out[s] = v
    return out


def rank(rows):
    mat = [[Fraction(x) for x in row] for row in rows]
    r = 0
    for c in range(len(mat[0])):
        pivot = next((i for i in range(r, len(mat)) if mat[i][c]), None)
        if pivot is None:
            continue
        mat[r], mat[pivot] = mat[pivot], mat[r]
        d = mat[r][c]
        mat[r] = [x/d for x in mat[r]]
        for i in range(len(mat)):
            if i != r:
                d = mat[i][c]
                mat[i] = [x-d*y for x, y in zip(mat[i], mat[r])]
        r += 1
    return r


def main():
    certificate = {3: 2, 4: 1, 5: 1, 6: 2, 8: 1, 9: 1}
    totals = [0]*4
    rows = []
    for n, weight in certificate.items():
        row = [current(n).get(s, 0) for s in range(4)]
        check(max(bits(n) | bits(step(n))) < 4, "certificate support")
        rows.append(row)
        totals = [a+weight*b for a, b in zip(totals, row)]
        print(f"edge {n}->{step(n)} weight={weight} digit_current={row}")
    check(totals == [0]*4, "positive current cancellation")
    check(rank(rows) == 4, "certificate has full column rank")
    print("positive weighted current=0; rank=4; all four digit weights forced zero")

    for k in range(2, 257):
        check(current(1 << k) == {k-1: 1, k: -1}, "power current")
        check(current((1 << k)-1) == {k: 1, k-1: -1}, "Mersenne current")
    print("power/Mersenne opposite currents: 255 exact pairs through 256 bits")

    # For an eventual rank, all sufficiently high weights are the same c.
    # These lifts then have exactly the original low current and zero sum
    # in their separated high current, without assuming any particular c.
    count = 0
    for n in range(1, 4096):
        for gap in (2, 7, 19):
            m = max(n.bit_length(), step(n).bit_length())+gap
            lifted = ((3 if n & 1 else 1) << m) + n
            expected = ((9 << (m-1)) if n & 1 else (1 << (m-1))) + step(n)
            check(step(lifted) == expected, "integer lift")
            c = current(lifted)
            check({s: v for s, v in c.items() if s < m-1} == current(n), "low current")
            check(sum(v for s, v in c.items() if s >= m-1) == 0, "high current mass")
            count += 1
    print(f"eventual-to-global lift: {count} exact ordinary edges/paddings")

    for m in range(5, 81):
        for n in certificate:
            lifted = ((3 if n & 1 else 1) << m)+n
            check(lifted > (1 << m), "large source")
            check(sum(v for s, v in current(lifted).items() if s >= m-1) == 0, "template mass")
    print("six template lifts: 456 exact cases M=5..80")

    # Hostile control: this cancellation is not a circulation of actual
    # integer vertices, nor one source's legal closed orbit.
    divergence = {}
    for n, weight in certificate.items():
        divergence[n] = divergence.get(n, 0)-weight
        y = step(n)
        divergence[y] = divergence.get(y, 0)+weight
    divergence = {n: v for n, v in divergence.items() if v}
    check(divergence == {5: 1, 6: -2, 9: -1, 2: 1, 14: 1}, "integer graph boundary")
    print(f"nonzero actual-vertex boundary={dict(sorted(divergence.items()))}")
    print("boundary multiset: 6+6+9 -> 2+5+14; digit inventory agrees")

    # Positive control: n itself is digit additive but increases at 9->14.
    check(sum((1 << s)*v for s, v in current(9).items()) == 5, "height hostile")
    # Nonlinear interactions are not constrained by the linear certificate.
    nonlinear_sum = sum(w*(step(n)**2-n**2) for n, w in certificate.items())
    check(nonlinear_sum == 72, "nonlinear control")
    print(f"weighted current of n^2={nonlinear_sum}, not zero")
    print("PASS: exact certificate, rank, lifts, scale templates and hostile controls")


if __name__ == "__main__":
    main()
