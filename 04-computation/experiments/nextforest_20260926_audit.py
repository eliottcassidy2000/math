"""Independent bounded-carry-defect descent certificate audit.

Does not import producer code. Carries are computed from residue floors;
large-input certificates are checked by exact affine arithmetic; the finite
exception set is generated independently from each top-bit/offset pair.
"""
from fractions import Fraction
from hashlib import sha256
import json


def require(ok, why):
    if not ok:
        raise RuntimeError(why)


def odd_step(n):
    v = 3*n+1
    a = (v & -v).bit_length()-1
    return v >> a, a


def defect(n):
    c = sum((1 << s)*((3*(n % (1 << s))+1) // (1 << s))
            for s in range(1, n.bit_length()+3))
    return 8*n-c


def descent(n):
    y = n
    for j in range(1, 10001):
        y, _ = odd_step(y)
        if y < n:
            return j
    raise RuntimeError("finite exception did not descend within cap")


def main():
    eq, bound_count = 0, 0
    for n in range(1, 65536, 2):
        m = 3*n+1
        d = m-(1 << (m.bit_length()-1))
        e = defect(n)
        require(e >= 2*d, "defect/top remainder inequality")
        no00 = "00" not in bin(n)[2:]
        require((e == 2*d) == no00, "equality binary language")
        require(0 <= e <= 6*n-6, "sharp carry bounds")
        eq += (e == 2*d)
        bound_count += 1
    print(f"independent residue-floor inequality checks={bound_count}, equalities={eq}")

    rows, exceptions = [], []
    max_cert = 1
    max_exception = 0
    for d in range(0, 2049, 2):
        # An integer source needs 2^t = 1-d (mod3).
        rem = (1-d) % 3
        if rem == 0:
            continue
        parity = 0 if rem == 1 else 1
        if d == 0:
            rows.append([0, parity, "one-step-to-one"])
            continue
        a = (d & -d).bit_length()-1
        core = d >> a
        q, j, total = core, 0, 0
        while q != 1 and j < 10000:
            q, v = odd_step(q)
            total += v
            j += 1
        require(q == 1, "finite core certification")
        while 3**(j+1) >= 2**(a+total):
            j += 1
            total += 2  # U(1)=1, exact valuation2
        denominator = 2**(a+total)
        gain = denominator-3**(j+1)
        cutoff = max(a+total+1, d.bit_length())
        while (1 << cutoff)*gain <= (4-d)*denominator:
            cutoff += 1
        require(gain > 0 and cutoff > a+total, "strict lift gate")
        require((1 << cutoff)*gain > (4-d)*denominator, "endpoint descent gate")
        max_cert = max(max_cert, j+1)
        rows.append([d, parity, a, core, j, total, cutoff])

        # Check all small top exponents; the proof covers every larger one.
        for t in range(cutoff):
            if t % 2 != parity or (1 << t) <= d:
                continue
            n = ((1 << t)+d-1)//3
            if n <= 1 or not n & 1:
                continue
            require(3*n+1 == (1 << t)+d, "finite source integrality")
            steps = descent(n)
            max_exception = max(max_exception, steps)
            exceptions.append([d, t, n, steps])

        # Endpoint identity tested at three representative large exponents.
        t0 = cutoff+((parity-cutoff) % 2)
        for t in (t0, t0+2, t0+20):
            n = ((1 << t)+d-1)//3
            y = n
            for _ in range(j+1):
                y, _ = odd_step(y)
            expected = 3**j * 2**(t-a-total)+1
            require(y == expected and y < n, "large lift actual replay")
    print(f"offset certificates={len(rows)}, finite exceptions={len(exceptions)}")
    print(f"maximum large-input odd horizon={max_cert}, finite-exception first descent={max_exception}")
    require(len(rows) == 684, "683 positive offsets plus the zero-offset case")
    require(max(max_cert, max_exception) <= 104, "declared uniform horizon")
    data = json.dumps([rows, exceptions], separators=(",", ":"), sort_keys=True).encode()
    print(f"independent certificate transcript sha256={sha256(data).hexdigest()}")

    family_times = []
    for k in range(3, 61):
        n = (4**k+17)//3
        require(3*n == 4**k+17, "27 family integrality")
        require(defect(n) == 36, "27 family defect")
        s = descent(n)
        if k <= 10:
            family_times.append([k, n, s])
        if k >= 6:
            y = n
            for _ in range(6):
                y, _ = odd_step(y)
            require(y == 243*2**(2*k-10)+5 < n, "six-step family descent")
    print(f"27-family first descents (k,n,steps)={family_times}")
    # 27's post-first-step value41 initially shadows the smaller core9;
    # the fourth core step consumes more binary precision than is available.
    large, small = 41, 9
    for _ in range(3):
        large, al = odd_step(large)
        small, ass = odd_step(small)
        require(al == ass, "initial core valuation agreement")
    require((large, small) == (71, 17), "27 precision boundary")
    require(odd_step(large)[1] == 1 and odd_step(small)[1] == 2, "finite precision hostile")
    print("PASS: all-height certificate E<=4096 -> a smaller odd iterate within104 U steps")
    print("Scope: first descent only; the later orbit can leave the bounded-defect class")


if __name__ == "__main__":
    main()
