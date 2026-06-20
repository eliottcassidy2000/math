"""
INDEPENDENT verification (attempt to REFUTE): the "product kernel reality" claim.

Setup
-----
Let omega = exp(-2 pi i / 7). For sector index j in 1..6 the Fourier coefficient
of the sector-j indicator is

    shat(n, j) = exp(-2 pi i n j / 7) * (1 - exp(-2 pi i n / 7)) / (2 pi i n),   n != 0.

Claim under test (two-far product kernel reality):
For all j, j' in 1..6 and all integers m, n:

    P = shat(m, j) * shat(n, j') + shat(-m, j) * shat(-n, j')

is REAL, i.e. Im(P) = 0.

Strategy
--------
Write our OWN cmath implementation from the definition (no repo imports),
then scan max |Im(P)| over j, j' in 1..6 and 1 <= m, n <= 10. Also do a
broader scan (including negative m,n and larger range) to try to break it.
"""

import cmath
import math


def shat(n, j):
    """Fourier coefficient of the sector-j indicator. n != 0 required."""
    if n == 0:
        raise ValueError("n must be nonzero")
    return (cmath.exp(-2j * math.pi * n * j / 7.0)
            * (1.0 - cmath.exp(-2j * math.pi * n / 7.0))
            / (2j * math.pi * n))


def product_kernel(m, n, j, jp):
    """P = shat(m,j)*shat(n,j') + shat(-m,j)*shat(-n,j')."""
    return shat(m, j) * shat(n, jp) + shat(-m, j) * shat(-n, jp)


def scan(m_range, n_range, j_range, jp_range, label):
    max_imag = 0.0
    argmax = None
    count = 0
    for j in j_range:
        for jp in jp_range:
            for m in m_range:
                if m == 0:
                    continue
                for n in n_range:
                    if n == 0:
                        continue
                    P = product_kernel(m, n, j, jp)
                    im = abs(P.imag)
                    count += 1
                    if im > max_imag:
                        max_imag = im
                        argmax = (j, jp, m, n, P)
    print(f"[{label}] cases={count}  max|Im P|={max_imag:.3e}")
    if argmax is not None:
        j, jp, m, n, P = argmax
        print(f"    argmax: j={j} j'={jp} m={m} n={n}  P={P!r}")
    return max_imag


def main():
    print("=" * 70)
    print("Product kernel reality — independent cmath check")
    print("=" * 70)

    # Primary scan requested: j,j' in 1..6, 1<=m,n<=10
    primary = scan(range(1, 11), range(1, 11),
                   range(1, 7), range(1, 7),
                   "PRIMARY j,j'=1..6, m,n=1..10")

    # Adversarial scan 1: include negative and larger m,n
    adv1 = scan(range(-20, 21), range(-20, 21),
                range(1, 7), range(1, 7),
                "ADV neg+large m,n=-20..20")

    # Adversarial scan 2: include j,j' = 0 and 7 (out of stated range, probe)
    adv2 = scan(range(1, 11), range(1, 11),
                range(0, 8), range(0, 8),
                "ADV j,j'=0..7 (probe outside range)")

    # Adversarial scan 3: cross check shat(-n,j) vs conj(shat(n,j))
    # Symmetry reason: shat(-n, j) should equal conj(shat(n, j)) because the
    # indicator is a real-valued function => Fourier coeffs satisfy c_{-k}=conj(c_k).
    print("-" * 70)
    print("Conjugation symmetry check: shat(-n,j) vs conj(shat(n,j))")
    sym_err = 0.0
    for j in range(1, 7):
        for n in range(1, 21):
            err = abs(shat(-n, j) - (shat(n, j)).conjugate())
            sym_err = max(sym_err, err)
    print(f"    max|shat(-n,j) - conj(shat(n,j))| = {sym_err:.3e}")

    print("-" * 70)
    overall = max(primary, adv1, adv2)
    print(f"OVERALL max|Im P| across all scans = {overall:.3e}")
    print(f"Conjugation symmetry max error      = {sym_err:.3e}")


if __name__ == "__main__":
    main()
