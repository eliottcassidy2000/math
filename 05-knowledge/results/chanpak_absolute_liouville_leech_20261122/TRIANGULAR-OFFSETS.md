# Exact perfect-number offsets from signed factors and all eight Leech labels

**PROVED.** For every integer \(p\ge13\) and
\[
           j\in\{-3,-15,-21,-45,-55,-91,-105\},
\]
\(E_p+j\) is an ordinary Hamilton-path count, where
\(E_p=2^{p-1}(2^p-1)\). In particular the prescribed exponent-31 target
\(E_{31}-3=2,305,843,008,139,952,125\) is hit exactly. The stronger
\(E_p-3\) formula works for \(p\ge8\). These are scalar witnesses, not
a single-arc correction or an assertion that the resulting ordinal hosts
are strong. Primality of \(2^p-1\) is not needed. A separate [compressed-port construction](MODULE-OFFSETS.md) gives a strong \(E_p-15\) host already for \(p\ge6\).

**PROVED.** There is also a rule for all eight pure-\(A\) Leech-side
Coxeter labels \(a\in\{2,3,4,5,7,9,13,25\}\). For \(A=a2^k,\ k\ge10\),
both
\[
          2A^2-A-3,\qquad 7(2A^2-A-3)
\]
are ordinary Hamilton-path counts. The new ingredient is a signed
\(-3\) reader at every one of the eight labels, not just the inherited
paired \(+1\) readers.

## 1. The factor pair is selected by the target

Set \(t=2^{p-1}\). Completing the square gives
\(8(E_p+j)=(4t-1)^2-(1-8j)\). The discriminant \(1-8j\) is an odd square
precisely for the integral offsets \(j=-r(r+1)/2,\ r\ge0\).
For odd \(j\) this requires \(r\equiv1,2\pmod4\). Matching the
integer roots of the factorization gives
\[
 E_p-\frac{r(r+1)}2=
 \begin{cases}
 (t-(r+1)/2)(2t+r),&r\equiv1\pmod4,\\
 (t+r/2)(2t-(r+1)),&r\equiv2\pmod4.
 \end{cases}
\]
Thus an exact target first selects \(r=(\sqrt{1-8j}-1)/2\). There is
no fitting of a residue or use of a near-perfect center. The replay
supplies both factors for \(r=2,5,6,9,10,13,14\), with each starting
power of two chosen by the displayed factorization. The already known
\(r=1\) is the lower perfect-number neighbor
[constructed earlier](../chanpak_perfect_wieferich_seams_20261122/PROOF.md).

For example at \(r=2\), the two exact factors are \(t+1\) and \(2t-3\).
At \(p=31\) these are \(1,073,741,825\) and \(2,147,483,645\).
A transitive spine with its extreme arc reversed gives the first;
the \((G,\beta)=(128,-3)\) row below, read at width 24, gives the second.
Ordinal joining forces a Hamilton path to traverse the two hosts in order,
which multiplies the two counts. This gives no arc between \(E_{31}-3\)
and \(E_{31}-1\).

## 2. Literal signed profiles

Use the [one-return scaffold in the Leech-label note](LEECH-LABELS.md).
For deletion rows \(D\), contacts \(I\), and a transitive module of
width \(w\ge1\), its read is \(G2^w+\beta\) with
\(\beta=H(D)-D_{II}\). Each core listed in the replay reaches \(I\)
from every vertex, so each *factor* host is strong.
The new dyadic signed factors available are:

| \(\beta\) | \(-15\) | \(-11\) | \(-7\) | \(-5\) | \(-3\) | \(3\) | \(5\) | \(7\) | \(9\) | \(13\) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| \(G\) | 2048 | 128 | 512 | 2048 | 128 | 32 | 32 | 32 | 16 | 128 |

For each desired factor \(2^d+\beta\), choose \(w=d-\log_2G\);
the uniform cutoff \(p\ge13\) ensures all required widths are positive.
For \(r=2\) alone \(p\ge8\) suffices (\(d=p\) for the \(-3\) factor).
The positive \(+1\) factor is inherited.

For the quadratic Leech rails, the \(-3\) readers are as follows;
the labels 2 and 4 are dyadic shifts of the odd part 1.

| Odd part of \(a\) | 1 | 3 | 5 | 7 | 9 | 13 | 25 |
|---|---:|---:|---:|---:|---:|---:|---:|
| \(G\) | 128 | 768 | 5120 | 448 | 9 | 416 | 200 |

The target \(2A-3\) selects \(w\) by \(2^w=2A/G\), an integer power
of two at least two for \(k\ge10\). Now \(A+1\) and \(7(A+1)\) are
the eight paired +1 labels of [the earlier six](../chanpak_six_signed_height_pairs_20261118/README.md)
and [the completion at 9 and 25](LEECH-LABELS.md).
Their two ordinal products with the newly constructed \(2A-3\)
give the claimed absolute values. The indexing is inherited from the
pure-\(A\) Niemeier/Leech catalogue; the proof requires the physical
profiles and does not make a lattice identification.

The companion `python3 check_triangular.py` records the complete rows
and contacts, computes each claimed profile by subset DP, and independently
checks an assembled strong factor host at width two. Exact integer tests
select the widths and multiply the factors, including exponent 31.
