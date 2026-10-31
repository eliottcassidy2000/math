# Exact norm–resolvent bridge

`THEORY.md` proves that for primitive rational γ=a/b,

    λγ=α/(1−γα) is integral  iff  F_k(b,a)=±1.

The same condition characterizes integral closed-loop/resolvent matrices and
integral U_k(γ), which is exactly λγ's constant power-basis coefficient.
At k=3 the cited complete Thue classification gives precisely
`γ=0,1/2,56/103`, with generators `β,2β⁴,−103β¹⁷`.

The packet separates same field from same order, integral from unimodular,
internal contraction from global stability, and exact feedback solvability
from dynamical conjugacy. Order indices are `b^(k(k−1))`; image indices are
`b^k`. Localized orders agree after inverting b.

Run `python3 -B verify_release.py` for integrity and relocated normal/optimized
replay. The checker is standalone; no repository or third-party packages are
needed. Its2,436 finite checks corroborate, not replace, the proofs.
