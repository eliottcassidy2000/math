# Span-three counting and density: what follows, and what does not

**Scope.** This audits the counting row and density discussion in `chanpak_waves_20261113/UNIFYING_FRAME.md`, and the calculation in `scripts/w39_frame_params.py`, as of repository commit `a2ec2ffce0d1a7ab3f57b976f1368570c35b7e0d`. It is not an obstruction to the span-three density conjecture itself.

## Correct word-counting implication

Let `A_1,...,A_b` be the permitted transfer matrices, and let `H(w)=h A_w v` with fixed initial column and readout row. Suppose at least `c b^n` permitted words of length `n` satisfy `H(w) <= C G^n`, where `b,G>1` and `c,C>0` are uniform in `n`. In particular, a uniform bound on *all* words works with `c=1` for a full alphabet. Count words with multiplicity, setting `W(X)=#{w : H(w)<=X}` (possibly infinite). For `n=floor(log(X/C)/log G)`, the length-`n` words alone give

```
W(X) >= c b^n >= (c/b) (X/C)^(log b / log G).
```

For a finite tuple let `a_n=max_{|w|=n} ||A_w||` under an induced matrix norm, and let `r=lim_n a_n^(1/n)` be the joint spectral radius. Submultiplicativity gives `r=inf_n a_n^(1/n)`. For every `G>r` choose an integer `k` such that `a_k<G^k` and write `n=qk+s`, `0<=s<k`. Then `a_n <= a_k^q a_s <= G^n max_{0<=j<k}(a_j/G^j)`. Thus a **certified upper** bound `r<=U` with `U>1` proves the lower logarithmic word-counting exponent `log b/log U`, taking `G` down to `U`; a bound with that exact power and a uniform constant requires the stronger estimate at `G=U`. Norm estimates concern the specified permitted language and alphabet: if restrictions reduce the number of available words, use the corresponding growth rate for their count.

Conversely, for each word, `rho(A_w)^(1/|w|) <= r`, by applying the definition to repetitions of that word. Wave 39 maximizes floating-point eigenvalue estimates over words of length at most six. Even granting a rigorously enclosed value for the displayed witness near `3.32260`, this is only a **lower** bound on `r`. Since `log b/log r` decreases in `r>1`, it does not yield the stated lower bound `log 8/log 3.32260 ≈ 1.732`. The sampled length-400 average near `2.2449` is also not a uniform upper bound or a proof of the positive-fraction hypothesis above. Neither observation establishes the reported counting exponent. A finite maximum of *operator norms*, rather than spectral radii, would give a rigorous route if exactly bounded and sufficiently small: `r<=a_k^(1/k)`.

## The counting object matters

A word count does not bound the number of distinct semigroup elements without a word-to-matrix multiplicity bound. A matrix count does not bound the number of distinct integers without a matrix-to-readout multiplicity bound. For a finite collection of `M` words each with readout at most `X`, if each integer has at most `B` preimages in that collection, one obtains only `M/B` distinct values. The sharper elementary second-moment bound is `M^2/(sum_t m_t^2)`, where `m_t` is the number of words with readout `t`. Neither is automatically a density-one statement. For example, take eight identical scalar matrices `[3]` and `h=v=[1]`: the uniform bound and JSR are exactly `3`, the number of words of size at most `X` is of order `X^(log 8/log 3)`, and every value `3^n` is odd, but there are only `O(log X)` distinct matrices or readouts. The proved 13-unit cone bound is a **lower** bound on `H` for bottleneck-free words in that language, so it does not supply the **upper** bound needed by the word-counting lemma or control either multiplicity.

## All local residues, even with fixed-modulus equidistribution, do not imply density

Put `s_n=n!+2n+1` for `n>=2` and `S={s_n:n>=2}`. These are strictly increasing odd integers. For each fixed modulus `q`, all `n>=q` have `n! ≡ 0 (mod q)`, so `s_n ≡ 2n+1 (mod q)`. If `q` is odd, the latter sequence cycles uniformly through all residues; if `q` is even, it cycles uniformly through all odd residues. Thus every residue class compatible with oddness occurs infinitely often, and each fixed modulus even has equidistribution among those classes as the number of sequence terms tends to infinity. The same limiting proportions hold when the set is truncated by an upper bound `X`, because `s_n` is strictly increasing. Nevertheless `n!>=2^(n-1)` for `n>=2`, so `#(S intersect [1,X])=O(log X)`, and its relative density among the odds is zero.

This is an implication counterexample, not a model of the ribbon semigroup. In particular, finite residue experiments, fixed-modulus equidistribution, a word-counting exponent, and an analogy to Zaremba do not replace the problem-specific uniform distribution, multiplicity, and analytic estimates needed for a local-to-global theorem. No density-one or cofinite conclusion follows from the current inputs alone.
