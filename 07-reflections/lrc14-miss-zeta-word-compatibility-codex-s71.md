# Miss-Zeta Word Compatibility

HYP-2721 left a gap: the abstract atom cone is too permissive.  It lets high
missed-count atoms fake low-moment invisibility cheaply.  The S71 compatibility
scout says the next quotient is not another scalar; it is the generated
miss-zeta word from HYP-2698.

The important move is to keep fixed-`x` product structure:

```text
Z_context,x(A) = product_i Z_i,x(A).
```

Only after this product word is OR-convolved with the candidate cluster should
we take missed-count atoms `q_t` and the boundary value `q_0`.  If we take the
atom cone first, we admit fake directions that no coherent LRC context seems
to generate.

The sparse-tail frontier is a good stress test because raw residual
coordinates genuinely fail there.  Nonconsecutive shapes beat consecutive
blocks on named sparse coordinates.  But after generated contexts are applied,
all `318` tested atom moves have positive `q0`, positive `U4/q0`, and
nontrivial `W1+W2` leakage.  This is exactly the behavior HYP-2721 needed:
the cheap abstract LP directions are not silently present in the generated
word cone.

The proof route should be:

1. Prove the singleton-context exclusion using the death-chain kernel
   `g_r(t)`.
2. Prove coherent context merging cannot weaken the relevant `W1+W2` or q0
   margin below the singleton frontier.
3. Route low-support relation packets through HYP-2719/HYP-2714.
4. Use THM-558 only for the Bonferroni-visible part; it is not enough alone.

The lesson is the same as the Vitali atom lesson in tournament language:
low-observer data is allowed only after the compatibility/completion channel
has been retained long enough to pay for the target functional.

## Pull Integration: Depth Law Versus Generated Word

The mac-mini S10 pull reframes the full actual-row crux as depth-law
convex-order occupancy.  That is the right global row-level object.  It should
not be collapsed with the generated-context compression object.

After adding `B2` to the S71 scout, the distinction is concrete:

```text
decorrelated generated contexts:
  B2/q0 nonpositive = 59/318
  min B2/q0 = -1135/282
  U4/q0 nonpositive = 0/318
```

So full-row even `B2` cleanliness is not a product-word compression theorem.
For generated words, the usable witnesses are low-factorial leakage and
`B4=U4`, plus relation-support handoff.  For full rows, the depth-law
convex-order program can still pursue even Krawtchouk extremality.  These are
two layers of the proof stack, not competing explanations.
