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

## Pull Integration: Relation Code Versus Generated Compatibility

The later KPS HYP-2724 pull gives the right name for the relation-support
handoff: view `Lambda(E)` as a relation code and measure low-support packets by
`dmin,A2,A3,A4`, with MDS/arc sets as easy general-position cases and AP-like
sets as anti-MDS hard cases.  A concurrent HYP-2725 note then separates this
relation-support axis from the factorial-support axis used by the q0 atom.

I tested whether this relation-code lens can replace the generated-word
compatibility filter.  It cannot.  The bridge scout
`lrc14_miss_zeta_relation_code_bridge_codex_s71.py` reuses the exact `318`
generated-context tests and attaches relation spectra to the `72` unique
sparse-tail challenger shapes.  Globally the relation spectra correlate with
q0/U4/B2 barriers, but this is partly a shape-size effect.  Within fixed size
the signs are mixed, and the size-3 death-chain frontier has the worst q0/U4
witnesses while the simple `|coef|<=2` relation spectrum is flat.  KPS's
pattern-linearity probe points the same way from another angle: low-support
pattern counts explain only `R^2=0.572814` of exact `corr(E)`, so the tail and
the generated-word compatibility layer cannot be discarded.  The later
HYP-2726 Delsarte/MacWilliams unification is compatible with this: it gives a
row-level LP umbrella for moment/Krawtchouk/relation-code certificates, but it
does not identify the decorrelated product-word cone where `B2` is dirty and
relation-proxy order flips generated-risk order.

So the proof order should be:

1. generated miss-zeta singleton/death-chain compatibility;
2. coherent context merging;
3. relation-code/Delsarte `A3/dmin/MDS` finite-packet classification;
4. q0/Vitali atom boundary evaluation.

The useful theorem is not `A3 implies compatibility`.  It is that
compatibility reduces the remaining packets to a low-support relation ledger
whose finite leading cases are organized by `A3,dmin,MDS`.
