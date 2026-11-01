# Three FRACTRAN programs and an exact square-shell fibre theorem

**Proved program semantics and all-order restricted-family identities;
finite exact replays are separately stated.**

| Program | Ordered fractions | Input | Readout |
|---|---:|---|---|
| Any tournament's H(T) | 4,581 | `131*5^N`, `N=2^v(2b+1)-1`, triangular adjacency bits b | At halt30689, exponent3=H(T). Invalid graph codes return0. |
| Full one-apex word family | 32 | `38*5^E`, `E=2^m+sum b_i2^i` | At halt107, exponent2=H and exponent3=all odd cycles; q-sidecar retained. |
| Square-shell sequence G(q) | 1,371 | `71*5^q`, q>=1 | At halt7937, exponent3=`q²+(q mod2)+2f(q)-1`. q0 returns0. |

The compact32-fraction list is printed in [THEORY.md](THEORY.md). Its input
is a marked apex adjoining a transitive chain, not an arbitrary formal
independence graph. All odd cycles share that apex, so it computes the entire
polynomial `1+cz`, with `H=1+2c`.

## Structural refinement

Let R_h count **pointed nontrivial strong one-apex isomorphism classes** with
H=h, over all orders. The inherited digital inverse gives an exact finite
fibre, with no assumed order cutoff. Converse pairs its nonfixed words.
The fixed words are `u + complement(reverse(u))`, where q is u's binary
value, and their H-values are

```
G(q)=3,9,11,25,31,45,51,81,91,117,127,165,179,...
q² < G(q) <= (q+1)².
```

Consequently **exactly one H-value in each shell `(q²,(q+1)²]` has R_h odd**,
and exactly N-1 such values lie at or below N². Upper-shell equality occurs
exactly at q=2,4,8,16,... . The word/pointing convention is essential:
R_5=2, but forgetting the apex identifies the two words100,110 as one
self-converse tournament. The inherited deletion-apex classification gives
an exact correction for the unpointed version; see the proof.

The square-shell conclusion is a converse-fixed-point theorem, not a claim
about arbitrary tournament H-fibres or the unrelated square-sum graphs.

## What it does and does not prove about 7 and21

The one-apex inverse excludes7 and21 at every word length, but it also
excludes13, which an unrestricted five-tournament attains. The extra gap is
a precise certificate that this grammar cannot prove unrestricted spectrum
completeness. The universal program counts every encoded tournament but
supplies no new universal non-reachability proof by itself. Existing
THM-1370 already excludes unrestricted7 and21 by strong floors and SCC
multiplicativity; that inheritance is explicit.

## Run

Python3 standard library only; no network or repository checkout required:

```sh
python3 check.py
python3 -O check.py
python3 fibres.py --bound 1000000
```

`build.py` deterministically regenerates both compiled lists; `one_apex.py`
regenerates its hand-written list. `VERIFIED.normal.json` and
`VERIFIED.optimized.json` are matching producer receipts. The independent
audit is in the sibling Session133 audit packet.

`FIBRES.json` records the exact million-bound restricted inverse:693,169
pointed strong classes,249,076 distinct H-values,999 odd fibres, maximum
fibre size14 first attained at176077. It is **not** a claim that one million
large FRACTRAN inputs were executed, or that the unrestricted spectrum is
complete. Small lists, sources, explicit witnesses and integer-only
reproduction commands accompany all claims.
