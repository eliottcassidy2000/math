# Anti-Cosets Everywhere S632

The request was exactly the right shove: stop treating `Anti(T)` as a cute
self-converse tournament object and look for the transporter pattern itself.

The general shape is:

```text
Aut(X)  = Iso(X, X)
Anti(X) = Iso(X, JX).
```

Once one anti-map exists, all the other anti-maps are obtained by composing
with automorphisms.  That sentence is almost tautological, but it changes what
we should compute.  The data object is not a preferred flip; it is the whole
transporter from a carrier to its opposite.

S632 checks three finite faces.

First, self-converse tournaments remain the clean theorem face.  Through
`n=7`, the generated SC classes have no coset failures: `|Anti|=|Aut|`, both
left and right coset checks pass, and `H=7,21` are still absent.  The long
anti-cycle types `(6,)` and `(6,1)` are not weird exceptions; they are the
warning label saying "do not force this into a vertex involution."

Second, the LRC `2n-1` shell machinery is already using transporters.  The
action of `<2,-1>` on folded shells has exact maps from a shell to its doubled
shell, with transporter size equal to stabilizer size in every checked even row
`n=6..20`.  At `n=14`, `C=27`, the whole thing snaps into the known gcd strata
`1,3,9`.  The new phrasing adds a useful caution: folding `j ~ -j` turns
reflection anti-data into an ordinary stabilizer.  That quotient is convenient,
but it has eaten the orientation side channel.

Third, the Eisenstein/unit-distance prefixes give the geometric version.
Centered hex shells `n=1,7,19` keep the full `D6` anti-coset.  Partial shell
prefixes beyond `n=7` often keep only one reflection or none.  The prefix
`n=21` in the triangular shell order has `47` unit edges and just one anti-map,
which reinforces HYP-2206: the exact/Moser `n=21` story is not "keep placing
Eisenstein points."  Its carrier is the `57=20+37` spine/bulk split.

The practical rule for future scripts is simple:

```text
name X;
name J;
compute Aut(X);
compute Anti(X)=Iso(X,JX);
record what quotienting destroys.
```

This should replace raw "merge complements" wherever the downstream proof cares
about observer roots, endpoint labels, shell orientation, unit spines, or
coimage cancellation choices.
