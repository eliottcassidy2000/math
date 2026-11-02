# Independent bounded audit: prime source-padding cycles

**Verdict: PASS; release-ready.** This is a proof review and replay of the
compact supplied verifier only. No new tournament/state probes, large
searches, Git mutations, or changes to frozen source files were performed.

## Pinned inputs and replay

- `PRIME_PADDING_ORBITS.md`: SHA256 `8b2f4a83e2a6d96d18e52d3d5a0d750b05af797a51f1c56c8efda5228d226887`.
- `verify_padding_orbits.py`: SHA256 `3ec2bfa83d6aa1a7bd9cb52d1aafac58cfadd19898106365669f93c4252f9d6e`.
- Producer receipts and independent replay receipts, normal and `-O`:
  SHA256 `46910b91c7e3785522b3ba4e745383987d5c3c79a3797b6eade7126e25ef1be5`.

Both executions completed and all four outputs agreed byte-for-byte.
The frozen main-proof hash and all three historical Git-blob pins in
`PADDING_SOURCE_PINS.json` were checked. The cited census and merged-node
source were read at the specified baseline commit.

## Proof audit

Multiplication by x is invertible on functions vanishing at0, and exchanges
the two parity spaces. Thus every odd power fixes only their intersection,
zero. For even d the permitted support consists of gcd(d,p-1) roots of
unity, paired under x -> -x. Each parity has gcd(d,p-1)/2 free coordinates;
subtracting the duplicated zero yields exactly2p^(gcd(d,p-1)/2)-1 fixed
functions. The actual transitive seed x forces operator order p-1.

Möbius inversion correctly gives least-period state counts, not cycle
counts; division by the period gives cycle counts. Independently, a nonzero
opposite-pair support S has2(p-1)^|S| functions, and the simultaneous
conditions at x and -x give least period lcm_{x in S}(2,ord_p(x)). This
confirms the support formula and excludes every odd period greater than1.
The all-prime assertion concerns the formal parity-function space. Actual
complete tables are restricted to the previously audited primes3,5,7,11.
At other primes the actual subset is a union of cycles, since both P and
its inverse P^(p-2) are positive source-padding constructions.

The reduced-coordinate padding formula is well-defined because
m_c divides m_(c-1). Multiplying by d_c gives exactly the stated weighted
bidiagonal matrix, using d_c/d_(c-1)=c(c-1); no modular division occurs.
An arbitrary formal tuple lifts to an integer falling-factorial polynomial
with zero constant coefficient. Its integer monomial expansion expresses
the weighted padded profiles as a linear combination of transitive clocks,
so the sharp clock identity extends to the full formal coordinate product.
The K1 state preserves sharpness. This does not imply formal-state
realizability by tournaments. The empty-coordinate M=1 case is harmless.

## Historical cardinality firewall

The pinned Session138 literal census records8 full profiles at order5;
the pinned THM-4857 records10 merged nodes. Ordered versus unordered path
covers differ by coordinatewise multiplication by the fixed nonzero
integers c!, an injective map on integer profiles, so their distinct-profile
counts agree. Path reversal also preserves each
profile, so F descends to the merged nodes. Therefore no injective map on
those10 nodes can factor through F alone. The ten padding cycles are a
different, all-orders contextual-state object; their numerical coincidence
does not establish a node correspondence.

## Compact replay coverage

The receipt verifies24 formal primes through97,2,068 root-of-unity counts,
and agreement with support/order counting at10 primes through31. It checks
131 matrix moduli,3,009 formal basis instances,202,222 weighted matrix
updates, and131 sharp K1 clocks. Parity-loss and doubled-zero controls are
retained. The actual p5 table has1/8/40 states of periods1/2/4; the p11 table
has1/20/322080 states of periods1/2/10. No scope correction was required.
