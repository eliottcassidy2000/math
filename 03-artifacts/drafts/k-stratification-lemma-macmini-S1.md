# The k-stratification lemma and the k-reduction for the 3/38 cell

**mac-mini-2026-07-06-S1 (HYP-4232).** Composes the S55b binding-pair derivation
with kps-S17's cluster-gcd ladder (HYP-4217).  Verification:
`04-computation/lrc_kstrat_verify_macmini_S1.py`.

## Setting

W a primitive 12-family with M(W) = 3/38 (a hypothetical attainer of the first
Farey cell of the second-value gap).  By grid attainment (THM-592/HYP-4108,
formal), the maximizer t* = m/q* (lowest terms) exists with margins
dist(v·t*) ≥ 3/38 for all v ∈ W, equality for the binders.

## Lemma (k-stratification; elementary, proved below)

Write d_q(x) for the distance of x to 0 in Z/q.  Then:

**(0) 38 | q*.**  The binder's grid distance q*·(3/38) = d_{q*}(v_b m) is an
integer; gcd(3, 38) = 1 forces 38 | q*.  Write q* = 38k.

**(i) Every binder is divisible by k.**  v_b·m ≡ ±3k (mod 38k) ⟹
v_b·m ≡ 0 (mod k); gcd(m, 38k) = 1 ⟹ gcd(m, k) = 1 ⟹ k | v_b.

**(ii) The scaling identity and the quotient structure.**  For k | v, v = kv′:
d_{38k}(k v′ m) = k·d_{38}(v′ m).  Hence the QUOTIENT family
W′ = {v/k : v ∈ W, k | v} satisfies: every quotient clears level 3 mod 38 at
dilation m (d₃₈(v′m) ≥ 3), the binders' quotients sit at EXACTLY 3, and the
binding quotient pair sums to exactly 38.  The k = 1 mod-38 structure recurses
onto the quotient subfamily.

**(iii) Non-multiples clear on the fine grid.**  k ∤ v ⟹ d_{38k}(vm) ≥ 3k with
NO scaling reduction (their residues are genuinely mod-38k data).

**(iv) The k-cluster.**  K = {v ∈ W : k | v} ⊇ the binding pair, |K| ≥ 2,
and k | gcd(K).

## The k-reduction (via the cluster-gcd ladder)

The attainer has M = 3/38 < 2/25, so it is a no-2/25-point family and
kps-S17's ladder applies: for every S with 1 ≤ |S| ≤ 6,
(25 − 4|S|)·gcd(W∖S) ≤ 50·Σ_{w∈S} |w|.

Take S = W∖K (the non-multiples of k):
- **|S| = 0** (all twelve divisible by k): k | gcd(W) = 1 (primitivity) ⟹ k = 1.
- **1 ≤ |S| ≤ 6**: gcd(W∖S) = gcd(K) ≥ k ⟹
      k ≤ 50·Σ_{w∈S}|w| / (25 − 4|S|).
  The k ≥ 2 strata are HEIGHT-BOUNDED relative to the non-multiples — no
  CRT-ray escape (this is exactly what the ladder was built to do; the S55
  periodicity theorem says nothing weaker could).
- **7 ≤ |S| ≤ 10** (the named residual): |K| ≤ 5, so the mod-38 witness
  structure of (ii) lives on ≤ 5 quotients — with the binding pair summing to
  38 among them, i.e. a ≤ 5-element level-3 mod-38 configuration containing a
  (v′, 38−v′) pair; the 38k-grid tight cover must then be carried mostly by
  the ≥ 7 non-multiples on the fine grid.  This stratum is left OPEN here,
  with its structure pinned (census-shaped; small quotient side).

## The k = 1 anchor (the finite base)

At k = 1 the binding pair sums to EXACTLY 38: shapes (v, 38−v), v = 1..18,
both ≤ 37; 9 both-odd pairs {(1,37),(3,35),(5,33),(7,31),(9,29),(11,27),
(13,25),(15,23),(17,21)} and 9 both-even {(2,36),...,(18,20)}.  Both-even
pairs halve to (u, 19−u) — a pair summing to 19, pushing the analysis onto
the mod-19 world (the E_p structure); both members even also feeds the
parity dividend.  Each shape is a concrete finite anchor: pair ∈ W, the
other ten runners constrained by the level-4 tight cover mod 38, the full
gap profile, and — decisively, per S55b — the REALIZATION-level exactness
that no residue system captures.

## Status

(0)–(iv) and the |S| ≤ 6 reduction: PROVED (elementary; verification script
confirms all identities and the lemma's predictions on synthetic families).
The |S| ≥ 7 stratum and the k=1 realization kill: OPEN, structure pinned.
Lean next (this session): the lemma core (0)+(i)+(ii) as pure integer
arithmetic — the scaling identity and the divisibility transport.
