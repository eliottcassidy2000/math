# Independent source and inheritance pins

Checked directly on 2026-10-21 UTC. These records support bounded
attributions only; no full-paper or LRC computational-proof audit is made.

## External primary source

Jaan Allikvere, *Fourteen and fifteen lonely runners*,
[arXiv:2609.02604v2](https://arxiv.org/abs/2609.02604v2).
The version record gives initial submission 2 September 2026 and v2
24 September 2026.

The independently opened [Section 3.1](https://arxiv.org/html/2609.02604v2#S3.SS1)
uses the same projected cube and ellipsoid, with
`q_i=1+2u_i-u_i^2`, `u_i=v_i/max(v)`. Lemma 3.1 proves inclusion;
Lemma 3.2 gives the ellipsoidal covolume through compressed covariance.
Neither local lemma claims maximality. The introduction indexes `LRC(k)`
by nonzero speeds, corresponding to `k+1` runners. These local formulas
and normalization, not the paper's global results, are the comparison
inputs. The web inspection establishes version and mathematical location;
no raw external-page SHA-256 is claimed.

## Repository read-only pins

The three inherited files below were read and hashed. Their bytes at the
producer's declared snapshot `de0463bc8ddf` were independently recovered
by read-only `git show` and match the observed working-tree bytes.
The separately observed repository HEAD during this audit was
`f55a3c10df77805f01cd977ee6d882672cc28cce`; this is a read observation,
not an audit publication or a promise about later repository freshness.

- `01-canon/theorems/THM-4014-lrc14-diagonal-polar-ellipsoid-fastest-coordinate-relation-compression.md`
  SHA-256 `85774a46b9a8d509cc36591c65687998b0661df0176254663e79665f38d9e802`.
  The exact pair criterion and corrected compressed inverse agree with
  the base packet's inheritance description; the `(1,2)` hostile is present.
- `01-canon/theorems/THM-4009-euclidean-covering-transference-short-relation-compression.md`
  SHA-256 `f027a1893ebe81cc14f0deb6ddab9d48bc58010aa644bbedaec7dc9fb2a15aa9`.
  The projected zonotope and speed-dependent inradius context is present.
- `01-canon/theorems/THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction.md`
  SHA-256 `53e703085c31e888c5ad8daa749c8663ba6e0df15d6bcc8770e629965c17aad8`.
  The projected-lattice dual context and superseded numerical-cap notice
  are present.

These inheritance references are not unproved premises of either audited
ellipsoid theorem: both packets give the finite-dimensional arguments they
need. No repository file was edited, no fetch was needed, and no Git write
was performed by this audit.
