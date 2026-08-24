# Imaginary-quadratic class-rank challenge: `ell=3` attack

> **SEARCH REPORT / NO CHALLENGE CERTIFICATE, 2026-08-24.**  The requested
> pairs `(3,r)`, `9 <= r <= 16`, remain **OPEN**.  No line is written to
> `classrank_ell3_candidates_20260824.txt`, because no rank-nine subgroup was
> found.  Full invariant factors returned by PARI `quadclassunit` are used as
> **GRH-assisted discovery data**.  A displayed subgroup is promoted only
> after the repository's exact binary-form verifier accepts it.

## Outcome

The attack reached, but did not pass, the known rank-eight frontier twice:

1. Elkies's rank-sixteen clean Mordell curve gives the already audited
   discriminant
   `-217541503961543485618350976479` and an exact rank-eight form subgroup.
2. The rank-seventeen Mordell record gives a different fundamental
   discriminant
   `-908800736629952526116772283648363`.  PARI discovers 3-rank eight there,
   and eight extracted forms pass the repository's exact order and
   independence verifier.  No ninth direction appears in PARI's
   GRH-assisted invariant-factor basis.

Thus the rank-seventeen record does not solve `(3,9)`.  The clean sufficient
target remains a square-free `D = 1 (mod 4)`, `3 not dividing D`, for which
`E_(16D): y^2=x^3+16D` has eighteen independent rational points.  Elkies's
Proposition 2 would then force the real and imaginary mirror 3-ranks to sum
to at least seventeen, and Scholz reflection would force the imaginary
mirror to have rank at least nine.  The current rational Mordell-curve record
is seventeen, not eighteen.  The record coefficient does not lie in this
clean `16D` family and its two mirror ranks remain eight and seven in
GRH-assisted computation.

## Exact rank-eight hostile at the rank-seventeen record

Elkies's record curve has

```text
k = -908800736629952526116772283648363
  = -2195745961 * 413891567044514092637683.
```

The following line is below the challenge threshold, but is
**FINITE-EXACT** as a rank-eight lower bound:

```text
3|8|908800736629952526116772283648363|8525249643490871,-7893554946833949,28477434296310971|10014450078965347,-3926333777860899,23072081499172303|12270135103487629,8646924841398227,20039918826223337|15133648406371307,-5992356835391825,15606102568678921|9984208154774713,1114533378924643,22787058001376981|3704086961184571,-2729890827134629,61840680993691331|10110405876787561,2702193319127479,22652468074183991|7116930434060473,1458695044572495,31998645213091339
```

Reproduction command:

```bash
python3 04-computation/imaginary_quadratic_class_rank_certificate_tool_20260824.py \
  search --pari-gp /opt/homebrew/bin/gp --ell 3 --target-rank 8 \
  --candidate 908800736629952526116772283648363
```

The search emitted

```text
cyc=[1759845629196,3,3,3,3,3,3,3]
ell_rank=8
EXTRACTED_CERTIFICATE_PURE_VERIFY=PASS
```

The leading invariant is itself divisible by three; counting only the seven
explicit `3` entries would incorrectly report rank seven.  Direct replay of
the emitted line, both normally and under `python3 -O`, returned

```text
PURE_PYTHON_VERIFY=PASS|ell=3|r=8|D=908800736629952526116772283648363|squarefree_core_factors=[(2195745961,1),(413891567044514092637683,1)]|mitm_halves=81,81
```

The pure verifier checks fundamentality, every form discriminant and
primitivity, exact order three, and the absence of a relation by intersecting
two 81-element half-spans.  Its lower-bound conclusion does not depend on
GRH.  The statement that the complete group has no ninth direction does.
The full frozen replay is in
[`classrank_ell3_rank17_record_rank8_control_20260824.out`](classrank_ell3_rank17_record_rank8_control_20260824.out).

## Public-data and construction audit

The following sources and routes were checked.

* The public Bagshaw--Jacobson--Scheidler--Rollick repository contains its
  algorithms, the small `p=5` discriminant file, and no release, branch,
  historical commit, issue attachment, or raw `p=3` form bucket containing
  the billions-scale computation described in the paper.
* GitHub code search for the two rank-eight discriminants, certificate
  grammar, `3-rank 9`, and high-rank quadratic-class-group phrases found no
  qualifying form line.  Several apparent `rank9` hits were coding-theory
  matrix ranks, not quadratic class groups.
* The public Mordell record ladder and Elkies's 2009 and 2016 announcements
  were checked.  The rank-fifteen minimum is not a clean fundamental-
  discriminant specialization.  The special rank-fourteen field in the 2009
  announcement has imaginary 3-rank seven.  The clean rank-sixteen and the
  later rank-seventeen records both stop at imaginary rank eight.
* Random fundamental-discriminant search was rejected as noncompetitive:
  the public audit already records that a biased search of over twenty
  billion `p=3` discriminants did not publish a rank-nine example, and the
  Cohen--Lenstra scale is roughly `3^(-81)` before constants.

The highest-value missing input remains the unpublished sorted `p=3` raw
bucket with every norm-equation form, not merely a discriminant list.  The
highest-value new computation is a Mestre--Nagao specialization search on a
high generic-rank `j=0` surface, scored first for the clean `16D` local
conditions and then subjected to exact 3-isogeny descent.  An arbitrary
rank-eighteen elliptic curve or a rank-eighteen geometric K3 surface does not
suffice.

## Submission interface

No `submit`, `submit_answer`, or problem-specific verifier tool is present in
the enabled tool inventory.  Epoch's public Open Problems page says verifier
access is external and the FAQ directs prospective submissions to
`math@epoch.ai`.  Consequently this task can prepare a certificate file but
cannot transmit one from the present environment.  Since no qualifying line
was found, there is in any event nothing honest to submit.

## Claim ledger

* **FINITE-EXACT:** the displayed rank-eight form subgroup at the
  rank-seventeen record discriminant.
* **FINITE-COMPUTED / GRH-ASSISTED:** the complete PARI invariant factors and
  the negative discovery statement that they contain no ninth 3-direction.
* **CITED / PROVED:** Elkies Proposition 2 and Scholz reflection give the
  clean rank-eighteen-to-imaginary-rank-nine implication under the displayed
  square-free and congruence hypotheses.
* **OPEN:** every requested pair `(3,r)`, `9 <= r <= 16`.
