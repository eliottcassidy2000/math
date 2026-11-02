# Extremal source check — 2026-11-02

**CITED PREPRINTS. Full-text statements checked; complete proofs not independently audited. These are separate context, not dependencies of the five new frameworks.**

## Divoux–Wang

Alexander Divoux and Shouda Wang, [The Maximum Number of Hamiltonian Cycles in Tournaments, v1](https://arxiv.org/html/2610.23755v1), submitted 2026-10-16. Theorem 1.4 states

    max_{|T|=n} H(T) = (c+o(1)) n!/2^(n-1),
    c=cosh(1)/cos(1)=2.8559578925651135... .

Section 10 supplies the path argument. Section 1.2 leaves eventual exact cyclic-tournament maximization conjectural. The authors credit GPT-5.6 for original proof arguments; the paper's authors are Divoux–Wang. No journal publication was listed on the inspected record.

For minimum realization order mu(h), our elementary Stirling inversion gives, along attained h tending to infinity,

    mu(h) >= L/W(L/(2e)) - O(1),    L=ln h,

where W is the positive Lambert-W branch. In particular mu(h)>=(1-o(1))ln h/ln ln h. Prior polynomial-factor bounds already give that scale; this is not a new growth-order claim. No effective finite cutoff was extracted. An extremal upper endpoint does not fill interior spectrum values.

## Wormald

N. C. Wormald, [Tournaments with many Hamilton cycles, v1](https://arxiv.org/html/2610.16245v1), posted 2026-10-12, is an older manuscript conditionally accepted in 2006 whose requested revisions were not made. Theorem 5 and the R101 table give a path lower bound exceeding `2.85588 n!/2^(n-1)` infinitely often. Its sharp limiting constant was conjectural, unlike the newer cited asymptotic result. Do not present this as a newly proved 2026 theorem.

Neither paper was listed in the inspected CORE-PAPERS.md. This note records source context without promoting either to independently audited canon or proving full-spectrum completeness.
