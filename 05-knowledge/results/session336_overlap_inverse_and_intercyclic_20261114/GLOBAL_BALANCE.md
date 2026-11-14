# Whole-tournament balance, not only frame balance

Use the hub-free intercyclic classification. For a uniformly chosen directed Hamiltonian path write $u\prec v$ if $u$ appears before $v$. These are statistics over **all** paths, not just a marked frame. Let $z$ be the unique spine digit 2 in each normalized representative.

| Tournament | $H$ | $\#(x\prec y)$ | $\#(x\prec z)$ | $\#(y\prec z)$ |
|---|---:|---:|---:|---:|
| $1^a21^b$, $A=2^a,B=2^b$ | $1+3AB$ | $H-A-B-1$ | $2AB$ | $AB+1$ |
| $1^a20$, $A=2^a$ | $9A-3$ | $6A-3$ | $6A-2$ | $3A$ |

On $1^a21^b$, the pair $(x,z)$ is always balanced in $[1/3,2/3]$, with probability $2AB/(1+3AB)$. The pair of hubs, however, can be arbitrarily biased: $\Pr(x\prec y)\to1$ when both lengths increase. In the second family, $\Pr(x\prec z)=2/3$ exactly at **every** length. The converse just reverses each probability. For the supplied counts 385 and 573, the two witnesses have probabilities $256/385$ and $382/573=2/3$.

**Proof.** Keep $x,y,z$ distinguished and fix their order in each candidate Hamiltonian path. At each homogeneous block all its allowable assignments into any fixed set of occupied run slots are the surjections onto those slots. Order constraints involving $z$ restrict the allowed slot set, not that surjection count, since $z$ is outside the varied block. So each marked total still lies in the span of $1,2^r,3^r$ in each positive block length. Direct literal path permutation counts at lengths 1,2,3 per independent axis give precisely the table, and the basis is invertible. This verifies the marked identities at all lengths, not only the sample orders. The marked count at count-385 is independent of which of its six splits is used, even though many other pair probabilities vary.

This supplies true whole-path balance in the requested families without turning the family of paths into the linear extensions of a poset. The analogue is more useful as a diagnostic than as an inverse: two balanced pairs and an increasingly imbalanced hub pair can coexist with *zero* ternary overlap. Balance does not tell us which number of paths or which overlap a different tournament will realize.
