# One-hard actual-H3 link and bad-triangle closure

**Status: THM-2903 PROVED + FINITE-EXACT + VERIFIED.**

The actual `H_3` geometry is sparse: the `61` non-direct one-hard roots
have `461` core vertices and `1,961` pair flags, with core sizes `2..23`.
Literal singleton links close `1,907` pairs.  The `54` unresolved pairs
form ten bad graphs, nine triangle-free.  The last graph has exactly two
triangles, both of whose literal two-label children close by exact global
pair caps.  Therefore all `61` roots close.

The full sorted root list is:

```text
(1,2,3,4,5,6,10)
(1,2,3,4,5,8,13)
(1,2,3,4,5,9,11)
(1,2,3,4,5,10,12)
(1,2,3,4,5,11,12)
(1,2,3,4,5,12,13)
(1,2,3,4,6,7,13)
(1,2,3,4,6,8,9)
(1,2,3,4,6,8,13)
(1,2,3,4,7,9,14)
(1,2,3,4,8,9,13)
(1,2,3,4,9,11,13)
(1,2,3,5,6,12,13)
(1,2,3,5,7,9,14)
(1,2,3,5,7,10,14)
(1,2,3,5,7,13,14)
(1,2,3,6,11,12,13)
(1,2,4,5,6,8,9)
(1,2,4,5,6,9,13)
(1,2,4,5,6,10,12)
(1,2,4,5,7,8,9)
(1,2,4,5,8,9,11)
(1,2,4,5,8,11,12)
(1,2,4,6,7,8,11)
(1,2,4,6,8,9,11)
(1,2,4,6,8,11,13)
(1,2,4,7,8,12,13)
(1,2,4,8,11,12,13)
(1,2,5,6,8,10,11)
(1,2,5,6,10,12,13)
(1,3,4,5,6,11,12)
(1,3,4,6,7,12,14)
(1,3,4,6,8,12,13)
(1,3,5,6,7,10,12)
(1,3,5,6,7,12,14)
(1,3,5,6,9,10,11)
(1,3,5,7,9,10,14)
(1,3,5,9,10,11,12)
(1,3,6,7,11,12,14)
(1,3,6,9,11,12,13)
(1,3,7,11,12,13,14)
(2,3,4,5,8,9,11)
(2,3,4,5,8,10,11)
(2,3,4,6,8,9,11)
(2,3,5,6,9,11,12)
(2,4,5,6,8,11,12)
(2,4,5,7,8,10,11)
(2,4,5,8,10,11,13)
(2,4,5,8,11,12,13)
(2,4,6,7,8,12,13)
(2,4,6,8,9,11,12)
(2,4,6,8,11,12,13)
(3,4,5,6,7,12,14)
(3,4,5,6,8,11,13)
(3,4,5,6,9,11,12)
(3,5,6,8,9,12,13)
(4,5,6,8,10,11,12)
(5,6,7,8,11,12,14)
(5,6,7,9,10,11,12)
(5,6,8,9,10,11,12)
(6,7,8,9,11,12,14)
```

Canonical digests:

```text
61-root list
f58c4143f329d215ff9bb7ec594d172e831fbccbab713a2398dc6cb53c60b8b7

51 roots with every pair child individually closed
e2a7f23111c3dc8b41a6518fa83571c3a334b1f601e8216ee5eb0a701539c92a

two-triangle ledger
635acf7fbc9a23257c9a746fc5f5ed1eb4db4538baff7ffe42632cfbd678ece3

complete proof
08dc4a539544a417ff884ef4631b10d6eb14fd63d7b62b06e76d92e3d4d9b162
```

The set overlaps THM-2902 in exactly its two named roots.  It adds `59`
new roots to the current proved union, making the union `76` and the
official residual `3356`.
