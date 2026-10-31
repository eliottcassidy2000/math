# The same triple, two different future operations

**Status: elementary proved identities, independently reviewed. Classical D(1) structure; no novelty claim.**

Suppose a,b are positive integers and ab+1=r^2, with r positive. The regular extension c=a+b+2r satisfies

    ac+1=(a+r)^2,    bc+1=(b+r)^2.

Thus (a,b,c) is a D(1) triple. Starting from (a,b)=(1,3), r=2 gives c=8. What comes next depends on which pair is retained.

## Hold the graph port1 fixed

The operation (1,b) -> (1,c) has b=k^2-1, r=k, and c=(k+1)^2-1. It gives

    (1,3,8), (1,8,15), (1,15,24), (1,24,35), ... .

The moving labels3,8,15,24,35,... are exactly all positive neighbors of vertex1 in the infinite square-sum graph. The next available closing port after8 is15, because1+15=16. This is a quadratic sequence, not the even-Fibonacci sequence.

## Slide the pair instead

The operation (a,b) -> (b,c), starting at (1,3), preserves r=b-a. Indeed the new square root is b+r=2b-a=c-b. Hence

    c=3b-a,

and the moving sequence is1,3,8,21,55,..., the even-index Fibonacci weights used by Chung–Graham numeration. This is the norm-one canonical shift orbit discussed in the numeration packet.

The first divergence is therefore **15 versus21**. Although21 is the next sliding term, it is not a neighbor of port1:1+21=22 is not square. The invariant quadratic form a^2-3ab+b^2 equals1 along the sliding orbit; the anchored next pair(1,8) has value41 instead and is not the canonical pair(1,P(1)).

These are two exact state-update rules with a common initial triple, not evidence that the future is determined by the three numbers without specifying the operation. Neither rule constructs a Hamiltonian cycle, and product-plus-one adjacency is not ordinary square-sum adjacency except in the anchored relations involving1. No classification of all D(1) pairs or triples is asserted.

The associated test verifies both universal identities on generated branches and checks the complete port1 neighbor set through order10,000. The proof is the algebra above, not extrapolation from the finite controls.
