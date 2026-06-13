# Triune Carrier Applications to LRC and Sibling Problems (S663)

The pi trinity became genuinely useful once it touched the LRC carry seam.
The translation is:

```text
sum      -> active wall and pair-sum packets
product  -> C=27 gcd shell and local obstruction product
fraction -> carry/owner boundary state
```

The exact seam is `v = r + 27k`.  Since `27 == -1 mod 14`, the carry `k`
decides the `n=14` clock residue.  Since `27` is odd, it also toggles parity.
So the carry vector is the LRC version of a continued-fraction boundary state:
not a scalar value, but the recursive data saying how a shadow is lifted.

The small S663 experiment makes that concrete.  AP, `Vstar`, and `2AP` are all
floor atoms with `M=1/14`.  They share the same product face:

```text
gcd shells = ((1,9), (3,3), (9,1))
mass = 27.
```

Then the script adds one or two `+27` carries.  All `273` perturbations are
strict, but if we group by only the `Res_27` additive/product shadow we get
three mixed groups:

```text
AP:    1 floor + 91 strict
Vstar: 1 floor + 91 strict
2AP:   1 floor + 91 strict
```

After adding the carry-continuant/fraction state, mixed groups drop to zero.
That is the whole point in one line: sum and product faces identify the correct
floor pocket, but they do not decide floor-vs-strict inside the lift.  The
fraction face does.

For the next LRC theorem, I would now write the target as:

```text
fixed odd-wall packets
+ fixed C=27 gcd/product shell
+ fixed carry-owner continuant
  => AP, Vstar, 2AP, or strict looseness.
```

This also clarified the sibling-problem transfer.

For OCF, the sum face is odd-cycle packet coefficients; the product face is
strong-component multiplication; the fraction face should be a deletion or
substitution boundary continuant for Hamiltonian paths.

For tournament decks, the sum face is deleted-card loss data; the product face
is component/scissors packets; the fraction face is exactly S660's paired
`(card, deleted score)` boundary state.

For unit distance, the sum face is spine plus bulk edge packets; the product
face is direction/norm support; the fraction face is point-deletion frontier or
ear ownership.

For Goldbach/Lemoine and pi/e trace-norm, the fraction face is reconstruction:
`q=O-E`, `p=2E-O`, or the branch sheet of `T^2-S*T+P`.

So the creative rule coming out of S663 is:

```text
When a scalar quotient almost works, ask what its fraction face is.
```

In this repo that usually means: who owns the boundary, which lift was chosen,
what deletion state was attached, and which recursive route produced the
visible packet.
