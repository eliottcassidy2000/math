# Goldbach/Lemoine Same-Pair Bridge S643

HYP-2218/S642 pins down the exact projection identity.  This S643 note keeps
the same algebra but shifts the object to the induced even/odd companion graph.

The fun version is exactly right:

```text
even: E = p + q
odd:  O = p + 2q
```

Even representations are unordered pairs.  Odd representations are ordered
pairs, because one prime is singled out and doubled.

The delightful part is that an even target and one odd companion reconstruct the
ordered pair:

```text
q = O - E
p = 2E - O.
```

So this is a linear Vieta carrier.  `E` is the trace-like shadow; `O` is the
tilted ordered shadow.  Together they recover the hidden pair exactly.

For an unordered pair `{p,q}`, the even node is `E=p+q`.  The odd companions are

```text
E+p, E+q.
```

They are symmetric around `3E/2`, since `(E+p)+(E+q)=3E`.  Their difference is
`q-p`.  When `p=q`, the two odd branches fold together:

```text
p=q:  E=2p, O=3p.
```

That is the branch locus.  In the language of the pi/e thread, it is the
discriminant-zero face.  In the doubled-prime thread, it is the diagonal where
`2p=p+p` and `3p=p+2p` live on consecutive additive/multiplicative rungs.

The prime `2` is special in exactly the way the LRC apex is special.  A Lemoine
representation with doubled prime `2` has the form

```text
O = p + 4.
```

But the natural pair `{p,2}` sums to `p+2`, which is odd, not an even Goldbach
target.  So `q=2` reps are apex/boundary reps, not shared-pair bridge reps.
The S643 scan found that up to `121`, only `7=3+2*2` is apex-only; all larger
tested odd targets have at least one odd-prime bridge representation.

This gives a clean graph:

```text
even node E=p+q
  -> odd node E+p  (double p)
  -> odd node E+q  (double q)
```

Duplicates have only one outgoing odd node.  The branch line is:

```text
3 -> (6,9)
5 -> (10,15)
7 -> (14,21)
11 -> (22,33)
13 -> (26,39)
...
```

That line is a prime detector in disguise: an even `E` lies on the duplicate
branch exactly when `E/2` is prime, and its odd companion is `3E/2`.

The new application I like most is to treat Goldbach/Lemoine as a companion
graph, not two separate conjectural covers.  A Goldbach pair does not merely
represent one even number; it casts one or two ordered odd shadows.  Conversely,
a non-apex Lemoine representation points back to an even parent by subtracting
the doubled prime.

That suggests three finite diagnostics:

- duplicate branch price: how much representation mass lives on `p=q`;
- apex price: how much Lemoine mass uses `q=2`;
- orientation price: how often the two odd branches collide with branches from
  other even nodes.

This fits HYP-2215 perfectly.  The scalar parity slogan is weak.  The retained
side channels are pair identity, orientation, duplicate branch, and the `2`
boundary.
