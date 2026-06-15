# A000568: how far number theory goes

I pushed the Burnside exponent as far into pure arithmetic as I could without
lying to myself about the runtime. The clean reformulation is real:

```text
gcd(a,b) = Σ_{d|a,b} φ(d)
```

turns the cross-term in the tournament count into a divisor-profile update.
That means a large part of A000568 is not really about tournaments anymore. It
is about odd partitions carrying a totient-weighted divisibility profile.

That is the good news.

The bad news is also useful: in Python, this does not speed things up. The odd
partitions of `n` already have very few distinct sizes, so the original pairwise
`gcd` work was not the bottleneck. At `n=100`, the average odd partition only has
about `5.7` distinct sizes. The pairwise cross-size scan is tiny. The real cost
is just the number of odd partitions themselves.

So the honest conclusion is:

- number theory completely explains the local interaction term;
- number theory alone does not yet compress the global state space.

That is a satisfying boundary. It tells me what to stop doing. Another pure-Python
rewrite of the same Burnside term is not the path. The arithmetic recurrence
belongs in the compiled CRT / GMP engines, where shaving local interaction work
actually matters.

The deeper open problem is now cleaner:

> can the odd partitions themselves be quotiented by a second arithmetic profile
> without losing the denominator `Π k^{m_k} m_k!`?

That is where a real speedup would live. The current divisor-profile only
compresses the exponent, not the whole weight.
