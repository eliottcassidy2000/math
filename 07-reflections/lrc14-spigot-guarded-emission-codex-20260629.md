# LRC14 Spigot Guarded Emission Reflection

The spigot prompt gave a better word for something the repo has been circling:
emission.

HYP-3523 reserves the concrete random031 certificate-spigot stream, and
HYP-3524 adds the hydrotope scout.  HYP-3525 is the guardrail layer around
both.  Recent LRC work already says "emit named debt" everywhere, but HYP-3525
makes that operational.  A proof token is a digit.  A quotient is the visible
head state.  The unresolved sidecars are the tail.  You are allowed to print
the digit only when the hidden tail cannot change it.

For random031 the pattern is concrete:

```text
private_firewall_status can emit through C/F/N/T or I/Q
h3490_route can emit only through R or route reconstruction
terminal_class can emit through flow_class/allowed_exit/sheet_pgf_bucket
owner_residual_pair can emit through transport+branch_boundary+residual+R
raw counts are checksums
```

This clarifies the relationship among HYP-3513, HYP-3520, HYP-3521, and
HYP-3522.  HYP-3513 is the route-carry guard.  HYP-3520 is the tail-error
canary for component compression.  HYP-3521 is the bounded preallocation of
terminal branches.  HYP-3522 is the owner predigit holdback that prevents
printing `(45,173)` before transport and bracket lift are attached.

The best proof target now looks like:

```text
GuardedEmission(target, visible_packet, hidden_tail):
  if target is constant over every legal hidden-tail completion,
  emit target;
  otherwise hold the first missing sidecar or emit named debt.
```

That sounds abstract, but it is exactly the random031 problem.  We do not need
more ways to admire the number `282`.  We need the theorem saying when `282`
may be printed as terminal class digits:

```text
230 ordinary
40 free-hole
12 bypass
```

and when the bypass digit may further emit:

```text
transport (23,93,113)
branch lift (147,169)
residual (45,173)
```

The spigot reframe also gives a good anti-pattern: a fake proof is a digit
stream without carry guards.  It looks productive, but one hidden carry can
rewrite the last emitted digit.  In repo language, that is raw `12`, raw `242`,
raw endpoint rank, or raw mirror closure pretending to be a route.

Two next experiments feel worth doing:

1. Write a small generic `GuardedEmission` checker over packet rows: given a
   target and a set of visible columns, report mixed fibers and the first
   missing sidecar.
2. Apply it to the random031 terminal stack first, then to older HYP-2963
   route/status banks.  This would turn the spigot metaphor into a reusable
   proof-interface validator.

The broader lesson is nice and ruthless: do not output the theorem digit until
the proof knows no future carry can change it.
