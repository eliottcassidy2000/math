# Frozen zero-heavy and physical-port audit

Independent audit of the session-18 packet whose manifest SHA-256 is
`a5490b17a42062db222676a44a1eb2ef1c3f0b55499b2389cc70691886f9112c`.
Read `AUDIT.md` for the exact proof review, scope, and limitations.

The checker uses Python's standard library only. No producer checker is
imported or executed, and no input, repository, or remote is modified.
The evolving later-corridor work is deliberately excluded.

Portable replay (paths may point anywhere):

```
python3 verify.py --source /path/to/frozen/session18-zero-heavy --repo /path/to/math
```

The repository need only contain the eleven dependency files listed by the
frozen source manifest, with their exact pinned contents. The checker does
not rely on the checkout's current HEAD. `verify.py` authenticates this
audit manifest, runs the independent checker under normal Python and `-O`,
and compares both outputs byte-for-byte with the included receipts.

For a single run:

```
python3 independent_audit.py --source /path/to/frozen/session18-zero-heavy --repo /path/to/math
```

The all-order claims are established by the reviewed mathematics, not by
extrapolating the finite receipts. No audit approval extends to unpinned
addenda or to unconditional six-colour forcing.
