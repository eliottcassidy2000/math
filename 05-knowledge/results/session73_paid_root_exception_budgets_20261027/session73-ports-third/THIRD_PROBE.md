# Final bounded probe from the single saved Paley11 near miss

**FINITE-EXACT negative result in the declared search domain only.** Neither
Paley11 two-complement nor three-complement forest reachability is decided.
This mini-packet does not change any frozen earlier artifact.

Input is ONLY the saved two-complement state with FAS13 in
`input_two_step_receipt.json`, copied byte-for-byte from the depth supplement.
Seed73073311 chooses32 optimal rank ties; for each, exact Hamilton-path DP
supplies4 maximum-backward-overlap path ties. Each resulting complement is
checked by exact FAS subset DP. If its FAS were at most10, at most32 final rank
ties with bounded slack would be checked for an actual forest; none qualified.

The fixed budget exhausted128 path draws, reaching93 distinct labelled target
tournaments. Minimum FAS was12; no forest candidate was found. This says
nothing about third complements from other two-step states, other paths from
this state, or larger depth. In particular it supplies NO three-step upper
bound. No further search was launched.

The best additional literal path is

    R=(4,7,5,10,1,3,0,8,2,6,9).

Its maximum backward overlap in the chosen source rank is5. The exact named
sequence from Paley11 has FAS20→15→13→12. All three Hamilton paths are legal
at their actual stages, with whole-complement clock3 and a separate active-
path remarking sidecar; no elementary relocation cost is claimed.

Reproduce:

    python3 search_third.py > receipt.json
    python3 -O search_third.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json
    python3 verify_third.py
    sha256sum -c SHA256SUMS

The search and independent literal/prepend-DP verification use exact integers
and explicit exceptions. They do not rerandomize any native source masses.
