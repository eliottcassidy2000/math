# Rerooting an exact return supplies two coprime physical parents

**PROVED for the constructed recontact families.** This is a new type of source pair: equal height and weighted start, unequal deletion count. It does not assert such a pair exists at every target.

## Exact final-vertex rerooting

Take any rooted state (h,e,t,d). Append letter1 (the new vertex v beats only the old root a), but mark v as the new root. The resulting profile is exactly

    (H,E,T,D)=(h+t-e, h-e, h+t, h).

The ordinary count is the inherited one-back-edge formula. Paths ending at v correspond to old Hamiltonian paths not ending at a, giving h-e. Deletion gives the old tournament. For the weighted start, the singleton path v contributes h; every longer path starts v,a and corresponds to an old optional path-cover term starting at a, contributing t. Thus T=h+t.

This supplies a literal smaller parent whenever t>e, rather than an arithmetic parent assumed to be physical.

## Recontact gives a deletion-axis choice

Use two sources with common deletion d and an exact recontact produced by a width-(r+1) edit with remote coefficient A and common high word V>0 of bit length s. Both final letters are1 and their common upper tail makes their pre-final end counts equal, say eta. (For s=1, use the common bit at position r+1.) Let their common final height be N. Mark the final vertex of each graph.

The new profiles are

    L': (N, dL-eta, N+eta, dL),
    R': (N, dR-eta, N+eta, dR),
    dL-dR=2^s A.

Indeed the pre-final weighted-start difference tR-tL is exactly2^s A, from the binary transfer and reversed signed edit. Since the final counts agree and pre-final end reads agree, hL-hR=tR-tL. These pre-final heights are precisely the new deletion counts dL,dR.

In particular for unit response A=1, both dL,dR are odd and

    gcd(dL,dR)=gcd(dR,2^s)=1.

They are actual strictly smaller-count parents in the present families (their pre-final t exceeds eta, which is0 or5). Thus every such return supplies two coprime smaller physical cores while preserving both absolute H and weighted-start T. For general odd A, their common divisor must divide A.

## Concrete literal example and the continuation it preserves

The first original-source return has N=12,233 and shared high address V=34, so s=6. Independently counting the actual15-vertex graphs after rerooting gives

    (12233,817,12233,817),
    (12233,753,12233,753).

The two deletion cores have817 and753 paths; gcd(817,753)=1 and817-753=64. A common letter1 at these NEW roots gives counts23,649 and23,713, respectively. This is a physical controlled difference, not just two factorizations of an integer.

More generally, with Delta=dL-dR and a common subsequent binary address W, their ordinary outputs differ by

    H_L'(W)-H_R'(W)=-Delta [K(W)+(W mod2)].

The H and T terms cancel; both E and D differ by Delta. A preliminary common0 resets E and changes the response to -Delta K(W). Thus information formerly carried by the weighted-start difference has been transferred to the deletion coordinate without changing H or T.

## What the Collatz-style interpretation does and does not buy

These are genuine descending physical parents, not merely a scalar shrinking map. Their coprimality offers a possible new interface for combining repairs; it does not by itself provide a CRT construction at the desired height. The unresolved universal step remains: choose an analogous pair from an arbitrary N, instead of first generating N through a recontact.

There is also a precise reset warning. A common word with only one1 after at least one0 has K(W)=W mod2=0, so its counts agree. But rerooting at THAT final vertex makes the whole new profile identical: (N+2^qT,N,N+2^qT,N). This discards the deletion-axis choice. A repeatable construction must preserve or replenish that data rather than assume an ordinary count-preserving reset preserves it.

Replay: `python3 finite_edit_recontact/coprime_parent_replay.py`. It independently recounts the15-vertex roots and both deletion cores, tests literal next-step counts, and checks the all-order transfer identities on the first ten constructed family levels using exact integers.
