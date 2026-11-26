# Carry windows and a neutral repair that changes a later arc response

2026-11-26. This session revisits three older mechanisms, not just three names: bounded carry in the [June certificate spigot](../../../07-reflections/lrc14-random031-certificate-spigot-codex-20260629.md), the distinction between residue and the full address in the [Mahler carry frontier](../../../00-navigation/MAHLER-THREE-HALVES-FRONTIER-2026-08-23.md), and physical state retention in the [Zeckendorf/shared-owner repair](../zeckendorf_target_selection_20261113/selection/PROOF.md). The diagnostic targets come from the [endpoint-changing grammar](../endpoint_context_recursion_20261113/README.md) and the [prime-lane packet](../prime_terminal_20261113_1610/README.md). **Chan–Pak remains open.**

## 1. An exact carry-window inverse for absolute targets

The existing square-root inverse asks for a particular **unit crossing**. This calculation asks the weaker but different question whether a prescribed **odd path count itself** occurs in the nested seam bank, not assuming a unit response on the binary arm. For a fixed left endpoint orientation and fixed inner profile, let the count at odd address $c=v+1$ be $A F(v)+B(v+1)-D$ as in the [nested seam](../overnight_chanpak_moore_20261126/NESTED-FORK-INVERSE.md), and let $K$ be the desired absolute count plus $D$.

Split $v=\ell q+r$ at $\ell=2^k$, with $k\ge1$ and even $0\le r<\ell$. Set $Z$ to the $k$-bit reverse-complement of $r$, $\epsilon=Z\bmod2$, $w=B\ell+AZ$ and $C=B(r+1)+A(F(r)-1)$. In the inverted-pair formula for $F$, comparisons crossing this cut contribute exactly $qZ+\epsilon[q\text{ odd}]$: the first term sums the weights of all distant cross pairs; the second remembers the extra adjacent-pair term at the join. Thus the requirement is the exact equality

    K = w q + C + A(F(q)+ε[q odd]).

The prefix $q$ lies between $\max(0,\lceil(K-C-2A)/(w+A)\rceil)$ and $\lfloor(K-C-A)/w\rfloor$ by $1\le F(q)\le q+1$. It must also satisfy $w q\equiv K-C\pmod A$. The useful constancy is that $w\equiv B\ell$ modulo $A$ for **every** suffix. Since $A$ is odd, $g=\gcd(A,B\ell)=\gcd(A,B)$. The prefix progresses with stride $A/g$, with a **single** modular inverse of $B\ell/g$ reused across the entire suffix scan. If $g$ does not divide $K$, there are no solutions. Otherwise the algorithm tests the original equality, so a surviving phase alone is never reported.

The window width before the congruence is $O(1+AK/(B^2\ell^2))$. Across all suffixes, including the cost of visiting them, this gives $O(\ell+gK/(B^2\ell))$ arithmetic visits, hence $O(1+\sqrt{gK}/B)$ by choosing a dyadic $\ell$ on the scale $\sqrt{gK}/B$ (and at least two). These are arithmetic-visit bounds, not bit-complexity claims. Unlike the unit inverse this is exponential in a worst-case target's bit length. If an ordered incoming prefix is allowed, $B$ becomes $2^jB$ and the offset changes once to $D_0$; no $j$ with $2^jB>K$ can contribute. This yields an exact finite inverse for that whole prefix family.

### Exact replays on old targets

Use the seven marked outer cores from [the phase study](../ordowski_cfh_response_frontiers_20261126/README.md), all uniquely terminal inner cores through order 12, both outer orientations and arbitrary ordered incoming prefix length. The exact replay considers 123,907 distinct profile/context entries and checks 13,824 known inverse controls. It recovered these historical targets with a marked physical arc:

| Prescribed count | Two physical orientations, before / after outer edit | Literal order check |
|--:|:--|--:|
| 289 | 291, 289 | 10, both strong |
| 889 | 889, 887 | 12, both strong |
| 9,011 | 9,011, 9,009 | 14, both strong, using no redundant high zero |
| 47,184,047 | 47,184,049, 47,184,047 | 28, exact nested count |

For the last row the inner core is a singleton with two ordered incoming vertices, and $c=156477$ in the $(23,24)$ outer frame. The arithmetic formula gives the count and the marked bridge establishes strength at that address. These values were already known attainable; what is new here is the exact inverse and the marked constructions in this architecture. The inverse found none of $13,27,1489,53084231,34379738367$ **in this bank only**. A separate exact *double-unit* scan of all nine targets in the same bank found none even as a corner of those particular commuting unit forks. In particular 53,084,231 still separates the older return language from the tested nested one; nothing here says it is globally absent.

## 2. An old neutral repair is distinguishable by the same internal arc

Take the two old skeletons $S=(46,12,56,48,35,2)$ and $T=(110,8,122,112,99,66,2)$, with root 4 and the uniform module $M$ at vertex 0. Use $h=H(M)$ and $p=Q_M(2)$, the latter allowing empty ordered slots. The old proof gives identical seam-visible profiles $(H,E,R)=(15p-h,5h,15p+5h)$.

Now reverse the **same retained arc on vertices 3 and 5**, in each skeleton. The differences are

    Δ_S(H,E,R) = (0,−2h,−2h),
    Δ_T(H,E,R) = (−2h,−2h,−2h).

This all-order identity is verified by two physical coefficient probes. In either edited or unedited skeleton the root remains the only outer in-neighbor of $M$. An individual Hamilton path can use at most two module runs; a root-starting path has at most one. In the ending-first two-cover, no path can use the arc from the root to $M$ (the first ends at the root and the second omits it), so there are at most two runs total, one per slot. Uniform incidences make the filling counts $h$ for one run and $p-2h$ for two. Hence both profile differences are linear in these two variables. A singleton and an ordered two-vertex module give $(h,p-2h)=(1,0),(1,2)$, so the literal subset counts determine the exact coefficients. The regularity premise is unchanged by the edit, which is why this is a proof for arbitrary $M$ and not interpolation from numerical data alone. The old Hamilton cycles never use the edited arc, so both orientations of both skeletons remain strong for any nonempty module.

For any fixed right rooted block of profile $(h_R,e,R')$, the difference of the two full seam responses is $2h h_R$. With a singleton on the right and an ordered $m$-vertex module, the two original seams have precisely the same absolute count $N=30\cdot2^m-1$, yet the same old edit leaves one at $N$ and sends the other to $N-2$. Thus these are **unbounded physically marked unit-ready counts obtained from a prior neutral repair**, and a literal demonstration that equality for exterior continuation does not certify equality for internal edits to the retained skeleton. The two starting representations differ in order; this is not an independent two-arc fork on one tournament.

## 3. A three-reading diagnostic for the state hidden by a count

Place the root respectively first, second and third in an ordered three-vertex left probe. The three counts against a fixed right block are

    N0=h_R−e+R',   N1=h_R+2R',   N2=h_R+4R'.

The second to third placement uses a single adjacent reversal. Together the three readings recover the whole profile: $R'=(N2-N1)/2$, $h_R=2N1-N2$, $e=h_R+R'-N0$. The last two deliberately fail to distinguish the actual right profiles $(23,8,24)$ and $(23,10,24)$; their readings are $(39,71,119)$ and $(37,71,119)$. Like the older split carry/sidecar construction, the missing coordinate is identified by giving an explicit collision and an explicit observation that breaks it, rather than treating the scalar as a complete state.

## Replay

- `node absolute_inverse.js 12`: absolute-target carry windows.
- `node target_fork.js 12`: paired unit-crossing screen.
- `node history_probes.js`: old repair coefficients, three diagnostic reads and literal strong target witnesses.
