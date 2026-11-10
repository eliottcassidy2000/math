---
id: THM-5092
title: "Four-test tail content and invisible reroot escape"
status: "PROVED; independently audited. Chan–Pak completeness remains OPEN."
session: Session237, 2026-11-09
source: 05-knowledge/results/session237_chanpak_synthesis_20261109/README.md
proof: 05-knowledge/results/session237_chanpak_synthesis_20261109/symmetry_grafts/PROOF.md
audit: 05-knowledge/results/session237_chanpak_synthesis_20261109/symmetry_grafts_audit/AUDIT.md
---

# THM-5092: Four-test tail content and invisible reroot escape

For the physical fixed-port state z=(h,R,e,d), the inherited maps are
S(z)=(h,2R,0,d) and B(z)=(h+R-e,2R-d,d,d). The gcd of every descendant
Hamilton-path count equals gcd(h,R,e,d), already certified by the four words
empty, B, SB, BB. Content divides every future content; an actual rooted
automorphism group divides it, but equality need not hold.

The strong rigid tournament T=[24,1,19,22,2] at root1 has state(9,9,3,3),
so every fixed-port continuation is divisible3. Rerooting at0 and taking B
gives19; taking converse at the old root and B gives29. Rigidity alone does
not remove arithmetic confinement.

Prepend a source to T to obtain X=[62,48,2,38,44,4], rooted at2. It has
order6 and state(9,9,3,3), identical to C3=>C3 rooted in its second component.
Only X has a content-one reroot; their converse states differ too. Therefore
even (order,h,R,e,d) cannot determine primitive rerootability or converse
state, despite determining every fixed-port S/B continuation count.

This is an information-loss and operation theorem, not a coverage theorem.
Actual subgroup orbit omissions are inherited from THM-4987/4990 and are
not applied to formal division by arithmetic content or endpoint counts.

[Full proof](../../05-knowledge/results/session237_chanpak_synthesis_20261109/symmetry_grafts/PROOF.md) · [Independent audit](../../05-knowledge/results/session237_chanpak_synthesis_20261109/symmetry_grafts_audit/AUDIT.md).
