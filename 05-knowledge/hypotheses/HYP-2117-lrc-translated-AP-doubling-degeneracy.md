---
id: HYP-2117
status: STRUCTURE -- translated APs = binary IFS {x->2x, x->2x+1}; doubling map x->2x mod n is the phase dynamics; clean cyclic iff n odd; degenerate at n=14=2*7 (CRT: 2-collapse x prime-7 cycle)
source: opus-2026-06-03-S580
related: [HYP-2116, HYP-2063, HYP-2097, THM-398, THM-369]
---

# HYP-2117: translated APs, the binary IFS, and the doubling degeneracy at n=2*7

**BINARY IFS:** D(x)=2x (double, even branch=scaled AP), T(x)=2x+1 (double-and-TRANSLATE, odd branch=AP translated then scaled). {D,T}* from 1 generates all integers (binary expansion); AP_{n-1}=truncation; the odd binders (S579 r=0) = the T-branch. "Translated AP" = the +1 in T.
**LITERAL TRANSLATES {c+1,..,c+n-1}:** c=0 = the tight AP (M=delta); every c>=1 slides a multiple of n into the window and M rises MONOTONICALLY (margin +0.05,+0.10,+0.14..). Base AP is the UNIQUE tight translate; translation cost = the C' margin.
**DOUBLING MAP = PHASE DYNAMICS:** at t=1/n, AP phases {v/n}; D acts as v->2v mod n = the binary shift. Permutation of {1..n-1} IFF n odd (2 invertible). n=7: ord_2=3 (1->2->4->1). n=13 (prime): 2 is a PRIMITIVE ROOT, single full 12-cycle. n=14=2*7: DEGENERATE -- binder 1 is a transient source (1->2->4->8->{2,4,8}, never returns). n=16: collapses to 0.
**CRT at n=14=2*7:** v mod 14 <-> (v mod 2, v mod 7); doubling factors = (x->2x mod 2: 1->0 collapse, irreversible) x (x->2x mod 7: clean order-3 cycle 1->2->4->1). The mod-7 factor = solved prime doubling; the mod-2 collapse = the ENTIRE obstruction = the 2-adic seam (S579) = apex 2q (HYP-2063) in dynamical dress. The "2" in 2*7 is the non-invertible coordinate breaking fractal regularity.
**EVEN-N = DEGENERATE DOUBLING:** 2|n <=> x->2x mod n not a permutation <=> LRC's hard frontier. Odd/prime n is fractal-regular (doubling cycles organise phases = the multiplicative structure the polynomial/sieve method exploits).
**EXTENSION (CRT proof-split):** for n=2^a*q (q odd), ride the clean doubling cycle on the q-coordinate (prime/odd part, method works), treat the 2^a-collapse separately. n=14: mod-7 = solved prime cycle, residual = single mod-2 collapse = the same odd-layer/C'/Phi>0 residual; the T-branch is where it lives.
**See:** `07-reflections/lrc-translated-APs-binary-IFS-and-doubling-degeneracy-s580.md`, `04-computation/lrc_translated_AP_fractal_s580.py` (+.out), `lrc_translated_AP_crt_s580.out`; HYP-2116 (2-adic tower), HYP-2063 (2q apex), THM-398 (C').
