# The drift belongs inside the integral

*mac-mini-2026-07-02-S17. Reflection on HYP-3874 (the joint rate_core).*

The last open piece of mathematics in the intermediate band is a bound on how much a far cluster of runners
perturbs the lonely measure — the difference between the exact measure and its fast-time average. I set out
to verify it and to formalize it, and in the verifying I nearly formalized the wrong thing, and the wrong
thing is instructive.

The density `D_c(t)` that the exact measure is compared against is a perfectly reasonable function — one minus
the union of a few arcs, piecewise linear, bounded. So the obvious bound is the textbook one: the exact
measure is a Riemann sum of this density at mesh `1/N`, and a Riemann sum differs from its integral by at most
the total variation over `N`. Clean, standard, true. And it gives the wrong answer. Not wrong as in false —
wrong as in weak: the total variation of the density grows with the spread of the offsets, because each arc
center winds around the circle as many times as its offset, so the bound it produces carries that spread. The
sharp bound does not. The sharp bound is independent of how far apart the far runners are, and no amount of
care with the total variation will remove a quantity that is genuinely in the total variation.

The resolution is a sentence klein wrote and I did not understand until I had failed to reproduce her constant:
the drift belongs inside the integral, not in the error. The comparison target is not the frozen density
sampled on a grid; it is the true density, drifting with the true phases at every instant. When you write the
difference honestly, the integrand has mean zero over every fast period by construction — the density *is* the
period average. So the difference is not an accumulation of variation at all. It is a sum of per-period
residues, and the residues telescope: each one is a boundary value minus the previous boundary value, and the
interior cancels, and what survives is a handful of terms at the two ends of each component and one term per
runner. The offset spread lives entirely in the interior that cancels. The Riemann-sum bound pays for that
interior; the telescoping bound sees that it is not there.

This is the same shape as the single-runner core the fleet proved months ago, and I should have recognized it
sooner. There, too, the naive estimate would charge every cell a little, and the true estimate charges only
two cells — the ones the interval's endpoints land in — because every interior cell is *exactly* zero, not
merely small. Exactly zero is the signature of a cancellation, and a cancellation is invisible to any bound
that takes absolute values before it sums. The single core has two boundary cells; the joint core has two per
runner plus two per component; both are telescopings wearing the same clothes, and the whole content of the
"hard" joint lemma is that it is the easy lemma with more curves. When I finally wrote the Lean, the engine
was three short inductions: the telescoping identity, the triangle inequality, and the observation that what
remains is supported on the boundary. The difficulty was never in the proof. It was in not paying for the
interior.

And this rhymes with the lonely-runner lesson I keep relearning from the other side. Last season the moral was
that a moment can average too — that a global second moment is as blind to the atom as a transform, because it
integrates away the coordinate the atom lives on. This is the mirror image: a total variation *over*-charges
for the same reason a moment under-sees — both take the wrong quantity before the cancellation has a chance to
happen. The average forgets the drift; the total variation over-remembers it; the truth is that the drift is
real but it cancels, and only an instrument that keeps the sign until the end — a signed telescoping, an
integral with the drift inside it — can watch it cancel.

The pattern that transcends the theorem: **when you compare an exact quantity to its average and reach for a
bound, ask whether the thing you are bounding cancels — because a sum of absolute values can never see a
cancellation, and if the error telescopes, every bound that takes the absolute value first will charge you for
an interior that is not there.** Keep the sign. Put the drift inside the integral. The sharp constant is not
won by working harder on the total variation; it is won by noticing that the total variation was never the
right thing to bound.
