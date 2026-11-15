# The leaf can remember its parent, but the two signs cannot share that memory

We use literal tournament growth: a newborn beats exactly its selected old
parent, and loses to everyone else. The reversed arcs in chronological order
form a tree. Write $\delta_q=H(L_qT)-2H(T)$. No scalar identity below is used
as a substitute for an actual tournament or continuation.

## Results

1. [Fan law and moving fan vertices](FAN_CHANNELS.md). The first leaf freezes
   its parent's defect forever under further attachments **at that parent**:
   $\delta_q'=\delta_q+2(d_q-K_q)$, where $d_q$ counts paths ending at $q$
   and $K_q=H(T-q)$. An initially negative unit defect need not stay negative
   unit: one chronological example reads $5,9,15$, not $5,9,17$.
   Each individual fan leaf has an exact ternary/binary formula. Eventually
   only the parent itself is a unit port along a stable +1 fan.
2. [Cross-port law](CROSS_PORTS.md). Two fixed, same-direction opposite unit
   ports cannot support arbitrary switching. Even allowing the port to move,
   chains restricted to repeat-stable roots can switch identities at most once.
   A separate [movable-root census](MOVABLE_PORTS.md) finds no positive-unit
   move retaining both signs at sources through order eight; this is not a
   theorem for unrestricted higher-order sources.
3. [Opposite-end construction](OPPOSITE_PENTAGON.md), with a second
   [strong seed](OPPOSITE_FANS.md). From the regular pentagon, either first move
   has count 33. Along either axis, continuing at the same end gives $-1$;
   switching end gives $+1$. This yields same-order pairs
   $16\cdot2^j+1,16\cdot2^j+3$ for every $j\ge2$.
   Once both ends have been used, the two designated defects become $-3$.
   Thus an exact width-two window repeats forever, but the switch does not reset.
4. [Changing the source](PARALLEL_FANS.md) yields thirteen consecutive odds
   $2^j+1,2^j+3,\ldots,2^j+25$ for every $j\ge10$, and eleven around
   $17\cdot2^j$ for every $j\ge6$. Each arm has a certified chronological
   core on at most nine vertices, and the target offset chooses the arm.
5. [All chronological parent choices](PARENT_COMPILER.md). The whole one-parent
   language admits an exact bounded inverse at each prescribed count via SCC
   stripping and Busch's bound. A literal replay through eleven vertices
   supplies every odd count through 365 except 7 and 21, including witnesses
   for all three local misses 23, 33 and 35. This is finite coverage, not an
   all-count theorem.

## Proof target

One-site stability and a complete menu at one height are different from a
menu that returns after either move. The pentagon separates them physically:
the bulk term survives switching, but the constant term changes from the
boundary atom to the mixed-contact atom. One direction for a genuine return
is to change cores at that transition while matching both bulk and offset;
another is to avoid returning to the same menu and instead grow the certified
offset bank enough that its absolute windows overlap. The finite bank here
has logarithmic capacity. Neither route is yet a proof of Chan–Pak.

Earlier chronological tree structure: [Session 199](../../../session199/hidden-structure/PROOF.md).
External problem and lower bound: [Chan–Pak](https://arxiv.org/abs/2308.10214),
[Busch](https://doi.org/10.37236/1141).
