# Two actual ordinary counts expose the binary cycle count without signed fitting

The [cycle digit calculation](CYCLE_DIGIT_GAP.md) extracted
\(2K\) through four signed ordinary counts while retaining the
cycle gadgets. A more direct physical choice is possible if
we delete the first hub from the binary cycle encoder.

Start with the unpadded binary word \(d\) for \(K>0\).
Delete \(u\) from that encoder, then consider deleting or
keeping the vertex \(a\). In both graphs the remaining
spine is transitive, with \(b\) before \(a\) if both are
present. Both \(b,a\) beat the remaining apex \(v\).
Thus the two actual one-apex contact words are \(d0\)
and \(d00\). They give an induced pair of strong tournaments:
the first bit is 1 and the last is 0, so every spine vertex
lies on a directed cycle through the apex.

The exact one-apex word recurrence from
[APEX_RENEWAL.md](APEX_RENEWAL.md) is
\(H(2p)=H(p)+p+(p\bmod2)\).
At the already even word \(p=d0\) its value is \(2K\),
so appending that zero adds **exactly** \(2K\) ordinary
Hamiltonian paths. No difference of differences is needed.
In particular
\[
 H(d0)=H(d)+K+(K\bmod2),\qquad
 H(d00)=H(d)+3K+(K\bmod2)\le4K+1.
\]
The last inequality is \(H(d)\le1+K-(K\bmod2)\)
from the nonnegative even-endpoint correction in the contact
word. When the binary digits contain at least two ones the
correction has size at least two, giving the improved
\(\le4K-1\). The two orders are
\(\lfloor\log_2 K\rfloor+3,\lfloor\log_2 K\rfloor+4\).
This makes a literal cycle-to-ordinary **gap** bridge in
the same initial graph. The unpadded choice matters for
the size claim; leading zeros would only add neutral
earlier spine vertices.

All positive even gaps, including induced pairs of linear
count cost, were already known in the repo
([Session 243](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session243_prime_gap_bridges_20261110/README.md)).
This proof gives the explicit bound \(4K+1\), simple
orders, and the correspondence with the binary cycle
encoder. It does not choose the baseline. In fact its lower
count \(H(d0)\ge K+1\) is already strictly higher than one;
subtracting it arithmetically is not a tournament that
realizes \(1+2K\). The missing Chan–Pak step remains an
absolute choice of the count at either end, not even-gap
availability.

Run positive_cycle_gap.py for literal cycle and path counts
on the same original graphs and their two induced deletions.
