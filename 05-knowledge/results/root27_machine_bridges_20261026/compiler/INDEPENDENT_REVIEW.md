# Independent review of the marked-program compiler

2026-10-26. **PASS in the stated scope.** The producer proof and code were
read independently; both producer modes passed with identical output.
The pre-padding receipt hash was
`9226b3d46a764c99db43f7f077fbe906f8e3b3fe98c3cdcff9986298d2b3ec94`.
After the independently audited padding addendum, both modes again pass
identically, with saved `normal.json` / `optimized.json` receipt hash
`b43c96480b8f3bc2755d973f6a3704434f78632e16203b4794b95a197b8d6ab7`.
No producer file was changed by this reviewer.

The unary header recovers s before target widths are needed. The fixed-width
records and zero padding therefore recover every table field. The slot
count is `C(n,2)-(n-1)=C(n-1,2)`. Injectivity and transported-frame
equivariance follow directly; neither requires unframed tournaments to
have unique presentations. The decoder can accept valid padded tables at
larger orders than the encoder chooses; the proof asserts injectivity of
the canonical encoder, not equality of its image with the entire decoder
domain, so this is not a gap.

The K operation flips the first header bit on every encoded image. It
therefore leaves this positive-state code language immediately. This proves
failure of closure for this code, not a general obstruction to simulation
with some other code. The two-state hostile retains targets, reads and
directions but changes a write; its two-step loop and three-step halt are
literal configuration facts. A control graph cannot replace the table.

The claimed runtime carrier is separate: ordered, colored tape-cell
tournaments require an unbounded vertex supply, the exact transition table,
head marker and tape data. They do not turn the finite program tournament
into an autonomous graph dynamical system. No graph move or complexity
reduction is certified merely by this serialization.

**Padding addendum, proof audit.** For any k>=1, remap every old halt target
s to s+k, leave all old live targets unchanged, and append k arbitrarily
specified states. A run from an old live initial state never enters these new
states: induction preserves the complete tape/head/old-state configuration
at every nonhalting native time. At halt only its label changes. Thus all
blank execution times and halting predicates agree while state count grows
arbitrarily. This is not a machine isomorphism. The added `pad_states`
implementation performs exactly this target remapping and gives appended
states self-loops. Its `prefix` comparison retains tape, head and native
time, renaming only the terminal marker. All 262 added finite controls
were replayed successfully, including all 64 one-state tables on four
tapes and the two hostile tables padded to 5, 15 and 25 states. These
checks support the stated induction; they do not replace its all-time proof.

## Independent operation controls

The following block uses producer encoding/decoding, but changes actual
tournament arcs from independently calculated serialized field positions.
It checks all write, direction and target replacements for deterministic
test tables at state counts crossing several target-width boundaries.
This is an operation check, not another replay of the producer's small-table
enumeration. Run from the repository root, normally or under `python -O`.

```python
import importlib.util
import json
from pathlib import Path

def need(ok, why):
    if not ok:
        raise RuntimeError(why)

path = Path('05-knowledge/results/root27_machine_bridges_20261026/compiler/check.py')
spec = importlib.util.spec_from_file_location('compiler_producer', path)
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)
checks = 0
for s in (3, 4, 5, 7, 8, 16):
    width = s.bit_length()
    for seed in range(4):
        machine = tuple((j % 2, 1 if (j + seed) % 2 else -1,
                         (j * j + seed) % (s + 1)) for j in range(2 * s))
        rows, frame = p.encode(machine)
        slots = [(i, j) for i in range(len(rows)) for j in range(i + 2, len(rows))]
        for record in range(2 * s):
            old = machine[record]
            for field, allowed in enumerate(((0, 1), (-1, 1), range(s + 1))):
                for value in allowed:
                    if value == old[field]:
                        continue
                    expected = list(machine)
                    changed = list(old)
                    changed[field] = value
                    expected[record] = tuple(changed)
                    base = s + 1 + record * (2 + width)
                    if field < 2:
                        offsets = [base + field]
                    else:
                        offsets = [base + 2 + j for j in range(width)
                                   if ((old[2] ^ value) >> (width - 1 - j)) & 1]
                    edited = list(rows)
                    for offset in offsets:
                        a, b = slots[offset]
                        edited[a] ^= 1 << b
                        edited[b] ^= 1 << a
                    need(p.decode(tuple(edited), frame) == tuple(expected), 'field edit commutes')
                    need(p.encode(tuple(expected))[0] == tuple(edited), 'exact edited image')
                    checks += 1
        permutation = tuple(reversed(range(len(rows))))
        rr, pp = p.relabel(rows, frame, permutation)
        kk, kp = p.relabel(p.k_move(rows, frame), frame, permutation)
        need(kk == p.k_move(rr, pp) and kp == pp, 'K transports with the frame')
print(json.dumps({'status': 'PASS', 'field_edit_squares': checks,
                  'frame_K_squares': 24}, sort_keys=True))
```
