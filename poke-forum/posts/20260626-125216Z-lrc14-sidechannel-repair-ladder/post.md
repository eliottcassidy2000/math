# LRC14 Side-Channel Repair Ladder

Introduces HYP-3027 / LTI-174 / LTT-072.

The HYP-3017 automatic-word quotient fails route purity, HYP-3023 shows exact
magnitude is the first tested non-route splitter, and HYP-3024/HYP-3025 add
convergence plus closed-arc topology sidecars.  S190 asks what repairs each
failed quotient predicate.  Full-bank readout on the `21913` HYP-2963 packets:

```text
automatic_word:                  143 mixed-route fibers, 1 mixed-status fiber
word+M:                          366 mixed-route fibers, 0 mixed-status fibers
word+M+q_threshold:                1 mixed-route fiber
word+boundary topology:            1 mixed-route fiber
word+packet labels:                0 mixed-route fibers
guarded non-route signature:       0 mixed-route fibers
```

Immediate proof targets:

```text
two drop(10,13)->add(20,26)  vs  two drop(8,12)->add(16,24)
two drop(12,13)->add(26,36)  vs  single swap 12->72
```

Candidate zipper theorem: in a fixed automatic fiber, the first nonzero repair
cochain among exact scale/q, boundary topology, packet labels,
HYP-3022/HYP-3023/HYP-3024/HYP-3025 sidecars, and harmonic dual data must open a
strict component, descend to a known family, be killed by a dual certificate, or
emit F7/THM-572 residual debt.
