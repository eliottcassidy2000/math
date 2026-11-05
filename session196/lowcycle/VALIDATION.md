# Reproduce the checks

Use ordinary Python with assertions enabled. These research verifiers use
`assert` for test conditions; **do not use `python -O`, `python -OO`, or
`PYTHONOPTIMIZE`**. Optimized execution removes those checks and is not a
verification run. The reported runs used unoptimized Python throughout.

From this directory, run:

```
python apex_transfer.py
python verify_language.py
python bounded_image.py 67
python bounded_image.py 501
```

`verify_language.py` requires SymPy for the exact rational-function identity.
The other listed programs use the Python standard library. All graph inputs
are constructed locally; there are no repository or network dependencies.

Expected results:

- transfer: 1,364 independent graph comparisons, 251 raw / 32 minimized states;
- language: 66-state complete equivalence product, all 32 quotient states
  distinguishable, stored machine agrees exactly with the rebuilt machine;
- bounded image 67: 3,388 states, 3,999 transitions, strong odd gaps 7,21,67;
- bounded image 501: 50,749 states, 58,358 transitions, strong odd gaps 7,21,67;

Independent review of the frozen primary packet is recorded separately, with
independent implementations and conclusions. The deletion extension is a
separate packet and is not needed to reproduce the primary results above.

`explore.py`, `explore_spectrum.py`, and their output, if retained in scratch,
are exploratory aids, not needed for these proofs. Publication can omit them.
