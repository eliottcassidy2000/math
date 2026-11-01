# Exact costed FRACTRAN word certificates

`accelerate.py` computes the exact maximum repetition count of an instruction
word, including every higher-priority guard. It can certify a repeating word
forever, but is not a general halting decider. Finite certificates preserve
both endpoint valuations and the exact number of original fraction steps.

Read [THEORY.md](THEORY.md). Run `python3 check.py` and `python3 -O check.py`.
Only the Python standard library is needed. No large literal integer is
silently substituted for a compact factorized input.

`periodic_points.py` supplies an exact finite-horizon cap and rational-density
algorithm for periodic **points** with any specified finite period bound.
It reproduces PRIMEGAME's inherited period-two point density; no arbitrary
program's global period bound or eventual-basin classification is assumed.
Run it normally and with `-O` and compare the frozen periodic receipts.
