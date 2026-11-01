# PRIMEGAME: all-input basin decidability

The fixed, never-halting fourteen-fraction PRIMEGAME has exactly three eventual
behaviors: a two-cycle, a six-cycle, or exact coalescence with a spectator-scaled
seed-2 orbit. A new positive-control budget and an exact onehot escape bit give
a total polynomial-time basin decider. The three basin densities exist and are
polynomial-time computable reals in precision bits; no closed values or general FRACTRAN/Collatz conclusion follow.

See [THEORY.md](THEORY.md), the runnable [decider](primegame_basin.py), and
[verification](verify.py). Input proof sources are pinned and retained in
[SOURCE_PINS.json](SOURCE_PINS.json).
