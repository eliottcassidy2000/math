# Divisor norm ball and conductor rotation gate

[The proof](DIVISOR_NORM_BALL.md) identifies the proper divisors of p^kqr with every nonzero point of norm at most k² in Z+kZω, transporting both actual reflections. Its six-point rotation-stable core is exactly the conductor-ideal intersection. An explicit reduced-form proof recovers the classical even-conductor idoneal list k=2,4,8 without any global completeness assumption.

Run `python3 -B verify_norm_ball.py` and `python3 -B -O verify_norm_ball.py`; both reproduce the 72,384-check receipt. The source-specific norm mismatch rules out identifying the k=2 ball with the native norm-one/norm-three cancellation alphabet. No multiplicative, ideal-class, or dynamical identification is asserted.
