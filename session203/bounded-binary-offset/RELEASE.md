# Release verification

Independent audit: PASS, frozen proof SHA-256
ccc0f66fbdfc2ed19b5767b8e12ffbbc02dfc7e81d1c3c7ac0307bf21447bf8e.

Only the status line was promoted. The author verifier passes normally and
under python -O; outputs are identical. Removing the intentional proof hash
change, promoted output is identical to the frozen author output. The
independent verifier also passes normally and under -O and remains byte for
byte identical to its audited receipt. Dependencies are included as text.
No unreviewed mathematical exploration is included in this packet.

Publication verification is recorded separately by the releasing worker.
