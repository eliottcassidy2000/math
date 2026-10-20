EFFECTIVE CUBIC ORDER BOUND FROM BLACK HAMILTON-PATH DEFECT
2026-10-20
Status: proved; independently audited; supplementary finite checks passed.

For a strong tournament T on n>=3 vertices let h be the number of its
Hamilton-path automorphism orbits, b the number of blue orbits, and Delta=h-b.
A Hamilton word is blue when its positional reversal is an anti-automorphism;
non-self-converse tournaments have b=0.

For any integer d>=0 with Delta<=d, set s=d+2 and let r be the least power of
two at least 4s(d+1). Then

    n <= 2sr < 16(d+2)^2(d+1),
    Delta > (n/16)^(1/3)-2.

This gives an elementary explicit cubic order bound and therefore general
bounded-black-defect finiteness without a Dixon automorphism-size bound or
a Thomassen Hamilton-cycle counting bound. It retains Moon's classical
strong-tournament cycle theorem and elementary group and degree arguments.
It applies to both self-converse and non-self-converse strong tournaments.

This deliberately coarse general theorem does not replace or alter the
sharper Delta=2 classification. It does not claim sharp constants, a sharp
extremal function, or a classification at other individual defects.

CONTENTS
PROOF.txt: complete theorem and proof, with primary reference.
AUDIT_RECEIPT.txt: independent audit verdict, scope, and provenance hashes.
VERIFICATION.txt: exact finite-check coverage and reproduction commands.
exhaustive_small.cpp / exhaustive_small.out: all labelled tournaments n=3..6.
check_rotation_frames.cpp / check_rotation_frames.out: larger cyclic frames.
provenance/ORIGINAL_CANDIDATE.txt: candidate mathematical text retained,
with an internal author label omitted; its historical "audit pending"
status is superseded by AUDIT_RECEIPT.txt.

The general theorem is proved analytically; the finite checks are supporting
validation, not premises. The source programs require only a C++17 compiler.
