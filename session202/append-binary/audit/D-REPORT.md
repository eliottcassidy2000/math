# Independent audit: actual source-count interface

**PASS, without mathematical correction.** The D-interface addendum extends the base append-only obstruction in the stated physical grammar. It does not extend the theorem to arbitrary graph operations.

Frozen proof: `D_INTERFACE.md`, SHA-256
`4d74882efe34ccaa95970a31ce50d942693f36d01b9bd5d142bf568be530b5b2`.
Frozen author verifier: `verify_d_interface.py`, SHA-256
`ebfb109a9cd88073f407c912fb9a7776458b950e1a34fa8b9b48c451c3004992`.

The physical inequalities U>=D and E<=(D-1)/2 imply F>=(3D+1)/2. Hence two D-changing Z/J letters, with any number of intervening O letters, give at least (9n+5)/4>2n+1. Zero D-changing letters cannot give either child when n>=3. Exactly one such letter remains. Prefix O letters are absorbed into the existential parent representative, and suffix O letters preserve the output D.

For a successful intermediate source s with D=n and F=2n±1, a final O is impossible: it gives E=0 and U>=2D. The source is not the seed because n>=3, so its actual last letter is Z or J. Its predecessor t has F(t)=n, and the required offset is precisely one of the base theorem's one-letter offsets. The converse is also physical: if t.c realizes a base one-letter witness, take s=t.c as the D=n source and append Z or J. Its new D is F(s), as required. Thus the same complete plus/minus parent sets apply.

The n=3 boundary is treated correctly and is not covered merely by invoking the base theorem's n>=5 growth bound: its triangle source has D=3,F=5, and the base one-letter slice classification excludes the plus child7. The singleton D=1 is explicitly outside the addendum's domain. Empty append words are harmless because n cannot equal either child for n>=3.

The independent `check_d.py` imports only this audit's independent `check.py`, never author code. It independently finds61 plus parents and20 minus parents through D=3333, with both branches exactly at5,9,17. This scan is complete because U=D+E±1<=(3D+1)/2<=5000 for these targets. It rechecks the affine D-growth inequality on all126,910 physical profiles with U<=5000 and literally counts parent and child graphs for all six signed controls. Normal/optimized outputs agree. Author normal/optimized receipts also replay exactly.

Reproduction:

    python3 check_d.py > check-d-normal.json
    python3 -O check_d.py > check-d-optimized.json
    cmp check-d-normal.json check-d-optimized.json
    cmp check-d-normal.json D-NORMAL.json

Finite controls support the separate all-order proof. The eventual odd-count threshold remains open.
