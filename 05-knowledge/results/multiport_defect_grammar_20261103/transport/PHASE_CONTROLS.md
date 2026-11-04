# Existing field and literal hostile controls survive the forest readouts

**PROVED readout compatibility; FINITE-EXACT small replay.** This note adds no
physical Frobenius action and changes neither `PROOF.md` nor its checker.

For the forest state, let d be its total number of components. Retain every
physical root and spare that will be observed, including39,50 and6 here.

* H is the total weight with d=1.
* `e_r` and `s_r` restrict that layer to a component ending/starting at r.
* `E_r(2)` and `S_r(2)` restrict layers d=1,2 to a component ending/starting at r.
* `J(r,u)` restricts d=2 to different components ending at r and starting at u.
* `H(T-r)` restricts to a singleton-r component and d=1 or2; d=1 supplies the
  singleton-T convention `H(empty)=1`.
* Ordered two-slot covers allowing empty slots equal `2H+2pc_T(2)`.

Consequently the inherited integer golden observer is the exact readout

```
(H,X,Y) = (H, (S_r(2)-e_r)/2, (S_r(2)+e_r)/2).
```

The halves are integral: the converse single-in-port birth has count
`H+S_r(2)-e_r`, and both Hamilton counts are odd by Rédei. The saved decoder
is in `zero_ready_strong_switch_20261103/history_phase/validate.py`, functions
`observe` and `chart`; its earlier conventions are documented in
`return6325_completion_20261103/TRANSPORT_AND_RECOVERY.md`.

The two current field decoders, with all coordinates modulo7, are

```
F49 plus charge: (a,b,c)=(X+4(H-Y), Y, H-Y),
                  Frobenius: (a,b,c) -> (a+4b,-b,c).

F343 vector codec: xi=H+Xt+Yt^2,  t^3+t^2+1=0,
                  Frobenius: (H,X,Y) -> (H+4X+6Y,2X,3X+4Y).

Inverse coordinate change: (H,X,Y)=(b+c,a-4c,b).
```

The inverse is an F7-vector-space identification, not a ring isomorphism
`F49×F7 ≅ F343`; their displayed Frobenius maps have orders2 and3. Neither
field replaces integer precision modulo49. After an actual graph change,
recompute the integer observer at its actual retained root and decode again.
An accidental field equality does not certify a tournament automorphism,
native operation, or ancestral clock step.

The small replay uses the exact source rows from
`zero_ready_strong_switch_20261103/history_phase/validation_new.json`:

| Physical root on the same170807 graph | Golden (H,X,Y) | F49/charge | F343 vector |
|---|---|---|---|
|39|(170807,1344477,1344530)|(2,5;2)|(0,1,5)|
|50|(170807,10747054,10747107)|(3,0;0)|(0,3,0)|

The root39 cubic orbit continues `(6,2,2),(5,4,0)`; the root50 orbit continues
`(5,6,2),(6,5,5)`. Joint responses with retained spare6 and between39/50
are also checked against their actual saved records.

The literal81/27 hostile is read from the two **actual12-vertex matrices** at
`return6325_completion_20261103/halving_response.json`, key `hostile.rows`.
Cut256 has sole old in-port8, so it is a legal uniform-bulk birth with8 active.
The exact readout is `H'=H+E_8(2)-s_8`:

| Matrix | H | E8(2) | s8 | Birth H' |
|---|---:|---:|---:|---:|
|first|462915|646002|47142|1061775|
|second|462915|963549|31185|1395279|

Fresh endpoint DP checks both the formula and the literal newborn graphs.
Thus the forest interface distinguishes this named-cut control despite equal
H. The81/27 automorphism labels remain inherited source assertions; this
small replay does not recompute group orders or the older full cut census.

Reproduce with `python check_phase_controls.py` in this directory; `python -O
check_phase_controls.py --output phase_controls_optimized.json` gives the same
receipt bytes. The two receipts pin the inherited data and this small checker.
The established forest proof/checker hashes are unchanged.
