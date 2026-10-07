# HN rotation-field tower (HYP-2277): spindle with arm sqrt(N) about an Eisenstein pivot needs rotation e^{i th},
# cos th = 1-1/(2N), field Q(sqrt(-(4N-1))).  Plane L_N = Q(sqrt-3, sqrt(-(4N-1))).
# Lower bound 4 if an Eisenstein spindle exists: N = 3n with n Loeschian (points at sq-distance N are
# forced equal-coloured in the unique 3-colouring of the triangular lattice iff 3 | N; N = norm of sqrt-3*y).
# Upper bound from residue colourings (s4_multiquad).
from s4_multiquad import best, sqfree
def loeschian(n):
    return any(a*a+a*b+b*b == n for a in range(0, int(n**.5)+2) for b in range(0, int(n**.5)+2))
lucky = {2, 3, 5, 11, 17, 41}
rows = []
for N in range(1, 61):
    m = sqfree(-(4*N-1))
    if m == -3:
        L = [-3]
    else:
        # independence mod squares with -3
        L = [-3, m] if sqfree(-3*m) != 1 else [-3]
    b, res = best(L, pmax=300)
    spindle = (N % 3 == 0) and loeschian(N//3)
    ub = b[2] if b else None
    rows.append((N, 4*N-1, m, spindle, ub, b[0] if b else None))
print(" N  4N-1  field         spindle(chi>=4)  residue-UB  via p   note")
for N, D, m, sp, ub, p in rows:
    note = []
    if N in lucky: note.append("lucky-Euler (Heegner 4N-1)")
    if sp and ub == 4: note.append("chi=4 EXACT")
    if (not sp) and ub == 3: note.append("chi=3 EXACT (no 4-chromatic UDG at all)")
    print(f"{N:2d}  {D:4d}  Q(sqrt-3,sqrt{m:5d})  {str(sp):5s}            {str(ub):4s}       {str(p):4s}  {'; '.join(note)}")
