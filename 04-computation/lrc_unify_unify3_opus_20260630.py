"""
UNIFICATION 3: the tournament metagraph H-gradient (transitive H=1 <-> regular H=max) maps to the LRC
two columns: TRANSITIVE (H=1, cusp, binds at 3-cycle) <-> MEASURE/doublet; REGULAR (H=max, half-turn =
additive interval) <-> EXISTENCE/comb (AP = additive interval). Both extremes, one gradient.
"""
import math
print("THE METAGRAPH H-GRADIENT <-> THE LRC TWO COLUMNS (one unified gradient):")
print()
print(f"   {'metagraph extreme':>30} | {'LRC column':>28} | {'shared object':>26}")
rows=[("TRANSITIVE (H=1, the cusp)","MEASURE (off-cusp, doublet)","the MINIMAL relation"),
      ("  binds at the 3-CYCLE (C_3)","  binds at the DOUBLET","  3-cycle <-> doublet (mac-mini)"),
      ("  = OCF origin (odd cycle)","  = 4cos^2(3pi/7), Q(sqrt-7)","  the apex obstruction atom"),
      ("REGULAR (H=max, half-turn)","EXISTENCE (cusp, the comb)","the ADDITIVE INTERVAL"),
      ("  = connection {1..(n-1)/2}","  = AP {1..n-1} (the comb)","  both = the interval/Dirichlet"),
      ("  = max Ham paths (A038375)","  = M=n/Phi_6(n), Q(sqrt-3)","  the max-packing / tightest"),]
for r in rows: print(f"   {r[0]:>30} | {r[1]:>28} | {r[2]:>26}")
print()
print("  => the H-gradient (1 -> max) IS the measure->existence axis. Transitive/regular = doublet/comb")
print("     = 3-cycle/interval = Q(-7)/Q(-3) = proved/open. The SAME gradient, tournament & LRC.")
print()
# the master synthesis: ONE object, two regimes, the whole map
print("="*78)
print("THE MASTER UNIFICATION -- LRC(14) is ONE object in TWO complementary regimes:")
print("="*78)
print("""
   ONE OBJECT: doublet = empty tooth = cusp form f_14 = curve 14a = 7-cycle C_7
               = atom 4cos^2(3pi/7) = the genus-1 mode = the minimal relation.

   TWO REGIMES (the two columns, the two Heegner fields, the two cusps):

     MEASURE / Q(sqrt-7) / d=7 (-)        EXISTENCE / Q(sqrt-3) / d=14 (S)
     -------------------------------       --------------------------------
     off-cusp (proper core)                cusp (core = Z_7)
     the DOUBLET (per-level atom)          the EMPTY TOOTH (comb witness)
     rho_j >= 4cos^2(3pi/7)=0.198          M >= n/Phi_6(n) = 14/183
     density (the bulk)                    witness (the point)
     PROVED (THM-590, finite)              OPEN (covering-min >= 1/n)
     metagraph: transitive/3-cycle         metagraph: regular/half-turn/interval
     razor-thin: measure -> 0              razor-thin: witness robust

   THE CONJECTURE = the EXISTENCE column (Q(sqrt-3), the covering-min, the Phi_6/Eisenstein side).
   THE PROVED HALF = the MEASURE column (Q(sqrt-7), the doublet, the apex cyclotomic gap).
   THEY MEET at 7 = Phi_6(3) (the apex = the Eisenstein norm at n=3) -- the two Heegner worlds touch.
""")
