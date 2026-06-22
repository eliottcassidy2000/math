"""
The 3 recursion modes are PARITY-STRATIFIED on the half-tiling (complement-quotient) (kps-S31r).
h(n) = floor((n-1)^2/4) = #cells of the half-tiling (THM-549, the T<->T^op fundamental domain).
  EISENSTEIN (EVEN n=2k): A+B-C, h(n)=2h(n-1)-h(n-2), char (x-1)^2, reductions {1,2}, PRONIC k(k-1).
  LEGENDRE  (ODD n=2k+1): A+B-C+D-E-F+G, char (x-1)^3(x+1), SQUARE k^2.
      Slot sizes are {1,2,3,4}: A,B reduce by 1; C,D by 2; E,F by 3; G by 4.
      The scalar projection is h(n)=2h(n-1)-2h(n-3)+h(n-4) because -C+D cancels at depth 2.
  MOBIUS    (ALWAYS): A+B+C-D-E-F+G, inclusion-exclusion over the 3 modes A,B,C (reduce by 1,2,3).
The SIZES of the subtournaments (letters) DIFFER per mode: Eisenstein slots {1,2};
Legendre slots {1,2,3,4} with net coefficients {1,3,4}; Mobius slots {1,2,3}.
For LRC(14): 14 EVEN => Eisenstein, half-tiling pronic = 7*6 EXPOSES the apex prime 7 = n/2.
"""
def h(n): return ((n-1)**2)//4
def even_rec(n): return 2*h(n-1)-h(n-2)            # A+B-C
def odd_rec(n):  return 2*h(n-1)-2*h(n-3)+h(n-4)   # A+B-C+D-E-F+G
print("n  parity  h(n)=floor((n-1)^2/4)   shape          even-rec A+B-C   odd-rec A+B-C+D-E-F+G")
for n in range(4,16):
    k=n//2
    shape = f"PRONIC {k}*{k-1}={k*(k-1)}" if n%2==0 else f"SQUARE {((n-1)//2)}^2={((n-1)//2)**2}"
    e = even_rec(n); o = odd_rec(n)
    tag = "EISENSTEIN" if n%2==0 else "LEGENDRE"
    print(f"{n:2d}  {('EVEN ' if n%2==0 else 'ODD  ')}{tag:10s} {h(n):3d}   {shape:16s} {e:3d}{'=h OK' if e==h(n) else 'X':6s} {o:3d}{'=h OK' if o==h(n) else ' X'}")
print("\nTHE THREE MODES' REDUCTION-SIZE-SETS (subtournament sizes, DIFFERENT per mode):")
print("  MOBIUS (always)  : reduce by {1,2,3} = Modes A,B,C (incl-excl skeleton, char of the 7-term lattice)")
print("  EISENSTEIN (even): reduce by {1,2}   = order-2 (x-1)^2, PRONIC, the complement-halving / binary Z[i]")
print("  LEGENDRE (odd)   : slots reduce by {1,2,3,4}: A,B n-1; C,D n-2; E,F n-3; G n-4")
print("                     net scalar offsets are {1,3,4} only because -C + D cancels at n-2")
print("\nLRC(14) COMPOSITION = parity factorization 14 = 2 * 7:")
print(f"  14 EVEN  -> EISENSTEIN: half-tiling = 7*6 = {h(14)} (pronic), k=7=14/2 EXPOSES the apex prime 7")
print(f"   7 ODD   -> LEGENDRE  : half-tiling = 3^2 = {h(7)} (square), the QR(7)={{1,2,4}} / NQR={{3,5,6}} character")
print(f"  => the even Eisenstein mode FACTORS OUT THE 2 (complement Z/2 = sector-halving x->-x, THM-280/549),")
print(f"     leaving the odd apex 7 to the Legendre mode. 14=2*7 is literally Eisenstein(even) o Legendre(odd).")
print(f"  Parity ALTERNATES down the tower: 14e 13o 12e 11o 10e 9o 8e 7o ... each step its own mode.")
