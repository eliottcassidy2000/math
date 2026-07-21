from math import comb, factorial
# MORE frames + the REPRESENTABILITY axis (basis order / forbidden set / image density).
# The key contrast: additive figurate image sets are SPARSE (density->0), the multiplicative Ham
# spectrum is DENSE (density 1/2, odds), the factorial arb spectrum is SPARSE (bands).
def is_triangular(m): 
    k=int(((8*m+1)**0.5-1)/2); return any((j*(j+1)//2==m) for j in (k-1,k,k+1))
def sum_of_3_triangulars(m):  # Gauss Eureka: every n>=0 is a sum of 3 triangular numbers
    ts=[t*(t+1)//2 for t in range(0,60) if t*(t+1)//2<=m]
    S=set(ts)
    return any((a+b in S or (m-a-b)>=0 and (m-a-b) in S) for a in ts for b in ts if a+b<=m)
# check Gauss on 0..40
gauss=all(sum_of_3_triangulars(m) for m in range(0,41))
print("Gauss Eureka (every n = sum of 3 triangulars), n<=40:", gauss)
# image density of triangular numbers up to X (additive frame image is sparse)
X=1000; tri=[t*(t+1)//2 for t in range(1,50) if t*(t+1)//2<=X]
print(f"triangular numbers <= {X}: {len(tri)} of {X}  (density {len(tri)/X:.3f} -> 0, ~sqrt)")
print("Ham-path spectrum <= 1000: all odds except {7,21} => 499 of 1000 (density 0.499 -> 1/2)")
print()
# NEW frames catalog with mode/representability
NEW=[
 ("centered triangular 1+3T_{n-1}", "ADD", [1+3*comb(n,2) for n in range(1,8)], "polynomial deg2"),
 ("q-triangular [n+1,2]_q at q=2", "q-ADD", [ (2**(n+1)-1)*(2**n-1)//((2**2-1)*(2**1-1)) for n in range(1,7)], "exp (q-def)"),
 ("doubly-triangular T_{T_n}", "ADD^2", [ (t:=n*(n+1)//2)*(t+1)//2 for n in range(1,7)], "polynomial deg4"),
 ("Gauss: n = sum of 3 triangulars", "ADD-basis", "every n>=0 (Gauss 1796)", "basis order 3"),
 ("Pollock: n = sum of <=5 tetrahedral", "DIM-basis", "conjectured (Pollock 1850)", "basis order ~5"),
 ("tournament scores sum to T_{n-1}", "PARTITION", "Landau: s_1<=..<=s_n, sum prefix >= C(k,2)", "constrained partitions of T_{n-1}"),
 ("Ham spectrum = odds\\{7,21}", "MUL", "monoid, 2 forbidden", "density 1/2"),
 ("arborescence Sum a_r", "FAC", "(n-1)!-bands, forbidden {4,5,8,11..}", "density->0"),
]
print(f"{'FRAME':40}{'MODE':12}{'representability/image':42}growth/order")
for name,mode,rep,g in NEW:
    print(f"{name:40}{mode:12}{str(rep):42}{g}")
print("\nDENSITY LADDER of the image/spectrum (how many integers <=X the frame realizes):")
print("  additive figurate (triangular, k-gonal): ~sqrt(X)  (density -> 0)")
print("  multiplicative Ham monoid: X/2 - 2       (density -> 1/2)   <-- DENSEST, and the {7,21} live here")
print("  factorial arborescence: ~log X bands     (density -> 0)")
print("  => the MULTIPLICATIVE frame is where a *finite* forbidden set ({7,21}) can even exist;")
print("     additive & factorial frames are too sparse for a finite forbidden set to be meaningful.")
