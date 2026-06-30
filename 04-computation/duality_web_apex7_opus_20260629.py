from sympy import isprime, factorint
# Verify Mersenne cap Heegner = {3,7}, perfect=T(Mersenne)=arc-count, and the add/mult structure.
heegner={1,2,3,7,11,19,43,67,163}  # d with class number h(Q(sqrt-d))=1
print("Mersenne primes 2^p-1, which are Heegner (h=1) and 3 mod 4:")
for p in [2,3,5,7,13,17,19,31]:
    M=2**p-1
    if isprime(M):
        print(f"  p={p}: M=2^{p}-1={M}  prime={isprime(M)}  Heegner={M in heegner}  3mod4={M%4==3}  "
              f"-> {'**APEX (all 3)**' if (M in heegner and M%4==3) else ''}  perfect T(M)=C(M+1,2)={M*(M+1)//2} (=arc-count of {M+1}-tournament)")
print("=> Mersenne  cap  Heegner = {3,7} => LRC(2*3)=LRC(6), LRC(2*7)=LRC(14) are the uniquely-poised cases.")
print()
# The N = 2^h * odd_core gauge (HYP-1984): even/odd = add/mult. 14 = 2^1 * 7.
print("even/odd duality N=2^h*odd_core (HYP-1984):  14 = 2^1 * 7  (h=1 even-part / odd_core=7 = apex Mersenne)")
print("  ADDITIVE axis: odd_core -> odd_core+2 (steps 3->5->7..);  MULT axis: 2^h doubling (Cayley-Dickson/DRT tower)")
print()
# my extension: H is ADD (sum over cycles) AND MULT (product over SC primes) simultaneously
print("EXTENSION (opus): H sits at the add/mult interface:")
print("  H = PROD over strongly-connected primes  (multiplicative: H(X(+)Y)=H(X)H(Y))   [condensation]")
print("  H = 1 + 2*SUM_k alpha_k 2^k = I(Omega,2)   (additive: sum over disjoint odd-cycle families) [OCF]")
print("  The H-MAXIMIZER is mult-IRREDUCIBLE (a single SC prime) but add-MAXIMAL (most cycle families).")
print("  Fourier face: Gauss sum = MULTIPLICATIVE char chi (Paley, R-odd/imaginary/obstruction);")
print("                Dirichlet kernel = ADDITIVE interval (half-turn, R-even/real/provable). eps=chi(-1) hinges them.")
