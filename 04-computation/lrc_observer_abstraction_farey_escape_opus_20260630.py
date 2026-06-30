"""
The OBSERVER abstraction. The observer sits at the origin; a runner of speed s RETURNS to it (blocks it) at
t=1/s and every k/s. A covering set blocks the observer at every harmonic point 1/q (q=2..n). The lonely
time = the observer's ESCAPE. CLAIM: the covering-min is the observer's escape between its TWO TIGHTEST
blockings 1/n and 1/(n-1) -- and the two families are the FAREY MEDIANT (drop-2) and the CONTINUED-FRACTION
convergent (construction), with the blockings (n-1, n) as the partial quotients.
"""
from fractions import Fraction
def cf(fr):  # continued fraction of a Fraction
    a=[]; p,q=fr.numerator,fr.denominator
    while q:
        a.append(p//q); p,q=q,p-(p//q)*q
    return a
print("OBSERVER's escape between its two tightest blockings 1/n and 1/(n-1):")
print(f"{'n':>3} {'1/n':>8} {'1/(n-1)':>9} {'mediant 2/(2n-1)':>16} {'CF n/Phi6 [0;n-1,n]':>20} {'cov-min':>10}")
for n in range(4,15):
    one_n=Fraction(1,n); one_nm1=Fraction(1,n-1)
    mediant=Fraction(1+1,(n-1)+n)  # Farey mediant of 1/(n-1),1/n = 2/(2n-1)
    Phi6=n*n-n+1; constr=Fraction(n,Phi6)
    c=cf(constr)
    covmin = mediant if n<=6 else constr
    fam = "drop-2(mediant)" if n<=6 else "constr(convergent)"
    print(f"{n:>3} {str(one_n):>8} {str(one_nm1):>9} {str(mediant):>16} {str(constr)+' '+str(c):>20} {str(covmin):>8} {fam}")
print()
print("=> BOTH covering-min families live STRICTLY BETWEEN the observer's two tightest blockings 1/n, 1/(n-1):")
for n in [6,7,14]:
    print(f"   n={n}: 1/n={1/n:.5f} < drop-2 2/(2n-1)={2/(2*n-1):.5f} < constr n/Phi6={n/(n*n-n+1):.5f} < 1/(n-1)={1/(n-1):.5f}")
print()
print("THE OBSERVER READING:")
print(" * the construction's escape n/Phi6 = [0; n-1, n] -- the CONTINUED FRACTION whose partial quotients")
print("   ARE the observer's two tightest blockings (q=n-1 and q=n). The escape IS the blockings, nested.")
print(" * the drop-2 escape 2/(2n-1) = the FAREY MEDIANT of 1/(n-1) and 1/n -- the simplest point between.")
print(" * the covering-min = which escape a covering set can REALIZE: the mediant (drop-2) for n<=6, the")
print("   convergent (construction) for n>=7. The transition is a realizability question about the observer.")
print(" * the empty tooth = the observer (deleted origin); the lonely set = where the observer is unblocked.")
print(" * margin: escape - 1/n = (n-1)/(n.Phi6) ~ 1/n^2 -- the observer escapes 1/n only by a Farey hair.")
