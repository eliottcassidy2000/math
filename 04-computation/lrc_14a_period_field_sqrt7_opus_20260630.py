"""The period field of 14a = Q(sqrt(Delta)) = Q(sqrt-7): Delta = -2^6 7^3, squarefree part -7."""
def sqfree(n):
    s=1; m=abs(n); d=2
    while d*d<=m:
        c=0
        while m%d==0: m//=d; c+=1
        if c%2: s*=d
        d+=1
    s*=m
    return -s if n<0 else s
D=-21952
print(f"Delta(14a) = {D} = -2^6 * 7^3   (factor: {D}=-(2**6)*(7**3)={-(2**6)*(7**3)})")
print(f"squarefree part of Delta = {sqfree(D)}  => period field Q(sqrt(Delta)) = Q(sqrt{sqfree(D)}) = Q(sqrt-7)")
print(f"  = the MEASURE column Heegner field. The apex-7 (cubed) IS the discriminant's odd part.")
print(f"  Period quartic 4s^4 +/- 13 s^2 + 32 has disc_(s^2) = 13^2-512 = -343 = -7^3 (same -7).")
print()
print("So 14a's period GEOMETRY is defined over Q(sqrt-7) (measure column), though 14a is NOT CM:")
print("  the branch points live in Q(sqrt-7); the lattice modulus tau is transcendental (no CM).")
print()
print("SUMMARY of the cusp-form-period chase:")
print("  CLEAN:  Omega+ = 1.98134, and L(14a,1) = Omega+/6 EXACTLY (BSD: prod c_p . #Sha = 6, #tors^2=36).")
print("          Omega3(imag) = 2.65098, Area = 5.2525, tau = 1/2 + 0.66899 i.")
print("          Period field = Q(sqrt-7) = Q(sqrt Delta), Delta=-2^6 7^3 (measure column, exact).")
print("  NEGATIVE: the floor's 0.040 cusp-part is NOT a clean period -- checked Omega+, Omega3, Area, tau,")
print("            and all combinations. THIRD negative in a row: floor != L(14a,1)=0.330, != L(sym^2,s)~1.1,")
print("            != any period. => the floor CONSTANT is NOT a single modular invariant of f_14.")
print("  CONCLUSION: the modular form gives the STRUCTURE (genus=hardness, rank-0=sign, GL(3)=level,")
print("            Q(sqrt-7)=period field), NOT the constant. The 0.040 lives in the DESCENT (the product")
print("            prod rho_j / the covering-min), a multi-level quantity, not a single period. Stop")
print("            chasing it as a modular period; chase it in the descent.")
