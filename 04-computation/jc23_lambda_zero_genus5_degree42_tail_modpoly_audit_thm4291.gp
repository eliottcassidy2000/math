\\ Independent PARI/GP modular-polynomial audit for THM-4291.

P = polmodular(7);
special = subst(P, x, 0);
A = 34848505552896000;
B = 11356800389480448000000;
expected = y^2 * (y^2 + A*y + B)^3;

if (special != expected, error("Phi_7(0,y) factor identity failed"));
if (A^2 - 4*B != 4 * 21 * 3802283679744000^2, error("nonzero-root discriminant identity failed"));

print("THM4291_PARI_MODPOLY_AUDIT_V1");
print("PHI7_AT_J0 ", factor(special));
print("NONZERO_ROOTS -17424252776448000 +/- 3802283679744000*sqrt(21)");
print("VERDICT PASS PARI_EXACT_MODULAR_POLYNOMIAL");
