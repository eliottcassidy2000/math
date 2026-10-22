\\ Closed form for Wormald's blow-up constant rho_T(T*) = sqrt(d/H):  rho = chi_S(d) / |chi_S(i d)|,
\\ chi_S(t) = det(t I - S), S = A - A^T (skew adjacency) of the regular base tournament T* on d vertices.
\\ Compared against the exact rationals in Wormald, "Tournaments with many Hamilton cycles", table p. 5.
circ(d, D) = matrix(d, d, i, j, if(i == j, 0, if(setsearch(Set(D), (j - i) % d), 1, 0)));
rho(A) = { my(d = #A, S = A - A~, P = charpoly(S, t)); subst(P, t, d) / abs(subst(P, t, I * d)); }
rot(d) = circ(d, vector((d - 1) / 2, r, r));
lex3(B) = { my(m = #B, d = 3 * m); matrix(d, d, i, j, my(bi = (i - 1) \ m, bj = (j - 1) \ m); if(bi == bj, B[(i - 1) % m + 1, (j - 1) % m + 1], if((bj - bi) % 3 == 1, 1, 0))); }
C3 = rot(3);
print("R3  ", rho(rot(3)), "   Wormald 2");
print("R5  ", rho(rot(5)), "   Wormald 44/19");
print("R7  ", rho(rot(7)), "   Wormald 3032/1231");
print("R9  ", rho(rot(9)), "   Wormald 437584/171769");
print("R'9 = C3[C3]  ", rho(lex3(C3)), "   Wormald 5488/2197");
print("R11 ", rho(rot(11)), "   Wormald 108859552/41833969");
print("R13 ", rho(rot(13)), "   Wormald 41637757888/15769379963");
print("9-circulants: ", [rho(circ(9, [1,3,4,5]))*1., rho(circ(9, [1,2,4,5]))*1., rho(circ(9, [1,2,3,5]))*1., rho(circ(9, [2,3,4,5]))*1., rho(circ(9, [1,3,5,6]))*1., rho(circ(9, [1,4,5,6]))*1., rho(circ(9, [2,3,5,6]))*1.]);
print("   Wormald 1.4582 1.6811 2.4557 1.5618 1.5529 1.9325 1.8117");
print("R21 ", rho(rot(21))*1., "  Wormald 2.7215188");
print("R41 ", rho(rot(41))*1., "  Wormald 2.7867038");
print("R81 ", rho(rot(81))*1., "  Wormald 2.820801899");
print("R101", rho(rot(101))*1., "  Wormald 2.827747088");
print("P7  ", rho(circ(7, [1,2,4])), "  P11 ", rho(circ(11, [1,3,4,5,9])), "  (Paley: ((d+1)/(d-1))^((d-1)/2) = ", (8/6)^3, ", ", (12/10)^5, ")");
print("cosh(1)/cos(1) = ", cosh(1)/cos(1));
quit
