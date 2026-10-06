\\ Complement to square11_octic_langlands_collatz_20261006.md (opus S16), whose open thread 4 asks for a
\\ polredabs-type small model of K = Q(T) and its class group.  mac-mini six-seven session, 2026-10-06.
\\ Run: gp -q sixseven_20261006_octic_class_group.gp < /dev/null
default(parisize, 400000000);
Ps = x^8-20*x^7+178*x^6-842*x^5+1923*x^4-496*x^3-6754*x^2+12420*x-6865;
Pu = 5*x^8-10*x^7-2*x^6+14*x^5+12*x^4-6*x^3+2*x^2+2*x-1;
R = polresultant(subst(Pu,x,u), s*(1+2*u-u^2)-(6*u+4), u);
print("Res_u(Pu(u), s(1+2u-u^2)-(6u+4)) = ", pollead(R), " * Ps(s):  ", R == pollead(R)*subst(Ps,x,s));
f = polredabs(Ps);
print("polredabs(Ps) = ", f, "    polredabs(Pu) = ", polredabs(Pu), "    same field: ", nfisisom(Ps, Pu) != 0);
K = bnfinit(f, 1);
print("signature ", K.sign, "   d_K = ", K.disc, " = ", factor(K.disc));
print("proper subfields (other than Q, K): ", #nfsubfields(f) - 2);
print("class group ", K.cyc, "  (class number ", K.no, ")   regulator ", K.reg, "   torsion units ", K.tu[1], "   unit rank ", #K.fu);
print("bnfcertify (unconditional, removes GRH): ", bnfcertify(K));
print("ramified primes: ", factor(K.disc)[,1]~, "   splitting of small primes (e,f): ");
foreach([2,3,5,7,11,13,17,19,23,29,31], p, print("   p = ", p, ": ", apply(P -> [P.e, P.f], idealprimedec(K, p))));
quit
