\\ audit G1 common setup: F = Q(zeta15, sqrt-11, w) in PARI, from g1_buildF.gp (integral basis verified there).
V = read("g1_F.bin"); T2 = V[1]; zF = V[2]; jF = V[3]; wF = V[4]; cX = V[5]; nf = V[6];
one = Mod(1, T2);
cj(y) = subst(lift(y), x, cX);                 \\ complex conjugation
Nm(y) = y * cj(y);                              \\ relative norm to F+
isU(y) = (Nm(y) == 1);                          \\ unit vector
z3 = zF^5; z5 = zF^3;
s3 = 2*z3 + 1;  s5 = 2*(z5 + z5^4) + 1;  s15 = s3 * s5;
if (s3^2 != -3 || s5^2 != 5 || s15^2 != -15 || jF^2 != -11, error("square roots"));
\\ cubic recipe (golden_contact_circuit_20261017i.md section 4), rebuilt from the formulas in that note
dd = 1 + z5; aa = cj(dd);
hh = wF + 3; uu = (hh - jF)/(hh + jF);
m2 = -(20 + 11*jF)/19; m1 = (20 - 11*jF)/19;
cc = (aa*m2 + 1 - aa^2*m1 + aa^3)/(Nm(dd) - 1);
bb = cc*dd + aa*m1 - aa^2;
ZZ = dd + uu; YY = -ZZ/cj(ZZ);
\\ recipe sets: ports, spindle M, C = M + first five ports, A = ports + ports + M
etM = 1 + z3; rho = (5 + jF)/6;
Mset = [0*one, one, etM, 1 + etM, rho, rho*etM, rho*(1 + etM)];
ports = [0*one, one, 1 + z5, 1 + z5 + z5^2, -z5^4, (3 - jF*s5)/8, 3/4*one, (3 + jF*s5)/8];
uniqv(L) = my(out = List(), seen = Map()); for (i = 1, #L, my(k = lift(L[i])); if (!mapisdefined(seen, k), mapput(seen, k, 1); listput(out, L[i]))); Vec(out);
Cset = uniqv(concat(Mset, ports[1..5]));
Aset = uniqv(concat(ports, concat(vector(#ports, i, vector(#Mset, k, ports[i] + Mset[k])))));
\\ complex embeddings of F: the 48 roots of T2
emb = polroots(T2);
ev(y, r) = subst(lift(y), x, r);
house(y) = vecmax(vector(#emb, k, abs(ev(y, emb[k]))));
\\ square test in F (complete): roots of Y^2 - y over F
YV = varhigher("YV");
sqrtF(y) = my(r = nfroots(nf, YV^2 - lift(y))); if (#r, Mod(r[1], T2), 0);
