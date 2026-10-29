\\ full Smith form of the audit's own M_3(P7) (210 x 210), to compare with THREE_RUN_GLUE (R9)
default(parisizemax, 4000000000)
read("mats/P7_k3_M.gp");
t = getabstime();
D = matsnf(M);
print("matsnf of M_3(P7) in ", getabstime() - t, " ms");
nu = select(x -> x != 1, D);
print("number of unit elementary divisors: ", #D - #nu);
\\ group equal divisors
v = vecsort(nu);
out = List(); i = 1;
while (i <= #v, j = i; while (j < #v && v[j+1] == v[i], j++); listput(out, [v[i], j - i + 1]); i = j + 1);
print("non-unit elementary divisors [value, multiplicity]: ", Vec(out));
print("order of the cokernel: ", factor(abs(matdet(M))));
quit;
