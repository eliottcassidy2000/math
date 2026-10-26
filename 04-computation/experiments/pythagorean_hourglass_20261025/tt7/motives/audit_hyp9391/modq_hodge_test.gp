\\ Audit of HYP-9391 (c): a Hodge-type test that needs no exact separation.
\\ If a new motive had no (5,0) Hodge part (e.g. q times a weight-3 motive, or q^2 times an abelian surface), its
\\ Frobenius traces would be divisible by q (Newton above Hodge), and then Y(q) mod q would equal the known part mod
\\ q, i.e. a FIXED small-integer combination of 1, a84, a86, a46, |b|^2 (the columns q, q^2, q^3, q a84 vanish mod q),
\\ optionally also (2/q), (3/q) and the class numbers h(-q), (2/q)h(-q), h(-8q).  We look for short integer vectors v = (e, c_1..c_k), e != 0, with
\\     e Y(q) - sum c_j col_j(q) == 0 (mod q)
\\ at 40 training primes (CRT + LLL; under the null the shortest such vector has size ~ (prod q)^(1/(k+1)), ~10^15+),
\\ and validate every short vector on the remaining primes.  Positive controls: the twisted fourth moment
\\ M4t = -2q^2 - 3q - 1 - 3q a84 - a86 - a46 (expect v = (1, -1, 0, -1, -1, 0, ...)) and a synthetic
\\ T = 7 - 3 a84 + 2|b|^2 + q Y09.  Data file (written by the audit) defines QS, Y02, Y05, Y09, Y18, M4T, A84, A86, A46, B2.
default(parisize, 400000000);
read(DATAFILE);
n = #QS;
chi2 = vector(n, i, kronecker(2, QS[i]));
chi3 = vector(n, i, kronecker(3, QS[i]));
one = vector(n, i, 1);
NT = 40;
test(name, T, cols, colnames) =
{
  my(k = #cols, Q = prod(i = 1, NT, QS[i]), W = 2^400, R, M, L, best = [], vals, ok, nv);
  \\ CRT lift of each column over the training primes
  R = vector(k + 1, j, my(c = if (j == 1, T, cols[j - 1])); lift(chinese(vector(NT, i, Mod(c[i], QS[i])))));
  M = matrix(k + 2, k + 2);
  for (j = 1, k + 1, M[j, j] = 1; M[k + 2, j] = W * R[j] * if (j == 1, 1, -1));
  M[k + 2, k + 2] = W * Q;
  L = qflll(M, 0);
  L = M * L;
  nv = 0;
  for (c = 1, k + 2,
    my(v = L[, c]);
    if (v[k + 2] != 0, next);
    if (v[1] == 0, next);
    nv++;
    \\ validate on the other primes
    ok = 1;
    for (i = NT + 1, n,
      my(s = v[1] * T[i] - sum(j = 1, k, v[j + 1] * cols[j][i]));
      if (s % QS[i] != 0, ok = 0; break));
    best = concat(best, [[vecmax(abs(v[1..k+1])), ok, v[1..k+1]~]]));
  best = vecsort(best, 1);
  print(name, ": ", nv, " kernel vectors with e != 0; shortest sup-norm = ", if (#best, best[1][1], "none"),
        "; validated on ", n - NT, " primes: ", if (#best, best[1][2], "-"), "; vector (e, ", colnames, ") = ", if (#best, best[1][3], "-"));
}
cols5 = [one, A84, A86, A46, B2]; nm5 = "c1, ca84, ca86, ca46, cb2";
cols7 = [one, A84, A86, A46, B2, chi2, chi3]; nm7 = "c1, ca84, ca86, ca46, cb2, c(2/q), c(3/q)";
print("training primes: ", QS[1], "..", QS[NT], "; validation primes: ", QS[NT + 1], "..", QS[n], " (", n - NT, ")");
print("-- positive controls --");
test("M4t (twisted 4th moment)", M4T, cols5, nm5);
SYN = vector(n, i, 7 - 3 * A84[i] + 2 * B2[i] + QS[i] * Y09[i]);
test("synthetic 7 - 3 a84 + 2|b|^2 + q Y09", SYN, cols5, nm5);
print("-- the triple products (5 columns) --");
test("Y02", Y02, cols5, nm5);
test("Y05", Y05, cols5, nm5);
test("Y09", Y09, cols5, nm5);
test("Y18", Y18, cols5, nm5);
test("Y05 + Y18", Y05 + Y18, cols5, nm5);
print("-- the triple products (7 columns, adding (2/q), (3/q)) --");
test("Y02", Y02, cols7, nm7);
test("Y05", Y05, cols7, nm7);
test("Y09", Y09, cols7, nm7);
test("Y18", Y18, cols7, nm7);
test("Y05 + Y18", Y05 + Y18, cols7, nm7);
print("-- the triple products (10 columns, adding class numbers h(-q), (2/q) h(-q), h(-8q)) --");
hq = vector(n, i, qfbclassno(-QS[i])); h8q = vector(n, i, qfbclassno(-8 * QS[i]));
hq2 = vector(n, i, chi2[i] * hq[i]);
cols10 = [one, A84, A86, A46, B2, chi2, chi3, hq, hq2, h8q]; nm10 = "c1, ca84, ca86, ca46, cb2, c(2/q), c(3/q), ch(-q), c(2/q)h(-q), ch(-8q)";
test("Y02", Y02, cols10, nm10);
test("Y05", Y05, cols10, nm10);
test("Y09", Y09, cols10, nm10);
test("Y18", Y18, cols10, nm10);
test("Y05 + Y18", Y05 + Y18, cols10, nm10);
HQSYN = vector(n, i, 3 * hq[i] - 5 + 2 * A84[i] + QS[i]^2 * Y02[i]);
test("synthetic 3 h(-q) - 5 + 2 a84 + q^2 Y02 (positive control)", HQSYN, cols10, nm10);
print("-- residues: share of primes with Y = 0 mod q --");
{
foreach([["Y02", Y02], ["Y05", Y05], ["Y09", Y09], ["Y18", Y18]], P,
  print("  ", P[1], ": ", sum(i = 1, n, P[2][i] % QS[i] == 0), " of ", n));
}
quit;
