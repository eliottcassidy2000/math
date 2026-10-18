\\ audit G1: which Q(sqrt(-D)) lie in F, for D = squarefree part of 4n^2 - 1, n = 2..40 (THM-4700(2) two-step criterion).
read("g1_common.gp");
sqf(m) = core(m);
Ds = Set(vector(39, k, sqf(4*(k+1)^2 - 1)));
Ds = setunion(Ds, [3, 11, 15, 55, 1, 5, 33, 165]);
{
foreach(Ds, D, my(r = nfroots(nf, YV^2 + D)); print("D = ", D, ": sqrt(-D) in F: ", #r > 0));
}
quit;
