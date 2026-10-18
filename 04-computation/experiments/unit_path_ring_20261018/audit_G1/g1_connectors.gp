\\ audit G1: house table of THM-4700(5d) and an exact 1- and 2-step connector census between coefficient copies.
\\ run: gp -q -s 800000000 g1_connectors.gp
read("g1_common.gp");
default(realprecision, 50);
emb = polroots(T2);
nemb = #emb;
\\ embedding values of the generators, so that each element is evaluated once per embedding
EV(y) = vector(nemb, k, ev(y, emb[k]));
VC = vector(#Cset, i, EV(Cset[i])); VM = vector(#Mset, i, EV(Mset[i]));
res = vector(2);
{
for (ci = 1, 2,
  my(k = [bb, cc][ci], Vk = EV(k), tab = List(), small = List());
  for (i = 1, #Mset - 1, for (i2 = i + 1, #Mset,
    my(best = [10^9, 0, 0]);
    for (xi = 1, #Cset, for (yi = 1, #Cset,
      my(hv = vecmax(vector(nemb, e, abs(VC[yi][e] - VC[xi][e] + Vk[e]*(VM[i2][e] - VM[i][e])))));
      if (hv < best[1], best = [hv, xi - 1, yi - 1]);
      if (hv <= 2 + 10^-30, listput(small, [i - 1, i2 - 1, xi - 1, yi - 1, hv]))));
    listput(tab, [i - 1, i2 - 1, best[1], best[2], best[3]])));
  res[ci] = [Vec(tab), Vec(small)];
  my(hs = vector(#tab, r, tab[r][3]));
  print("k = ", ["b", "c"][ci], ": min over (x, y) of house(y - x + k (m' - m)) ranges ", vecmin(hs), " .. ", vecmax(hs),
        " ; ceil range ", ceil(vecmin(hs)), " .. ", ceil(vecmax(hs)));
  print("   per pair (m, m', house_min, x, y):");
  for (r = 1, #tab, print("   ", tab[r][1], " ", tab[r][2], " ", precision(tab[r][3], 12), " ", tab[r][4], " ", tab[r][5]));
  print("   displacements with house <= 2 (the only candidates for 1- or 2-step connectors): ", #small);
  for (r = 1, #small,
    my(s = small[r], t = Cset[s[4] + 1] - Cset[s[3] + 1] + k*(Mset[s[2] + 1] - Mset[s[1] + 1]), N = Nm(t), q, verdict);
    q = sqrtF(N*(4 - N)/3);
    verdict = if (N == 1, "ONE STEP", if (N == 4, "2U (two equal steps)", if (q != 0 && cj(q) == q, "TWO STEPS", "no 1- or 2-step walk")));
    print("     pair (", s[1], ", ", s[2], ") x = ", s[3], " y = ", s[4], " house ", precision(s[5], 8), " : N(t)(4 - N(t))/3 square in F: ", q != 0, " -> ", verdict)));
}
quit;
