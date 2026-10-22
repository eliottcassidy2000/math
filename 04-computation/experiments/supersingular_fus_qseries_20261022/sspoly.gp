\\ Shared: the supersingular polynomial ss_p(j) in F_p[x] (monic), p >= 5, via the Hasse
\\ invariant E_{p-1} == 1 (mod p): see the header of supersingular_census.gp (route A).
ssA(p) = {
  my(r = (p - 1) % 12, m = (p - 1) \ 12, de, ep, N, E4, E6, D, J, G, f = 0, c, i);
  de = (r == 4 || r == 10); ep = (r == 6 || r == 10);
  N = m + 3;
  E4 = Mod(1, p) * (1 + 240 * sum(n = 1, N, sigma(n, 3) * 'q^n) + O('q^(N + 1)));
  E6 = Mod(1, p) * (1 - 504 * sum(n = 1, N, sigma(n, 5) * 'q^n) + O('q^(N + 1)));
  D = Mod(1, p) * ('q * eta('q + O('q^(N + 1)))^24);
  J = E4^3 / D;
  G = 1 / (D^m * E4^de * E6^ep);
  forstep (i = m, 0, -1,
    c = polcoeff(G, -i, 'q);
    f += lift(c) * 'x^i;
    G -= c * J^i);
  if (valuation(G, 'q) < 1, error("route A peel failed at p = ", p));
  f = Mod(1, p) * f * 'x^de * ('x - 1728)^ep;
  f / pollead(f);
}
