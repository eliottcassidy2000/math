\\ rho_T of the ternary towers S_k = Cay(Z_(3^k), {x : unit part of x = 1 mod 3}) by the determinant formula,
\\ against the level product prod_(i<k) ((3^(2i+1)+1)/(3^(2i+1)-1))^(3^i)  (lexicographic spectrum law).
rho(A) = { my(d = #A, S = A - A~, P = charpoly(S, t)); subst(P, t, d) / abs(subst(P, t, I * d)); }
tern(k) = { my(n = 3^k); matrix(n, n, i, j, my(x = (j - i) % n); if(x == 0, 0, while(x % 3 == 0, x /= 3); if(x % 3 == 1, 1, 0))); }
lev(k) = prod(i = 0, k - 1, ((3^(2*i+1) + 1) / (3^(2*i+1) - 1))^(3^i));
for(k = 1, 3, my(r = rho(tern(k))); print("k=", k, "  rho_T(S_k) = ", r, "  level product = ", lev(k), "  equal: ", r == lev(k)));
quit
