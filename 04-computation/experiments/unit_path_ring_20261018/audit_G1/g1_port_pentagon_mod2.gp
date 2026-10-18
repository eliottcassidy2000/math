\\ audit G1: what the golden port pentagon {0, 1, 1 + z5, 1 + z5 + z5^2, -z5^4} becomes at the prime above 2 (residue F_16),
\\ versus the residue pentagon mu_5 (THM-4701 remark "Six-colour programme", results note section 3).
P5 = polcyclo(5); nf5 = nfinit(P5); z = Mod(x, P5);
ports = [0, 1, 1 + z, 1 + z + z^2, -z^4];
cj5(y) = subst(lift(y), x, Mod(x^4, P5));
print("port pentagon: consecutive differences are unit vectors: ", prod(k = 1, 5, my(d = ports[k % 5 + 1] - ports[k]); d*cj5(d) == 1));
print("diagonal squared lengths equal phi^2: ", prod(k = 1, 5, my(d = ports[(k + 1) % 5 + 1] - ports[k]); d*cj5(d) == (1 + z + z^4)^2));
pr = idealprimedec(nf5, 2)[1]; m = nfmodprinit(nf5, pr);
red(y) = nfmodpr(nf5, lift(y + 0*z), m);
rp = vector(5, k, red(ports[k]));
adj(a, b) = (a != b) && (a - b)^5 == 1;      \\ G_4: x ~ y iff N(x - y) = (x - y)^5 = 1
print("residues of the port pentagon in F_16: ", rp);
print("adjacent pairs among them in G_4 (Clebsch): ", sum(a = 1, 5, sum(b = a + 1, 5, adj(rp[a], rp[b]))), " of 10  (an independent set would give 0)");
print("is the residue set equal to mu_5? ", Set(rp) == Set(vector(5, k, red(z^(k-1)))));
mu = vector(5, k, red(z^(k-1)));
print("mu_5 = residues of the side vectors z^k: adjacent pairs inside mu_5: ", sum(a = 1, 5, sum(b = a + 1, 5, adj(mu[a], mu[b]))), " (independent)");
print("diagonal vectors of the port pentagon are the consecutive two-sums z^k + z^(k+1): ", prod(k = 1, 5, ports[(k + 1) % 5 + 1] - ports[k] == z^(k-1) + z^k));
quit;
