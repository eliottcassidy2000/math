\\ The weight-3 level-32 newform orbit with character (-4/.) that evaluates the Paley K6-P6 terminal.
default(parisizemax, 10^9);   \\ mfcoefs to 30000 needs more than gp's default 8 MB stack (audit C)
mf = mfinit([32, 3, -4], 0);
L = mfeigenbasis(mf);
print("newform orbits in S_3^new(32, (-4/.)): ", #L, "; dimension ", mfdim(mf));
{
for(i = 1, #L,
  f = L[i]; co = mfcoefs(f, 30);
  K = if(type(co[2]) == "t_POLMOD", co[2].mod, x);
  print("orbit ", i, ": coefficient field polynomial ", K, ", first coefficients ", co[2..16]);
  print("   CM? ", mfisCM(f))
);
}
\\ exact check of the conjectured identity at all odd primes up to B, using Tr(a_p^2)/2
B = 30000;
co = mfcoefs(L[1], B);
P = primes([5, B]);
v = vector(#P, j, my(c = co[P[j] + 1]); if(type(c) == "t_POLMOD", trace(c^2)/2, c^2));
write("g_trace_a2_half.txt", v);
print("written Tr(a_p^2)/2 for ", #P, " primes 5 <= p <= ", B);
quit
