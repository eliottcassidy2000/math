\\ Wilson quotients w_p = ((p-1)! + 1)/p mod p for all primes 3 <= p < N (exact, Mod arithmetic mod p^2).
\\ Prints the primes with w_p == 0 (Wilson primes) and w_p == -1 (Giuga sum S(p) == -1 mod p^2, by Lerch).
N = 200000;
W0 = List(); Wm1 = List(); Wp1 = List();
forprime(p = 3, N, f = Mod(1, p^2); for(k = 2, p - 1, f *= k); w = (lift(f) + 1) / p % p; \
  if(w == 0, listput(W0, p)); if(w == p - 1, listput(Wm1, p)); if(w == 1, listput(Wp1, p)));
print("N = ", N);
print("w_p == 0  (Wilson primes):          ", Vec(W0));
print("w_p == -1 (S(p) == -1 mod p^2):      ", Vec(Wm1));
print("w_p == +1:                           ", Vec(Wp1));
quit
