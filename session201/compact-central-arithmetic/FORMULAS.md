# The twelve formulas in ordinary arithmetic

This is an expansion of `MENUS.json`, not a new hypothesis. Put
A=g(u), B=g(v), p=u mod2, q=v mod2, with u,v>=2.
Both index orders are allowed; A,B,p,q are not free parameters.

The OPEN occupancy lemma asks whether every odd n>2134117 equals at
least one of the following twelve expressions.

- empty: A + B + p*q - p - q + u*v - u - v + 1
- J: A*B - p*q + u*v
- JJ: A*B - A*q + A*v - B*p + B*u + p*q + u*v
- JJJ: A*B + 2*A*v + 2*B*u - p*v - q*u + 2*u*v
- JOZJ: 3*A*B - A*q + 3*A*v - B*p + 3*B*u + 4*u*v
- JJZJ: 3*A*B - A*q + 3*A*v - B*p + 3*B*u + p*q - p*v - q*u + 3*u*v
- JJJJ: 3*A*B - A*q + 3*A*v - B*p + 3*B*u - p*v - q*u + 4*u*v
- JOOZJ: 5*A*B - A*q + 5*A*v - B*p + 7*B*u + 8*u*v
- JJZZJ: 5*A*B - A*q + 5*A*v - B*p + 5*B*u + p*q - p*v - q*u + 5*u*v
- JOZJJ: 5*A*B - A*q + 7*A*v - 2*B*p + 6*B*u + p*q - 3*p*v - q*u + 7*u*v
- JJZJJ: 5*A*B - A*q + 5*A*v - 2*B*p + 6*B*u - 2*p*v - 2*q*u + 6*u*v
- JJJJJ: 5*A*B - 2*A*q + 6*A*v - 2*B*p + 6*B*u + p*q - 2*p*v - 2*q*u + 7*u*v
