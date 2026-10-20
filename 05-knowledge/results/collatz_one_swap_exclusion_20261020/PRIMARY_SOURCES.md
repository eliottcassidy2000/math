# Imported analytic input

The only external theorem used by the final closure is the explicit real Matveev bound in Bugeaud, Mignotte and Siksek, *Classical and modular approaches to exponential Diophantine equations I. Fibonacci and Lucas perfect powers*, Theorem 9.4, printed page 16 (PDF page 16 in one-based numbering):
https://arxiv.org/pdf/math/0403046

The authors state the real-field bound with factor 1.4*30^(n+3)*n^4.5*D^2*(1+logD)*(1+logB), and derive it from Corollary 2.3 of Matveev. The current application uses positive rationals 2 and 3, integer coefficients q and -p, D=1, n=2, and height upper bounds 0.7 and 1.1. This exact form was source-checked on 20 October 2026.

Original theorem: E. M. Matveev, *An explicit lower bound for a homogeneous rational linear form in the logarithms of algebraic numbers. II*, Izvestiya: Mathematics 64(6) (2000), 1217–1269, DOI 10.1070/IM2000v064n06ABEH000314.
https://www.mathnet.ru/eng/im314

The final proof does not use any remembered Rhin constant, Hercher cycle bound, Barina verification limit, or conjectured irrationality exponent. Earlier exploratory certificates in this staging directory that mention those routes are not part of the final theorem or publication inventory.
