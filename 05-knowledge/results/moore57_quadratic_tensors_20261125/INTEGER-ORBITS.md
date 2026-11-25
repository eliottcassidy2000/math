# What the first seven integer parameters see

Take \(f_z(x)=x^2-z\). The Mandelbrot parameter in the usual convention is \(c=-z\). It is the orbit starting at the critical point \(0\) that determines membership. For real \(z\ge0\) that orbit is bounded exactly for \(0\le z\le2\). Indeed \(b=(1+\sqrt{1+4z})/2\) is fixed. If \(z\le2\), then \(z\le b\) and \([-b,b]\) maps into itself. If \(z>2\), the second critical iterate \(z^2-z>z>b\) and the subsequent iterates increase to infinity. Individual periodic orbits can still exist at parameters outside this interval.

| \(z\) | All integer cycles | Orbit beginning at zero |
|---:|---|---|
| 0 | \(0\) and \(1\) fixed | \(0\) fixed; also \(-1\to1\) |
| 1 | \(0\leftrightarrow-1\) | on that cycle; \(1\to0\) |
| 2 | \(2\) and \(-1\) fixed | \(0\to-2\to2\); also \(1\to-1\) |
| 3 | \(1\leftrightarrow-2\) | \(0\to-3\to6\to33\to\cdots\) |
| 4, 5 | none | escapes |
| 6 | \(3\) and \(-2\) fixed | escapes |

For every integer \(z\), rational periodic points are integers: in lowest terms, a nonunit denominator squares at each step. A fixed integer \(a\) requires \(z=a(a-1)\); its partner is \(1-a\). A genuine two-cycle is \(a\leftrightarrow-1-a\) with \(z=a(a+1)+1\). These give the observed missing intervals in parameter space. They exhaust the integer cycles: on a nonfixed orbit, taking consecutive differences and multiplying \(f_z(x)-f_z(y)=(x-y)(x+y)\) round the orbit gives \(\prod_i(x_i+x_{i+1})=1\). Every factor is an integer unit, so all successive differences have the same positive absolute value. A finite arithmetic lattice with steps of one fixed length is a path and cannot contain a simple cycle of length greater than two.

At the rational parameter \(29/16\) the three-cycle is \(-7/4\to5/4\to-1/4\to-7/4\). Its adjacent sums \(-1/2,1,-2\) still multiply to one, illustrating exactly where the integer proof loses its force. The global period-four and period-five exclusions need the arithmetic of their respective dynatomic curves. Poonen's conditional polynomial portrait bound is nine rational points including infinity. The count of twelve portrait types, the twelve-point result for a different rational-map family, and the repository's twelve points across two parameters have distinct meanings. [Poonen](https://arxiv.org/abs/math/9512217) · [Flynn–Poonen–Schaefer](https://people.maths.ox.ac.uk/flynn/arts/art11.pdf) · [Portraits, the other twelve-point bound, and the Moore transport test](DYNAMICS.md)
