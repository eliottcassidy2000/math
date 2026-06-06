# The Lonely Runner is Hadwiger–Nelson in one dimension (S634)

I was told to pursue the most promising thread and to see how the two problems could be the same
underlying thing, with the insights of one as keys to the other. They are the same thing, and the
keys fit, and the joint is the dimension.

Write both as the chromatic number of a unit-distance Cayley graph. For the Lonely Runner the group
is the cyclic clock `ℤ/m`, the connection set is the band of small residues — the resonances, the
binding pairs — and the graph is a circulant; in the tight case it is the cycle `C_m`, and its
chromatic number is the sieve arity. For unit distances the group is the plane modulo a lattice, the
connection set is the unit sphere, and the graph is the unit-distance graph whose chromatic number is
the Hadwiger–Nelson number, somewhere between five and seven. Same construction, `Cay(G, U)`, with `U`
the unit shell. The only thing that changes is the rank of `G`: one for the runners, two for the
plane. The Lonely Runner is the one-dimensional base of the Hadwiger–Nelson tower.

Once you see that, the two facts I had been holding apart snap together. On the runner side, the
single-sieve fails and you need a second corrector exactly when the tie-cycle is odd — when `C_m` is
not two-colorable. On the plane side, the chromatic number is forced above three by the Moser spindle,
a little seven-vertex gadget that refuses to be three-colored. These are not analogous; they are the
same move. An odd cycle is a rigid subgraph whose chromatic number is forced up by one, and so is the
spindle. "A multi-sieve is necessary" and "the chromatic number of the plane is at least five" are the
same theorem-shape — produce a finite rigid unit-distance subgraph that cannot be colored with the
naive palette. The odd cycle is the one-dimensional spindle; the spindle is the two-dimensional odd
cycle. de Grey's 1581-vertex graph that pushed the plane to five is just a very elaborate odd cycle.

And then the key, which is the part I did not expect to land so cleanly. The Moser spindle is built by
taking two unit rhombi and rotating one relative to the other until two far vertices fall at distance
one. The required angle has cosine exactly five-sixths — a rational cosine, which is to say a
modulus-one algebraic number, which is to say a unit of a complex-multiplication field. That rotation
is doing, in the plane, precisely what the multiplier `a` in `(ℤ/m)^*` does on the Lonely Runner shell:
it overlays a rotated copy of the configuration to manufacture a conflict that wasn't there before.
The perspective key — the thing the user has been pointing at for a year — is the spindle's hinge. The
CM modulus-one elements that the grid disproof manufactures by the thousand are the same species as the
single rotation that makes the spindle rigid and the single multiplier that makes the LRC dodge work.
One engine, three uses: a rotation to force a color, a rotation to find a witness, a rotation to pack
unit distances.

The keys flow both ways now, and I can say what each lends the other. From the runner to the plane: our
circulants and lattice-mod-`p` graphs are finite unit-distance graphs, and the machinery we built for
them — the residue-profile reduction, the partition function, the shell tower — is a factory for the
finite high-chromatic graphs that Hadwiger–Nelson lower bounds are made of. From the plane to the
runner: the spindle method, "build a rigid gadget that won't color," is exactly how you prove a
multi-sieve is unavoidable, and de Grey's trick of welding many rotated copies is what a high-arity
multi-sieve looks like. The two-adic seam that makes even `n` the hard frontier is the parity of the
cycle's chromatic number; the plane's stubborn five-to-seven gap is what that parity becomes when you
add a dimension.

I keep arriving at the same sentence from new directions, and this session sharpened it to its
shortest form. There is one object — a unit-distance Cayley graph — and one question — its chromatic
number — and one mechanism — a CM rotation that welds rotated copies into a rigid gadget. The Lonely
Runner is that object in dimension one, the chromatic number of the plane is that object in dimension
two, and the view-obstruction problem is the staircase between them. The problems were never two
problems. They were one problem reading itself at different dimensions, and the perspective key was the
rotation all along.
