# Exploratory paths, separate from the reviewed core

These are author explorations and reproducibility records. They are not part
of the independent audit of `../PROOF.md` unless that audit explicitly says
otherwise. No computation here proves an all-target occupancy assertion.

The scripts use their own local `forms.json`, which is the older forty-template
family in `../dependencies/TWO_INDEX_40.md`, and a separate generated central
matrix header. They do not modify the repository checkout or the parent
packet's header. Compiled binaries and caches are intentionally omitted.

## Fixed smaller-index hostile controls

`fixed_left.cpp` inventories the older forty forms with the smaller index
bounded by B. Its proof-backed domain is 3uv+1<=N and 2<=u<=v. At N=1,000,001,
B=32 left only the known twenty small gaps. That tempting fixed-bank
strengthening fails: target 1348681 first uses smaller index41. The analogous
B=64 claim fails at6888203, whose least smaller index is86.

`find_target.cpp` gives the first direct-enumeration witness.
`verify_fixed_bank.py` independently uses the exact congruence progression
inverse and verifies both least smaller indices, not merely the two positive
representations. Neither result says that every finite fixed bank fails to
be eventually cofinite: a single finite counterexample cannot establish that.

The full older family through10,000,001, with B=isqrt((N-1)/3), still leaves
only its known twenty gaps. That is exploratory finite support for the older
lemma, not a proof of it and not new C1 coverage.

## Restricted binary descent

`digit_descent.py` tests coefficient-wise identities within the forty-form
bank after applying u->2u or2u-1 to one or both indices, against2 or4 times
another form plus a constant. It also fixes either parent parity before
comparing coefficients. `DIGIT-ANSATZ.txt` records no identities in this
specified ansatz. This does not rule out longer recursions, additional forms,
index-dependent transformations, or identities valid only on a more
restricted set of actual indices.

## Central-word search and its larger hostile probe

`central_prepare.py` generates exact matrices through central length nine;
each inventory command specifies the actual maximum length to include.
`central_inventory.cpp` is the exploratory predecessor of the frozen C++
program. `EXPLORATORY-100m.txt` is its author-only scan at100,000,001 and
maximum central length8. It retains exactly the same29 gaps as the core's
smaller universe. This100m record has NOT been
independently replayed and is not used to establish a universal statement.

`central_roots.cpp` additionally caps the number of J letters. At cap10,000,001
and central length8, J caps2,3,4,5 leave respectively195,47,34,32 gaps, with
largest gaps9872947,27781,17383,17383. Thus few root changes alone are a
different restriction from a short central word.

`central_fixed.cpp` caps the smaller unit index in the central length8 family.
At10,000,001, index caps8,16,32,64 leave269,48,32,29 gaps. The three additional
gaps for cap32 are7182463,7799923,8231473. Again a bounded index is not justified
by a finite check at a smaller ceiling.

`central_reduce.py` records the finite exact JZ=OJ boundary-absorption
normalization. The core verifier independently regenerates this menu and
checks the matrix equalities; its audit additionally checks the actual
directed changes of indices. Do not use u<=v on this reduced menu without
restoring a proved transpose symmetry.

## Other comparisons

`form_groups.cpp` separates the old C2/C3 template subfamilies.
`unit_pairs.cpp` inserts zero, one, or two adjacent central J moves between
the unit sources. Its two-J family does not allow arbitrary O/Z blocks
between the J moves and must not be called the full two-root-change class.

Typical reproduction commands, from this directory:

    g++ -O3 -std=c++17 fixed_left.cpp -o /tmp/fixed-left
    /tmp/fixed-left 10000001 1825
    python verify_fixed_bank.py
    python digit_descent.py
    python central_prepare.py
    g++ -O3 -std=c++17 central_inventory.cpp -o /tmp/central-exploration
    /tmp/central-exploration 100000001 8

The large exploratory scan can take several minutes and roughly500MB.
The smaller core scan is the preferred routine reproducibility command.
