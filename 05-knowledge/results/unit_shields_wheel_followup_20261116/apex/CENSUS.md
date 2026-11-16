# Spine 14: second apex with at most three spine contacts

This extends [Klein's November 15 sparse-contact census](https://github.com/eliottcassidy2000/math/commit/e8ac9fe3f37f20a99a750831a0c077b21d5f520a)
from spine 13 to 14. It is not the unrestricted census, which
the earlier [S93 packet](https://github.com/eliottcassidy2000/math/blob/d8ab558d0d58e11ce69a09e7155575cd3ce0f203/04-computation/experiments/s93_chanpak_wave_20261111/twoapex/census_L17_outcome.txt)
had already taken farther. We enumerate literal tournaments
with ordered transitive spine, a labelled apex A with arbitrary
out-neighbors, and a labelled B which beats at most three spine
sites. We allow both orientations of their mutual arc and
count presentations, not isomorphism types.

| Spine | Tapes per orientation | Distinct counts, both orientations | First odd miss other than 7,21 | Highest count |
|---|---:|---:|---:|---:|
| 13 | 3,096,576 | 401,873 | 256,521 | 1,904,361 |
| 14 | 7,700,480 | 1,138,644 | 729,421 | 5,889,173 |

Spine 14 with A→B alone already gives every odd below 388,061
except 7 and 21; with B→A alone the first further miss is 400,199.
The union covers every odd below 2^19 other than 7 and 21 and
misses 139 odd values below 2^20, including these two. With
|S_B|≤2 instead the first further miss for the combined
orientations is 164,333; below 2^18, exactly 70 attained counts
need the third contact. With |S_B|≤1 the first further miss is
7,477. These are all finite statements.

The old first miss 256,521 has a new spine-14 certificate:
A beats mask 16,147 and B beats mask 19, A→B, where bit i
refers to spine i in increasing order. Independent literal
Hamilton path subset DP reproduces that count. See
[census.txt](census.txt) for all coverage tables and witnesses.

## An independent 12-state transfer for each apex order

Fix the order A-before-B in a Hamilton path and scan the
increasing spine. Each site must be assigned to exactly one
of the three increasing runs, before A, between A and B, or
after B. These assignments have a one-to-one correspondence
to paths satisfying four endpoint tests. The first run is
empty or ends at a site that beats A. The middle run, unless
empty, begins at a site beaten by A and ends at a site that
beats B. The last run is empty or begins at a site beaten by
B. The middle run may be empty exactly if A→B is present.

The first run has two states: legal last endpoint (can also
be empty) and illegal last endpoint. Assignment of a site
to this run resets the endpoint status. The middle run has
three states (empty, nonempty with legal endpoint, nonempty
with illegal endpoint); new assignment is allowed unless the
run is empty and the site is not beaten by A. The last run has two
states (empty or already activated), and a site first activates
it iff B beats the site. Their product gives 12 states. The
12-state transfer with A,B swapped handles the other order.
All counters at length n are bounded by 2·3^n, because
there are at most 3^n colorings per order.

The exact census enumerates all pairs of masks with the
support bound; each recursion from a given B mask extends
both A-out choices. It independently checks every tape
through spine 5 against literal subset DP, and verifies
the displayed large witnesses by the same literal recurrence.
Prepending a 00 spine letter is an actual source, so all
earlier-length counts already occur at exactly the listed
length without using a B contact.

Reproduce both lengths:

    c++ -std=c++17 -O3 census.cpp -o /tmp/apex-census
    /tmp/apex-census

The first-miss census is not a cutoff theorem for Chan–Pak.
