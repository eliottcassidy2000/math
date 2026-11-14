# What the full digit bank does not preserve

The [full 99-offset bank](FULL_99_DIGIT.md) gives all the output digits, but the output must be usable as a next reservoir before the result becomes an induction. The immediate optional-two-slot response supplies a decisive check: **the 99 selected templates all have different uniform next-state responses, even allowing different amounts of singleton padding after the output**. This is a scoped statement about the published bank, not all charts at each digit.

For each fixed star with the nonsingleton SCC retained, compress a two-slot cover of the full output to itineraries, with each maximal run in U replaced by its mark. Each slot can visit U once without a return. Every extra visit has to traverse the apex: leaving U for W cannot return to U directly, and there are no arcs from U back to I. Since the apex belongs to only one of the two slots, there is at most one extra visit. Hence the exact response is a linear combination of F_1(U),F_2(U),F_3(U), the ordered nonempty path-cover counts. In optional-slot coordinates it has the unique form

    Q_output(2) = A H(U)+B Q_U(2)+C Q_U(3).

Singleton SCCs removed from the star are fixed outside U; the same marked order and run bound apply after the deletion. For the transitive inputs of lengths 1,2,3 the displayed expression is A+B·2^m+C·3^m. This three-point probe is exact for arbitrary U because the itinerary proof already limits the response to these three basis functions. It does not assume an empirical polynomial fit. For the phase-r digit template, its response on the original U before trimming is therefore the coefficient vector

    (A, B/2^r, C/3^r).

Appending a singleton SCC to the output leaves its H unchanged and doubles Q(2), so neutral output padding only multiplies the whole vector by a power of two. Normalize each vector by dividing by the two-adic part of its first nonzero rational coefficient. The literal replay finds 99 distinct normalized rational vectors, one for each c. There is no equality of uniform Q(2) transfer functions even after such padding, since the admissible transitive inputs (all lengths at least five) already separate the three basis functions. Two outputs can still have equal Q(2) for a special chosen U; this proof neither rules that out nor classifies unselected templates.

The new concrete search obligation is thus stronger than finding a complete digit alphabet: either choose the cut and its reservoir so that the desired output also has the next exact Q(2), or discover a second alphabet whose output responses can be synchronized. A uniform claim of continuation using just the 99 existing choices and neutral padding is not valid.
