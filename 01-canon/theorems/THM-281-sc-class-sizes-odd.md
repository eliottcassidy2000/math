# THM-281: Self-Complementary Class Sizes Are Always Odd

**Status:** PROVED  
**Session:** opus-2026-04-03-S27  
**Depends on:** THM-280 (grid reflection = complement)  

## Statement

For any tournament T on n vertices, if the isomorphism class [T] is self-complementary (T ≅ T^op), then the tiling fiber |[T]| is odd.

Equivalently: the number of distinct labeled tournaments isomorphic to T, restricted to those with the fixed base path n → (n-1) → ... → 1, is always odd when T is self-complementary.

## Proof

By THM-280, the grid reflection R acts as an involution on the tiling space that induces T → T^op at the class level. For a self-complementary class [T], the reflection R maps the fiber (the set of all tilings representing [T]) to itself.

Since R is an involution (R² = id), it partitions the fiber into:
- **Fixed points**: tilings b where R(b) = b (grid-symmetric tilings)
- **Free orbits**: pairs {b, R(b)} where R(b) ≠ b

Therefore: |fiber| = |fixed points| + 2 × |free orbits|.

The fixed points are exactly the grid-symmetric tilings in this class. By THM-048 (gs-class-size-odd), the number of grid-symmetric tilings in a self-complementary class is always odd.

Therefore: |fiber| = (odd) + 2k = odd. □

## Consequence

In the merged metagraph G_n/Z_2:
- **Backbone nodes** (SC classes): tiling count is always **odd**
- **Paired nodes** (NS classes collapsed): combined tiling count is always **even** (2 × class_size, since complement pairs have equal size by |Aut(T)| = |Aut(T^op)|)

Therefore: 2^m = (sum of odd terms from SC classes) + (sum of even terms from NS pairs).

The **parity of the tiling count determines the merge type**: odd = backbone, even = paired. This is visible directly from the numbers.

## Verified

Computationally verified for all classes at n = 3, 4, 5, 6, 7.
