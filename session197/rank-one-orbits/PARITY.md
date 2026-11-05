# Parity of the strong rank-one class count

This short consequence was suggested by the rank-one packet's author and independently checked during review. It is separate from that frozen packet.

Let C_t count ordinary unmarked strong rank-one self-converse tournament classes on2t+1 vertices, with the anti-fixed-root rank-one convention of the reviewed correspondence. For every t>=1,

    C_t is even exactly when3 divides t.

Indeed each class contributes the positive odd integer c(T)=H(S)/|Aut(S,a)| to its c-mass. Hence the ordinary count and the weighted count have the same parity. The reviewed weighted series is

    C_mass(x)=x M'(x)/M(x),
    M(x)=1+sum_(t>=1)2^((t-1)(t-2)/2)x^t.

Modulo2, M(x)=1+x+x². Therefore, as a formal power series over the two-element field,

    C_mass(x)=x/(1+x+x²)
             =(x+x²)/(1+x³)
             =sum_(k>=0)(x^(3k+1)+x^(3k+2)).

The coefficient at x^t is zero exactly for t divisible by3, proving the claim. The constant term is also zero, consistently with the absence of an empty rank-one class; the theorem itself concerns t>=1.

The independent physical inventories through t=5 give1,1,2,7,33. The independent checker also calculates the integer weighted recurrence through t=50 and verifies the displayed parity. Those finite controls are not the all-order proof.
