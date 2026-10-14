Q = 4*t^2*(4*(a-(t^2-1)/(4*t))^2+1);
E = ellinit(ellfromeqn(y^2 - subst(Q, t, x)));
J = E.j; print("j(E_a) = ", J);
print("equals 16*(16a^4+16a^2+1)^3/(a^2*(a^2+1)): ", J == 16*(16*a^4+16*a^2+1)^3/(a^2*(a^2+1)));
quit
