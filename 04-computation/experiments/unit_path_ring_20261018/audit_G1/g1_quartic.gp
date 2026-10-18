\\ audit G1: the walk-length quartic of THM-4700 (remarks): N(x - u_l) (l^2 + 3) = A_x(l), u_l = (l + s3)/(l - s3).
read("g1_common.gp");
A(xx, l) = Nm(xx - 1)*l^2 + 2*s3*(xx - cj(xx))*l + 3*Nm(xx + 1);
ok = 1;
{
foreach([bb, cc, aa, rho, 1/3 + zF, wF/7 + jF], xx, foreach([0, 1, 2/3, -5, 7/11], l,
  my(ul = (l + s3)/(l - s3));
  ok = ok && isU(ul) && Nm(xx - ul)*(l^2 + 3) == A(xx, l) && cj(A(xx, l)) == A(xx, l)));
}
print("u_l unit and N(x - u_l)(l^2 + 3) = A_x(l), A_x(l) real, for 6 x and 5 l: ", ok);
\\ the Cayley parametrisation covers U: for a unit z != 1, l = s3 (z + 1)/(z - 1) is real and u_l = z
z = cj((1 + 5*s3)/2)/((1 + 5*s3)/2); l = s3*(z + 1)/(z - 1);
print("inverse map: l real ", cj(l) == l, " ; u_l = z ", (l + s3)/(l - s3) == z);
quit;
