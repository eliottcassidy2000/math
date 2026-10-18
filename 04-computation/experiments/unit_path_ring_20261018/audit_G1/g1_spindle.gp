read("g1_common.gp");
e = 0; for (a = 1, #Mset, for (b = a + 1, #Mset, if (isU(Mset[a] - Mset[b]), e++)));
print("unit edges in M: ", e, " (Moser spindle: 7 vertices, 11 edges)");
eC = 0; for (a = 1, #Cset, for (b = a + 1, #Cset, if (isU(Cset[a] - Cset[b]), eC++)));
print("unit edges in C: ", eC);
quit;
