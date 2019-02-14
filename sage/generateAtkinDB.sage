AtkinDB = AtkinModularPolynomialDatabase()
l = 2
while (l < 500):
    f = open("Atkin/phi_a" + str(l) + ".txt", 'w')
    L = list(AtkinDB[l])
    for (c,P) in L:
        s = str(P.degrees()).replace("(","[").replace(")","]").replace(", ",",")
        f.write(s + " " + str(c) + "\n")
    f.close()
    l = next_prime(l)