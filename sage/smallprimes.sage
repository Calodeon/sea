
l = 2

def findNonSquareInFiniteField(l):
    for i in range(1,l):
        if (jacobi_symbol(i,l) == -1) :
            return i

for i in range(0,10):
    l = next_prime(l)
    d = findNonSquareInFiniteField(l)
    x = polygen(GF(l))
    F.<a> = GF(l^2, modulus=x^2-d)
    g = F.multiplicative_generator()

    print str(l) + " " + str(d) + " " + str(g.vector())