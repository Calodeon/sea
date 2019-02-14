p = 2037035976334486086268445688409378161051468393665936250636140449354381299763336706183398631

F = GF(p)

R.<X,Y> = F['X','Y']

phi(X,Y) = 141359947154721358697753474691071362751004672000 + 53274330803424425450420160273356509151232000 * (X + Y) - 264073457076620596259715790247978782949376 * (X*Y) + 6692500042627997708487149415015068467200 * (X^2 + Y^2) + 36554736583949629295706472332656640000 * (X^2 * Y + Y^2 * X) + 5110941777552418083110765199360000 * X^2 * Y^2 + 280244777828439527804321565297868800 * (X^3 + Y^3) - 192457934618928299655108231168000 * (X^3 * Y + Y^3 * X) + 26898488858380731577417728000 * (X^3 * Y^2 + Y^3 * X^2) - 441206965512914835246100 * (X^3 * Y^3) + 1284733132841424456253440 * (X^4 + Y^4) + 128541798906828816384000 * (X^4 * Y + Y^4 * X) + 383083609779811215375 * (X^4 * Y^2 + Y^4 * X^2) + 107878928185336800 * (X^4 * Y^3 + X^3 * Y^4) + 1665999364600 * (X^4 * Y^4) + 1963211489280 * (X^5 + Y^5) - 246683410950 * (X^5 * Y + Y^5 * X) + 2028551200 * (X^5 * Y^2 + X^2 * Y^5) - 4550940 * (X^5 * Y^3 + X^3 * Y^5) + 3720 * (X^5 * Y^4 + X^4 * Y^5) + Y^6 + X^6 - X^5 * Y^5

phif(X,Y) = X^6 + Y^6 - X^5 * Y^5 + 4 * X*Y

FF(X) = (X^24 - 16)^3 / X^24

phiX = phi.derivative(X)
phiY = phi.derivative(Y)
phiXY = phiX.derivative(Y)
phiXX = phiX.derivative(X)
phiYY = phiY.derivative(Y)

phifX = phif.derivative(X)
phifY = phif.derivative(Y)
phifXY = phifX.derivative(Y)
phifXX = phifX.derivative(X)
phifYY = phifY.derivative(Y)

FP = FF.derivative()
FPP = FP.derivative()

A1 = F(164)
B1 = F(277)

EC1 = EllipticCurve(F, [0,0,0,A1,B1])

N1 = 2037035976334486086268445688409378161051468394649346661970376171306709029254574151895545000

j1 = F(1089313181734818163915301575565179061318882213618098684970663300430780893677503404292632181)
f1 = F(1178508179430829697899397473751835138435193504605125973124453448032052007006846870848434177)

j2 = F(38703108668702376420673779515777016027210246044506308289731263484777885353192969744583853)
f2 = F(617886390098565140587217812201880803383402853975116690813788162765371356263851914927743859)

phij = phi(j1,j2)
phiXj = phiX(j1,j2)
phiYj = phiY(j1,j2)
phiXYj = phiXY(j1,j2)
phiXXj = phiXX(j1,j2)
phiYYj = phiYY(j1,j2)

phiff = phif(f1,f2)
phifXf = phifX(f1,f2)
phifYf = phifY(f1,f2)
phifXYf = phifXY(f1,f2)
phifXXf = phifXX(f1,f2)
phifYYf = phifYY(f1,f2)

l = 5

m = F(18) * B1 / A1
j1P = m*j1
k = j1P / (F(1728) - j1)
j2P = - j1P * phiXj / (F(l) * phiYj)
m2 = j2P / j2
k2 = j2P/(F(1728) - j2)
A2 = F(l*l*l*l)*m2*k2 / F(48)
B2 = F(l**6) * m2 * m2 * k2 / F(864)

EC2 = EllipticCurve(F, [0,0,0,A2,B2])

r = -(j1P*j1P * phiXXj + F(2)*F(l)*j1P*j2P*phiXYj + F(l*l)*j2P*j2P*phiYYj)/(j1P * phiXj)

f1P = j1P / FP(f1)
f2P = j2P / FP(f2)
rf = -(f1P*f1P * phifXXf + F(2)*F(l)*f1P*f2P*phifXYf + F(l*l)*f2P*f2P*phifYYf)/(f1P * phifXf)

err1 = FPP(f1)/FP(f1)*f1P
err2 = F(l)*FPP(f2)/FP(f2)*f2P
r2 = rf + err1 - err2



p1 = F(l)*(r/F(2) + (k-l*k2)/F(4) + (l*m2 - m)/F(3))

d = (l-1)/2

t = []
t.append(F(d))
t.append(p1/F(2))
t.append((F(1 - 10*d) * A1 - A2)/F(30))
t.append((F(1 - 28*d)*B1 - F(42)*t[1]*A1 - B2)/F(70))

#c = []
#c.append(F(0))
#c.append(F(6) * t[2] + F(2) * A1 * t[0])
#c.append(F(10) * t[3] + F(6) * A1 * t[1] + F(4) * B1 * t[0])

s = []
s.append(F(1))
s.append(t[1])
s.append((t[1]*s[1] - t[2])/F(2))
s.append((t[1]*s[2] - t[2]*s[1] + t[3])/F(3))

R.<x> = F['x']

psi = x^2 - F(s[1]) * x^1 + F(s[2])












