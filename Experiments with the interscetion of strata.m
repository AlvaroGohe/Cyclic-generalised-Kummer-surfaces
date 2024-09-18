F<f0, f1, f2, f3, f4, f5, f6, g0, g1, g2, g3> :=PolynomialRing(GF(2),11);
GenC := HyperellipticCurve(Polynomial([f0, f1, f2, f3, f4, f5, f6]),Polynomial([g0, g1, g2, g3]));
GenC;
K<f0, f1, f2, f3, f4, f5, f6, g0, g1, g2, g3> :=BaseField(GenC);
IgusaInvariants(GenC);
F<f0, f1, f2, f3, f4, f5> :=PolynomialRing(GF(2),6);
GenC := HyperellipticCurve(Polynomial([f0, f1, f2, f3, f4, f5, 1]),Polynomial([0, 0, f1*0, 1]));
GenC;
K<f0, f1, f2, f3, f4, f5> :=BaseField(GenC);
IgusaInvariants(GenC);
F<f0, f1, f2, f3, f4, f5> :=PolynomialRing(GF(2),6);
GenC := HyperellipticCurve(Polynomial([f0, f1, f2, f3, f4, f5, 1]),Polynomial([0, f1*0, 1, 1]));
GenC;
K<f0, f1, f2, f3, f4, f5> :=BaseField(GenC);
IgusaInvariants(GenC);

F<f0, f1, f2, f3, f4, f5> := PolynomialRing(GF(3),6);
GenC := HyperellipticCurve(Polynomial([f0, (f2*f4)/f5, f2, f3, f4, f5, 1]));
GenC;
K<f0, f1, f2, f3, f4, f5> :=BaseField(GenC);
IgusaInvariants(GenC);

F<f1, f2, f3, f4, f5, f6> := PolynomialRing(GF(5),6);
GenC := HyperellipticCurve(Polynomial([(2*(3*f2^2*f4^2 + f1*f3*f4^2 + f2^2*f3*f5 + 2*f1*f3^2*f5 + 3*f1*f2*f4*f5 + f2^3*f6))/(3*f4^3 + 4*f3^2*f6 + f2*f4*f6), f1, f2, f3, f4, f5, f6]));
GenC;
K<f1, f2, f3, f4, f5, f6> :=BaseField(GenC);
IgusaInvariants(GenC);