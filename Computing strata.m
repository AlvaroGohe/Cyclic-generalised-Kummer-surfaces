K<f0,f1,f2,g0,g1>:= PolynomialRing(Rationals(),5);
GenC := C3Curve(K,[f0,f1,f2,g0,g1]);
K<f0,f1,f2,g0,g1>:= BaseRing(GenC);
II := [Numerator(i): i in IgusaInvariants(GenC)];
K!II;

Q0 := Rationals();
P := ProjectiveSpace(Q0,3);
K<s1,s2,s3,s4> := CoordinateRing(P);
C := HyperellipticCurve(Polynomial([s1,0,s2,0,s3,0,s4]));
K<s1,s2,s3,s4>:= BaseRing(C);
II := [Numerator(i): i in IgusaInvariants(C)];
[Evaluate(i,[-1,s2,-s1,1]): i  in II];
IP<I2,I4,I6,I10> := ProjectiveSpace(Q0,[1,2,3,5]);
phi := map<P->IP | II[[1,2,3,5]]>;
J := ideal< CoordinateRing(IP) | &cat[ DefiningEquations(Image(phi,P,d)) : d in  [1..15]] >;
MinimalBasis(J);
V4sp := Scheme(IP,MinimalBasis(J));
IP<I2,I4,I6,I8,I10> := ProjectiveSpace(Q0,[1,2,3,4,5]);
phi := map<P->IP | II>;
J := ideal< CoordinateRing(IP) | &cat[ DefiningEquations(Image(phi,P,d)) : d in  [1..15]] >;
MinimalBasis(J);
V4 := Scheme(IP,MinimalBasis(J));
E1 := EllipticCurve(Polynomial([s1*s4^2,s2*s4,s3,1]));
j1 := jInvariant(E1);
P2<w,t,z> :=ProjectiveSpace(K,2);
C2 := Curve(P2, w^2*z^2-(s4*t^4+s3*t^3*z+s2*t^2*z^2+s1*t*z^3));
E2 := EllipticCurve(C2,C2![0,0,1]);
j2 := jInvariant(E2);
E2 := EllipticCurve(Polynomial([s4*s1^2,s3*s1,s2,1]));
j2 := jInvariant(E2);
fac := s2^2*s3^2 - 4*s1*s3^3 - 4*s2^3*s4 + 18*s1*s2*s3*s4 - 27*s1^2*s4^2;
polcoeffs := [Numerator(j1*j2),-Numerator(j1+j2)*fac,Denominator(j1*j2)];
[Factorisation(Numerator(Evaluate(coef,[-1,s2,-s1,1]))): coef in polcoeffs];
Evaluate(Numerator(j1*j2),[-1,s2,-s1,1]);
Factorisation(Polynomial([Numerator(j1*j2),-Numerator(j1+j2)*fac,Denominator(j1*j2)]));
P3<s1,s2,s3,s4> := ProjectiveSpace(Q0,3);
P112<h1,h2,h3> := ProjectiveSpace(Q0, [1,1,2]);
phi := map <P3->P112 | [-32*(s2*s3 + 15*s1*s4), 64*(s2*s3 + 14*s1*s4), -16*(-65*s2^2*s3^2 + 4*s1*s3^3 + 4*s2^3*s4 - 1938*s1*s2*s3*s4 - 13605*s1^2*s4^2)]>;
[DefiningEquations(phi(Scheme(P3,Numerator(coeff))))[1]: coeff in polcoeffs];
// DefiningEquations(phi(Scheme(P3, 65536*s2^6*s3^6 - 589824*s1*s2^4*s3^7 + 1769472*s1^2*s2^2*s3^8 - 1769472*s1^3*s3^9 - 589824*s2^7*s3^4*s4 +         5308416*s1*s2^5*s3^5*s4 - 15925248*s1^2*s2^3*s3^6*s4 + 15925248*s1^3*s2*s3^7*s4 + 1769472*s2^8*s3^2*s4^2 -         15925248*s1*s2^6*s3^3*s4^2 + 47775744*s1^2*s2^4*s3^4*s4^2 - 47775744*s1^3*s2^2*s3^5*s4^2 - 1769472*s2^9*s4^3 +         15925248*s1*s2^7*s3*s4^3 - 47775744*s1^2*s2^5*s3^2*s4^3 + 47775744*s1^3*s2^3*s3^3*s4^3)));
// DefiningEquations(phi(Scheme(P3, -256*s1^2*s2^2*s3^8 + 1024*s1^3*s3^9 + 3328*s1^2*s2^3*s3^6*s4 - 13824*s1^3*s2*s3^7*s4 - 256*s2^8*s3^2*s4^2 +         3328*s1*s2^6*s3^3*s4^2 - 32256*s1^2*s2^4*s3^4*s4^2 + 103680*s1^3*s2^2*s3^5*s4^2 - 20736*s1^4*s3^6*s4^2 + 1024*s2^9*s4^3 -         13824*s1*s2^7*s3*s4^3 + 103680*s1^2*s2^5*s3^2*s4^3 - 304128*s1^3*s2^3*s3^3*s4^3 + 62208*s1^4*s2*s3^4*s4^3 -         20736*s1^2*s2^6*s4^4 + 62208*s1^3*s2^4*s3*s4^4 + 373248*s1^4*s2^2*s3^2*s4^4 - 186624*s1^5*s3^3*s4^4 -         186624*s1^4*s2^3*s4^5)));
// DefiningEquations(phi(Scheme(P3, s1^2*s2^4*s3^4*s4^2 - 8*s1^3*s2^2*s3^5*s4^2 + 16*s1^4*s3^6*s4^2 - 8*s1^2*s2^5*s3^2*s4^3 + 68*s1^3*s2^3*s3^3*s4^3 -         144*s1^4*s2*s3^4*s4^3 + 16*s1^2*s2^6*s4^4 - 144*s1^3*s2^4*s3*s4^4 + 270*s1^4*s2^2*s3^2*s4^4 + 216*s1^5*s3^3*s4^4 +         216*s1^4*s2^3*s4^5 - 972*s1^5*s2*s3*s4^5 + 729*s1^6*s4^6)));
P := ProjectiveSpace(Q0,6);
K<a0,a1,a2,a3,a4,a5,a6> := CoordinateRing(P);
C := HyperellipticCurve(Polynomial([a0,a1,a2,a3,a4,a5,a6]));
K<a0,a1,a2,a3,a4,a5,a6>:= BaseRing(C);
II := [Numerator(i): i in IgusaInvariants(C)];



Q0 := Rationals();
// Q0 := GF(2);
P := ProjectiveSpace(Q0,[2,2,2,2,1,1]);
K<f0,f1,f2,f3,g0,g1>:= CoordinateRing(P);
C := HyperellipticCurve(Polynomial([-f0, -f1, -f2, -f3, -f2, -f1, -f0]),Polynomial([g0, g1, g1, g0]));
C := HyperellipticCurve(Polynomial([0, -f1, 0,f3,0,-f1, 0]),Polynomial([g0, 0, 0, g0]));
K<f0,f1,f2,f3,g0,g1>:= BaseRing(C);
Factorisation(Numerator(Discriminant(C)));
II := [Numerator(i): i in IgusaInvariants(C)];
IP<I2,I4,I6,I8,I10> := ProjectiveSpace(Q0,[1,2,3,4,5]);
phi := map<P->IP | II>;
J := ideal< CoordinateRing(IP) | &cat[ DefiningEquations(Image(phi,P,d)) : d in  [1..15]] >;
MinimalBasis(J);
V4 := Scheme(IP,MinimalBasis(J));
IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(V4)));
[Factorisation(i): i in II];



Q0 := GF(2);
P := ProjectiveSpace(Q0,[2,2,2,2,1,1]);
K<f0,f1,f2,f5,g0,g1>:= CoordinateRing(P);
C := HyperellipticCurve(Polynomial([f0, f1, f2, f5, f1 + f2 + f5, f5, f5]),Polynomial([g0, g1, g1]));
K<f0,f1,f2,f5,g0,g1>:= BaseRing(C);
II := [Numerator(i): i in IgusaInvariants(C)];
IP<I2,I4,I6,I8,I10> := ProjectiveSpace(Q0,[1,2,3,4,5]);
phi := map<P->IP | II>;
J := ideal< CoordinateRing(IP) | &cat[ DefiningEquations(Image(phi,P,d)) : d in  [1..15]] >;
MinimalBasis(J);
V4 := Scheme(IP,MinimalBasis(J));
IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(V4)));
P114<s,t,n> := ProjectiveSpace(Q0,[1,1,4]);
phi := map<V4->P114 | [I2^2,I4,I2^4*I8+I2^2*I4^3+I4^4]>;
IrreducibleComponents(ReducedSubscheme(BaseScheme(phi)));
Pullback(phi, P2![0,0,1]);
P2<s,t,m> := ProjectiveSpace(Q0,2);
P1 := ProjectiveSpace(Q0,[4,5]);
proj := map<P2->P1| [(t+m)^4,t*m^4]>;
Image(proj);
inj := map<P2->V4 | [s, s*t, s*t^2, m^4+s*t^3+t^4, m^4*t]>;
[IrreducibleComponents(Pullback(inj,Scheme(V4,V4.i))): i in [1..4]];
IrreducibleComponents(ReducedSubscheme(BaseScheme(phi)));
//P := ProjectiveSpace(Q0,[2,2,2,2,2,2,2,1]);
//K<f0,f1,f2,f3,f4,f5,f6,g0>:= CoordinateRing(P);
//C := HyperellipticCurve(Polynomial([f0,f1,f2,f3,f4,f5,f6]),Polynomial([0*f0,g0]));
///K<f0,f1,f2,f3,f4,f5,f6,g0>:= BaseRing(C);
//II := [Numerator(i): i in IgusaInvariants(C)];
//IP<I2,I4,I6,I8,I10> := ProjectiveSpace(Q0,[1,2,3,4,5]);
//phi := map<P->IP | II>;
//J := ideal< CoordinateRing(IP) | &cat[ DefiningEquations(Image(phi,P,d)) : d in  [1..15]] >;
//MinimalBasis(J);
W1 := Scheme(IP,[I2,I4]);
Intersection(W1,V4);

Q0 := Rationals();
Q0 := GF(3);
P := ProjectiveSpace(Q0,[1,1,1]);
K<f1,f3,f5>:= CoordinateRing(P);
C := HyperellipticCurve(Polynomial([0,f1,0,f3,0,f5]));
K<f1,f3,f5>:= BaseRing(C);
II := [Numerator(i): i in IgusaInvariants(C)];
IP<I2,I4,I6,I8,I10> := ProjectiveSpace(Q0,[1,2,3,4,5]);
phi := map<P->IP | II>;
J := ideal< CoordinateRing(IP) | &cat[ DefiningEquations(Image(phi,P,d)) : d in  [1..15]] >;
MinimalBasis(J);
D4 := Scheme(IP,MinimalBasis(J));
Ipproj<J2,J4,J6> := ProjectiveSpace(Q0,[1,2,3]);
proj := map<D4 -> Ipproj|[I2,I4,I6]>;
D4c := Image(proj);

Q0 := Rationals();
Q0 := GF(3);
P := ProjectiveSpace(Q0,[1,1,1]);
K<f1,f3,f5>:= CoordinateRing(P);
C := HyperellipticCurve(Polynomial([0,f1,0,f3,0,f5]));
K<f1,f3,f5>:= BaseRing(C);
II := [Numerator(i): i in IgusaInvariants(C)];
IP<I2,I4,I6,I8,I10> := ProjectiveSpace(Q0,[1,2,3,4,5]);
phi := map<P->IP | II>;
J := ideal< CoordinateRing(IP) | &cat[ DefiningEquations(Image(phi,P,d)) : d in  [1..15]] >;
MinimalBasis(J);
D4 := Scheme(IP,MinimalBasis(J));
Ipproj<J2,J4,J6> := ProjectiveSpace(Q0,[1,2,3]);
proj := map<D4 -> Ipproj|[I2,I4,I6]>;
D4c := Image(proj);




K<t>:= PolynomialRing(Rationals());
C := HyperellipticCurve(Polynomial([t,0,0,1,0,0,1]));
K<t> := BaseRing(C);
II := IgusaInvariants(C);
GenC := C3Curve(K,[t,0,0,1,1]);
K<t> := BaseRing(C);
II := IgusaInvariants(GenC);

K<t>:= PolynomialRing(GF(2));
C := HyperellipticCurve(Polynomial([0,0,0,t,0,1]),1);
K<t> := BaseRing(C);
II := IgusaInvariants(C);

K<t>:= PolynomialRing(Rationals());
C := HyperellipticCurve(Polynomial([t,0,0,0,0,0,2]),Polynomial([1,t*0,0,1]));
K<t> := BaseRing(C);
Discriminant(C);
II := IgusaInvariants(C);

K<a,t>:= PolynomialRing(Rationals(),2);
C := HyperellipticCurve(Polynomial([t,0,0,0,0,0,a]),Polynomial([1,t*0,0,1]));
K<a,t> := BaseRing(C);
II := IgusaInvariants(C);

P2<f0,f1,f2,t>:= ProjectiveSpace(Rationals(),3);
S := Scheme(P2,720*f0^2 + 720*f0*f1 + 220*f1^2 + 80*f1*f2 + 16*f2^2-t^2);
HasRationalPoint(S);
RationalPoints(S);


K<z> := ext<GF(2)|Polynomial([1,1,1,1,1])>;
Factorisation(Polynomial([1,0,0,0,0,1]));
KK<t> := PolynomialRing(K);
Roots(t^6+t+K.1);
PolK<x> :=PolynomialRing(KK);
Factorisation(x^(16)+(t^6+t)^2*x^8+(t^6+t)*x^2+x);
F<b> := ext<K | Polynomial([1,z*(1+z),0,0,0,1])>;
PolF<x> := PolynomialRing(F);
Factorisation(x^(16)+(1+z)^2*z^2*x^8+(1+z)*z*x^2+x);
c1s := [r[1]: r in Roots(x^(16)+(1+z)^2*z^2*x^8+(1+z)*z*x^2+x)];
c2s := [Roots(Polynomial([c1+c1^2*z*(1+z),0,1]))[1,1]: c1 in c1s];
c3s1 := [Roots(Polynomial([c1s[i]^2*c2s[i],1,1]))[1,1]: i in [1..16]];
c3s2 := [Roots(Polynomial([c1s[i]^2*c2s[i],1,1]))[2,1]: i in [1..16]];
C<x,y,s> := HyperellipticCurve(Polynomial([b*0,0,0,(1+z),0,1]),1);
[map< C->C | [x+z^2*c1s[i]*s,y+z*c1s[i]*x^2*s+z^3*c2s[i]*x*s^2+c3s1[i]*s^3,s]>: i in [1..2]];


C := HyperellipticCurve(Polynomial([0,0,0,1,0,1]),1);
Factorisation(x^(16)+z^2*x^8+z*x^2+x);
r5 := [r[1]: r in Roots(x^5 + (z^2 + z)*x + 1)];
[[<r+s in r10>: r in r5]: s in r5];
SequenceToSet(&cat[[r+s: r in r5]: s in r5]);
#SequenceToSet(&cat[[r+s: r in r5]: s in r5]);
r10 := [r[1]: r in Roots(x^10 + z*(1+z)*x^6 + x^5 + 1)];
SequenceToSet(&cat[[r+s: r in r10]: s in r10]);
[[<r+s in r5,r+s in r10>: r in r10]: s in r10];
[[<r+s in r5>: r in r10]: s in r10];
[[<r+s in r5>: r in r5]: s in r10];
[[<r+s in r5>: r in r5]: s in r10];
Factorisation(x^10 + z*(1+z)*x^6 + x^5 + 1);
ra := [r[1]: r in Roots(x^5 + z^2*x^3 + (z + 1)*x^2 + (z^2 + z)*x + z^3 + z^2 + 1)];
rb := [r[1]: r in Roots(x^5 + z^2*x^3 + (z + 1)*x^2 + (z^3 + z^2 + z + 1)*x + z^3 + z^2)];



Pol<x> := PolynomialRing(GF(3));
C := HyperellipticCurve(x*(x^2-1)*(x^2+1));
HasseWittInvariant(C);