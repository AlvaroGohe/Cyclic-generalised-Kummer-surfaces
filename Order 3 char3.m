



load "Functions.m";
load "WeightedJacobian.m";
A6<x1,x2,x3,y1,y2,y3>:=AffineSpace(Rationals(),6);

function InvariantEmbedding(J)
 C :=  Curve(J);
 k := BaseRing(C);
 f,g := HyperellipticPolynomials(C);
 f0 := Coefficient(f,0);
 f1 := Coefficient(f,1);
 f2 := Coefficient(f,2);
 f3 := Coefficient(f,3);
 f4 := Coefficient(f,4);  
 f5 := Coefficient(f,5);
 f6 := Coefficient(f,6);
 g0 := Coefficient(g,0);
 g1 := Coefficient(g,1);
 g2 := Coefficient(g,2);
 g3 := Coefficient(g,3);
 if ([10*f0 + 5*f1 + 2*f2 + f3, -15*f0 - 5*f1 - f2 + f4, 6*f0 + f1 + f5, -f0 + f6, 3*g0 + g1 + g2, -g0 + g3] ne [0,0,0,0,0,0]) then 
 error "This function only works for Jacobians of curves of the form y^2+(g0 x^3 - 3 g0 x^2 - g1 x^2 + g1 x + g0)y =f0 x^6 - 6 f0 x^5 - f1 x^5 + 15 f0 x^4 + 5 f1 x^4 + f2 x^4 - 10 f0 x^3 - 5 f1 x^3 - 2 f2 x^3 + f2 x^2 + f1 x + f0. Sorry :()";
 end if;
 WJ<k1, k2, k3, k4, v1, v2, v3, v4, v5, v6> := WeightedJacobian(J);
 WP<c1, c2, cvw, d1, d2, e1, e2, e3, e4> := ProjectiveSpace(k,[1,1,2,2,2,3,3,3,3]);
 phi := map < WJ -> WP | [-6*f0*k1 - f1*k1 + 2*f0*k2 + f1*k3 + k4, 2*k1 - k2 + 2*k3, -k1^2 + k1*k2 - 3*k1*k3 + k2*k3 - k3^2, f1*k1^2 - 10*f0*k1*k3 - 5*f1*k1*k3 - 2*f2*k1*k3 - 6*f0*k3^2 - f1*k3^2 +   k2*k4 + 2*v4, -(g0*k1^2) + g0*k1*k2 + 3*g0*k1*k3 + 2*g1*k1*k3 - g0*k2*k3 + g0*k3^2 - 2*v2 + 2*v3, (k1 - k2)*(k1 - k3)*(k2 - k3),  k1^3 - 2*k1^2*k2 + k1*k2^2 + 4*k1^2*k3 - 3*k1*k2*k3 + 3*k1*k3^2 - k2*k3^2 + k3^3, -(g0^2*k1^3) + 2*g0*g1*k1^3 - f1*k1^2*k2 - 3*g0^2*k1^2*k2 - 2*g0*g1*k1^2*k2 +   2*f1*k1^2*k3 + 5*g0^2*k1^2*k3 - 3*g0*g1*k1^2*k3 - 2*g1^2*k1^2*k3 + 10*f0*k1*k2*k3 + 5*f1*k1*k2*k3 + 2*f2*k1*k2*k3 + 3*g0^2*k1*k2*k3 + g0*g1*k1*k2*k3 - g0^2*k2^2*k3 -   20*f0*k1*k3^2 - 10*f1*k1*k3^2 - 4*f2*k1*k3^2 - 4*g0^2*k1*k3^2 - g0*g1*k1*k3^2 + 6*f0*k2*k3^2 + f1*k2*k3^2 - 12*f0*k3^3 - 2*f1*k3^3 + g0^2*k3^3 - k2^2*k4 + 2*k2*k3*k4 +   2*g0*k1*v1 - 2*g0*k3*v1 - 2*g0*k1*v2 + 2*g1*k1*v2 + 2*g0*k3*v2 - 2*g1*k1*v3 - 2*g0*k2*v3 + 2*g0*k3*v3 - 2*k2*v4 + 4*k3*v4,  -2*f1*k1^3 - g0^2*k1^3 + 2*g0*g1*k1^3 + 2*f1*k1^2*k2 - 3*g0^2*k1^2*k2 - 5*g0*g1*k1^2*k2 + 8*g0^2*k1*k2^2 + 3*g0*g1*k1*k2^2 - 2*g0^2*k2^3 + 20*f0*k1^2*k3 + 10*f1*k1^2*k3 +   4*f2*k1^2*k3 + 3*g0^2*k1^2*k3 - 3*g0*g1*k1^2*k3 - 2*g1^2*k1^2*k3 - 20*f0*k1*k2*k3 - 10*f1*k1*k2*k3 - 4*f2*k1*k2*k3 - 11*g0^2*k1*k2*k3 + 2*g0*g1*k1*k2*k3 +   2*g1^2*k1*k2*k3 + 3*g0^2*k2^2*k3 - g0*g1*k2^2*k3 + 12*f0*k1*k3^2 + 2*f1*k1*k3^2 + 2*g0^2*k1*k3^2 + g0*g1*k1*k3^2 - 12*f0*k2*k3^2 - 2*f1*k2*k3^2 + g0*g1*k2*k3^2 +   g0^2*k3^3 - 2*k1*k2*k4 + 2*k2^2*k4 + 2*g0*k1*v1 - 4*g0*k2*v1 + 2*g0*k3*v1 - 2*g0*k1*v2 + 2*g1*k1*v2 + 4*g0*k2*v2 - 2*g1*k2*v2 - 2*g0*k3*v2 - 4*g0*k1*v3 - 2*g1*k1*v3 +   2*g0*k2*v3 + 2*g1*k2*v3 + 2*g0*k3*v3 - 4*k1*v4 + 4*k2*v4]>;
 I := ideal< CoordinateRing(WP) | &cat[ DefiningEquations(Image(phi,WJ,d)) : d in  [3..6]] >;
 X3 := Scheme(WP,MinimalBasis(I));
 phi := map < WJ -> X3 | [-6*f0*k1 - f1*k1 + 2*f0*k2 + f1*k3 + k4, 2*k1 - k2 + 2*k3, -k1^2 + k1*k2 - 3*k1*k3 + k2*k3 - k3^2, f1*k1^2 - 10*f0*k1*k3 - 5*f1*k1*k3 - 2*f2*k1*k3 - 6*f0*k3^2 - f1*k3^2 +   k2*k4 + 2*v4, -(g0*k1^2) + g0*k1*k2 + 3*g0*k1*k3 + 2*g1*k1*k3 - g0*k2*k3 + g0*k3^2 - 2*v2 + 2*v3, (k1 - k2)*(k1 - k3)*(k2 - k3),  k1^3 - 2*k1^2*k2 + k1*k2^2 + 4*k1^2*k3 - 3*k1*k2*k3 + 3*k1*k3^2 - k2*k3^2 + k3^3, -(g0^2*k1^3) + 2*g0*g1*k1^3 - f1*k1^2*k2 - 3*g0^2*k1^2*k2 - 2*g0*g1*k1^2*k2 +   2*f1*k1^2*k3 + 5*g0^2*k1^2*k3 - 3*g0*g1*k1^2*k3 - 2*g1^2*k1^2*k3 + 10*f0*k1*k2*k3 + 5*f1*k1*k2*k3 + 2*f2*k1*k2*k3 + 3*g0^2*k1*k2*k3 + g0*g1*k1*k2*k3 - g0^2*k2^2*k3 -   20*f0*k1*k3^2 - 10*f1*k1*k3^2 - 4*f2*k1*k3^2 - 4*g0^2*k1*k3^2 - g0*g1*k1*k3^2 + 6*f0*k2*k3^2 + f1*k2*k3^2 - 12*f0*k3^3 - 2*f1*k3^3 + g0^2*k3^3 - k2^2*k4 + 2*k2*k3*k4 +   2*g0*k1*v1 - 2*g0*k3*v1 - 2*g0*k1*v2 + 2*g1*k1*v2 + 2*g0*k3*v2 - 2*g1*k1*v3 - 2*g0*k2*v3 + 2*g0*k3*v3 - 2*k2*v4 + 4*k3*v4,  -2*f1*k1^3 - g0^2*k1^3 + 2*g0*g1*k1^3 + 2*f1*k1^2*k2 - 3*g0^2*k1^2*k2 - 5*g0*g1*k1^2*k2 + 8*g0^2*k1*k2^2 + 3*g0*g1*k1*k2^2 - 2*g0^2*k2^3 + 20*f0*k1^2*k3 + 10*f1*k1^2*k3 +   4*f2*k1^2*k3 + 3*g0^2*k1^2*k3 - 3*g0*g1*k1^2*k3 - 2*g1^2*k1^2*k3 - 20*f0*k1*k2*k3 - 10*f1*k1*k2*k3 - 4*f2*k1*k2*k3 - 11*g0^2*k1*k2*k3 + 2*g0*g1*k1*k2*k3 +   2*g1^2*k1*k2*k3 + 3*g0^2*k2^2*k3 - g0*g1*k2^2*k3 + 12*f0*k1*k3^2 + 2*f1*k1*k3^2 + 2*g0^2*k1*k3^2 + g0*g1*k1*k3^2 - 12*f0*k2*k3^2 - 2*f1*k2*k3^2 + g0*g1*k2*k3^2 +   g0^2*k3^3 - 2*k1*k2*k4 + 2*k2^2*k4 + 2*g0*k1*v1 - 4*g0*k2*v1 + 2*g0*k3*v1 - 2*g0*k1*v2 + 2*g1*k1*v2 + 4*g0*k2*v2 - 2*g1*k2*v2 - 2*g0*k3*v2 - 4*g0*k1*v3 - 2*g1*k1*v3 +   2*g0*k2*v3 + 2*g1*k2*v3 + 2*g0*k3*v3 - 4*k1*v4 + 4*k2*v4]>;
 return X3, phi;
end function;

function newInvariantEmbedding(J)
 C :=  Curve(J);
 k := BaseRing(C);
 f,g := HyperellipticPolynomials(C);
 f0 := Coefficient(f,0);
 f1 := Coefficient(f,1);
 f2 := Coefficient(f,2);
 f3 := Coefficient(f,3);
 f4 := Coefficient(f,4);  
 f5 := Coefficient(f,5);
 f6 := Coefficient(f,6);
 g0 := Coefficient(g,0);
 g1 := Coefficient(g,1);
 g2 := Coefficient(g,2);
 g3 := Coefficient(g,3);
 if ([10*f0 + 5*f1 + 2*f2 + f3, -15*f0 - 5*f1 - f2 + f4, 6*f0 + f1 + f5, -f0 + f6, 3*g0 + g1 + g2, -g0 + g3] ne [0,0,0,0,0,0]) then 
 error "This function only works for Jacobians of curves of the form y^2+(g0 x^3 - 3 g0 x^2 - g1 x^2 + g1 x + g0)y =f0 x^6 - 6 f0 x^5 - f1 x^5 + 15 f0 x^4 + 5 f1 x^4 + f2 x^4 - 10 f0 x^3 - 5 f1 x^3 - 2 f2 x^3 + f2 x^2 + f1 x + f0. Sorry :()";
 end if;
 WJ<k1, k2, k3, k4, v1, v2, v3, v4, v5, v6> := WeightedJacobian(J);
 WP<c1, c2, cvw, d1, d2, e1, e2, e3, e4, q> := ProjectiveSpace(k,[1,1,2,2,2,3,3,3,3,4]);
 phi := map < WJ -> WP | [-6*f0*k1 - f1*k1 + 2*f0*k2 + f1*k3 + k4, 2*k1 - k2 + 2*k3, -k1^2 + k1*k2 - 3*k1*k3 + k2*k3 - k3^2, f1*k1^2 - 10*f0*k1*k3 - 5*f1*k1*k3 - 2*f2*k1*k3 - 6*f0*k3^2 - f1*k3^2 +   k2*k4 + 2*v4, -(g0*k1^2) + g0*k1*k2 + 3*g0*k1*k3 + 2*g1*k1*k3 - g0*k2*k3 + g0*k3^2 - 2*v2 + 2*v3, (k1 - k2)*(k1 - k3)*(k2 - k3),  k1^3 - 2*k1^2*k2 + k1*k2^2 + 4*k1^2*k3 - 3*k1*k2*k3 + 3*k1*k3^2 - k2*k3^2 + k3^3, -(g0^2*k1^3) + 2*g0*g1*k1^3 - f1*k1^2*k2 - 3*g0^2*k1^2*k2 - 2*g0*g1*k1^2*k2 +   2*f1*k1^2*k3 + 5*g0^2*k1^2*k3 - 3*g0*g1*k1^2*k3 - 2*g1^2*k1^2*k3 + 10*f0*k1*k2*k3 + 5*f1*k1*k2*k3 + 2*f2*k1*k2*k3 + 3*g0^2*k1*k2*k3 + g0*g1*k1*k2*k3 - g0^2*k2^2*k3 -   20*f0*k1*k3^2 - 10*f1*k1*k3^2 - 4*f2*k1*k3^2 - 4*g0^2*k1*k3^2 - g0*g1*k1*k3^2 + 6*f0*k2*k3^2 + f1*k2*k3^2 - 12*f0*k3^3 - 2*f1*k3^3 + g0^2*k3^3 - k2^2*k4 + 2*k2*k3*k4 +   2*g0*k1*v1 - 2*g0*k3*v1 - 2*g0*k1*v2 + 2*g1*k1*v2 + 2*g0*k3*v2 - 2*g1*k1*v3 - 2*g0*k2*v3 + 2*g0*k3*v3 - 2*k2*v4 + 4*k3*v4,  -2*f1*k1^3 - g0^2*k1^3 + 2*g0*g1*k1^3 + 2*f1*k1^2*k2 - 3*g0^2*k1^2*k2 - 5*g0*g1*k1^2*k2 + 8*g0^2*k1*k2^2 + 3*g0*g1*k1*k2^2 - 2*g0^2*k2^3 + 20*f0*k1^2*k3 + 10*f1*k1^2*k3 +   4*f2*k1^2*k3 + 3*g0^2*k1^2*k3 - 3*g0*g1*k1^2*k3 - 2*g1^2*k1^2*k3 - 20*f0*k1*k2*k3 - 10*f1*k1*k2*k3 - 4*f2*k1*k2*k3 - 11*g0^2*k1*k2*k3 + 2*g0*g1*k1*k2*k3 +   2*g1^2*k1*k2*k3 + 3*g0^2*k2^2*k3 - g0*g1*k2^2*k3 + 12*f0*k1*k3^2 + 2*f1*k1*k3^2 + 2*g0^2*k1*k3^2 + g0*g1*k1*k3^2 - 12*f0*k2*k3^2 - 2*f1*k2*k3^2 + g0*g1*k2*k3^2 +   g0^2*k3^3 - 2*k1*k2*k4 + 2*k2^2*k4 + 2*g0*k1*v1 - 4*g0*k2*v1 + 2*g0*k3*v1 - 2*g0*k1*v2 + 2*g1*k1*v2 + 4*g0*k2*v2 - 2*g1*k2*v2 - 2*g0*k3*v2 - 4*g0*k1*v3 - 2*g1*k1*v3 +   2*g0*k2*v3 + 2*g1*k2*v3 + 2*g0*k3*v3 - 4*k1*v4 + 4*k2*v4,k1^4 + 2*k1^3*k2 + 2*k1*k2^3 + 2*k1^3*k3 + k1^2*k2*k3 + k1*k2^2*k3 + k1^2*k3^2 + k1*k2*k3^2 + k2^2*k3^2 + 2*k1*k3^3 + 2*k3^4 + k1^2*k2*k4 + k2^2*k3*k4 + 2*k2*k3^2*k4 + k1^2*v1 + k1*k2*v1 + k2*k3*v1 + 2*k1^2*v2 + 2*k1*k2*v2 + 2*k2^2*v2 + k1*k3*v2 + 2*k2*k3*v2 + 2*k1^2*v3 + k3^2*v3 + 2*k1^2*v4 + 2*k2*k3*v4 + k3^2*v4]>;
 I := ideal< CoordinateRing(WP) | &cat[ DefiningEquations(Image(phi,WJ,d)) : d in  [3..6]] >;
 X3 := Scheme(WP,MinimalBasis(I));
 return X3, phi;
end function;




function Deg2InvariantEmbedding()

end function

function Deg3InvariantEmbedding()

end function



function C3Curve(field,data)
 if ([data[4],4*data[1]+data[4]^2] eq [0,0]) then 
 error "This is not valid data";
 end if;
P<f0,f1,f2,g0,g1>:= PolynomialRing(field,5);
PolP<x>:= PolynomialRing(P,1);
C := HyperellipticCurve(Polynomial(Evaluate([f0, f1, f2, -10*f0 - 5*f1 - 2*f2, 15*f0 + 5*f1 + f2, -6*f0 - f1, f0],data)),Polynomial(Evaluate([g0, g1, -3*g0 - g1, g0],data)));
return C;
end function;

Q0<s> := GF(9);
// Q0 := GF(7);
C3 := C3Curve(Q0, [1,0,-1,1,0]);
G, gs := AutomorphismGroup(C3);
aut3 := [gs(g): g in G | Order(g) eq 3][8];
fixedpts := [pt: pt in Points(C3) | aut3(pt) eq pt];
J := Jacobian(C3);
fpJ := [J!0, fixedpts[1]-fixedpts[2], fixedpts[2]-fixedpts[1]];
Jac := GeneralJacobianSurface(J);
fp := [Gamma(Jac,pt): pt in fpJ];
X3, emb := InvariantEmbedding(J);
// nX3 := newInvariantEmbedding(J);
WJ := Domain(emb);
pi := Inverse(JacIso(WJ, Jac));
[pi(pt): pt in fp];
fpX3 := [emb(pi(pt)): pt in fp];
WJ<k1, k2, k3, k4, v1, v2, v3, v4, v5, v6> := Domain(emb);
HSX3 := HilbertSeries(Ideal(DefiningPolynomials(X3)));
// HSnX3 := HilbertSeries(Ideal(DefiningPolynomials(nX3)));
Ser<t> := PowerSeriesRing(Rationals(),8);
Ser!HSX3;
//P<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv, q> := Ambient(nX3);
//nquot := map<nX3->X3 | [c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv]>;
P<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := Ambient(X3);
// Pullback(nquot, X3![0,0,0,0,0,0,1,0,0]);
P5<c11, c12, c22, b1, b2, kvw> := ProjectiveSpace(Q0,5);
P15<v1, v2, v3, v4, v5, v6, k11, k12, k13, k14, k22, k23, k24, k33, k34, k44> := Ambient(Jac);
nm := map<Jac->P5 | [k14 + k33 + k44,  k13 + k33 + k44,    k12 + k22 + 2*k24 + k33 + 2*k34 + 2*k44,    k11 + k23 + 2*k24 + k33 + 2*k34 + 2*k44,    v4 + k44,    v2 + 2*v3 + k22 + 2*k23]>;
Scheme(P5,MinimalBasis(Image(nm)));         
ideal< CoordinateRing(P5) | &cat[ DefiningEquations(Image(nm,Jac,d)) : d in  [1..2]] >;
proj2 := map<X3->P5 |[c1^2, c1*c2, c2^2, d1, d2, cvw]>;
K3 := Scheme(P5,MinimalBasis(Image(proj2)));         
proj2 := map<X3->K3 |[c1^2, c1*c2, c2^2, d1, d2, cvw]>;
P<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv,q> := Ambient(nX3);
nproj2 := map<nX3->K3 |[c1^2, c1*c2, c2^2, d1, d2, cvw]>;
IrreducibleComponents(ReducedSubscheme(BaseScheme(proj2)));
 [#[eqs: eqs in DefiningEquations(X3) |  (Degree(eqs) eq i)]: i in [1..7]];
spK3 := SingularPoints(K3);
[<pt,ADEtype(pt)>: pt in spK3[1..3]];
_,F,seq,dat := IsHypersurfaceSingularity(spK3[4],10);
P<a,b,c> := Parent(F);
_,_,typ := NormalFormOfHypersurfaceSingularity(F);
TjurinaNumberAnalyticHypersurface(dat);
IC := [IrreducibleComponents(ReducedSubscheme(Pullback(proj2,pt))): pt in spK3];
nIC := [IrreducibleComponents(ReducedSubscheme(Pullback(nproj2,pt))): pt in spK3];
IC;
l1:= IrreducibleComponents(Scheme(K3,[c11,c12,kvw]))[1];
l2:= IrreducibleComponents(Scheme(K3,[c11,c12,kvw]))[2];
B1, phi1 := Blowup(K3, spK3[1]: Ordinary := false); //Blowing up the D5 singularity
AS<[t]> := AmbientSpace(B1);
P4<x7,x8,x9,x10,x11> := ProjectiveSpace(Q0,4);
pr := map<B1->P4 | [t[7],t[8],t[9],t[10],t[11]]>;
SingularPoints(Image(pr));
IrreducibleComponents(ReducedSubscheme(SingularSubscheme(Image(pr))));
E := IrreducibleComponents(ReducedSubscheme(Pullback(phi1,spK3[1])))[1];
sPt1 := B1![0,0,0,s^2,1,0,0,0,1,0,0];
IsSingular(sPt1);
sPt2 := B1![0,0,0,s^2,1,0,1,0,0,0,0];
IsSingular(sPt2);
[<pt,ADEtype(pt)>: pt in [sPt1, sPt2]];
L1 := IrreducibleComponents(ReducedSubscheme(Scheme(B1,[t[1],t[2],t[4]+s^6*t[5],t[6]] )))[2];
Intersection(L1,E); // The intersection of l1 with the exceptional divisor is on the first point which is consistent with E_6
B2, phi2 := Blowup(B1,sPt1: Ordinary := false);
AS<[t]> := AmbientSpace(B2);
// E := IrreducibleComponents(ReducedSubscheme(Pullback(phi2,sPt1)))[1];


A<x1,x2,x3,x4,x5> := AffinePatch(K3, spK3[1]);
tau := map<K3->A | [c11/b2,c12/b2,c22/b2,b1/b2,kvw/b2]>;
P := Scheme(A,[x1,x2,x3,x4+s^6,x5]);
l1 := Scheme(A,[x1,x2,x4+s^6,x5]);
LB := LocalBlowUp(A,P);



[[Dimension(V): V in IC[i]]:i in [1..4]];
IsSingular(X3,X3![0,0,0,0,0,0,1,0,0]);
IsSingular(X3,X3![1,0,0,0,0,0,0,0,0]);
P12<c111, c112, c122, c222, c1d1, c1d2, c1vw, c2d1, c2d2, c2vw, kvvv, kwww, kwdv> := ProjectiveSpace(Q0, 12);
proj3 := map<X3->P12 | [c1^3, c1^2*c2, c1*c2^2, c2^3, c1*d1, c1*d2, c1*cvw, c2*d1, c2*d2, c2*cvw, cvvv, cwww, cwdv]>;
J := ideal< CoordinateRing(P12) | &cat[ DefiningEquations(Image(proj3,X3,d)) : d in  [1..3]] >;
MinimalBasis(J);
X33 := Scheme(P12, MinimalBasis(Image(proj3,X3,2)));
proj3 := map<X3->X33 | [c1^3, c1^2*c2, c1*c2^2, c2^3, c1*d1, c1*d2, c1*cvw, c2*d1, c2*d2, c2*cvw, cvvv, cwww, cwdv]>;
P1 := X33!proj3(X3![0,0,0,0,0,0,1,0,0]);
P2 := X33!proj3(X3![1,0,0,0,0,0,0,0,0]);
proj31 := map<X33->K3 | [c111, c112, c122, c1d1, c1d2, c1vw]>;
proj32 := map<X33->K3 | [c112, c122, c222, c2d1, c2d2, c2vw]>;
IrreducibleComponents(ReducedSubscheme(Intersection(BaseScheme(proj31),BaseScheme(proj32))));
IC3 := [IrreducibleComponents(ReducedSubscheme(Intersection(Pullback(proj31,pt),Pullback(proj32,pt)))): pt in spK3];
[[Dimension(V): V in IC3[i]]:i in [1..4]];
IrreducibleComponents(SingularSubscheme(IC3[3,1]));
_,F,seq,dat := IsHypersurfaceSingularity(P2,10);
P<a,b,c> := Parent(F);
_,_,typ := NormalFormOfHypersurfaceSingularity(F);
TjurinaNumberAnalyticHypersurface(dat);
IrreducibleComponents(SingularSubscheme(IC3[3,1]));
IsHypersurfaceSingularity(P1,10);
IrreducibleComponents(ReducedSubscheme(Intersection(Scheme(X33,[c111,c112,c1d1+s^6*c1d2,c1vw]),Scheme(X33,[c112,c122,c2d1+s^6*c2d2,c2vw]))));

F;
P := spK3[4];
P := K3![0,0,1,0,0,0];
B,phi := Blowup(K3,P: Ordinary := false);
AS<[t]> := AmbientSpace(B);
IsReduced(Pullback(phi,P));
IrreducibleComponents(Pullback(phi,P));



terms := func< f,j | [[t:t in Terms(f) | Degree(t) eq i] : i in [2..j]] >;
f := Evaluate(F,[a-b,b,c]);
terms(F,3);
K<r,s,t,u,v,w> := RationalFunctionField(GF(3^7),6);
R<a,b,c> := PolynomialRing(K,3);
F := 2*a^6 + 2*b^6 + 2*c^6 + a^4 + 2*a^2*b^2 + b^4 + 2*a^2*c^2 + 2*b^2*c^2 + c^4 + c^3 + b^2;
A<a,b,c> := AffineSpace(Q0,3);
F := 2*a^10 + a^8*b^2 + 2*a^6*b^4 + 2*a^4*b^6 + a^2*b^8 + 2*b^10 + a^8*c^2 + 
    a^6*b^2*c^2 + a^2*b^6*c^2 + b^8*c^2 + 2*a^6*c^4 + 2*b^6*c^4 + 2*a^4*c^6 + 
    a^2*b^2*c^6 + 2*b^4*c^6 + a^2*c^8 + b^2*c^8 + 2*c^10 + 2*a^8 + 2*a^6*b^2 + 
    2*a^2*b^6 + 2*b^8 + 2*a^6*c^2 + 2*b^6*c^2 + 2*a^2*c^6 + 2*b^2*c^6 + 2*c^8 + 
    2*a^6 + 2*b^6 + 2*c^6 + a^4 + 2*a^2*b^2 + b^4 + 2*a^2*c^2 + 2*b^2*c^2 + c^4 
    + c^3 + b^2;
sing := Scheme(A, F);
B,phi := Blowup(sing,sing![0,0,0]: Ordinary := false);
AS<[t]> := AmbientSpace(B);
P := sing![0,0,0];
IrreducibleComponents(Pullback(phi,P));
IrreducibleComponents(ReducedSubscheme(SingularSubscheme(B)));
P2 := B![0,0,0,0,0,1];
_,F := IsHypersurfaceSingularity(P2,10);
_,_,typ := NormalFormOfHypersurfaceSingularity(F);
P<a,b,c> := Parent(F);
F;
B2,phi2 := Blowup(B,P2: Ordinary := false);
AS<[t]> := AmbientSpace(B2);


// An example with good reduction at 3
 //Q0 := Rationals();
 //Q0 := GF(169);
 Q0<s> := GF(9);
 Q0pol<x> := PolynomialRing(Q0);
 P5<c1,c2,c3,c4,c5,c6> := ProjectiveSpace(Q0,5);
 K3 := Scheme(P5,[-258*c1^2 + 1024*c1*c2 - 1046*c2^2 - 232*c4*c5 - 616*c5^2 + 98*c3*c6 + 326*c4*c6 - 53*c6^2, -211*c1^2 + 752*c1*c2 - 793*c2^2 - 124*c4*c5 + 124*c5^2 + 67*c3*c6 + 145*c4*c6 - 
  212*c5*c6, -c4^2 + c3*c5]);
 SSK3 := SingularSubscheme(K3);  
 spK3 := SingularPoints(K3);
 spK3ss := PointsOverSplittingField(SSK3);
 [ADEtype(pt): pt in spK3ss]; 
 C3 := HyperellipticCurve(-(-1 + x + x^2 + 9*x^3 - 15*x^4 + 7*x^5 - x^6), 1 + x + x^3);
 Jac3 := GeneralJacobianSurface(C3);
 P15<v1, v2, v3, v4, v5, v6, k11, k12, k13, k14, k22, k23, k24, k33, k34, k44> := Ambient(Jac3);
 WJ := WeightedJacobian(Jacobian(C3));
 PW<k1, k2, k3, k4, v1, v2, v3, v4, v5, v6> := Ambient(WJ);
 //P5<c11, c12, c22, b1, b2, kvw>:= ProjectiveSpace(Q0,5);
 //proj2 := map<WJ->P5 |[36*k1^2 - 18*k1*k2 + (9*k2^2)/4 - 12*k1*k4 + 3*k2*k4 + k4^2, -6*k1^2 + (9*k1*k2)/2 - (3*k2^2)/4 - 6*k1*k3 + (3*k2*k3)/2 + k1*k4 - (k2*k4)/2 + k3*k4,  k1^2 - k1*k2 + k2^2/4 + 2*k1*k3 - k2*k3 + k3^2, -3*k1^2 - 9*k1*k3 - 2*k2*k3 - 7*k3^2 + k2*k4 - 4*v2 + 2*v4, -k1^2 + k1*k2 + k1*k3 - k2*k3 + k3^2 - 2*v2 + 2*v3,  k2^2/4 - k1*k3]>;
 //J := ideal< CoordinateRing(P5) | &cat[ DefiningEquations(Image(proj2,WJ,d)) : d in  [1..2]] >;
 //MinimalBasis(J);
 //K3 := Scheme(P5,MinimalBasis(J));

 PW<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := ProjectiveSpace(Q0, [1,1,2,2,2,3,3,3,3]);
 Wpi := map<WJ->PW | [k2, k4, k2^2 + k1*k3 + 2*v1, k1*k3 + k2*k4 + 2*v4, k1*k3, k1^3, k3^3, k1*(k2*k3 + 2*v2), k3*(k3^2 + 2*v3)]>;
 J := ideal< CoordinateRing(PW) | &cat[ DefiningEquations(Image(Wpi,WJ,d)) : d in  [1..6]] >;
 MinimalBasis(J);
 X3 := Scheme(PW,MinimalBasis(J));
 P<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := Ambient(X3);
 proj2 := map<X3->P5 |[c1^2, c1*c2, c2^2, d1, d2, cvw]>;
 K3 := Scheme(P5,MinimalBasis(Image(proj2)));




 Q0t<t> :=FunctionField(Q0);
 K3t<c11, c12, c22, b1, b2, kvw> := BaseExtend(K3,Q0t);
 P1t<x1,x2> := ProjectiveSpace(Q0t,1);
 Fib1t := map<K3t->P1t |[c11,c12]>;
 Fib2t := map<K3t->P1t |[c12,c22]>;
 coordt :=[1,t];
 G1t := Intersection(Pullback(Fib1t,P1t!coordt),Pullback(Fib2t,P1t!coordt));
 [<Dimension(v),ArithmeticGenus(v)>: v in IrreducibleComponents(ReducedSubscheme(G1t))];
 Cut := IrreducibleComponents(G1t)[1];
 linest:= IrreducibleComponents(Scheme(K3t,[kvw]));
 //linest:= IrreducibleComponents(Scheme(K3t,[b1-b2]));
 ptO := Cut!Points(IrreducibleComponents(Intersection(linest[1],Cut))[1])[1];
 D := Divisor(ptO);
 phi3t:= DivisorMap(3*D);
 E1, mapE1 := EllipticCurve(Image(phi3t), phi3t(ptO));
 E, mapE2 := MinimalModel(E1);
 mapE:=mapE1*mapE2;
 invmapE := Inverse(mapE);
 BadPlaces(E);
 LocalInformation(E);
 KodairaSymbols(E);
 G, mapG := TorsionSubgroup(E);
 P2 := mapG(G.1);
 Inverse(mapE)(P2);
 




 Q0<r> := ext<Rationals() | Polynomial([1,-1,1])>;
 // Q0 := Rationals();
 Q0<r> := GF(997^2);
 Q0pol<x> := PolynomialRing(Q0);
 C3 := HyperellipticCurve(-(-5/18 + (13*x)/6 - 3*x^2 + (71*x^3)/18 - (7*x^4)/3 + (3*x^5)/2 - (5*x^6)/18),1 + x + x^3);
 Jac3 := GeneralJacobianSurface(C3);
 P15<v1, v2, v3, v4, v5, v6, k11, k12, k13, k14, k22, k23, k24, k33, k34, k44> := Ambient(Jac3);
 P5<c1,c2,c3,c4,c5,c6> := ProjectiveSpace(Q0,5);
 K3 := Scheme(P5,[11*(2304*c1^2 - 3168*c1*c2 + 4752*c2^2 + 1344*c4*c5 - 6032*c5^2 - 1728*c3*c6 + 4464*c4*c6 - 195*c6^2), 
 -2168*c1^2 + 2916*c1*c2 - 4374*c2^2 - 2088*c4*c5 + 4524*c5^2 + 2016*c3*c6 - 4038*c4*c6 - 715*c5*c6, -c4^2 + c3*c5]);
 SSK3 := SingularSubscheme(K3);  
 spK3 := SingularPoints(K3);
 [ADEtype(pt): pt in spK3]; 
 WJ := WeightedJacobian(Jacobian(C3));
 P<k1, k2, k3, k4, v1, v2, v3, v4, v5, v6> := Ambient(WJ);
 pi5 := map <Jac3->K3 |[(-25*k11)/6 - (71*k13)/18 - 2*k23 + k24 - (3*k33)/2 - 4*v2 + 2*v4, -k11 + k12 + k13 - k23 + k33 - 2*v2 + 2*v3,  (4*k11)/9 - (19*k12)/27 + (4*k14)/3 + (361*k22)/1296 - (19*k24)/18 + k44, (2*k11)/3 - (31*k12)/36 + (2*k13)/3 + k14 + (19*k22)/72 - (19*k23)/36 - k24/2 + k34,  k11 - k12 + 2*k13 + k22/4 - k23 + k33, -4*k13 + k22]>;
 Wpi5 := map <WJ->K3|[(-25*k1^2)/6 - (71*k1*k3)/18 - 2*k2*k3 - (3*k3^2)/2 + k2*k4 - 4*v2 + 2*v4, -k1^2 + k1*k2 + k1*k3 - k2*k3 + k3^2 - 2*v2 + 2*v3,  (4*k1^2)/9 - (19*k1*k2)/27 + (361*k2^2)/1296 + (4*k1*k4)/3 - (19*k2*k4)/18 + k4^2, (2*k1^2)/3 - (31*k1*k2)/36 + (19*k2^2)/72 + (2*k1*k3)/3 - (19*k2*k3)/36 + k1*k4 -   (k2*k4)/2 + k3*k4, k1^2 - k1*k2 + k2^2/4 + 2*k1*k3 - k2*k3 + k3^2, k1^2 - k1*k2 + k2^2 - k1*k3 - k2*k3 + k3^2]>;
 LinesA1s := IrreducibleComponents(Scheme(K3,[(-36*c3)/169 + c5, (-6*c3)/13 + c4]));
 Intersection(LinesA1s[1],LinesA1s[2]);
 IC:=IrreducibleComponents(Pullback(pi5,spK3[1]));
 [Dimension(v): v in IC];
 [Degree(v): v in IC];
 rIC:=[ReducedSubscheme(c): c in IC];
 ArithmeticGenus(rIC[1]);
 IC2:=IrreducibleComponents(Pullback(pi5,spK3[2]));
 [Dimension(v): v in IC2];
 [Degree(v): v in IC2];
 rIC2:=[Points(c): c in IC2];
 IC3:=IrreducibleComponents(Pullback(pi5,spK3[3]));
 [Dimension(v): v in IC3];
 [Degree(v): v in IC3];
 rIC3:=[Points(c): c in IC3];
 rIC3;
 genpt := Scheme(K3,[c1+c2,c3+c4]);
 Dimension(genpt);
 Degree(genpt);
 ICgen:=IrreducibleComponents(Pullback(pi5,genpt));
 [Dimension(v): v in ICgen];
 [Degree(v): v in ICgen];
 PW<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := ProjectiveSpace(Q0, [1,1,2,2,2,3,3,3,3]);
 Wpi := map<Jac3->PW | [(2*k1)/3 - (19*k2)/36 + k4, k1 - k2/2 + k3,(-k1^2 - k2*k3 - 2*v2 - 3*(k1^2 + k2*k3 + 2*v2) + 2*((-13*k1^2)/6 + f3*k1*k3 + f5*k3^2 + k2*k4 + 2*v4))/2, -k1^2 + k1*k2 - k2*k3 + k3*(k1 + k3) - 2*v2 + 2*v3]>;
 Kum3<k1,k2,k3,k4> := GeneralKummerSurface(C3);
 P3<y1,y2,y3,y4> := ProjectiveSpace(Q0,3);
 pi3 := map<Kum3->P3|[((-8*k1)/3 + (19*k2)/9 - 4*k4)^2/16, -1/8*((2*k1 - k2 + 2*k3)*((-8*k1)/3 + (19*k2)/9 - 4*k4)), (2*k1 - k2 + 2*k3)^2/4,  k1^2 - k1*k2 + k2^2 - k1*k3 - k2*k3 + k3^2,k1^2 - k1*k2 + k2^2 - k1*k3 - k2*k3 + k3^2,(k1 + k2 - 2*k3)*(2*k1 - k2 - k3)*(k1 - 2*k2 + k3),(k1 - k2)*(k1 - k3)*(k2 - k3)]>;
 K6 := Scheme(P3,MinimalBasis(Image(pi3))); //Does not provide an embedding
 P7<x1,x2,x3,x4,x5,x6,x7,x8>:= ProjectiveSpace(Q0,7);
 pi7 := map<Kum3->P7|[6912*k1^3 - 16416*k1^2*k2 + 6859*k2^3 + 51984*k1^2*k3 - 41154*k1*k2*k3 + 31104*k1^2*k4 - 49248*k1*k2*k4 + 77976*k1*k3*k4 + 46656*k1*k4^2 - 36936*k2*k4^2 + 23328*k4^3,  576*k1^3 - 1200*k1^2*k2 + 361*k2^3 + 2400*k1^2*k3 - 1634*k1*k2*k3 + 1728*k1^2*k4 - 2232*k1*k2*k4 + 4464*k1*k3*k4 - 1368*k2*k3*k4 + 1296*k1*k4^2 - 648*k2*k4^2 +   1296*k3*k4^2, 12*k1^3 - 31*k1^2*k2 + 19*k1*k2^2 + 43*k1^2*k3 - 38*k1*k2*k3 + 18*k1^2*k4 - 18*k1*k2*k4 + 54*k1*k3*k4 - 18*k2*k3*k4 + 18*k3^2*k4,  2*k1^3 - 6*k1^2*k2 + 6*k1*k2^2 - k2^3 + 6*k1^2*k3 - 6*k1*k2*k3 + 2*k3^3, -2*k1^2*k2 + 4*k1*k2^2 - k2^3 + 4*k1^2*k3 - 6*k1*k2*k3 + 2*k2*k3^2,  24*k1*k2^2 - 19*k2^3 - 96*k1^2*k3 + 76*k1*k2*k3 + 36*k2^2*k4 - 144*k1*k3*k4, 2*k1*k2^2 - k2^3 - 4*k1*k2*k3 + 2*k2^2*k3, k1^2*k3 - k1*k2*k3 + k1*k3^2]>;
 X6 := Scheme(P7,MinimalBasis(Image(pi7)));
 PW<c1,c2,cvw,cvvv,cwww> := ProjectiveSpace(Q0, [1,1,2,3,3]);
 Wpi := map<Kum3->PW | [(2*k1)/3 - (19*k2)/36 + k4, k1 - k2/2 + k3,k1^2 - k1*k2 + k2^2 - k1*k3 - k2*k3 + k3^2,(k1 + k2 - 2*k3)*(2*k1 - k2 - k3)*(k1 - 2*k2 + k3),(k1 - k2)*(k1 - k3)*(k2 - k3)]>;
 Image(Wpi);
 Wpi7 := map<WJ->X6|[6912*k1^3 - 16416*k1^2*k2 + 6859*k2^3 + 51984*k1^2*k3 - 41154*k1*k2*k3 + 31104*k1^2*k4 - 49248*k1*k2*k4 + 77976*k1*k3*k4 + 46656*k1*k4^2 - 36936*k2*k4^2 + 23328*k4^3,  576*k1^3 - 1200*k1^2*k2 + 361*k2^3 + 2400*k1^2*k3 - 1634*k1*k2*k3 + 1728*k1^2*k4 - 2232*k1*k2*k4 + 4464*k1*k3*k4 - 1368*k2*k3*k4 + 1296*k1*k4^2 - 648*k2*k4^2 +   1296*k3*k4^2, 12*k1^3 - 31*k1^2*k2 + 19*k1*k2^2 + 43*k1^2*k3 - 38*k1*k2*k3 + 18*k1^2*k4 - 18*k1*k2*k4 + 54*k1*k3*k4 - 18*k2*k3*k4 + 18*k3^2*k4,  2*k1^3 - 6*k1^2*k2 + 6*k1*k2^2 - k2^3 + 6*k1^2*k3 - 6*k1*k2*k3 + 2*k3^3, -2*k1^2*k2 + 4*k1*k2^2 - k2^3 + 4*k1^2*k3 - 6*k1*k2*k3 + 2*k2*k3^2,  24*k1*k2^2 - 19*k2^3 - 96*k1^2*k3 + 76*k1*k2*k3 + 36*k2^2*k4 - 144*k1*k3*k4, 2*k1*k2^2 - k2^3 - 4*k1*k2*k3 + 2*k2^2*k3, k1^2*k3 - k1*k2*k3 + k1*k3^2]>;
 P13<y1,y2,y3,y4,y5,y6,x1,x2,x3,x4,x5,x6,x7,x8>:= ProjectiveSpace(Q0,13);
 X3 := Scheme(P13,[-3168*x1^2 - 8568*x1*x2 + 8208*x2^2 + 2592*x1*x3 - 1944*x2*x3 - 2349*x3^2 + 1482*x1*x4 + 6936*x2*x4 - 18*x1*x5 - 264*x2*x5 + 54*x3*x5 - 38*x4*x5 + x5^2,  -19*x1^2 + 104*x1*x2 + 384*x2^2 + 19*x1*x3 - 152*x2*x3 + x1*x5 + 16*x2*x5 - 3*x3*x5 + 2*x4*x5 - 4*x1*x6 + 32*x2*x6,  -894*x1^2 - 5940*x1*x2 - 17604*x2^2 + 3522*x1*x3 + 11826*x2*x3 - 2052*x3^2 - 2788*x1*x4 - 8874*x2*x4 + 3079*x3*x4 - 867*x4^2 + 26*x1*x5 + 288*x2*x5 - 59*x3*x5 + 33*x4*x5 -   102*x1*x6 + 108*x3*x6 - 76*x4*x6 + 2*x5*x6, 6383*x1^2 + 35216*x1*x2 + 62208*x2^2 - 16680*x1*x3 - 45360*x2*x3 + 10557*x3^2 + 9670*x1*x4 + 28560*x2*x4 - 12316*x3*x4 +   3468*x4^2 - 86*x1*x5 - 488*x2*x5 + 144*x3*x5 - 94*x4*x5 + 408*x1*x6 - 432*x3*x6 + 304*x4*x6 - 2*x1*x7 + 16*x2*x7,  -15371*x1^2 - 77688*x1*x2 - 155016*x2^2 + 41221*x1*x3 + 114372*x2*x3 - 25940*x3^2 - 24194*x1*x4 - 74868*x2*x4 + 30790*x3*x4 - 8670*x4^2 + 267*x1*x5 + 2496*x2*x5 -   573*x3*x5 + 378*x4*x5 - 1420*x1*x6 + 1384*x3*x6 - 760*x4*x6 - 32*x6^2 + 6*x1*x7 - 6*x3*x7 + 4*x4*x7,  -25835*x1^2 - 230976*x1*x2 - 1204344*x2^2 + 265933*x1*x3 + 886140*x2*x3 - 109082*x3^2 - 243556*x1*x4 - 826812*x2*x4 + 211528*x3*x4 - 50286*x4^2 + 859*x1*x5 + 39096*x2*x5 -   5487*x3*x5 + 3726*x4*x5 - 25792*x1*x6 + 23380*x3*x6 - 11344*x4*x6 - 608*x6^2 + 12*x1*x7 - 6*x3*x7 + 2*x5*x7,  146219*x1^2 + 556560*x1*x2 + 892440*x2^2 - 364429*x1*x3 - 812268*x2*x3 + 198344*x3^2 + 187240*x1*x4 + 563244*x2*x4 - 211528*x3*x4 + 50286*x4^2 - 897*x1*x5 - 20400*x2*x5 +   3435*x3*x5 - 2282*x4*x5 + 25792*x1*x6 - 23380*x3*x6 + 11344*x4*x6 + 608*x6^2 - 12*x1*x7 + 6*x3*x7 - 2*x1*x8 + 16*x2*x8,  -328943*x1^2 - 1370592*x1*x2 - 2749536*x2^2 + 903309*x1*x3 + 2251152*x2*x3 - 489298*x3^2 - 528348*x1*x4 - 1629552*x2*x4 + 583164*x3*x4 - 145656*x4^2 + 3615*x1*x5 +   75264*x2*x5 - 13017*x3*x5 + 8674*x4*x5 - 58332*x1*x6 + 52376*x3*x6 - 26640*x4*x6 - 1216*x6^2 - 26*x1*x7 + 26*x3*x7 - 8*x6*x7 + 6*x1*x8 - 6*x3*x8 + 4*x4*x8,  -6206519*x1^2 - 32958864*x1*x2 - 86081040*x2^2 + 20093541*x1*x3 + 60693624*x2*x3 - 11477848*x3^2 - 14155302*x1*x4 - 44605416*x2*x4 + 15867630*x3*x4 - 4421700*x4^2 +   79791*x1*x5 + 1582224*x2*x5 - 278319*x3*x5 + 186046*x4*x5 - 820416*x1*x6 + 768164*x3*x6 - 422280*x4*x6 - 16912*x6^2 - 1358*x1*x7 + 926*x3*x7 - 152*x6*x7 + 12*x1*x8 -   6*x3*x8 + 2*x5*x8, 5087270*x1^2 + 31365000*x1*x2 + 63065232*x2^2 - 9831077*x1*x3 - 35531784*x2*x3 + 8848763*x3^2 + 6010784*x1*x4 + 15954840*x2*x4 - 10591394*x3*x4 +   3752376*x4^2 - 109155*x1*x5 + 394920*x2*x5 + 54015*x3*x5 - 32160*x4*x5 - 1007416*x1*x6 + 821420*x3*x6 - 246760*x4*x6 - 36592*x6^2 + 970*x1*x7 - 1454*x3*x7 + 152*x6*x7 -   4*x7^2 + 88*x1*x8 - 70*x3*x8 + 16*x6*x8, 67*x1*y1 + 40*x2*y1 - 30*x3*y1 + x5*y1 - 57*x1*y2 + 132*x2*y2 + 57*x3*y2 + 324*x2*y3 - 456*x2*y4,  x1*y1 + 4*x2*y1 - 3*x3*y1 + 2*x4*y1 + 12*x2*y2 - 3*x1*y3 + 12*x2*y3 + 3*x3*y3 - 24*x2*y4, x5*y3 - x1*y5 + 8*x2*y5,  9*x1*y1 + 36*x2*y1 + 27*x3*y1 - 66*x1*y2 - 120*x2*y2 + 27*x3*y2 + x5*y2 + 39*x1*y3 + 336*x2*y3 - 216*x2*y4 - x1*y5 + 8*x2*y5,  -108*x2*y1 + 19*x1*y2 + 412*x2*y2 - 57*x3*y2 + 38*x4*y2 - 67*x1*y3 + 380*x2*y3 + 30*x3*y3 - 216*x2*y4 - x1*y5 + 8*x2*y5,  -19*x1*y1 - 32*x2*y1 + 19*x3*y1 + 19*x1*y2 + 32*x2*y2 - 19*x3*y2 + 66*x1*y3 + 228*x2*y3 - 27*x3*y3 - 38*x1*y4 - 64*x2*y4 + 38*x3*y4 - x1*y5 + 8*x2*y5,  -779*x1*y1 + 220*x2*y1 - 570*x3*y1 + 2109*x1*y2 - 2880*x2*y2 + 1083*x3*y2 + 45*x1*y3 - 7524*x2*y3 + 2601*x3*y3 - 4218*x1*y4 + 2112*x2*y4 - 24*x1*y5 - 264*x2*y5 +   114*x3*y5 - 57*x1*y6, 1501*x1*y1 - 3868*x2*y1 - 152*x3*y1 + 152*x6*y1 - 2109*x1*y2 + 11544*x2*y2 - 1083*x3*y2 - 45*x1*y3 + 7524*x2*y3 - 2601*x3*y3 + 4218*x1*y4 -   2112*x2*y4 - 90*x1*y5 + 720*x2*y5 + 57*x1*y6 - 456*x2*y6, -70*x1*y1 + 92*x2*y1 + 57*x3*y1 + 10*x1*y2 - 404*x2*y2 - 30*x3*y2 + 39*x1*y3 + 12*x2*y3 + 240*x2*y4 + 2*x5*y4 -   x1*y5 + 8*x2*y5 - x1*y6 + 8*x2*y6, -2014*x1*y1 - 1852*x2*y1 + 1083*x3*y1 + 2147*x1*y2 + 788*x2*y2 - 1083*x3*y2 + 1592*x1*y3 + 9044*x2*y3 + 1836*x3*y3 - 1026*x4*y3 -   2850*x1*y4 + 4200*x2*y4 + 1444*x4*y4 - 46*x1*y5 + 368*x2*y5 - 19*x1*y6 + 152*x2*y6, 76627*x1*y1 + 29476*x2*y1 + 11058*x3*y1 - 178239*x1*y2 + 117888*x2*y2 - 29241*x3*y2 +   38229*x1*y3 + 489060*x2*y3 - 169065*x3*y3 + 228*x7*y3 + 274170*x1*y4 - 137280*x2*y4 - 5964*x1*y5 - 54432*x2*y5 - 912*x6*y5 + 5871*x1*y6 + 3648*x2*y6 - 2166*x3*y6,  -13224*x1*y1 + 76080*x2*y1 + 1083*x3*y1 + 5263*x1*y2 - 227548*x2*y2 + 23465*x3*y2 + 76158*x1*y3 - 28728*x2*y3 - 30456*x3*y3 + 32946*x4*y3 - 37962*x1*y4 + 19008*x2*y4 -   2040*x1*y5 - 1008*x2*y5 + 1444*x4*y5 + 209*x1*y6 + 1216*x2*y6 - 722*x3*y6, 15675*x1*y1 + 12612*x2*y1 - 3534*x3*y1 - 32699*x1*y2 + 2816*x2*y2 + 3971*x3*y2 - 405*x1*y3 +   67716*x2*y3 - 23409*x3*y3 + 37962*x1*y4 - 19008*x2*y4 - 2292*x1*y5 - 6288*x2*y5 + 38*x5*y5 + 1235*x1*y6 + 1216*x2*y6 - 722*x3*y6,  15675*x1*y1 + 12612*x2*y1 - 3534*x3*y1 - 32699*x1*y2 + 2816*x2*y2 + 3971*x3*y2 - 405*x1*y3 + 67716*x2*y3 - 23409*x3*y3 + 2*x8*y3 + 37962*x1*y4 - 19008*x2*y4 - 1570*x1*y5 -   14952*x2*y5 - 2*x7*y5 + 1235*x1*y6 + 1216*x2*y6 - 722*x3*y6, 1577*x1*y1 + 71900*x2*y1 + 11913*x3*y1 - 34808*x1*y2 - 172828*x2*y2 + 2888*x3*y2 + 89021*x1*y3 +   148884*x2*y3 - 93593*x3*y3 + 32946*x4*y3 + 2888*x6*y3 + 42180*x1*y4 - 21120*x2*y4 - 2306*x1*y5 - 7544*x2*y5 + 1292*x1*y6 + 1216*x2*y6 - 722*x3*y6,  1235*x1*y1 - 2044*x2*y1 + 114*x3*y1 - 2109*x1*y2 + 7592*x2*y2 - 1083*x3*y2 + 152*x6*y2 - 45*x1*y3 + 7524*x2*y3 - 2601*x3*y3 + 4218*x1*y4 - 2112*x2*y4 - 242*x1*y5 +   568*x2*y5 + 95*x1*y6 - 152*x2*y6 - 38*x3*y6, -7*x1*y1 + 28*x2*y1 + 7*x3*y1 - 84*x2*y2 + 97*x1*y3 - 8*x2*y3 - 21*x3*y3 - 40*x2*y4 + 8*x6*y4 - 4*x1*y5 - 4*x2*y5 - 8*x2*y6 +   2*x3*y6 - 2*x4*y6, 1591231*x1*y1 - 1101220*x2*y1 + 1788033*x3*y1 + 1444*x7*y1 - 5986140*x1*y2 + 1722900*x2*y2 - 1329924*x3*y2 + 4757079*x1*y3 + 30349308*x2*y3 -   10571481*x3*y3 + 2668626*x4*y3 + 9667656*x1*y4 - 2206848*x2*y4 - 361320*x1*y5 + 785208*x2*y5 + 230280*x1*y6 + 81168*x2*y6 - 181944*x3*y6 + 82308*x4*y6,  7363621*x1*y1 + 11774548*x2*y1 + 3777504*x3*y1 - 24105813*x1*y2 - 18033528*x2*y2 - 1374327*x3*y2 + 4332*x7*y2 + 16992243*x1*y3 + 86058828*x2*y3 - 35255673*x3*y3 +   8302392*x4*y3 + 34591818*x1*y4 - 13577664*x2*y4 - 1828842*x1*y5 + 57888*x2*y5 - 17328*x6*y5 + 917985*x1*y6 + 635664*x2*y6 - 567492*x3*y6 + 116964*x4*y6,  -318744*x1*y1 + 146952*x2*y1 - 243675*x3*y1 + 691087*x1*y2 - 2783468*x2*y2 + 765681*x3*y2 - 1049836*x1*y3 - 8379304*x2*y3 + 3573246*x3*y3 - 889542*x4*y3 - 3302694*x1*y4 +   1653696*x2*y4 + 30456*x1*y5 - 243648*x2*y5 - 29469*x1*y6 - 58824*x2*y6 + 38988*x3*y6 - 27436*x4*y6 + 722*x5*y6,  1582529*x1*y1 + 1908140*x2*y1 - 1837851*x3*y1 - 969228*x1*y2 - 3686844*x2*y2 + 2291628*x3*y2 - 998799*x1*y3 - 20768748*x2*y3 + 4627053*x3*y3 - 1186056*x4*y3 -   3060672*x1*y4 - 1636872*x2*y4 + 4332*x7*y4 - 105654*x1*y5 - 1801620*x2*y5 - 8664*x6*y5 + 5472*x1*y6 - 472644*x2*y6 + 129960*x3*y6 - 106134*x4*y6 - 8664*x6*y6,  113138103*x1*y1 + 176093628*x2*y1 + 60316602*x3*y1 - 377860999*x1*y2 - 254791856*x2*y2 - 17962277*x3*y2 + 1444*x8*y2 + 231525441*x1*y3 + 1317343644*x2*y3 -   493510077*x3*y3 + 107931096*x4*y3 + 506754738*x1*y4 - 267305616*x2*y4 - 26994572*x1*y5 + 4038024*x2*y5 - 467856*x6*y5 - 1444*x7*y5 + 13431955*x1*y6 + 10501832*x2*y6 -   8104450*x3*y6 + 1520532*x4*y6 + 77976*x6*y6, 94075821*x1*y1 - 1654956*x2*y1 + 61584795*x3*y1 + 1444*x8*y1 - 301281252*x1*y2 + 101953068*x2*y2 - 39356220*x3*y2 +   184746321*x1*y3 + 1407336156*x2*y3 - 457357887*x3*y3 + 104076414*x4*y3 + 450128088*x1*y4 - 254027376*x2*y4 - 18778704*x1*y5 + 32719800*x2*y5 - 233928*x6*y5 +   11532468*x1*y6 + 7939416*x2*y6 - 8659668*x3*y6 + 3210012*x4*y6 + 164616*x6*y6, 9092735*x1*y1 + 44652428*x2*y1 - 110841801*x3*y1 + 108756456*x1*y2 - 308292324*x2*y2 +   143779080*x3*y2 - 191931945*x1*y3 - 1973201604*x2*y3 + 588177945*x3*y3 - 151518654*x4*y3 - 515836092*x1*y4 + 214618368*x2*y4 + 8664*x8*y4 + 2128710*x1*y5 -   99164400*x2*y5 - 372552*x6*y5 - 4332*x7*y5 - 5248788*x1*y6 - 26940480*x2*y6 + 10342650*x3*y6 - 7455372*x4*y6 - 259920*x6*y6 - 4332*x7*y6,  12*y1*y3 - 19*y2*y3 - y1*y5 + y2*y5 - 2*y4*y5 + y3*y6]);
 pi13 := map<WJ ->X3 | [(-38*k1^3)/3 - 3*k1^2*k2 - k1*k2^2 - (214*k1^2*k3)/9 - k1*k2*k3 + 3*k2^2*k3 - 8*k1*k3^2 - 3*k2*k3^2 + 4*k1*k2*k4 - 4*k1*v1 + 4*k3*v1 - 4*k1*v2 + 2*k2*v2 - 8*k3*v2 -   8*k1*v3 + 2*k2*v3 + 8*k1*v4, (-22*k1^2*k2)/3 - k1*k2^2 - 2*k2^3 - 2*k1^2*k3 - (98*k1*k2*k3)/9 + 3*k2^2*k3 - 6*k2*k3^2 + 2*k2^2*k4 - 4*k2*v1 + 4*k3*v1 + 2*k2*v2 - 4*k3*v2 -   4*k1*v3 - 2*k2*v3 + 4*k2*v4, -2*k1^3 + 3*k1^2*k2 - k1*k2^2 - k1*k2*k3 + k2^2*k3 + 4*k1*k3^2 - 3*k2*k3^2 + 2*k3^3 - 4*k1*v2 + 2*k2*v2 - 4*k3*v2 + 4*k1*v3 - 2*k2*v3 +   4*k3*v3, -2*k1^3 + 2*k1^2*k2 - k2^3 - (16*k1^2*k3)/3 - 5*k1*k2*k3 + 2*k2^2*k3 - (62*k1*k3^2)/9 - 4*k2*k3^2 - 3*k3^3 + 2*k2*k3*k4 + 2*k1*v1 - 2*k2*v1 - 6*k1*v2 + 4*k2*v2 -   6*k3*v2 + 2*k1*v3 - 2*k2*v3 + 4*k3*v4, -24*k1^3 + 43*k1^2*k2 - 19*k1*k2^2 + 24*k1^2*k3 - 43*k1*k2*k3 + 19*k2^2*k3 + 24*k1*k3^2 - 19*k2*k3^2 - 36*k1^2*k4 + 36*k1*k2*k4 +   36*k1*k3*k4 - 36*k2*k3*k4 + 36*k3^2*k4 - 48*k1*v2 + 38*k2*v2 - 72*k4*v2 + 48*k1*v3 - 38*k2*v3 + 72*k4*v3,  -48*k1^3 + 55*k1^2*k2 - 7*k1*k2^2 - 38*k2^3 + 58*k1^2*k3 - 141*k1*k2*k3 + 97*k2^2*k3 + 24*k1*k3^2 - 21*k2*k3^2 - 300*k1^2*k4 - 284*k1*k3*k4 - 144*k2*k3*k4 - 108*k3^2*k4 +   72*k2*k4^2 + 48*k1*v1 - 76*k2*v1 + 28*k3*v1 - 144*k1*v2 + 166*k2*v2 + 20*k3*v2 - 288*k4*v2 + 20*k1*v3 - 62*k2*v3 + 144*k4*v4, 2*k1*k2^2 - k2^3 - 4*k1*k2*k3 + 2*k2^2*k3,  k1^2*k3 - k1*k2*k3 + k1*k3^2, -2*k1^2*k2 + 4*k1*k2^2 - k2^3 + 4*k1^2*k3 - 6*k1*k2*k3 + 2*k2*k3^2, 2*k1^3 - 6*k1^2*k2 + 6*k1*k2^2 - k2^3 + 6*k1^2*k3 - 6*k1*k2*k3 + 2*k3^3,  24*k1*k2^2 - 19*k2^3 - 96*k1^2*k3 + 76*k1*k2*k3 + 36*k2^2*k4 - 144*k1*k3*k4, 12*k1^3 - 31*k1^2*k2 + 19*k1*k2^2 + 43*k1^2*k3 - 38*k1*k2*k3 + 18*k1^2*k4 - 18*k1*k2*k4 +   54*k1*k3*k4 - 18*k2*k3*k4 + 18*k3^2*k4, 576*k1^3 - 1200*k1^2*k2 + 361*k2^3 + 2400*k1^2*k3 - 1634*k1*k2*k3 + 1728*k1^2*k4 - 2232*k1*k2*k4 + 4464*k1*k3*k4 - 1368*k2*k3*k4 +   1296*k1*k4^2 - 648*k2*k4^2 + 1296*k3*k4^2, 6912*k1^3 - 16416*k1^2*k2 + 6859*k2^3 + 51984*k1^2*k3 - 41154*k1*k2*k3 + 31104*k1^2*k4 - 49248*k1*k2*k4 + 77976*k1*k3*k4 +   46656*k1*k4^2 - 36936*k2*k4^2 + 23328*k4^3]>;
 Jpi7 := map<Jac3->X6|[6912*k11^2 - 16416*k11*k12 + 51984*k11*k13 + 31104*k11*k14 + 6859*k12*k22 - 41154*k11*k23 - 49248*k11*k24 + 77976*k11*k34 + 46656*k11*k44 - 36936*k12*k44 + 23328*k14*k44,  576*k11^2 - 1200*k11*k12 + 2400*k11*k13 + 1728*k11*k14 + 361*k12*k22 - 1634*k11*k23 - 2232*k11*k24 + 4464*k11*k34 - 1368*k12*k34 + 1296*k11*k44 - 648*k12*k44 +   1296*k13*k44, 12*k11^2 - 31*k11*k12 + 43*k11*k13 + 18*k11*k14 + 19*k11*k22 - 38*k11*k23 - 18*k11*k24 + 54*k11*k34 - 18*k12*k34 + 18*k13*k34,  2*k11^2 - 6*k11*k12 + 6*k11*k13 + 6*k11*k22 - k12*k22 - 6*k11*k23 + 2*k13*k33, -2*k11*k12 + 4*k11*k13 + 4*k11*k22 - k12*k22 - 6*k11*k23 + 2*k12*k33,  -96*k11*k13 + 24*k11*k22 - 19*k12*k22 + 76*k11*k23 + 36*k12*k24 - 144*k11*k34, 2*k11*k22 - k12*k22 - 4*k11*k23 + 2*k12*k23, k11*k13 - k11*k23 + k11*k33]>;
 SSX6 :=Scheme(X6,MinimalBasis(SingularSubscheme(X6)));
 spX6 := SingularPoints(X6);
 pi7 := map<Kum3->X6|[6912*k1^3 - 16416*k1^2*k2 + 6859*k2^3 + 51984*k1^2*k3 - 41154*k1*k2*k3 + 31104*k1^2*k4 - 49248*k1*k2*k4 + 77976*k1*k3*k4 + 46656*k1*k4^2 - 36936*k2*k4^2 + 23328*k4^3,  576*k1^3 - 1200*k1^2*k2 + 361*k2^3 + 2400*k1^2*k3 - 1634*k1*k2*k3 + 1728*k1^2*k4 - 2232*k1*k2*k4 + 4464*k1*k3*k4 - 1368*k2*k3*k4 + 1296*k1*k4^2 - 648*k2*k4^2 +   1296*k3*k4^2, 12*k1^3 - 31*k1^2*k2 + 19*k1*k2^2 + 43*k1^2*k3 - 38*k1*k2*k3 + 18*k1^2*k4 - 18*k1*k2*k4 + 54*k1*k3*k4 - 18*k2*k3*k4 + 18*k3^2*k4,  2*k1^3 - 6*k1^2*k2 + 6*k1*k2^2 - k2^3 + 6*k1^2*k3 - 6*k1*k2*k3 + 2*k3^3, -2*k1^2*k2 + 4*k1*k2^2 - k2^3 + 4*k1^2*k3 - 6*k1*k2*k3 + 2*k2*k3^2,  24*k1*k2^2 - 19*k2^3 - 96*k1^2*k3 + 76*k1*k2*k3 + 36*k2^2*k4 - 144*k1*k3*k4, 2*k1*k2^2 - k2^3 - 4*k1*k2*k3 + 2*k2^2*k3, k1^2*k3 - k1*k2*k3 + k1*k3^2]>;
 P4<k1,k2,k3,k4,k5>:=ProjectiveSpace(Q0,4);
 phi := map<K3->P4|[c1, c2, c3 - 169/144*c6, c4 - 13/24*c6, c5 - 1/4*c6]>;
 // X := Image(phi);
 X := Surface(P4,[k2^3*k3 - 13/3*k2^3*k4 + 32/3*k1*k3*k4^2 - 32/27*k1*k4^3 + 32/9*k2*k4^3 + 169/36*k2^3*k5 - 32/3*k1*k3^2*k5 - 40/9*k1*k3*k4*k5 - 4*k2*k3*k4*k5 - 208/9*k1*k4^2*k5 - 52/9*k2*k4^2*k5 + 1070/27*k1*k3*k5^2 + 50/9*k2*k3*k5^2 + 208/27*k1*k4*k5^2 + 65/9*k2*k4*k5^2 - 53911/1458*k1*k5^3 - 4901/486*k2*k5^3, k1*k2*k3 - 3/2*k2^2*k3 - 13/3*k1*k2*k4 + 13/2*k2^2*k4 - 24*k3*k4^2 + 8/3*k4^3 + 169/36*k1*k2*k5 - 169/24*k2^2*k5 + 24*k3^2*k5 + 10*k3*k4*k5 + 52*k4^2*k5 - 535/6*k3*k5^2 - 52/3*k4*k5^2 + 53911/648*k5^3, k1^2 - 27/19*k1*k2 + 81/38*k2^2 - 36/19*k4^2 + 36/19*k3*k5 - 261/76*k5^2]);
 KodairaEnriquesType(X);
 [ADEtype(pt): pt in SingularPoints(X)];
 SSX := SingularPoints(X);
 pi4 := map <Jac3->X |[(-75*k11 - 71*k13 - 36*k23 + 18*k24 - 27*k33 - 72*v2 + 36*v4)/18, -k11 + k12 + k13 - k23 + k33 - 2*v2 + 2*v3,  (144*k11 - 228*k12 + 1521*k13 + 432*k14 - 290*k22 - 342*k24 + 324*k44)/324, (24*k11 - 31*k12 + 102*k13 + 36*k14 - 10*k22 - 19*k23 - 18*k24 + 36*k34)/36,  k11 - k12 + 3*k13 - k23 + k33]>;











 P4<k1,k2,k3,k4,k5>:=ProjectiveSpace(Q0,4);
 phi := map<K3->P4|[c1, c2, 4*c3-c6, 4*c4-c6,4*c5-c6]>;
 X := Image(phi);
 [ADEtype(pt): pt in SingularPoints(X)];



 K := GF(3);
 // Q0 := Rationals();
 P5<c1,c2,c3,c4,c5,c6> := ProjectiveSpace(K,5);
 K3 := Scheme(P5,[-258*c1^2 + 1024*c1*c2 - 1046*c2^2 - 232*c4*c5 - 616*c5^2 + 98*c3*c6 + 326*c4*c6 - 53*c6^2, -211*c1^2 + 752*c1*c2 - 793*c2^2 - 124*c4*c5 + 124*c5^2 + 67*c3*c6 + 145*c4*c6 - 
  212*c5*c6, -c4^2 + c3*c5]);
 SSK3 := SingularSubscheme(K3);  
 spK3 := SingularPoints(K3);
 [ADEtype(pt): pt in spK3[1..3]];
 [TjurinaP5(pt): pt in spK3]; 
 P4<k1,k2,k3,k4,k5>:=ProjectiveSpace(K,4);
 X := Surface(P4,MinimalBasis([18256*k2^3*k3 - 36512*k2^3*k4 - 416*k1*k3*k4^2 + 320*k2*k3*k4^2 - 3328*k1*k4^3 + 7968*k2*k4^3 + 18256*k2^3*k5 + 416*k1*k3^2*k5 - 320*k2*k3^2*k5 + 3848*k1*k3*k4*k5 -   9044*k2*k3*k4*k5 + 2704*k1*k4^2*k5 - 6136*k2*k4^2*k5 + 1222*k1*k3*k5^2 - 6179*k2*k3*k5^2 - 9412*k1*k4*k5^2 + 27858*k2*k4*k5^2 + 4966*k1*k5^3 - 14467*k2*k5^3,  416*k1*k2*k3 - 304*k2^2*k3 - 832*k1*k2*k4 + 608*k2^2*k4 + 16*k3*k4^2 + 128*k4^3 + 416*k1*k2*k5 - 304*k2^2*k5 - 16*k3^2*k5 - 148*k3*k4*k5 - 104*k4^2*k5 - 47*k3*k5^2 +   362*k4*k5^2 - 191*k5^3, 80*k1^2 - 224*k1*k2 + 256*k2^2 - 4*k4^2 + 4*k3*k5 - 31*k5^2]));
 spX := SingularPoints(X);
 [TjurinaP5(pt): pt in spX]; 

 Q0<r> := GF(121);
 Q0pol<x> := PolynomialRing(Q0);
 C3 := HyperellipticCurve(1 + x^6, x^3);
 Jac3 := GeneralJacobianSurface(C3);
 P15<v1, v2, v3, v4, v5, v6, k11, k12, k13, k14, k22, k23, k24, k33, k34, k44> := Ambient(Jac3);
 WJ := WeightedJacobian(Jacobian(C3));
 P<k1, k2, k3, k4, v1, v2, v3, v4, v5, v6> := Ambient(WJ);
 PW<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := ProjectiveSpace(Q0, [1,1,2,2,2,3,3,3,3]);
 Wpi := map<WJ->PW | [k2, k4, k2^2 + k1*k3 + 2*v1, k2*k4 + 2*v4, k1*k3, k1^3, k3^3, k1*(k2*k3 + 2*v2), k3*(k3^2 + 2*v3)]>;
 J := ideal< CoordinateRing(PW) | &cat[ DefiningEquations(Image(Wpi,WJ,d)) : d in  [1..6]] >;
 MinimalBasis(J);
 X3 := Scheme(PW,MinimalBasis(J));
 Wpi3 := map<WJ->X3 | [k2, k4, k2^2 + k1*k3 + 2*v1, k2*k4 + 2*v4, k1*k3, k1^3, k3^3, k1*(k2*k3 + 2*v2), k3*(k3^2 + 2*v3)]>;
 P<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := Ambient(X3);
 P13<c111, c112, c122, c222, c1d1, c1d2, c1vw, c2d1, c2d2, c2vw, kvvv, kwww, kvdw, kwdv> := ProjectiveSpace(Q0, 13);
 proj3 := map<X3->P13 | [c1^3, c1^2*c2, c1*c2^2, c2^3, c1*d1, c1*d2, c1*cvw, c2*d1, c2*d2, c2*cvw, cvvv, cwww, cvdw, cwdv]>;
 ProjX3 := Scheme(P13, MinimalBasis(Image(proj3)));
 proj3 := map<X3->ProjX3  | [c1^3, c1^2*c2, c1*c2^2, c2^3, c1*d1, c1*d2, c1*cvw, c2*d1, c2*d2, c2*cvw, cvvv, cwww, cvdw, cwdv]>;
 P5<c11, c12, c22, b1, b2, kvw> := ProjectiveSpace(Q0,5);
 proj2 := map<X3->P5 | [c1^2,c1*c2,c2^2,d1,d2,cvw]>;
 //proj2 := map<Jac3->P5 |[k22, k24, k44, k13 + k22 + 2*v1, k24 + 2*v4, k13]>;
 Proj2X3 := Scheme(P5,MinimalBasis(Image(proj2)));
 sPtX3 := [Pullback(proj2, pt): pt in SingularPoints(Proj2X3)];
 [IrreducibleComponents(ReducedSubscheme(pt)): pt in sPtX3];
 notsing := [X3![0,0,0,0,0,-6,0,1,0], X3![0,0,0,0,0,-5,0,1,0], X3![0,0,0,0,0,0,-3,0,1], X3![0,0,0,0,0,0,-8,0,1]];
 singpts := [X3![10,7,4,1,0,0,0,0,0], X3![10,4,7,1,0,0,0,0,0], X3![r^(90),r^(54),7,1,0,0,0,0,0], X3![3*r^(54),r^(54),4,1,0,0,0,0,0], X3![0,1,0,0,0,0,0,0,0]];
 X3![10,7,4,1,0,0,0,0,0] eq X3![1,4,4,1,0,0,0,0,0];
 [ADEtype(proj3(pt)): pt in singpts];
 // [ADEtype(proj3(pt)): pt in notsing];
 HSWJ := HilbertSeries(Ideal(DefiningPolynomials(WJ)));
 HSX3 := HilbertSeries(Ideal(DefiningPolynomials(X3)));
 FindFirstGenerators(HSWJ);
 HilbertNumerator(HSWJ,[1,1,1,1,2,2,2,2,2,2]);
 Ser<t> := PowerSeriesRing(Rationals(),8);
 Ser!HSX3;
 Ptrunc2 := ProjectiveSpace(Q0, [1,1,2,2,2]);
 trunc2 := map<X3->Ptrunc2|[c1, c2, d1, d2, cvw]>;
 truncim := Scheme(Ptrunc2,MinimalBasis(Image(trunc2)));
 HStrunc := HilbertSeries(Ideal(DefiningPolynomials(truncim)));
 HSX3ww := (HSWJ-HSX3)/2;


 PW5:= ProjectiveSpace(Q0, [1,1,2,3,3]);
 Wpi := map<X3->PW5| [c1,c2,cvw,cvvv,cwww]>;
 P<k1,k2,kvw,kvvv,kwww>:= Ambient(PW5);
 J := ideal< CoordinateRing(PW5) | &cat[ DefiningEquations(Image(Wpi,X3,d)) : d in  [1..6]] >;
 X6 := Scheme(PW5,MinimalBasis(J));
 pi2 := map<X3->X6| [c1,c2,cvw,cvvv,cwww]>;
 SX6 := JacobianSubrankScheme(X6);
 SXred := ReducedSubscheme(SX);





K := GF(7);
C := C3Curve(K,[1,0,-1,1,0]);
J:= Jacobian(C);
X3 := InvariantEmbedding(J);

Qs<s> := FunctionField(Rationals());
P2<X,Y,Z> := ProjectiveSpace(Qs,2);
C := Curve(P2, Y^2*Z^2-(X^3-3*s^2*X*Z^2+Z^3)*(X+2*s*Z));
HC := HyperellipticCurve();
E := EllipticCurve(C);


K := GF(7);
K := Rationals();
C := C3Curve(K,[1,0,-1,1,0]);
J:= Jacobian(C);
//X3 := InvariantEmbedding(J);
//[#[eqs: eqs in DefiningEquations(X3) |  (Degree(eqs) eq i)]: i in [1..6]];
Kum := GeneralKummerSurface(J);
K := BaseRing(Kum);
P<k1,k2,k3,k4>:= PolynomialRing(K,4);
I := ideal<P | DefiningEquations(Kum)>;
A := MatrixRing(P,4)![[1, -1, 1, 0], [2, -1, 0, 0], [1, 0, 0, 0], [-4, -2, 6, 1]];
R := InvariantRing(I,A: PolynomialRing := P, LinearlyReductive := true);
[InvariantsOfDegree(R,d): d in [1..5]];
MB :=MinimalBasis(&cat[InvariantsOfDegree(R,d): d in [1..5]]);
WP<c1,c2,d1,e1,e2> := ProjectiveSpace(K, [1,1,2,3,3]);
pi6 := map<Kum->WP | MB>;
MinimalBasis(ideal< CoordinateRing(WP) | &cat[ DefiningEquations(Image(pi6,Kum,d)) : d in  [1..2]] >);

[
    x2 + x3 + 6*x4,
    x1 + 5*x3 + 3*x4,
    x1*x3 + 4*x2*x3 + 3*x2*x4 + 2*x3^2 + 3*x3*x4 + 2*x4^2,
    x1*x2*x4 + 5*x1*x3^2 + 6*x1*x3*x4 + 3*x1*x4^2 + 5*x2^2*x3 + 2*x2^2*x4 + 6*x2*x3^2 + x2*x3*x4 + 3*x3^2*x4 + x3*x4^2 + 3*x4^3,
    x1^2*x4 + 2*x1*x3^2 + 5*x1*x3*x4 + 3*x1*x4^2 + 2*x2^2*x3 + x2*x3^2 + x2*x3*x4 + 3*x2*x4^2 + x3^2*x4 + 4*x3*x4^2 + 2*x4^3
]

Q0<s> := GF(9);
//K := Rationals();
C := C3Curve(Q0,[1,0,-1,1,0]);
J:= Jacobian(C);
WJ := WeightedJacobianDiag(J);
Jac<v1, v2, v3, v4, v5, v6, k11, k12, k13, k14, k22, k23, k24, k33, k34, k44> := GeneralJacobianSurface(J); 
P<k1,k2,k3,k4,b1,b2,b3,b4,b5,b6> := PolynomialRing(K,[1,1,1,1,2,2,2,2,2,2]);
I := ideal<P | DefiningEquations(WJ)>;
A := MatrixRing(P,10)![[1, -1, 1, 0, 0, 0, 0, 0, 0, 0], [2, -1, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-4, -2, 6, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -2, 3, 0, 2, 0, 0], [0, 0, 0, 0, -1, 2, -1, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 132/5, -9/5, -453/5, -43, 5, -1/5], [0, 0, 0, 0, 990, -14, -2824, -1291, 155, -6]];
R := InvariantRing(I,A: PolynomialRing := P, LinearlyReductive := true);
//FundamentalInvariants(R);
mb:=MinimalBasis(&cat[InvariantsOfDegree(R,d): d in [1..6]]);
mbs:=&cat[InvariantsOfDegree(R,d): d in [1..3]];


function sDegree(mon)
P<k1,k2,k3,k4,b1,b2,b3,b4,b5,b6> := Parent(mon);
v := [k1,k2,k3,k4,b1,b2,b3,b4,b5,b6];
w := [1,1,1,1,2,2,2,2,2,2];
ws := &+[Degree(LeadingMonomial(mon),v[i])*w[i]: i in [1..10]];
return ws;
end function;

P<k1,k2,k3,k4,b1,b2,b3,b4,b5,b6> := Parent(mb[1]);
mb2 := [eqs: eqs in mb | sDegree(eqs) lt 3];
mbs3 := [eqs: eqs in mbs | sDegree(eqs) eq 3];
mb3 := [eqs: eqs in mb | sDegree(eqs) lt 4];
mb4 := [eqs: eqs in mb | sDegree(eqs) lt 5];
[sDegree(eqs): eqs in mb3];
[sDegree(eqs): eqs in mb4];
Ps2<c11,c12,c22,d1,d2,d3> := ProjectiveSpace(Q0,5);
phi2 := map<WJ->Ps2 | [mb2[1]^2,mb2[1]*mb2[2],mb[2]^2,mb[3],mb[4],mb[5]]>;
I := ideal< CoordinateRing(Ps2) | &cat[ DefiningEquations(Image(phi2,WJ,d)) : d in  [1..3]] >;
nX2 := Scheme(Ps2,MinimalBasis(I));
Ps3<c1,c2,d1,d2,d3,e1,e2,e3> := ProjectiveSpace(Q0,[sDegree(eqs): eqs in mb3]);
BaseRing(Ps3);
phi3 := map<WJ->Ps3 | mb3>;
I := ideal< CoordinateRing(Ps3) | &cat[ DefiningEquations(Image(phi3,WJ,d)) : d in  [1..7]] >;
nX3 := Scheme(Ps3,MinimalBasis(I));
[Dimension(V): V in IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(nX3)))];
npi := map<nX3->nX2 | [c1*c1,c1*c2,c2*c2,d1,d2,d3]>;
[IrreducibleComponents(Pullback(npi,pt)): pt in SingularPoints(nX2)];
[Dimension(v): v in &cat [IrreducibleComponents(Pullback(npi,pt)): pt in SingularPoints(nX2)]];
IrreducibleComponents(BaseScheme(npi));
Ps4<c1,c2,d1,d2,d3,e1,e2,e3,f1,f2,f3> := ProjectiveSpace(Q0,[sDegree(eqs): eqs in mb4]);
phi4 := map<WJ->Ps4 | mb4>;
I := ideal< CoordinateRing(Ps4) | &cat[ DefiningEquations(Image(phi4,WJ,d)) : d in  [1..6]] >;
nX4 := Scheme(Ps4,MinimalBasis(I));
HSnX3 := HilbertSeries(Ideal(DefiningPolynomials(nX3)));
HSnX4 := HilbertSeries(Ideal(DefiningPolynomials(nX4)));
Ser!HSnX3;






ids := [k2*b1 + 6*k3*b1 + 6*k2*b3 + k3*b3 + 6*k3*b4 + 2*k4*b4 + 5*k2*b5 + 2*k3*b5 + 2*k1*b6 + 3*k2*b6 + 2*k3*b6,
    k2*b1 + 6*k3*b1 + k3*b2 + 4*k4*b2 + 6*k3*b3 + 3*k4*b3 + 2*k4*b4 + 5*k1*b5 + 5*k2*b5 + 4*k3*b5 + k2*b6 + 6*k3*b6,
    k2*b1 + 6*k3*b1 + 2*k3*b2 + 6*k4*b2 + k1*b3 + 5*k2*b3 + 6*k3*b3 + k4*b3 + k3*b4 + 6*k4*b4,
    k1*b1 + k2*b1 + 5*k3*b1 + k3*b2 + 2*k4*b2 + 4*k2*b3 + 2*k3*b3 + 5*k4*b3 + 4*k3*b4 + 5*k4*b4];
Basis(ids,I);%29
MinimalBasis([NormalForm(f,I): f in ids]);

K := GF(3);
//K := Rationals();
C := C3Curve(K,[1,0,-1,1,0]);
J:= Jacobian(C);
WJ := WeightedJacobianDiag(J);
Jac<v1, v2, v3, v4, v5, v6, k11, k12, k13, k14, k22, k23, k24, k33, k34, k44> := GeneralJacobianSurface(J); 
P<k1,k2,k3,k4,b1,b2,b3,b4,b5,b6> := PolynomialRing(K,[1,1,1,1,2,2,2,2,2,2]);
I := ideal<P | DefiningEquations(WJ)>;
A := MatrixRing(P,10)![[1, -1, 1, 0, 0, 0, 0, 0, 0, 0], [2, -1, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-4, -2, 6, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -2, 3, 0, 2, 0, 0], [0, 0, 0, 0, -1, 2, -1, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 132/5, -9/5, -453/5, -43, 5, -1/5], [0, 0, 0, 0, 990, -14, -2824, -1291, 155, -6]];
R := InvariantRing(I,A: PolynomialRing := P, LinearlyReductive := true);
MinimalBasis(&cat[InvariantsOfDegree(R,d): d in [1..4]]);

K := GF(3);
//K := Rationals();
C := C3Curve(K,[1,0,-1,1,0]);
J:= Jacobian(C);
Jac<v1, v2, v3, v4, v5, v6, k11, k12, k13, k14, k22, k23, k24, k33, k34, k44> := GeneralJacobianSurface(J); 
P<v1, v2, v3, v4, v5, v6, k11, k12, k13, k14, k22, k23, k24, k33, k34, k44> := PolynomialRing(K,16);
I := ideal<P | DefiningEquations(Jac)>;
A := MatrixRing(P,16)![[-2, 3, 0, 2, 0, 0, 2, 1, 0, -11, 0, -6, 0, 1, 0, 0], [-1, 2, -1, 0, 0, 0, 1, 1, -1, -3, 2, -1, 0, 0, 0, 0], [-1, 1, 0, 0, 0, 0, -1, 3, -1, -3, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 11, -4, -1, -6, 3, -3, -1, 1, 0, 0], [132/5, -9/5, -453/5, -43, 5, -1/5, -118/5, -497/5, 11/5, 1296/5, -17, 582/5, 12/5, -157/5, 14/5, 0], [990, -14, -2824, -1291, 155, -6, -513, -2940, -211, 7485, 48, 3343, -11, -900, 86, 0], [0, 0, 0, 0, 0, 0, 1, -2, 1, 2, -2, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, -3, 1, 2, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 4, -4, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -4, 2, 2, 2, -8, 6, 1, -1, 1, 0], [0, 0, 0, 0, 0, 0, -8, 0, 2, 12, -6, 0, 2, -1, 0, 0], [0, 0, 0, 0, 0, 0, -4, -2, 0, 6, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 16, 16, 4, -48, -24, 36, -8, -4, 12, 1]];
R := InvariantRing(I,A: PolynomialRing := P, LinearlyReductive := true);
F := InvariantField(I,A: LinearlyReductive := true);
Group(F);
MB := MinimalBasis(&cat[InvariantsOfDegree(R,d): d in [1..2]]);
WP<[x]> := ProjectiveSpace(K,3);
pi := map<Jac->WP | MB[1..4]>;
MinimalBasis(ideal< CoordinateRing(WP) | &cat[ DefiningEquations(Image(pi,Jac,d)) : d in  [1..2]] >);







K<r> := GF(2);
C := C3Curve(K,[1,0,-1,1,0]);
J:= Jacobian(C);
X3 := InvariantEmbedding(J);
[#[eqs: eqs in DefiningEquations(X3) |  (Degree(eqs) eq i)]: i in [1..6]]; // We observe that there are more relations than in characteristic zero
WP<c1, c2, cvw, d1, d2, e1, e2, e3, e4> := Ambient(X3);
P5<k11, k12, k22, kvw, t1, t2> := ProjectiveSpace(K,5);
proj2 := map< X3->P5 | [c1^2, c1*c2, c2^2, cvw, d1, d2]>;
X2 := Scheme(P5,MinimalBasis(Image(proj2,X3,2)));
proj2 := map< X3->X2 | [c1^2, c1*c2, c2^2, cvw, d1, d2]>;
sptsX2 := SingularPoints(X2);
BaseScheme(proj2);
[Dimension(Pullback(proj2,pt)): pt in sptsX2];
 HSX3 := HilbertSeries(Ideal(DefiningPolynomials(X3)));
 Ser<t> := PowerSeriesRing(Rationals(),8);
 Ser!HSX3;
lines := IrreducibleComponents(Scheme(X2,[k11,k12,t2]));
[Pullback(proj2,l): l in lines];
v := [1,0,-1,1,0];

Kum := GeneralKummerSurface(J);
K := BaseRing(Kum);
P<k1,k2,k3,k4>:= PolynomialRing(K,4);
I := ideal<P | DefiningEquations(Kum)>;
A := MatrixRing(P,4)![[1, -1, 1, 0], [2, -1, 0, 0], [1, 0, 0, 0], [-4, -2, 6, 1]];
R := InvariantRing(I,A: PolynomialRing := P, LinearlyReductive := true);
R := InvariantRing(I,A: PolynomialRing := P, Reductive := true);
[InvariantsOfDegree(R,d): d in [1..5]];
MinimalBasis(&cat[InvariantsOfDegree(R,d): d in [1..5]]);
Wp<x1, y1, x2, y3>:= ProjectiveSpace(K,[1,1,2,3]);
Wp<x1, y1, x2, x3,y3>:= ProjectiveSpace(K,[1,1,2,3,3]);
phi := map<Kum->Wp| MinimalBasis(&cat[InvariantsOfDegree(R,d): d in [1..5]])>;
X63 := Image(phi);
IrreducibleComponents(JacobianSubrankScheme(X63));
IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(X63)));
[Dimension(V): V in IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(X63)))];



J := Jacobian(C3);
WJ<k1, k2, k3, k4, v1, v2, v3, v4, v5, v6> := WeightedJacobian(J);
WP<c1, c2, d1, d2, d3, e1, e2, e3> := ProjectiveSpace(Q0,[1,1,2,2,2,3,3,3]);
WP4<c1, c2, d1, d2, d3, e1, e2, e3, f1> := ProjectiveSpace(Q0,[1,1,2,2,2,3,3,3,4]);
phi := map < WJ -> WP | [k1 + k2 + k3, 2*k2 + k4, k2^2 + 2*k1*k3, k2^2 + k2*k4 + 2*v4, 2*k2*k3 + k3^2 + 2*v2 + v3, k1^2*k2 + k1*k2^2 + 2*k1*k3^2 + 2*k2*k3^2 + 2*k2^2*k4 + k1*k3*k4,  k1^2*k2 + 2*k1*k2^2 + k1*k3^2 + k1^2*k4 + k2^2*k4 + k3^2*k4 + k1*k4^2 + k2*k4^2 + k3*k4^2, 2*k1^2*k2 + 2*k2^3 + k3^3 + 2*k1^2*k4 + 2*k2*k3*k4 + k3^2*k4 + 2*k1*k4^2 +   2*k3*k4^2 + k3*v2 + 2*k3*v3 + 2*k3*v4 + 2*k4*v4 + k2*v5 + 2*k3*v5 + k1*v6 + 2*k2*v6]>;
phi4 := map < WJ -> WP4 | [k1 + k2 + k3, 2*k2 + k4, k2^2 + 2*k1*k3, k2^2 + k2*k4 + 2*v4, 2*k2*k3 + k3^2 + 2*v2 + v3, k1^2*k2 + k1*k2^2 + 2*k1*k3^2 + 2*k2*k3^2 + 2*k2^2*k4 + k1*k3*k4,  k1^2*k2 + 2*k1*k2^2 + k1*k3^2 + k1^2*k4 + k2^2*k4 + k3^2*k4 + k1*k4^2 + k2*k4^2 + k3*k4^2, 2*k1^2*k2 + 2*k2^3 + k3^3 + 2*k1^2*k4 + 2*k2*k3*k4 + k3^2*k4 + 2*k1*k4^2 +   2*k3*k4^2 + k3*v2 + 2*k3*v3 + 2*k3*v4 + 2*k4*v4 + k2*v5 + 2*k3*v5 + k1*v6 + 2*k2*v6,k1^4 + 2*k1^3*k2 + 2*k1*k2^3 + 2*k1^3*k3 + k1^2*k2*k3 + k1*k2^2*k3 + k1^2*k3^2 + k1*k2*k3^2 + k2^2*k3^2 + 2*k1*k3^3 + 2*k3^4 + k1^2*k2*k4 + k2^2*k3*k4 + 2*k2*k3^2*k4 + 
 k1^2*v1 + k1*k2*v1 + k2*k3*v1 + 2*k1^2*v2 + 2*k1*k2*v2 + 2*k2^2*v2 + k1*k3*v2 + 2*k2*k3*v2 + 2*k1^2*v3 + k3^2*v3 + 2*k1^2*v4 + 2*k2*k3*v4 + k3^2*v4]>;
I := ideal< CoordinateRing(WP) | &cat[ DefiningEquations(Image(phi,WJ,d)) : d in  [3..6]] >;
I4 := ideal< CoordinateRing(WP4) | &cat[ DefiningEquations(Image(phi4,WJ,d)) : d in  [3..8]] >;
fX3 := Scheme(WP,MinimalBasis(I));
fX4 := Scheme(WP4,MinimalBasis(I4));
HSfX3 := HilbertSeries(Ideal(DefiningPolynomials(fX3)));
HSfX4 := HilbertSeries(Ideal(DefiningPolynomials(fX4)));
Ser<t> := PowerSeriesRing(Rationals(),8);
Ser!HSfX3;
Ser!HSfX4;
P5<k11, k12, k22, kvw, t1, t2> := ProjectiveSpace(Q0,5);
proj2 := map< fX3->P5 | [c1^2, c1*c2, c2^2, d1, d2, d3]>;
fX2 := Scheme(P5,MinimalBasis(Image(proj2,fX3,2)));
proj2 := map< fX3->fX2 | [c1^2, c1*c2, c2^2, d1, d2, d3]>;
sptsX2 := SingularPoints(fX2);
IrreducibleComponents(ReducedSubscheme(BaseScheme(proj2)));
[ReducedSubscheme(Pullback(proj2,pt)): pt in sptsX2];
proj24 := map< fX4->fX2 | [c1^2, c1*c2, c2^2, d1, d2, d3]>;
sptsX2 := SingularPoints(fX2);
IrreducibleComponents(ReducedSubscheme(BaseScheme(proj24)));
[IrreducibleComponents(ReducedSubscheme(Pullback(proj24,pt))): pt in sptsX2];

 Jac<v1, v2, v3, v4, v5, v6, k11, k12, k13, k14, k22, k23, k24, k33, k34, k44> := GeneralJacobianSurface(J);
 PJ<c11,c12,c22,d1,d2,d3,e1,e2> := ProjectiveSpace(Q0,[1,1,1,1,1,1,2,2]);
 phi := map<Jac->PJ | [k11 + 2*k12 + 2*k13 + k22 + 2*k23 + k33, 2*k12 + k14 + 2*k22 + 2*k23 + k24 + k34, 4*k22 + 4*k24 + k44, 2*k13 + k22, k22 + k24 + 2*v4, 2*k23 + k33 + 2*v2 + v3,  k11^2 + 2*k11*k12 + 2*k11*k13 + k11*k22 + k12*k22 + k11*k24 + k12*k24 + 2*k11*k33 + 2*k12*k33 + 2*k13*k33 + 2*k22*k33 + 2*k33^2 + 2*k11*k34 + 2*k13*k34 + 2*k22*k34 +   2*k11*v1 + 2*k12*v1 + 2*k23*v1 + 2*k11*v2 + k12*v2 + 2*k13*v2 + k22*v2 + 2*k23*v2 + k33*v2 + 2*k23*v3 + k33*v3 + 2*k11*v4 + 2*k23*v4 + k13*v5 + 2*k22*v5 + k13*v6 +   2*k22*v6, 2*(k11*k12 + 2*k11*k13 + 2*k12*k22 + 2*k22^2 + 2*k12*k23 + 2*k22*k23 + k11*k24 + k22*k24 + 2*k11*k33 + 2*k22*k33 + k23*k33 + k12*k34 + 2*k13*k34 + k22*k34 +    2*k12*k44 + k13*k44 + k24*k44 + 2*k11*v1 + k14*v1 + k22*v1 + 2*k24*v1 + k22*v2 + 2*k24*v2 + 2*k33*v2 + k34*v2 + 2*k11*v3 + k14*v3 + k22*v3 + 2*k24*v3 + 2*k11*v4 +    k14*v4 + 2*k44*v4)]>;
 I := ideal< CoordinateRing(PJ) | &cat[ DefiningEquations(Image(phi,Jac,d)) : d in  [1..3]] >;
 X4 := Scheme(PJ,MinimalBasis(I));
 P5<w11, w12, w22, t1, t2, t3> := ProjectiveSpace(Q0,5);
 proj2 := map< X4->P5 | [c11, c12, c22, d1, d2, d3]>;
 fX2 := Scheme(P5,MinimalBasis(Image(proj2,X4,2)));
 proj2 := map< X4->fX2 | [c11, c12, c22, d1, d2, d3]>;
 sptsX2 := SingularPoints(fX2);
 IrreducibleComponents(ReducedSubscheme(BaseScheme(proj2)));
 [Dimension(V): V in &cat[IrreducibleComponents(ReducedSubscheme(Pullback(proj2,pt))): pt in sptsX2]];


eqs :=[2*k11^2 + 4*k11*k12 + 4*k11*k14 + 2*k11*k22 + 4*k12*k22 + 2*k11*k23 + 2*k22*k23 + 12*k11*k24 + 8*k12*k24 + 6*k11*k33 + 4*k12*k33 + 4*k22*k33 + 2*k23*k33 + 4*k33^2 +   4*k11*k34 + 12*k12*k34 + 8*k13*k34 + 4*k22*k34 + 12*k23*k34 + 8*k33*k34 + k11*(4*v2 + 2*v3) + k22*(4*v2 + 2*v3) + k33*(4*v2 + 2*v3) + k14*(8*v2 + 4*v3) +   k24*(8*v2 + 4*v3) + k34*(8*v2 + 4*v3), k11^2 + k11*k12 + k11*k22 + k22^2 + k11*k23 + k12*k23 + k22*k23 + 2*k22*k24 + 2*k22*k33 + 2*k23*k33 + 2*k33^2 + 2*k11*k34 +   2*k12*k34 + 2*k23*k44 + k12*(v1 + v2) + k24*(v1 + v2 + v3) + k11*(2*v2 + v3) + k33*(2*v2 + v3) + k14*(2*v1 + 2*v3) + k22*(2*v1 + 2*v2 + 2*v3) + k34*(2*v2 + v4) +   k23*(2*v2 + 2*v3 + 2*v4), 4*k11*k22 + 2*k12*k22 + 4*k22*k23 + 16*k11*k24 + 8*k12*k24 + 2*k22*k33 + 16*k22*k34 + 8*k23*k34 + 16*k11*k44 + 8*k12*k44 + 16*k23*k44 +   8*k33*k44 + k22*(2*v2 + 4*v3) + k24*(8*v2 + 16*v3) + k44*(8*v2 + 16*v3), k11^2 + 2*k11*k12 + 2*k11*k13 + 2*k12*k22 + k11*k23 + k12*k23 + k11*k24 + k11*k33 + k12*k33 +   2*k13*k33 + k22*k33 + 2*k33^2 + k22*k34 + 2*k23*k34 + k13*v2 + 2*k22*v2 + k12*(v1 + 2*v2) + k33*(v3 + v4) + k23*(v1 + 2*v2 + 2*v4) + k11*(v1 + 2*v2 + 2*v3 + 2*v4),  2*k11^2 + 8*k11*k12 + 8*k11*k22 + 6*k11*k23 + 12*k12*k23 + 6*k11*k33 + 12*k12*k33 + 4*k22*k33 + 10*k23*k33 + 4*k33^2 + k11*(4*v2 + 2*v3) + k33*(4*v2 + 2*v3) +   k12*(8*v2 + 4*v3) + k23*(8*v2 + 4*v3), k11*k13 + 2*k11*k22 + k12*k23 + k11*k24 + k11*k33 + k12*k33 + 2*k12*k34 + k12*v2 + 2*k22*v2 + 2*k33*v2 + k23*(v1 + v3) +   k13*(v2 + v4) + k11*(2*v1 + 2*v3 + 2*v4), k11^2 + 2*k11*k12 + 2*k12*k22 + k12*k23 + 2*k11*k24 + 2*k12*k24 + k11*k33 + k22*k33 + 2*k33^2 + k13*v2 + 2*k22*v2 +   k23*(v1 + 2*v2) + k33*v3 + k12*(v1 + 2*v2 + v4) + k11*(v1 + 2*v2 + 2*v3 + v4), 4*k11*k12 + 4*k11*k13 + 4*k11*k22 + 2*k12*k22 + 2*k22^2 + 4*k12*k23 + 2*k22*k23 +   4*k22*k24 + 2*k12*k33 + 2*k22*k33 + 4*k23*k33 + 2*k13*k34 + 4*k22*k34 + 2*k12*k44 + 2*k22*k44 + 2*k12*v2 + 2*k13*v2 + 4*k34*v2 + k23*(2*v1 + 2*v3) +   k22*(4*v1 + 2*v2 + 4*v3) + k24*(2*v1 + 2*v2 + 2*v3 + 4*v4) + k14*(4*v1 + 4*v3 + 4*v4), 4*k11*k13 + 2*k11*k22 + k12*k22 + 2*k11*k23 + 2*k22*k23 + 4*k12*k33 + 2*k13*k33 +   k22*k33 + k22*(v2 + 2*v3) + k13*(2*v2 + 4*v3), k11*k13 + 2*k11*k22 + k11*k24 + 2*k22*k24 + 2*k11*k33 + k12*k33 + k12*v2 + k13*v2 + 2*k33*v2 + k23*(v1 + v3) +   k22*(2*v2 + v4) + k11*(2*v1 + 2*v3 + 2*v4), k11^2 + 2*k11*k12 + 2*k12*k22 + k22^2 + k11*k23 + k12*k23 + 2*k11*k24 + 2*k12*k33 + k22*k33 + 2*k33^2 + 2*k22*k34 +   k13*(v1 + v2) + k12*(v1 + 2*v2) + k22*(2*v1 + 2*v2) + k33*v3 + k23*(v1 + 2*v2 + v4) + k11*(v1 + 2*v2 + 2*v3 + v4),  2*k11^2 + 4*k11*k12 + 4*k11*k14 + 2*k12*k22 + 2*k11*k23 + 2*k12*k23 + 2*k22*k23 + 4*k11*k33 + 2*k22*k33 + 2*k23*k33 + 4*k23*k34 + 2*k33*k34 + 4*k13*v2 + 2*k22*v3 +   4*k33*v5 + k23*(4*v1 + 2*v2 + 4*v3 + 2*v5) + k11*(4*v2 + 2*v3 + 2*v4 + 2*v6) + k12*(2*v1 + 4*v6), k11^2 + 2*k11*k22 + 2*k12*k22 + 2*k11*k23 + 2*k12*k23 + k11*k24 +   2*k12*k24 + 2*k11*k33 + 2*k12*k33 + k23*k33 + 2*k33^2 + k11*k34 + k13*k34 + k12*k44 + 2*k13*v2 + k34*(v1 + v2) + k33*(v2 + v3) + k22*(2*v2 + v3) + k23*(v1 + 2*v2 + v3) +   k24*(2*v2 + 2*v3) + k11*(2*v1 + 2*v2 + 2*v4) + k14*(2*v1 + v3 + 2*v4), k11^2 + 2*k11*k12 + 2*k11*k13 + k11*k22 + k12*k22 + k11*k24 + k12*k24 + 2*k11*k33 + 2*k12*k33 +   2*k13*k33 + 2*k22*k33 + 2*k33^2 + 2*k11*k34 + 2*k13*k34 + 2*k22*k34 + k12*(2*v1 + v2) + k33*(v2 + v3) + k11*(2*v1 + 2*v2 + 2*v4) + k23*(2*v1 + 2*v2 + 2*v3 + 2*v4) +   k13*(2*v2 + v5 + v6) + k22*(v2 + 2*v5 + 2*v6), 2*k11^2 + 2*k11*k12 + 4*k11*k14 + 2*k11*k22 + 2*k12*k22 + 2*k11*k23 + 4*k12*k23 + 4*k22*k23 + 2*k11*k24 + 2*k12*k33 +   2*k22*k33 + 2*k23*k33 + 4*k13*k34 + 2*k33*k34 + k12*(4*v2 + 2*v5) + k13*(4*v2 + 2*v5) + k33*(2*v2 + 4*v5) + k22*(4*v2 + 4*v3 + 4*v5) +   k11*(2*v1 + 4*v2 + 4*v3 + 4*v4 + 2*v6) + k23*(4*v1 + 4*v3 + 2*v4 + 4*v6), 2*k11^2 + 2*k11*k12 + 4*k11*k13 + 4*k11*k14 + 2*k11*k22 + 2*k12*k22 + 2*k11*k23 + 4*k12*k23 +   2*k11*k24 + 4*k12*k33 + 4*k22*k33 + 4*k33^2 + 2*k11*k34 + 2*k12*k34 + 2*k23*k34 + k13*(4*v2 + 2*v5) + k22*(2*v2 + 4*v5) + k12*(4*v1 + 2*v2 + 4*v5) +   k23*(4*v1 + 2*v2 + 2*v4 + 4*v5) + k11*(4*v1 + 4*v3 + 4*v4 + 4*v5 + 2*v6) + k33*(4*v3 + 4*v6), k11^2 + 2*k11*k12 + k11*k13 + 2*k11*k22 + 2*k11*k23 + k12*k23 + k22*k23 +   k11*k24 + k11*k33 + k12*k33 + 2*k13*k33 + 2*k33^2 + 2*k12*v2 + k13*v2 + 2*k22*v2 + 2*k23*v2 + k33*(v1 + v3) + k11*(2*v1 + 2*v2 + 2*v3 + 2*v4),  k11^2 + 2*k11*k12 + k11*k14 + k12*k22 + k11*k23 + k12*k23 + 2*k22*k23 + k11*k24 + 2*k12*k24 + k11*k33 + k12*k33 + 2*k22*k33 + 2*k12*k34 + 2*k13*k34 + k22*k34 + k33*k34 +   k11*k44 + k23*k44 + 2*k33*k44 + k12*v2 + k13*v2 + k23*(v1 + v3) + k34*v5 + k33*(2*v2 + 2*v5) + k11*(2*v1 + 2*v2 + v6) + k24*(v2 + 2*v3 + 2*v5 + v6) +   k14*(v2 + 2*v3 + 2*v4 + 2*v6) + k22*(2*v3 + v5 + 2*v6), k11^2 + k11*k14 + 2*k11*k22 + k12*k22 + k22^2 + k11*k24 + 2*k22*k24 + k11*k33 + 2*k12*k33 + 2*k23*k33 + 2*k33^2 +   k13*k34 + k22*k34 + 2*k23*k34 + 2*k33*k34 + k11*k44 + k12*k44 + k13*k44 + k23*k44 + k33*v3 + k11*(v1 + 2*v2 + 2*v3 + v4) + k12*(v1 + v5) +   k23*(2*v1 + 2*v2 + v3 + 2*v4 + v5) + k13*(2*v2 + 2*v5) + k24*(v1 + v2 + v3 + 2*v5) + k22*(2*v1 + v2 + v3 + 2*v5) + k34*(2*v2 + v6) +   k14*(2*v1 + v2 + v3 + 2*v4 + v5 + 2*v6), 2*k11*k12 + 4*k11*k13 + 4*k12*k22 + 4*k22^2 + 4*k12*k23 + 4*k22*k23 + 2*k11*k24 + 2*k22*k24 + 4*k11*k33 + 4*k22*k33 + 2*k23*k33 +   2*k12*k34 + 4*k13*k34 + 2*k22*k34 + 4*k12*k44 + 2*k13*k44 + 2*k24*k44 + 4*k33*v2 + 2*k34*v2 + k22*(2*v1 + 2*v2 + 2*v3) + k24*(4*v1 + 4*v2 + 4*v3) + 4*k44*v4 +   k14*(2*v1 + 2*v3 + 2*v4) + k11*(4*v1 + 4*v3 + 4*v4)];
Pr<[x]> := ProjectiveSpace(Q0,#eqs-1);
phi := map<Jac->Pr | eqs>;
 I := ideal< CoordinateRing(Pr) | DefiningEquations(Image(phi,Jac,1)) >;
neweqs := eqs[[3,12,13,14,15,16,17,18,19,20]];
Pr<[x]> := ProjectiveSpace(Q0,#neweqs-1);
phi := map<Jac->Pr | neweqs>;
I := ideal< CoordinateRing(Pr) | &cat[ DefiningEquations(Image(phi,Jac,d)) : d in  [1..3]] >;
X4 := Scheme(Pr,MinimalBasis(I));
HSX4 := HilbertSeries(Ideal(DefiningPolynomials(X4)));
Ser<t> := PowerSeriesRing(Rationals(),8);
Ser!HSX4;
P5<w11, w12, w22, t1, t2, t3> := ProjectiveSpace(Q0,5);
proj2 := map< Jac->P5 | [k11 + 2*k12 + 2*k13 + k22 + 2*k23 + k33, 2*k12 + k14 + 2*k22 + 2*k23 + k24 + k34, 4*k22 + 4*k24 + k44, 2*k13 + k22, k22 + k24 + 2*v4, 2*k23 + k33 + 2*v2 + v3]>;
I := ideal< CoordinateRing(P5) | &cat[ DefiningEquations(Image(proj2,Jac,d)) : d in  [1..2]] >;
X2 := Scheme(P5,MinimalBasis(I));