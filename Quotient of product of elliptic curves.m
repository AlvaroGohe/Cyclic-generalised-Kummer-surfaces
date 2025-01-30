load "Functions.m";

function IntMatDes(pt)
X := Scheme(pt);
A,P := AffinePatch(X,pt);
Q<[xx]> := Ambient(A);
coords := Eltseq(P);
P := Scheme(A,[xx[i]-coords[i]: i in [1..#coords]]);
dsd := ResolveSingByBlowUp(A,P);
return IntersectionMatrix(dsd);
end function;

function DesInfo(pt)
X := Scheme(pt);
A,P := AffinePatch(X,pt);
Q<[xx]> := Ambient(A);
coords := Eltseq(P);
P := Scheme(A,[xx[i]-coords[i]: i in [1..#coords]]);
dsd := ResolveSingByBlowUp(A,P);
S := Surface(Ambient(A),DefiningEquations(A));
return dsd, S;
end function;

function CoordPoint(ptscheme)
k := BaseRing(ptscheme);

// Extract the defining equations of the scheme
eqns := DefiningEquations(ptscheme);

// Solve equations in affine space
num_vars := #Gradings(AmbientSpace(ptscheme))[1]; // Number of variables
A := AffineSpace(k, num_vars);
weights := Gradings(AmbientSpace(ptscheme))[1];

// Try normalizing each coordinate until one is non-zero
weighted_pts := [];
for i in [1..num_vars] do
    // Add normalization condition: A.i = 1
    normalize_eqn := A.i - 1;
    eqns_affine := [Evaluate(f, [A.j^weights[j] : j in [1..num_vars]]) : f in eqns] cat [normalize_eqn];
    
    // Solve for affine points
    affine_pts := Points(Scheme(A, eqns_affine));
    
    if #affine_pts gt 0 then
        // Convert affine solutions to weighted projective coordinates
        weighted_pts := [ [pt[j]^weights[j] : j in [1..num_vars]] : pt in affine_pts ];
        break;
    end if;
end for;

// Output the weighted projective coordinates
return weighted_pts[1];
end function;


K := Rationals();
K<s> := GF(7^3);
K<s> := GF(3^12);
//K<s> := GF(2^2);
l := 1;
P2<X,Y,Z> := ProjectiveSpace(K,2);
cub := Scheme(P2, X^3+Y^3+Z^3-l*X*Y*Z);
dcub := [JacobianMatrix(cub)[1,i]: i in [1..3]];
P2<x,y,z> := ProjectiveSpace(K,2);
phi := map<cub->P2 | dcub>;
PW<x,y,z,w> := ProjectiveSpace(K,[1,1,1,3]);
phi := map<cub->PW | dcub cat [0]>;
ram := DefiningEquations(Image(phi))[2];
X3 := Scheme(PW,w^2-ram);
spts3 := IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(X3)));
spts3;
pts3 := [X3![0,0,1,0],X3![0,1,0,0],X3![1,0,0,0]];
[ADEtype(pt): pt in pts3];
[IntMatDes(pt): pt in pts3];
//P10<[xx]> := ProjectiveSpace(K,10);
//tau := map<X3->P10 | [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3, w]>;
//X3proj := Scheme(P10,DefiningEquations(Image(tau)));
//tau := map<X3->X3proj | [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3, w]>;
//pts3proj := [tau(pt): pt in pts3];
//A,P := AffinePatch(X3proj,pts3proj[1]);
//dsds := ResolveSingularSurface(X3proj);

PW6<a,b,c,d> := ProjectiveSpace(K,[1,1,2,4]);
phi6 := map<cub->PW6| [dcub[3],dcub[1]+dcub[2],dcub[1]*dcub[2],0]>;
ram6 := DefiningEquations(Image(phi6))[2];
X6 := Scheme(PW6,d^2-(b^2 - 4*c)*ram6);
spts6 := IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(X6)));
spts6;
pts :=[X6![2,0,-1,0], X6![2,2,-1,0], X6![-2,2,-1,0],X6![-5,1,1,0], X6![s^277,-1,2,0],X6![s^235,-1,2,0],X6![s^229,-1,2,0],X6![1,-1,2,0]];
pts := [X6![s^12,-1,1,0],X6![s^10,-1,1,0],X6![s^4,-1,1,0]];
otherpts := [X6![1,0,0,0],X6![0,1,0,0]];
[ADEtype(pt): pt in pts];
IntMatDes(otherpts[1]);
IntMatDes(otherpts[2]);

K := Rationals();
// K<s> := GF(7);
// K<s> := GF(3^12);
// K<s> := GF(2^2);
// vec := [1, -1, 1, -1, 0];
vec := [0, 1, 1, 1, 0];
a1 := vec[1];
a2 := vec[2];
a3 := vec[3];
a4 := vec[4];
a6 := vec[5];
P2<X,Y,Z> := ProjectiveSpace(K,2);
ec := Scheme(P2, -X^3 - a2*X^2*Z + a1*X*Y*Z + Y^2*Z - a4*X*Z^2 + a3*Y*Z^2 - a6*Z^3);
dec := [JacobianMatrix(ec)[1,i]: i in [1..3]];
P2<x,y,z> := ProjectiveSpace(K,2);
phi := map<ec->P2 | dec>;
PW<x,y,z,w> := ProjectiveSpace(K,[1,1,1,3]);
phi := map<ec->PW | dec cat [0]>;
ram := DefiningEquations(Image(phi))[2];
phi := map<ec->PW | dec cat [0]>;
X3 := Scheme(PW,w^2-ram);
spts3 := IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(X3)));
spts3;
pts3 := [X3![1,-1,0,0],X3![1,1,1,0],X3![0,0,1,0]];
[ADEtype(pt): pt in pts3];
[IntMatDes(pt): pt in pts3];
P22<X1,Y1,Z1,X2,Y2,Z2> := ProductProjectiveSpace(K,[2,2]);
ecec := Scheme(P22, [-X1^3 - a2*X1^2*Z1 + a1*X1*Y1*Z1 + Y1^2*Z1 - a4*X1*Z1^2 + a3*Y1*Z1^2 - a6*Z1^3,-X2^3 - a2*X2^2*Z2 + a1*X2*Y2*Z2 + Y2^2*Z2 - a4*X2*Z2^2 + a3*Y2*Z2^2 - a6*Z2^3]);
phi22 := map<ecec->P2 | [-(Y2*Z1) + Y1*Z2, X2*Z1 - X1*Z2, -(X2*Y1) + X1*Y2]>;
Image(phi22);
IC := IrreducibleComponents(Pullback(phi22,P2![1,1,-1]));
IC;


// KK<x0,y0,z0> := PolynomialRing(K,3);
// P22K<X1,Y1,Z1,X2,Y2,Z2> := BaseExtend(P22,KK);
// P2K<x,y,z> := BaseExtend(P2,KK);
// ecec := Scheme(P22K, [-X1^3 - a2*X1^2*Z1 + a1*X1*Y1*Z1 + Y1^2*Z1 - a4*X1*Z1^2 + a3*Y1*Z1^2 - a6*Z1^3,-X2^3 - a2*X2^2*Z2 + a1*X2*Y2*Z2 + Y2^2*Z2 - a4*X2*Z2^2 + a3*Y2*Z2^2 - a6*Z2^3]);
// phi22 := map<ecec->P2K | [-(Y2*Z1) + Y1*Z2, X2*Z1 - X1*Z2, -(X2*Y1) + X1*Y2]>;
// IrreducibleComponents(Pullback(phi22,P2K![x0,y0,z0]));

Q := Rationals();
// K<s> := GF(23);
// K<s> := GF(3^12);
// K<s> := GF(2^2);
// vec := [1, -1, 1, -1, 0];
K<a1,a2,a3,a4,a6> := FunctionField(Q,5);
P2<X,Y,Z> := ProjectiveSpace(K,2);
ec := Scheme(P2, -X^3 - a2*X^2*Z + a1*X*Y*Z + Y^2*Z - a4*X*Z^2 + a3*Y*Z^2 - a6*Z^3);
dec := [JacobianMatrix(ec)[1,i]: i in [1..3]];
P2<x,y,z> := ProjectiveSpace(K,2);
phi := map<ec->P2 | dec>;
PW<x,y,z,w> := ProjectiveSpace(K,[1,1,1,3]);
phi := map<ec->PW | dec cat [0]>;
// DefiningEquations(Image(phi));
// ram := DefiningEquations(Image(phi))[2]; This takes a long time to evaluate
phi := map<ec->PW | dec cat [0]>;
X3 := Scheme(PW,w^2-ram);
PW6<a,b,c,d> := ProjectiveSpace(K,[1,1,2,4]);
phi6 := map<ec->PW6| [-2*dec[1] + a1*dec[2], -(a3*dec[1]) + a1*a3*dec[2] - a1*dec[3],dec[1]*(dec[1] - a1*dec[2]),0]>;
// ram6 := DefiningEquations(Image(phi6))[2]; This takes a long time to evaluate


function KummerEll2flawed(field, vec1, vec2)
a1 := vec1[1];
a2 := vec1[2];
a3 := vec1[3];
a4 := vec1[4];
a6 := vec1[5];
b1 := vec2[1];
b2 := vec2[2];
b3 := vec2[3];
b4 := vec2[4];
b6 := vec2[5];
P22<X1,Y1,Z1,X2,Y2,Z2> := ProductProjectiveSpace(field, [2,2]);
E1xE2 := Scheme(P22, [-X1^3 - a2*X1^2*Z1 + a1*X1*Y1*Z1 + Y1^2*Z1 - a4*X1*Z1^2 + a3*Y1*Z1^2 - a6*Z1^3, -X2^3 - b2*X2^2*Z2 + b1*X2*Y2*Z2 + Y2^2*Z2 - b4*X2*Z2^2 + b3*Y2*Z2^2 - b6*Z2^3]);
P4<x11,x12,x21,x22,w> := ProjectiveSpace(field, 4);
phi := map<E1xE2->P4 | [X1*X2, X1*Z2, X2*Z1, Z1*Z2, 2*Y1*Y2+(a1*X1+a3*Z1)*Y2+(b1*X2+b3*Z2)*Y1]>;
K2 := Image(phi);
phi := map<E1xE2->K2 | [X1*X2, X1*Z2, X2*Z1, Z1*Z2, 2*Y1*Y2+(a1*X1+a3*Z1)*Y2+(b1*X2+b3*Z2)*Y1]>;
return K2, phi;
end function;


function KummerEll2(field, vec1, vec2)
a1 := vec1[1];
a2 := vec1[2];
a3 := vec1[3];
a4 := vec1[4];
a6 := vec1[5];
b1 := vec2[1];
b2 := vec2[2];
b3 := vec2[3];
b4 := vec2[4];
b6 := vec2[5];
P22<X1,Y1,Z1,X2,Y2,Z2> := ProductProjectiveSpace(field, [2,2]);
E1xE2 := Scheme(P22, [-X1^3 - a2*X1^2*Z1 + a1*X1*Y1*Z1 + Y1^2*Z1 - a4*X1*Z1^2 + a3*Y1*Z1^2 - a6*Z1^3, -X2^3 - b2*X2^2*Z2 + b1*X2*Y2*Z2 + Y2^2*Z2 - b4*X2*Z2^2 + b3*Y2*Z2^2 - b6*Z2^3]);
R<x1,z1,x2,z2,w11,w12,w22> := PolynomialRing(field,7);
B := [ideal< R | x1, z1, w11, w12>, ideal< R | x2, z2, w12, w22>];
Zwts := [[1,1,0,0,2,1,0],[0,0,1,1,0,1,2]];
Qwts := [];
C := CoxRing(R,B,Zwts,Qwts);
X := ToricVariety(C);
phi := map<E1xE2->X | [X1, Z1, X2, Z2, (a1*X1*Y1 + Y1^2 + a3*Y1*Z1), 2*Y1*Y2+(a1*X1+a3*Z1)*Y2+(b1*X2+b3*Z2)*Y1, (b1*X2*Y2 + Y2^2 + b3*Y2*Z2)]>;
K2 := Scheme(X,MinimalBasis(Image(phi)));
phi := map<E1xE2->K2 | [X1, Z1, X2, Z2, (a1*X1*Y1 + Y1^2 + a3*Y1*Z1), 2*Y1*Y2+(a1*X1+a3*Z1)*Y2+(b1*X2+b3*Z2)*Y1, (b1*X2*Y2 + Y2^2 + b3*Y2*Z2)]>;
return K2, phi;
end function;

function KummerEll3(field, vec)
a1 := vec[1];
a2 := vec[2];
a3 := vec[3];
a4 := vec[4];
a6 := vec[5];
PW<x,y,z,w> := ProjectiveSpace(field,[1,1,1,3]);
K3 := Scheme(PW,w^2 - a6*x^6 + a1*a3^2*x^5*y + 3*a1*a6*x^5*y + a2*a3^2*x^4*y^2 +  2*a1*a3*a4*x^4*y^2 - 3*a1^2*a6*x^4*y^2 + 3*a2*a6*x^4*y^2 + a3^3*x^3*y^3 +  2*a2*a3*a4*x^3*y^3 + a1*a4^2*x^3*y^3 + a1^3*a6*x^3*y^3 - 6*a1*a2*a6*x^3*y^3 +  5*a3*a6*x^3*y^3 + 3*a3^2*a4*x^2*y^4 + a2*a4^2*x^2*y^4 +  3*a1^2*a2*a6*x^2*y^4 - 3*a2^2*a6*x^2*y^4 - 4*a1*a3*a6*x^2*y^4 +  5*a4*a6*x^2*y^4 + 3*a3*a4^2*x*y^5 + 3*a1*a2^2*a6*x*y^5 - 4*a2*a3*a6*x*y^5 -  4*a1*a4*a6*x*y^5 + a4^3*y^6 + a2^3*a6*y^6 - 4*a2*a4*a6*y^6 + 7*a6^2*y^6 +  a1*a3*x^5*z + a4*x^5*z + a2*a3*x^4*y*z - 2*a1*a4*x^4*y*z +  a1^3*a3*x^3*y^2*z - 2*a1*a2*a3*x^3*y^2*z - a3^2*x^3*y^2*z +  3*a1^2*a4*x^3*y^2*z - 2*a2*a4*x^3*y^2*z - 9*a6*x^3*y^2*z +  3*a1^2*a2*a3*x^2*y^3*z - 2*a2^2*a3*x^2*y^3*z - a1*a3^2*x^2*y^3*z +  4*a1*a2*a4*x^2*y^3*z - 7*a3*a4*x^2*y^3*z + 14*a1*a6*x^2*y^3*z +  3*a1*a2^2*a3*x*y^4*z - 4*a2*a3^2*x*y^4*z + a2^2*a4*x*y^4*z +  2*a1*a3*a4*x*y^4*z - 6*a4^2*x*y^4*z - 4*a1^2*a6*x*y^4*z + 9*a2*a6*x*y^4*z +  a2^3*a3*y^5*z - 4*a2*a3*a4*y^5*z + 3*a1*a4^2*y^5*z - 4*a1*a2*a6*y^5*z +  14*a3*a6*y^5*z - a2*x^4*z^2 + a1^3*x^3*y*z^2 + 2*a1*a2*x^3*y*z^2 -  a3*x^3*y*z^2 + 2*a2^2*x^2*y^2*z^2 + 7*a1*a3*x^2*y^2*z^2 + 8*a4*x^2*y^2*z^2 -  2*a1*a2^2*x*y^3*z^2 - a1^2*a3*x*y^3*z^2 + 14*a2*a3*x*y^3*z^2 -  7*a1*a4*x*y^3*z^2 - a2^3*y^4*z^2 - 4*a1*a2*a3*y^4*z^2 + 7*a3^2*y^4*z^2 +  3*a1^2*a4*y^4*z^2 + 5*a2*a4*y^4*z^2 - 13*a6*y^4*z^2 + x^3*z^3 -  a1*x^2*y*z^3 - a1^2*x*y^2*z^3 - 9*a2*x*y^2*z^3 + a1^3*y^3*z^3 +  5*a1*a2*y^3*z^3 - 13*a3*y^3*z^3 + 7*y^2*z^4 +  w*(a3*x^3 + a1*a3*x^2*y + a4*x^2*y + a2*a3*x*y^2 + a1*a4*x*y^2 + a2*a4*y^3 +    a6*y^3 + a1*x^2*z + a1^2*x*y*z + a1*a2*y^2*z + a3*y^2*z + y*z^2));
return K3;
end function;

function KummerEll6(field, vec)
a1 := vec[1];
a2 := vec[2];
a3 := vec[3];
a4 := vec[4];
a6 := vec[5];
PW6<a,b,c,d> := ProjectiveSpace(field,[1,1,2,4]);
K6 := Scheme(PW6, a^8*a1^3*a3^3 + 3*a^8*a1^2*a3^2*a4 + 3*a^8*a1*a3*a4^2 + a^8*a4^3 +  a^8*a2^3*a6 - 4*a^8*a1*a2*a3*a6 - 4*a^8*a2*a4*a6 + 7*a^8*a6^2 +  a^7*a1^2*a2^2*a3*b + a^7*a2^3*a3*b - 3*a^7*a1^3*a3^2*b - 4*a^7*a1*a2*a3^2*b +  a^7*a1*a2^2*a4*b - 6*a^7*a1^2*a3*a4*b - 4*a^7*a2*a3*a4*b - 3*a^7*a1*a4^2*b +  5*a^7*a1*a2*a6*b + 14*a^7*a3*a6*b - a^6*a2^3*b^2 + 3*a^6*a1^3*a3*b^2 +  10*a^6*a1*a2*a3*b^2 + 7*a^6*a3^2*b^2 + 3*a^6*a1^2*a4*b^2 + 5*a^6*a2*a4*b^2 -  13*a^6*a6*b^2 - a^5*a1^3*b^3 - 4*a^5*a1*a2*b^3 - 13*a^5*a3*b^3 + 7*a^4*b^4 +  a^6*a1^4*a2*a3^2*c + 3*a^6*a1^2*a2^2*a3^2*c - a^6*a2^3*a3^2*c -  10*a^6*a1^3*a3^3*c + 5*a^6*a1*a2*a3^3*c + 2*a^6*a1^3*a2*a3*a4*c +  7*a^6*a1*a2^2*a3*a4*c - 36*a^6*a1^2*a3^2*a4*c + 5*a^6*a2*a3^2*a4*c +  a^6*a1^2*a2*a4^2*c + 4*a^6*a2^2*a4^2*c - 42*a^6*a1*a3*a4^2*c -  16*a^6*a4^3*c - 3*a^6*a1^2*a2^2*a6*c - 16*a^6*a2^3*a6*c +  5*a^6*a1^3*a3*a6*c + 63*a^6*a1*a2*a3*a6*c - 13*a^6*a3^2*a6*c +  5*a^6*a1^2*a4*a6*c + 72*a^6*a2*a4*a6*c - 108*a^6*a6^2*c -  8*a^5*a1^2*a2^2*a3*b*c - 12*a^5*a2^3*a3*b*c + 29*a^5*a1^3*a3^2*b*c +  41*a^5*a1*a2*a3^2*b*c - 13*a^5*a3^3*b*c - 6*a^5*a1*a2^2*a4*b*c +  65*a^5*a1^2*a3*a4*b*c + 54*a^5*a2*a3*a4*b*c + 36*a^5*a1*a4^2*b*c -  4*a^5*a1^3*a6*b*c - 54*a^5*a1*a2*a6*b*c - 162*a^5*a3*a6*b*c +  a^4*a1^4*a2*b^2*c + 5*a^4*a1^2*a2^2*b^2*c + 12*a^4*a2^3*b^2*c -  26*a^4*a1^3*a3*b^2*c - 81*a^4*a1*a2*a3*b^2*c - 40*a^4*a3^2*b^2*c -  28*a^4*a1^2*a4*b^2*c - 54*a^4*a2*a4*b^2*c + 162*a^4*a6*b^2*c +  9*a^3*a1^3*b^3*c + 36*a^3*a1*a2*b^3*c + 108*a^3*a3*b^3*c - 54*a^2*b^4*c -  4*a^4*a1^4*a2*a3^2*c^2 - 10*a^4*a1^2*a2^2*a3^2*c^2 + 12*a^4*a2^3*a3^2*c^2 +  32*a^4*a1^3*a3^3*c^2 - 45*a^4*a1*a2*a3^3*c^2 + 7*a^4*a3^4*c^2 -  10*a^4*a1^3*a2*a3*a4*c^2 - 36*a^4*a1*a2^2*a3*a4*c^2 +  152*a^4*a1^2*a3^2*a4*c^2 - 54*a^4*a2*a3^2*a4*c^2 - 6*a^4*a1^2*a2*a4^2*c^2 -  24*a^4*a2^2*a4^2*c^2 + 216*a^4*a1*a3*a4^2*c^2 + 96*a^4*a4^3*c^2 +  3*a^4*a1^4*a2*a6*c^2 + 36*a^4*a1^2*a2^2*a6*c^2 + 96*a^4*a2^3*a6*c^2 -  45*a^4*a1^3*a3*a6*c^2 - 324*a^4*a1*a2*a3*a6*c^2 + 162*a^4*a3^2*a6*c^2 -  54*a^4*a1^2*a4*a6*c^2 - 432*a^4*a2*a4*a6*c^2 + 648*a^4*a6^2*c^2 +  a^3*a1^6*a3*b*c^2 + 9*a^3*a1^4*a2*a3*b*c^2 + 40*a^3*a1^2*a2^2*a3*b*c^2 +  48*a^3*a2^3*a3*b*c^2 - 82*a^3*a1^3*a3^2*b*c^2 - 108*a^3*a1*a2*a3^2*b*c^2 +  108*a^3*a3^3*b*c^2 + a^3*a1^5*a4*b*c^2 + 8*a^3*a1^3*a2*a4*b*c^2 +  24*a^3*a1*a2^2*a4*b*c^2 - 228*a^3*a1^2*a3*a4*b*c^2 - 216*a^3*a2*a3*a4*b*c^2 -  144*a^3*a1*a4^2*b*c^2 + 36*a^3*a1^3*a6*b*c^2 + 216*a^3*a1*a2*a6*b*c^2 +  648*a^3*a3*a6*b*c^2 - 5*a^2*a1^4*a2*b^2*c^2 - 28*a^2*a1^2*a2^2*b^2*c^2 -  48*a^2*a2^3*b^2*c^2 + 69*a^2*a1^3*a3*b^2*c^2 + 216*a^2*a1*a2*a3*b^2*c^2 +  84*a^2*a1^2*a4*b^2*c^2 + 216*a^2*a2*a4*b^2*c^2 - 648*a^2*a6*b^2*c^2 -  18*a*a1^3*b^3*c^2 - 72*a*a1*a2*b^3*c^2 - 216*a*a3*b^3*c^2 + 108*b^4*c^2 +  a^2*a1^6*a3^2*c^3 + 7*a^2*a1^4*a2*a3^2*c^3 - 48*a^2*a2^3*a3^2*c^3 -  29*a^2*a1^3*a3^3*c^3 + 144*a^2*a1*a2*a3^3*c^3 - 54*a^2*a3^4*c^3 +  3*a^2*a1^5*a3*a4*c^3 + 32*a^2*a1^3*a2*a3*a4*c^3 + 80*a^2*a1*a2^2*a3*a4*c^3 -  252*a^2*a1^2*a3^2*a4*c^3 + 216*a^2*a2*a3^2*a4*c^3 + 2*a^2*a1^4*a4^2*c^3 +  24*a^2*a1^2*a2*a4^2*c^3 + 64*a^2*a2^2*a4^2*c^3 - 480*a^2*a1*a3*a4^2*c^3 -  256*a^2*a4^3*c^3 - a^2*a1^6*a6*c^3 - 24*a^2*a1^4*a2*a6*c^3 -  144*a^2*a1^2*a2^2*a6*c^3 - 256*a^2*a2^3*a6*c^3 + 144*a^2*a1^3*a3*a6*c^3 +  720*a^2*a1*a2*a3*a6*c^3 - 648*a^2*a3^2*a6*c^3 + 216*a^2*a1^2*a4*a6*c^3 +  1152*a^2*a2*a4*a6*c^3 - 1728*a^2*a6^2*c^3 - 2*a*a1^6*a3*b*c^3 -  20*a*a1^4*a2*a3*b*c^3 - 64*a*a1^2*a2^2*a3*b*c^3 - 64*a*a2^3*a3*b*c^3 +  66*a*a1^3*a3^2*b*c^3 + 72*a*a1*a2*a3^2*b*c^3 - 216*a*a3^3*b*c^3 -  2*a*a1^5*a4*b*c^3 - 16*a*a1^3*a2*a4*b*c^3 - 32*a*a1*a2^2*a4*b*c^3 +  264*a*a1^2*a3*a4*b*c^3 + 288*a*a2*a3*a4*b*c^3 + 192*a*a1*a4^2*b*c^3 -  72*a*a1^3*a6*b*c^3 - 288*a*a1*a2*a6*b*c^3 - 864*a*a3*a6*b*c^3 +  a1^6*b^2*c^3 + 12*a1^4*a2*b^2*c^3 + 48*a1^2*a2^2*b^2*c^3 + 64*a2^3*b^2*c^3 -  36*a1^3*a3*b^2*c^3 - 144*a1*a2*a3*b^2*c^3 + 216*a3^2*b^2*c^3 -  72*a1^2*a4*b^2*c^3 - 288*a2*a4*b^2*c^3 + 864*a6*b^2*c^3 +  4*a1^4*a2*a3^2*c^4 + 32*a1^2*a2^2*a3^2*c^4 + 64*a2^3*a3^2*c^4 -  4*a1^3*a3^3*c^4 - 144*a1*a2*a3^3*c^4 + 108*a3^4*c^4 - 4*a1^5*a3*a4*c^4 -  32*a1^3*a2*a3*a4*c^4 - 64*a1*a2^2*a3*a4*c^4 + 120*a1^2*a3^2*a4*c^4 -  288*a2*a3^2*a4*c^4 - 4*a1^4*a4^2*c^4 - 32*a1^2*a2*a4^2*c^4 -  64*a2^2*a4^2*c^4 + 384*a1*a3*a4^2*c^4 + 256*a4^3*c^4 + 4*a1^6*a6*c^4 +  48*a1^4*a2*a6*c^4 + 192*a1^2*a2^2*a6*c^4 + 256*a2^3*a6*c^4 -  144*a1^3*a3*a6*c^4 - 576*a1*a2*a3*a6*c^4 + 864*a3^2*a6*c^4 -  288*a1^2*a4*a6*c^4 - 1152*a2*a4*a6*c^4 + 1728*a6^2*c^4 +  a*(a^3*a1*a2*a3 + a^3*a2*a4 + a^3*a6 + a^2*a1*a2*b + a^2*a3*b + a*b^2 +    a*a1^3*a3*c + a*a3^2*c + a*a1^2*a4*c + a1^3*b*c)*d + d^2 );
// K6 := Scheme(PW6,d^2 - a3^2*x^6 - 4*a6*x^6 + 2*a1*a3^2*x^5*y - 2*a3*a4*x^5*y + 12*a1*a6*x^5*y - a1^2*a3^2*x^4*y^2 + 2*a2*a3^2*x^4*y^2 + 4*a1*a3*a4*x^4*y^2 - a4^2*x^4*y^2 -  12*a1^2*a6*x^4*y^2 + 12*a2*a6*x^4*y^2 - 2*a1*a2*a3^2*x^3*y^3 + 4*a3^3*x^3*y^3 - 2*a1^2*a3*a4*x^3*y^3 + 4*a2*a3*a4*x^3*y^3 + 2*a1*a4^2*x^3*y^3 + 4*a1^3*a6*x^3*y^3 -  24*a1*a2*a6*x^3*y^3 + 18*a3*a6*x^3*y^3 - a2^2*a3^2*x^2*y^4 - 4*a1*a2*a3*a4*x^2*y^4 + 12*a3^2*a4*x^2*y^4 - a1^2*a4^2*x^2*y^4 + 2*a2*a4^2*x^2*y^4 + 12*a1^2*a2*a6*x^2*y^4 -  12*a2^2*a6*x^2*y^4 - 18*a1*a3*a6*x^2*y^4 + 18*a4*a6*x^2*y^4 - 2*a2^2*a3*a4*x*y^5 - 2*a1*a2*a4^2*x*y^5 + 12*a3*a4^2*x*y^5 + 12*a1*a2^2*a6*x*y^5 - 18*a2*a3*a6*x*y^5 -  18*a1*a4*a6*x*y^5 - a2^2*a4^2*y^6 + 4*a4^3*y^6 + 4*a2^3*a6*y^6 - 18*a2*a4*a6*y^6 + 27*a6^2*y^6 + 2*a1*a3*x^5*z + 4*a4*x^5*z - 4*a1^2*a3*x^4*y*z + 4*a2*a3*x^4*y*z -  10*a1*a4*x^4*y*z + 2*a1^3*a3*x^3*y^2*z - 12*a1*a2*a3*x^3*y^2*z - 6*a3^2*x^3*y^2*z + 8*a1^2*a4*x^3*y^2*z - 8*a2*a4*x^3*y^2*z - 36*a6*x^3*y^2*z + 8*a1^2*a2*a3*x^2*y^3*z -  8*a2^2*a3*x^2*y^3*z - 6*a1*a3^2*x^2*y^3*z - 2*a1^3*a4*x^2*y^3*z + 12*a1*a2*a4*x^2*y^3*z - 30*a3*a4*x^2*y^3*z + 54*a1*a6*x^2*y^3*z + 10*a1*a2^2*a3*x*y^4*z -  18*a2*a3^2*x*y^4*z - 4*a1^2*a2*a4*x*y^4*z + 4*a2^2*a4*x*y^4*z + 6*a1*a3*a4*x*y^4*z - 24*a4^2*x*y^4*z - 18*a1^2*a6*x*y^4*z + 36*a2*a6*x*y^4*z + 4*a2^3*a3*y^5*z -  2*a1*a2^2*a4*y^5*z - 18*a2*a3*a4*y^5*z + 12*a1*a4^2*y^5*z - 18*a1*a2*a6*y^5*z + 54*a3*a6*y^5*z - a1^2*x^4*z^2 - 4*a2*x^4*z^2 + 2*a1^3*x^3*y*z^2 + 8*a1*a2*x^3*y*z^2 -  6*a3*x^3*y*z^2 - a1^4*x^2*y^2*z^2 - 2*a1^2*a2*x^2*y^2*z^2 + 8*a2^2*x^2*y^2*z^2 + 24*a1*a3*x^2*y^2*z^2 + 30*a4*x^2*y^2*z^2 - 2*a1^3*a2*x*y^3*z^2 - 8*a1*a2^2*x*y^3*z^2 -  6*a1^2*a3*x*y^3*z^2 + 54*a2*a3*x*y^3*z^2 - 30*a1*a4*x*y^3*z^2 - a1^2*a2^2*y^4*z^2 - 4*a2^3*y^4*z^2 - 18*a1*a2*a3*y^4*z^2 + 27*a3^2*y^4*z^2 + 12*a1^2*a4*y^4*z^2 +  18*a2*a4*y^4*z^2 - 54*a6*y^4*z^2 + 4*x^3*z^3 - 6*a1*x^2*y*z^3 - 6*a1^2*x*y^2*z^3 - 36*a2*x*y^2*z^3 + 4*a1^3*y^3*z^3 + 18*a1*a2*y^3*z^3 - 54*a3*y^3*z^3 + 27*y^2*z^4);
return K6;
end function;

Q := GF(7,2);
E1 := EllipticCurve([1, 1, 1, -135, -660]);
E2 := EllipticCurve([1, -1, 1, -6, -4]);
vec1 := [1, 1, 1, -135, -660];
vec2 := [1, -1, 1, -6, -4];
K2, phi := KummerEll2(Q, vec1, vec2);
Am<x1,z1,x2,z2,w11,w12,w22> := Ambient(K2);
IrreducibleComponents(Scheme(K2,[x1,z1,x2,z2]));
IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(K2)));
BaseScheme(phi);
E1xE2 := Domain(phi);
pts2 := SingularPoints(K2);
gen1 := Generators(E1);
tt1 := [2*gen1[1],gen1[1],gen1[2],gen1[1]+gen1[2]];
gen2 := Generators(E2);
tt2 := [2*gen2[1],gen2[1],gen2[2],gen2[1]+gen2[2]];
E1xE2!(Coordinates(tt1[1]) cat Coordinates(tt2[2])); 
[[phi(Coordinates(tt1[1]) cat Coordinates(tt1[2]); 
)]]
[ADEtype(pt): pt in pts2];


DirectProduct(Ambient(E1),Ambient(E2));


K := Rationals();
K<s> := GF(7^12);
K<s> := GF(3^48);
K<s> := GF(2^24);
// vec := [1, -1, 1, -1, 0];
//vec := [1, 0, 1, 0, 0]; // Bad reduction at 2, ordinary at 3
vec := [1, -1, 1, -1, -14+s*0]; // 2-rank 1, 3-rank 0, good reduction at 7
vec := [1, 0, 1, 0*s, 1]; // 2-rank 0,
E := EllipticCurve(vec);
ClassGroupPRank(E);
X3 := KummerEll3(K,vec);
spts3 := IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(X3)));
spts3;
pts3 := [X3!CoordPoint(pt): pt in spts3];
[ADEtype(pt): pt in pts3];
[IntMatDes(pt): pt in pts3];
X := Scheme(pts3[1]);
A,P := AffinePatch(X,pts3[1]);
Q<[xx]> := Ambient(A);
coords := Eltseq(P);
P := Scheme(A,[xx[i]-coords[i]: i in [1..#coords]]);
dsd := ResolveSingByBlowUp(A,P);
IntersectionMatrix(dsd);
// Does this work?

X6 := KummerEll6(K,vec);
P<a,b,c,d> := Ambient(X6);
ambspts := IrreducibleComponents(Scheme(X6,[a,b]));
ambpts := [X6!CoordPoint(s): s in ambspts];
spts6 := IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(X6)));
spts6;
pts6 := [X6!CoordPoint(s): s in spts6];
[ADEtype(pt): pt in pts6];
[IntMatDes(pt): pt in pts6];
[ADEtype(pt): pt in pts6[[1..3]]];
IntMatDes(pts6[4]);
dsd, S := DesInfo(pts6[4]);
C := BlowUpDivisor(S, dsd, 3);
G,sngs,multi_ints,min := MinimalDualResolutionGraph(dsd);
vG := VertexSet(G);
multi_ints;
[SelfIntersection(G,i): i in [1..3]];
[ArithmeticGenus(G,i): i in [1..3]];
[CanonicalIntersection(G,i): i in [1..3]];
[Genus(G,i): i in [1..3]];

K := GF(2);
P3<k1,k2,k3,k4> := ProjectiveSpace(K,3);
X := Scheme(P3,k1^4 + k2^4 + k2^3*k3 + k1^2*k3^2 + k1*k2*k3^2 + k3^3*k4 + k2^2*k4^2);
pt := SingularPoints(X)[1];
A,P := AffinePatch(X,pt);
Q<[xx]> := Ambient(A);
coords := Eltseq(P);
P := IrreducibleComponents(Scheme(A,[xx[i]-coords[i]: i in [1..#coords]]))[1];
A;
P;

R<x1,z1,x2,z2,w11,w12,w22> := PolynomialRing(Rationals(),7);
B := [];
Zwts := [[1,1,0,0,2,1,0],[0,0,1,1,0,1,2]];
Qwts := [];
C := CoxRing(R,B,Zwts,Qwts);