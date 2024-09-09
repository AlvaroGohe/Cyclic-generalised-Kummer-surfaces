function ADEtype(pt)
_,F := IsHypersurfaceSingularity(pt,3);
_,_,typ := NormalFormOfHypersurfaceSingularity(F);
return typ;
end function;

load "Functions.m";

A6<x1,x2,x3,y1,y2,y3>:=AffineSpace(Rationals(),6);
X := Scheme(A6,[27*y1^2 - 81*y2 - 18*y2^2 - y2^3 + 54*y1*y3 + 27*y3^2, 3 + 3*x3 - y2, 243*x2^3 + 27*x2^3*y2 + 72*y2^3 + 8*y2^4 - 216*y1*y2*y3 - 216*y2*y3^2, 
 27*x2^3*y1 + 8*y1*y2^3 + 27*x2^3*y3 - 72*y2^2*y3, -314928*x2^3 + 729*x2^6 - 93312*y2^3 - 64*y2^6 + 279936*y1*y2*y3 - 31104*y1*y2^2*y3 + 3456*y1*y2^3*y3 + 279936*y2*y3^2 - 
  31104*y2^2*y3^2 + 1728*y2^3*y3^2, -3*x2^2 + 4*x1*y2, 324*x1*x2 + 27*x2^3 + 72*y2^2 + 8*y2^3 - 216*y1*y3 - 216*y3^2, 6*x1^2*y1 + x2*y1*y2 + 6*x1^2*y3 - 9*x2*y3, 
 108*x1^3 + 27*x2^3 + 4*y2^3 - 108*y3^2]);

 Q0<r> := ext<Rationals() | Polynomial([1,-1,1])>;
 // Q0 := Rationals();
 P5<c1,c2,c3,c4,c5,c6> := ProjectiveSpace(Q0,5);
 K3 := Scheme(P5,[11*(2304*c1^2 - 3168*c1*c2 + 4752*c2^2 + 1344*c4*c5 - 6032*c5^2 - 1728*c3*c6 + 4464*c4*c6 - 195*c6^2), 
 -2168*c1^2 + 2916*c1*c2 - 4374*c2^2 - 2088*c4*c5 + 4524*c5^2 + 2016*c3*c6 - 4038*c4*c6 - 715*c5*c6, -c4^2 + c3*c5]);
 SSK3 := SingularSubscheme(K3);  
 spK3 := SingularPoints(K3);
 [ADEtype(pt): pt in spK3]; 
 LinesA1s := IrreducibleComponents(Scheme(K3,[(-36*c3)/169 + c5, (-6*c3)/13 + c4]));
 Intersection(LinesA1s[1],LinesA1s[2]);
 P4<k1,k2,k3,k4,k5>:=ProjectiveSpace(Q0,4);
 phi := map<K3->P4|[c1, c2, c3 - 169/144*c6, c4 - 13/24*c6, c5 - 1/4*c6]>;
 // X := Image(phi);
 X := Surface(P4,[k2^3*k3 - 13/3*k2^3*k4 + 32/3*k1*k3*k4^2 - 32/27*k1*k4^3 + 32/9*k2*k4^3 + 169/36*k2^3*k5 - 32/3*k1*k3^2*k5 - 40/9*k1*k3*k4*k5 - 4*k2*k3*k4*k5 - 208/9*k1*k4^2*k5 - 52/9*k2*k4^2*k5 + 1070/27*k1*k3*k5^2 + 50/9*k2*k3*k5^2 + 208/27*k1*k4*k5^2 + 65/9*k2*k4*k5^2 - 53911/1458*k1*k5^3 - 4901/486*k2*k5^3, k1*k2*k3 - 3/2*k2^2*k3 - 13/3*k1*k2*k4 + 13/2*k2^2*k4 - 24*k3*k4^2 + 8/3*k4^3 + 169/36*k1*k2*k5 - 169/24*k2^2*k5 + 24*k3^2*k5 + 10*k3*k4*k5 + 52*k4^2*k5 - 535/6*k3*k5^2 - 52/3*k4*k5^2 + 53911/648*k5^3, k1^2 - 27/19*k1*k2 + 81/38*k2^2 - 36/19*k4^2 + 36/19*k3*k5 - 261/76*k5^2]);
 KodairaEnriquesType(X);
 [ADEtype(pt): pt in SingularPoints(X)];

// An example with good reduction at 3
 Q0 := Rationals();
 P5<c1,c2,c3,c4,c5,c6> := ProjectiveSpace(Q0,5);
 K3 := Scheme(P5,[-258*c1^2 + 1024*c1*c2 - 1046*c2^2 - 232*c4*c5 - 616*c5^2 + 98*c3*c6 + 326*c4*c6 - 53*c6^2, -211*c1^2 + 752*c1*c2 - 793*c2^2 - 124*c4*c5 + 124*c5^2 + 67*c3*c6 + 145*c4*c6 - 
  212*c5*c6, -c4^2 + c3*c5]);
 SSK3 := SingularSubscheme(K3);  
 spK3 := SingularPoints(K3);
 spK3ss := PointsOverSplittingField(SSK3);
 [ADEtype(pt): pt in spK3]; 
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