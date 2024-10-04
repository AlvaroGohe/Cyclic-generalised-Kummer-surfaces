function ADEtype(pt)
_,F := IsHypersurfaceSingularity(pt,3);
P<a,b,c>:= Parent(F);
_,_,typ := NormalFormOfHypersurfaceSingularity(F);
return typ, F;
end function;




load "Functions.m";
load "WeightedJacobian.m";
A6<x1,x2,x3,y1,y2,y3>:=AffineSpace(Rationals(),6);
X := Scheme(A6,[27*y1^2 - 81*y2 - 18*y2^2 - y2^3 + 54*y1*y3 + 27*y3^2, 3 + 3*x3 - y2, 243*x2^3 + 27*x2^3*y2 + 72*y2^3 + 8*y2^4 - 216*y1*y2*y3 - 216*y2*y3^2, 
 27*x2^3*y1 + 8*y1*y2^3 + 27*x2^3*y3 - 72*y2^2*y3, -314928*x2^3 + 729*x2^6 - 93312*y2^3 - 64*y2^6 + 279936*y1*y2*y3 - 31104*y1*y2^2*y3 + 3456*y1*y2^3*y3 + 279936*y2*y3^2 - 
  31104*y2^2*y3^2 + 1728*y2^3*y3^2, -3*x2^2 + 4*x1*y2, 324*x1*x2 + 27*x2^3 + 72*y2^2 + 8*y2^3 - 216*y1*y3 - 216*y3^2, 6*x1^2*y1 + x2*y1*y2 + 6*x1^2*y3 - 9*x2*y3, 
 108*x1^3 + 27*x2^3 + 4*y2^3 - 108*y3^2]);

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




// An example with good reduction at 3
 Q0 := Rationals();
 Q0pol<x> := PolynomialRing(Q0);
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

 Q0 := Rationals();
 Q0 := GF(5);
 Q0pol<x> := PolynomialRing(Q0);
 C3 := HyperellipticCurve(x^3, 1 + x^3);
 Jac3 := GeneralJacobianSurface(C3);
 P15<v1, v2, v3, v4, v5, v6, k11, k12, k13, k14, k22, k23, k24, k33, k34, k44> := Ambient(Jac3);
 P5<c11, c12, c22, b1, b2, kvw> := ProjectiveSpace(Q0,5);
 proj2 := map<Jac3->P5 |[k22, k24, k44, k13 + k22 + 2*v1, k13 + k24 + 2*v4, k13]>;
 Proj2X3 := Scheme(P5,MinimalBasis(Image(proj2)));
 P4<x1,x2,x3,x4,x5> :=  ProjectiveSpace(Q0,4);
 proj4 := map<Jac3->P4|[k13 - k22, k13 + k24 + 2*v4, k13 + k22 + 2*v1, -k22 + k44, k22 + k24]>;
 X34 := Scheme(P4,MinimalBasis(Image(proj4)));
 WJ := WeightedJacobian(Jacobian(C3));
 P<k1, k2, k3, k4, v1, v2, v3, v4, v5, v6> := Ambient(WJ);
 PW<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := ProjectiveSpace(Q0, [1,1,2,2,2,3,3,3,3]);
 Wpi := map<WJ->PW | [k2, k4, k2^2 + k1*k3 + 2*v1, k1*k3 + k2*k4 + 2*v4, k1*k3, k1^3, k3^3, k1*(k1^2 + k2*k3 + 2*v2), k3*(k1*k2 + k3^2 + 2*v3)]>;
 J := ideal< CoordinateRing(PW) | &cat[ DefiningEquations(Image(Wpi,WJ,d)) : d in  [1..6]] >;
 MinimalBasis(J);
 X3 := Scheme(PW,MinimalBasis(J));
 Wpi := map<WJ->X3 | [k2, k4, k2^2 + k1*k3 + 2*v1, k1*k3 + k2*k4 + 2*v4, k1*k3, k1^3, k3^3, k1*(k1^2 + k2*k3 + 2*v2), k3*(k1*k2 + k3^2 + 2*v3)]>;
 P1<x1,x2> := ProjectiveSpace(Q0,1);
 fib := map<X3->P1|[c1,c2]>;
 G1 := Pullback(fib,P1![1,-2]);
 P13<c111, c112, c122, c222, c1d1, c1d2, c1vw, c2d1, c2d2, c2vw, kvvv, kwww, kvdw, kwdv> := ProjectiveSpace(Q0, 13);
 proj3 := map<X3->P13 | [c1^3, c1^2*c2, c1*c2^2, c2^3, c1*d1, c1*d2, c1*cvw, c2*d1, c2*d2, c2*cvw, cvvv, cwww, cvdw, cwdv]>;
 Proj3X3 := Scheme(P13, MinimalBasis(Image(proj3)));
 proj3 := map<X3->Proj3X3 | [c1^3, c1^2*c2, c1*c2^2, c2^3, c1*d1, c1*d2, c1*cvw, c2*d1, c2*d2, c2*cvw, cvvv, cwww, cvdw, cwdv]>;
 Fib3 := map<Proj3X3->P1 |[c111,c112]>;
 G1 := Pullback(Fib3,P1![1,-2]);
 proj2 := map<X3->P5 |[c1^2, c1*c2, c2^2, d1, d2, cvw]>;
 // Proj2X3 := Scheme(P5,MinimalBasis(Image(proj2)));
 Proj2X3 := Scheme(P5,[c11*c22 + 4*b2^2 + 3*c12*kvw + c22*kvw + kvw^2, c12^2 + 4*b2^2 + 3*c12*kvw + c22*kvw + kvw^2, c11^2 + 3*c11*c12 + 4*b1^2 + 2*b1*b2 + c12*kvw + 3*kvw^2]);
 sPtX3 := [Pullback(proj2, pt): pt in SingularPoints(Proj2X3)];
 setPtX3:= SetToSequence(SequenceToSet(&cat [IrreducibleComponents(ReducedSubscheme(pt)): pt in sPtX3]));
 ArithmeticGenus(setPtX3[6]);
 sppts := [X3![0,0,0,0,0,4,0,1,0], X3![1,0,-1,0,0,0,0,0,0], X3![1,0,1,0,0,0,0,0,0], X3![2,2,1,1,0,0,0,0,0], X3![1,1,1,1,0,0,0,0,0], X3![0,0,0,0,0,0,1,0,-1], X3![0,0,0,0,0,0,1,0,1], X3![0,0,0,0,0,1,0,-1,0],X3![0,1,0,0,0,0,0,0,0] ];
 notsing := [X3![0,0,0,0,0,4,0,1,0],X3![0,0,0,0,0,0,1,0,-1], X3![0,0,0,0,0,0,1,0,1], X3![0,0,0,0,0,1,0,-1,0]];
 singpts := [X3![1,0,-1,0,0,0,0,0,0], X3![1,0,1,0,0,0,0,0,0], X3![2,2,1,1,0,0,0,0,0], X3![1,1,1,1,0,0,0,0,0], X3![0,1,0,0,0,0,0,0,0]];
 [IsSingular(Proj3X3!proj3(pt)): pt in singpts];
 // [ADEtype(Proj3X3!proj3(pt)): pt in singpts]; These are all A2 but take a long time to evaluate
 [IsSingular(Proj3X3!proj3(pt)): pt in notsing];
 // Wjnotsing := [Pullback(Wpi,pt): pt in pts];
 schsingpts := &cat [IrreducibleComponents(ReducedSubscheme(Pullback(Wpi,pt))):pt in singpts];
 singpts3WJ := [WJ![0,1,0,0,-1,0,0,0,0,0], WJ![0,1,0,0,0,0,0,0,0,0], WJ![0,1,0,1,-1,0,0,-1,0,0], WJ![0,1,0,1,0,0,0,0,0,0], WJ![0,0,0,1,0,0,0,0,0,0]];
 schnotsing := &cat [IrreducibleComponents(Scheme(WJ, [k3^3, k3*(k1*k2 + k3^2 + 2*v3), 4*k1^3+k1*(k1^2 + k2*k3 + 2*v2),k2,k4,k2^2 + k1*k3 + 2*v1])), IrreducibleComponents(Scheme(WJ, [k1^3, k1*(k1^2 + k2*k3 + 2*v2), k3^3-k3*(k1*k2 + k3^2 + 2*v3),k2,k4,k2^2 + k1*k3 + 2*v1])), IrreducibleComponents(Scheme(WJ, [k1^3, k1*(k1^2 + k2*k3 + 2*v2), k3^3+k3*(k1*k2 + k3^2 + 2*v3),k2,k4,k2^2 + k1*k3 + 2*v1])), IrreducibleComponents(Scheme(WJ, [k3^3, k3*(k1*k2 + k3^2 + 2*v3), k1^3+k1*(k1^2 + k2*k3 + 2*v2),k2,k4,k2^2 + k1*k3 + 2*v1]))];
 notsingpts3WJ :=[WJ![1,0,0,0,0,0,0,0,-1,0], WJ![0,0,1,0,0,0,0,0,0,-1], WJ![0,0,1,0,0,0,-1,0,0,0], WJ![1,0,0,0,0,-1,0,0,0,0]];
 pts3WJ := singpts3WJ cat notsingpts3WJ;
 IC2 := IrreducibleComponents(Scheme(X3,[c1, c2, d1, d2, cvw]));
 [Dimension(V): V in IC2];
 [Degree(V): V in IC2];
 HSWJ := HilbertSeries(Ideal(DefiningPolynomials(WJ)));
 HSX3 := HilbertSeries(Ideal(DefiningPolynomials(X3)));
 FindFirstGenerators(HSWJ);
 HilbertNumerator(HSWJ,[1,1,1,1,2,2,2,2,2,2]);
 Ser<t> := PowerSeriesRing(Rationals(),8);
 Ser!HSX3;

 // What happens to the Hilbert series of ALL negative elements?
 PW2<w1, w2, w3, w4, w5, w6, w11, w12, w13, w14, w15, w16, w21, w22, w23, w24, w25, w26, w31, w32, w33, w34> := ProjectiveSpace(Q0, [2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]);
 projsp := map< WJ -> PW2 | [v1, v2, v3, v4, v5, v6, k1*v1, k1*v2, k1*v3, k1*v4, k1*v5, k1*v6, k2*v1, k2*v2, k2*v3, k2*v4, k2*v5, k2*v6, k3*v1, k3*v2, k3*v3, k3*v4] >;
 // ideal< CoordinateRing(PW2) | &cat[ DefiningEquations(Image(projsp,WJ,d)) : d in  [1..4]]>;


 PWw<cw, c1v, c2v, dv1, dv2, c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2, c22dw2> := ProjectiveSpace(Q0, [1,2,2,2,2,3,3,3,3,3,3,3,3,3,4]);
 Wpiw := map<WJ->PWw | [k3, k1*k2, k1*k4, k3^2 + 2*v3, 2*v6, k2^3, k2^2*k4, k2*k4^2, k4^3, k1^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k2*k4 + 2*v4), k1*(k2*k3 + 2*v2), k1*(k3*k4 + 2*v5),  k4^2*(k3*k4 + 2*v5)]>;
 Jw := ideal< CoordinateRing(PWw) | &cat[ DefiningEquations(Image(Wpiw,WJ,d)) : d in  [1..8]] >;
 ideal< CoordinateRing(PWw) | &cat[ DefiningEquations(Image(Wpiw,WJ,d)) : d in  [1..4]] >;
 MinimalBasis(Jw);
 X3w := Scheme(PWw,MinimalBasis(Jw));
 Wpiw := map<WJ->X3w | [k3, k1*k2, k1*k4, k3^2 + 2*v3, 2*v6, k2^3, k2^2*k4, k2*k4^2, k4^3, k1^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k2*k4 + 2*v4), k1*(k2*k3 + 2*v2), k1*(k3*k4 + 2*v5),  k4^2*(k3*k4 + 2*v5)]>;
 HSX3w  := HilbertSeries(Ideal(DefiningPolynomials(X3w)));
 [Wpiw(pts3WJ[7]), Wpiw(pts3WJ[8])];
 P4<qww, q1v, q2v, tv1, tv2> := ProjectiveSpace(Q0,4);
 proj2w := map<X3w->P4 |[cw^2, c1v, c2v, dv1, dv2]>;
 Proj2X3w := Scheme(P4,MinimalBasis(ideal< CoordinateRing(P4) | &cat[ DefiningEquations(Image(proj2w,X3w,d)) : d in  [1..4]] >));
 proj2w := map<X3w->Proj2X3w |[cw^2, c1v, c2v, dv1, dv2]>;
 sPtX3w := [Pullback(proj2w, pt): pt in SingularPoints(Proj2X3w)];
 setPtX3w:= SetToSequence(SequenceToSet(&cat [IrreducibleComponents(ReducedSubscheme(pt)): pt in sPtX3w]));
 spptsw := [PWw![0, 0, 0, 0, 0, 4, 0, 0, 0, 0, -1, 0, 0, 0, 0], PWw![0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0], PWw![0, 0, 0, 0, 0, 4, 4, 4, 4, 0, 4, -1, 0, 0, 0], PWw![0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 4, -1, 0, 0, 0], PWw![0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, -1, 0, 0], PWw![0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, -1, 0], PWw![0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], PWw![3, 0, 0, 3, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], PWw![1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ];
 ICw := IrreducibleComponents(Scheme(X3w,[cw, c1v, c2v, dv1, dv2]));
 IrreducibleComponents(Scheme(X3w,[cw, c1v, c2v, dv1, dv2, c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2]));
 rpt := IrreducibleComponents(ReducedSubscheme(Scheme(Proj2X3w,[qww-q1v, q2v+3*tv1])));
 ptsP2X3w:=Points(Proj2X3w); 
 IrreducibleComponents(Pullback(proj2w,ptsP2X3w[2]));
 P13<kwww, k1wv, k2wv, kwdv1, kwdv2, k111, k112, k122, k222, kvvv, k1d1, k1d2, kvdw1, kvdw2> := ProjectiveSpace(Q0,13);
 proj3w := map<X3w->P13 |[cw^3, cw*c1v, cw*c2v, cw*dv1, cw*dv2,c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2]>;
 Proj3X3w := Scheme(P13,MinimalBasis(ideal< CoordinateRing(P13) | &cat[ DefiningEquations(Image(proj3w,X3w,d)) : d in  [1..2]] >));


 PWv<cv, c1w, c2w, dw1, dw2, c111, c112, c122, c222, cwww, c1d1, c1d2, cwdv1, cwdv2, c22dv2>  := ProjectiveSpace(Q0, [1,2,2,2,2,3,3,3,3,3,3,3,3,3,4]);
 Wpiv := map<WJ->PWv | [k1, k2*k3, k3*k4, k2*k3 + 2*v2, k3*k4 + 2*v5, k2^3, k2^2*k4, k2*k4^2, k4^3, k3^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k2*k4 + 2*v4), k3*(k3^2 + 2*v3), 2*k3*v6, 2*k4^2*v6]>;
 Jv := ideal< CoordinateRing(PWv) | &cat[ DefiningEquations(Image(Wpiv,WJ,d)) : d in  [1..8]] >;
 ideal< CoordinateRing(PWv) | &cat[ DefiningEquations(Image(Wpiv,WJ,d)) : d in  [1..4]] >;
 MinimalBasis(Jv);
 X3v := Scheme(PWv,MinimalBasis(Jv));
 Wpiv := map<WJ->PWv | [k1, k2*k3, k3*k4, k2*k3 + 2*v2, k3*k4 + 2*v5, k2^3, k2^2*k4, k2*k4^2, k4^3, k3^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k2*k4 + 2*v4), k3*(k3^2 + 2*v3), 2*k3*v6, 2*k4^2*v6]>;
 HSX3v  := HilbertSeries(Ideal(DefiningPolynomials(X3v)));
 [Wpiw(pts3WJ[6]), Wpiw(pts3WJ[9])];
 
 HSP12 := 1/((1-t)^4*(1-t^2)^6);



// This is how I found the invariant of degree 4 that I was missing
 pols := [k1^2*k2^2, k1^2*k2*k4, k1^2*k4^2, k2^3*k3, k2^2*k3*k4, k2*k3*k4^2, k3*k4^3, k1^3*k3, k1*k2*k3^2, k1*k3^2*k4, k3^4, k1^2*(k2^2 + k1*k3 + 2*v1), k2*k3*(k2^2 + k1*k3 + 2*v1),  k3*k4*(k2^2 + k1*k3 + 2*v1), k1^2*(k2*k4 + 2*v4), k2*k3*(k2*k4 + 2*v4), k3*k4*(k2*k4 + 2*v4), k1*k2*(k3^2 + 2*v3), k1*k4*(k3^2 + 2*v3), k3^2*(k3^2 + 2*v3), (k3^2 + 2*v3)^2,  2*k1*k2*v6, 2*k1*k4*v6, 2*k3^2*v6, 2*(k3^2 + 2*v3)*v6, 4*v6^2, k2^2*(k2*k3 + 2*v2), k2*k4*(k2*k3 + 2*v2), k4^2*(k2*k3 + 2*v2), k1*k3*(k2*k3 + 2*v2),  (k2^2 + k1*k3 + 2*v1)*(k2*k3 + 2*v2), (k2*k3 + 2*v2)*(k2*k4 + 2*v4), k2^2*(k3*k4 + 2*v5), k2*k4*(k3*k4 + 2*v5), k4^2*(k3*k4 + 2*v5), k1*k3*(k3*k4 + 2*v5),  (k2^2 + k1*k3 + 2*v1)*(k3*k4 + 2*v5), (k2*k4 + 2*v4)*(k3*k4 + 2*v5)];
 schpols := [Scheme(WJ,pol): pol in pols];
 schpolsim := [Wpi(sc): sc in schpols];
 PWw<cw, c1v, c2v, dv1, dv2, c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2, a> := ProjectiveSpace(Q0, [1,2,2,2,2,3,3,3,3,3,3,3,3,3,4]);
 colmaps:=[map<WJ->PWw | [k3, k1*k2, k1*k4, k3^2 + 2*v3, 2*v6, k2^3, k2^2*k4, k2*k4^2, k4^3, k1^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k2*k4 + 2*v4), k1*(k2*k3 + 2*v2), k1*(k3*k4 + 2*v5),  pol]>: pol in pols];
 [<#Basis(ideal< CoordinateRing(PWw) | &cat[ DefiningEquations(Image(map,WJ,d)) : d in  [1..4]] >),DefiningPolynomials(map)[15]>: map in colmaps];
 Wpi := map<WJ->PWw | [k3, k1*k2, k1*k4, k3^2 + 2*v3, 2*v6, k2^3, k2^2*k4, k2*k4^2, k4^3, k1^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k2*k4 + 2*v4), k1*(k2*k3 + 2*v2), k1*(k3*k4 + 2*v5),  k1^2*(k2^2 + k1*k3 + 2*v1)]>;
 ideal< CoordinateRing(PWw) | &cat[ DefiningEquations(Image(Wpi,WJ,d)) : d in  [1..4]] >;
 Jw := ideal< CoordinateRing(PWw) | &cat[ DefiningEquations(Image(Wpi,WJ,d)) : d in  [1..4]] >;
 MinimalBasis(Jw);
 X3w := Scheme(PWw,MinimalBasis(Jw));
 HSX3w  := HilbertSeries(Ideal(DefiningPolynomials(HSX3w)));

HSKum  := HilbertSeries(Ideal(DefiningPolynomials(GeneralKummerSurface(J3))));
 [
8*d1^2 + 4*d2^2 + c2*cvvv + 7*c2*cwww,
2*c1*c2*d1 + c1^2*d2 + 10*d2*cvw + 7*c2*cvdw,
c1^2*d1 + 10*c1*c2*d2 + 10*d1*cvw + c2*cwdv,
c1^2*c2^2 + 7*c2^2*cvw + 10*d1^2 + 6*c2*cwww,
c1^3*c2 + 8*c1*c2*cvw + 10*d1*d2,
c1^4 + 9*c1^2*cvw + 10*d2^2 + cvw^2 + c2*cwww,
c1*c2*cvvv + 7*c1*c2*cwww + 5*d1*cvdw + 4*d2*cwdv,
9*c1*d2*cvw + c1^2*cwdv + 10*d1*cwww + 10*cvw*cwdv,
9*c2*d1*cvw + 5*c1*d2*cvw + c1*c2*cvdw + 9*d1*cvvv,
2*c1*d1*cvw + 7*c2*d2*cvw + c1*c2*cwdv + 6*d2*cwww,
10*c1*d1*cvw + c1^2*cvdw + 9*d2*cvvv + 10*cvw*cvdw,
4*c1*d1*d2 + 7*c2*d2^2 + c2*cvw^2 + c1^2*cwww + 4*c2^2*cwww + 10*cvw*cwww + 2*d1*cwdv,
5*c1*d1*d2 + 6*c2*d2^2 + 4*c2*cvw^2 + c1^2*cvvv + 5*c2^2*cwww + 10*cvw*cvvv + 5*d2*cvdw,
c1^3*cvw + 10*c1*cvw^2 + 6*c1*c2*cwww + 5*d2*cwdv,
c2*cvw*cvvv + 7*c2*cvw*cwww + cvvv^2 + 6*cwww^2 + 8*cvdw^2 + cwdv^2,
cvw^3 + 10*cvvv*cwww,
9*d2*cvw^2 + c1*cvw*cwdv + 10*cwww*cvdw,
10*d1*cvw^2 + c1*cvw*cvdw + 10*cvvv*cwdv,
c1*c2^2*cwww + 2*d1*d2*cvw + 3*c1*cvw*cvvv + c1*cvw*cwww + 10*c2*d2*cwdv + 2*cvdw*cwdv,
2*c2^2*cvw^2 + 7*c2*cvw*cwww + c1*d1*cvdw + 10*c2*d2*cvdw + 9*cvvv^2 + 5*cvvv*cwww + 8*cwww^2 + 6*cvdw^2 + 5*cwdv^2
]