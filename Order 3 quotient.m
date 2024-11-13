function ADEtype(pt)
_,F := IsHypersurfaceSingularity(pt,3);
P<a,b,c>:= Parent(F);
_,_,typ := NormalFormOfHypersurfaceSingularity(F);
return typ, F;
end function;

function P5model(X3) // Takes the model of the quotient as a subvariety of the weighted projective space with weights 1^2, 2^3, 3^4 and returns the model in P5 that we get when we consider the functions of degree 2, and the quotient map
P<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := Ambient(X3);
P5<c11, c12, c22, b1, b2, kvw> := ProjectiveSpace(Q0,5);
proj2 := map<X3->P5 |[c1^2, c1*c2, c2^2, d1, d2, cvw]>;
Proj2X3 := Scheme(P5,MinimalBasis(Image(proj2,X3,2)));
return Proj2X3, proj2;
end function;





load "Functions.m";
load "WeightedJacobian.m";

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
 PWalt<cvdw2, cwdv2,c1, c2, d1, d2, cvw, cvvv, cwww, cvdw1, cwdv1> := ProjectiveSpace(Q0, [3,3,1,1,2,2,2,3,3,3,3]);
 Wpialt := map<WJ->PWalt | [k1*(k3*k4 + 2*v5), k3*(k3^2 + 2*v6),k2, k4, k2^2 + k1*k3 + 2*v1, k1*k3 + k2*k4 + 2*v4, k1*k3, k1^3, k3^3, k1*(k2*k3 + 2*v2), k3*(k3^2 + 2*v3)]>;
 Jalt := ideal< CoordinateRing(PWalt) | &cat[ DefiningEquations(Image(Wpialt,WJ,d)) : d in  [1..6]] >;
 X3alt := Scheme(PWalt,MinimalBasis(Jalt));
 PW<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := ProjectiveSpace(Q0, [1,1,2,2,2,3,3,3,3]);
 Wpi := map<WJ->PW | [k2, k4, k2^2 + k1*k3 + 2*v1, k1*k3 + k2*k4 + 2*v4, k1*k3, k1^3, k3^3, k1*(k2*k3 + 2*v2), k3*(k3^2 + 2*v3)]>;
 J := ideal< CoordinateRing(PW) | &cat[ DefiningEquations(Image(Wpi,WJ,d)) : d in  [1..6]] >;
 MinimalBasis(J);
 X3 := Scheme(PW,MinimalBasis(J));
 [#[eqs: eqs in DefiningEquations(X3) |  (Degree(eqs) eq i)]: i in [1..6]];
 Wpi := map<WJ->X3 | [k2, k4, k2^2 + k1*k3 + 2*v1, k1*k3 + k2*k4 + 2*v4, k1*k3, k1^3, k3^3, k1*(k2*k3 + 2*v2), k3*(k3^2 + 2*v3)]>;
 //cd:=Scheme(WJ,k1*(k3*k4 + 2*v5));
 //Wpirest := map<cd->X3 | [k2, k4, k2^2 + k1*k3 + 2*v1, k1*k3 + k2*k4 + 2*v4, k1*k3, k1^3, k3^3, k1*(k2*k3 + 2*v2), k3*(k3^2 + 2*v3)]>;
 //cvdw2 := Wpi();
 //cwdv2 := Wpi(k1*(k3*k4 + 2*v5));
 // cwdv2 := Scheme(X3, &cat[ DefiningEquations(Image(Wpi,Scheme(WJ,k1*(k3*k4 + 2*v5)),d)) : d in  [1..2]]);
 //P1<x1,x2> := ProjectiveSpace(Q0,1);
 //fib := map<X3->P1|[c1,c2]>;
 //coord := [2,-1];
 //Cu := Pullback(fib,P1!coord);
 //ptO := Cu![0,0,0,0,0,1,0,1,0]; 
 //D := Divisor(ptO);
 P13<c111, c112, c122, c222, c1d1, c1d2, c1vw, c2d1, c2d2, c2vw, kvvv, kwww, kvdw, kwdv> := ProjectiveSpace(Q0, 13);
 proj3 := map<X3->P13 | [c1^3, c1^2*c2, c1*c2^2, c2^3, c1*d1, c1*d2, c1*cvw, c2*d1, c2*d2, c2*cvw, cvvv, cwww, cvdw, cwdv]>;
 //Proj3X3 := Scheme(P13, MinimalBasis(Image(proj3)));
 Proj3X3 := Scheme(P13, MinimalBasis(Image(proj3,X3,2)));
 
 proj3 := map<X3->Proj3X3 | [c1^3, c1^2*c2, c1*c2^2, c2^3, c1*d1, c1*d2, c1*cvw, c2*d1, c2*d2, c2*cvw, cvvv, cwww, cvdw, cwdv]>;
 coord := [2,-1];
 Fib3 := map<Proj3X3->P1 |[c111,c112]>;
 Fib32 := map<Proj3X3->P1 |[c122,c222]>;
 G1 := Intersection(Pullback(Fib3,P1!coord),Pullback(Fib32,P1!coord));
 Cu := Scheme(P13,DefiningEquations(IrreducibleComponents(G1)[1]));
 ptO := Cu!Points(Cu)[1]; 
 D := Divisor(ptO);
 Basis(3*D);
 BP2<x,y,z> := ProjectiveSpace(Q0, 2);
 phi3:= DivisorMap(3*D, BP2);
 E := EllipticCurve(Image(phi3), phi3(ptO));
 TorsionSubgroup(E);
 lines :=IrreducibleComponents(Scheme(Proj3X3,[kvdw]))[2..3];
 Q0t<t> :=FunctionField(Q0);
 P3X3t<c111, c112, c122, c222, c1d1, c1d2, c1vw, c2d1, c2d2, c2vw, kvvv, kwww, kvdw, kwdv> := BaseExtend(Proj3X3,Q0t);
 P1t<x1,x2> := ProjectiveSpace(Q0t,1);
 Fib1t := map<P3X3t->P1t |[c111,c112]>;
 Fib2t := map<P3X3t->P1t |[c122,c222]>;
 coordt :=[1,t];
 G1t := Intersection(Pullback(Fib1t,P1t!coordt),Pullback(Fib2t,P1t!coordt));
 // Cut := IrreducibleComponents(G1t)[1];
 // linest :=IrreducibleComponents(Scheme(P3X3t,[kvdw]))[2..3];
//  ptO := Cut!Points(IrreducibleComponents(Intersection(linest[2],Cut))[1])[1];
 // D := Divisor(ptO);
 // phi3t:= DivisorMap(3*D);
 // phi5t:= DivisorMap(5*D);
 //P2t<x1,x2,x3> := ProjectiveSpace(Q0t,2);
 //E := EllipticCurve(Image(phi3t), phi3t(ptO));
 P<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := Ambient(X3);
 isoalt := map<X3->X3alt | [-1*(4*c2*d1 + c1*d2 + cvvv + cvdw),-1*(4*c2*d1 + 4*c1*d2 + 2*c2*d2 + c2*cvw + cwww + cwdv), c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv]>;
 IrreducibleComponents(BaseScheme(isoalt));
 proj3alt := map<X3alt->Proj3X3 | [c1^3, c1^2*c2, c1*c2^2, c2^3, c1*d1, c1*d2, c1*cvw, c2*d1, c2*d2, c2*cvw, cvvv, cwww, cvdw, cwdv]>;

 P<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := Ambient(X3);
 proj2 := map<X3->P5 |[c1^2, c1*c2, c2^2, d1, d2, cvw]>;
 // Proj2X3 := Scheme(P5,MinimalBasis(Image(proj2)));
 Proj2X3 := Scheme(P5,[c11*c22 + 4*b2^2 + 3*c12*kvw + c22*kvw + kvw^2, c12^2 + 4*b2^2 + 3*c12*kvw + c22*kvw + kvw^2, c11^2 + 3*c11*c12 + 4*b1^2 + 2*b1*b2 + c12*kvw + 3*kvw^2]);
 ptsP2X3 := Points(Proj2X3);
 sPtX3 := [Pullback(proj2, pt): pt in SingularPoints(Proj2X3)];
 setPtX3:= SetToSequence(SequenceToSet(&cat [IrreducibleComponents(ReducedSubscheme(pt)): pt in sPtX3]));
 ArithmeticGenus(setPtX3[6]);
 ICpb := IrreducibleComponents(BaseScheme(proj2));
 [[Dimension(Intersection(V,W)): V in ICpb]: W in ICpb];
 [IrreducibleComponents(Intersection(V,ICpb[1])): V in ICpb];
 [<Dimension(V), ArithmeticGenus(V)>: V in ICpb];
 [[Dimension(Intersection(V,W)): W in setPtX3[1..5] cat setPtX3[7..10]]:  V in ICpb];
 //
 // sppts := [X3![0,0,0,0,0,4,0,1,0], X3![1,0,-1,0,0,0,0,0,0], X3![1,0,1,0,0,0,0,0,0], X3![2,2,1,1,0,0,0,0,0], X3![1,1,1,1,0,0,0,0,0], X3![0,0,0,0,0,0,1,0,-1], X3![0,0,0,0,0,0,1,0,1], X3![0,0,0,0,0,1,0,-1,0],X3![0,1,0,0,0,0,0,0,0] ];
 //notsing := [X3![0,0,0,0,0,4,0,1,0],X3![0,0,0,0,0,0,1,0,-1], X3![0,0,0,0,0,0,1,0,1], X3![0,0,0,0,0,1,0,-1,0]];
 //singpts := [X3![1,0,-1,0,0,0,0,0,0], X3![1,0,1,0,0,0,0,0,0], X3![2,2,1,1,0,0,0,0,0], X3![1,1,1,1,0,0,0,0,0], X3![0,1,0,0,0,0,0,0,0]];
 //[IsSingular(Proj3X3!proj3(pt)): pt in singpts];
 // [ADEtype(Proj3X3!proj3(pt)): pt in singpts]; These are all A2 but take a long time to evaluate
 //[IsSingular(Proj3X3!proj3(pt)): pt in notsing];
 //ptsP3X3 := Points(Proj3X3);
 //[Degree(Pullback(proj3,pt)): pt in ptsP3X3];



 Q0t<t> :=FunctionField(Q0);
 P2X3t<c11, c12, c22, b1, b2, kvw> := BaseExtend(Proj2X3,Q0t);
 P1t<x1,x2> := ProjectiveSpace(Q0t,1);
 Fib1t := map<P2X3t->P1t |[c11,c12]>;
 Fib2t := map<P2X3t->P1t |[c12,c22]>;
 coordt :=[1,t];
 G1t := Intersection(Pullback(Fib1t,P1t!coordt),Pullback(Fib2t,P1t!coordt));
 Cut := IrreducibleComponents(G1t)[1];
 linest:= IrreducibleComponents(Scheme(P2X3t,[kvw]));
 ptO := Cut!Points(IrreducibleComponents(Intersection(linest[2],Cut))[1])[1];
 D := Divisor(ptO);
 phi3t:= DivisorMap(3*D);
 E1, mapE1 := EllipticCurve(Image(phi3t), phi3t(ptO));
 E, mapE2 := MinimalModel(E1);
 mapE:=mapE1*mapE2;
 invmapE := Inverse(mapE);
 BadPlaces(E);
 LocalInformation(E);
 G, mapG := TorsionSubgroup(E);
 P2 := mapG(G.1);
 Inverse(mapE)(P2);

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
 phi := JacIso(WJ,Jac3);
 J3 := Jacobian(C3);
 [InverseGamma(J3, phi(pt)): pt in pts3WJ];
 [3*InverseGamma(J3, phi(pt)): pt in pts3WJ];

 



 // What happens to the Hilbert series of ALL negative elements?
 PW2<w1, w2, w3, w4, w5, w6, w11, w12, w13, w14, w15, w16, w21, w22, w23, w24, w25, w26, w31, w32, w33, w34> := ProjectiveSpace(Q0, [2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]);
 projsp := map< WJ -> PW2 | [v1, v2, v3, v4, v5, v6, k1*v1, k1*v2, k1*v3, k1*v4, k1*v5, k1*v6, k2*v1, k2*v2, k2*v3, k2*v4, k2*v5, k2*v6, k3*v1, k3*v2, k3*v3, k3*v4] >;
 // ideal< CoordinateRing(PW2) | &cat[ DefiningEquations(Image(projsp,WJ,d)) : d in  [1..4]]>;
 
 P<k1, k2, k3, k4, v1, v2, v3, v4, v5, v6> := Ambient(WJ);
 PWw<cw, c1v, c2v, dv1, dv2, c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2, c22dw2> := ProjectiveSpace(Q0, [1,2,2,2,2,3,3,3,3,3,3,3,3,3,4]);
 Wpiw := map<WJ->PWw | [k3, k1*k2, k1*k4, k3^2 + 2*v3, k3^2 + 2*v6, k2^3, k2^2*k4, k2*k4^2, k4^3, k1^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k1*k3 + k2*k4 + 2*v4), k1*(k2*k3 + 2*v2), k1*(k3*k4 + 2*v5), k4^2*(k3*k4 + 2*v5)]>;
 Jw := ideal< CoordinateRing(PWw) | &cat[ DefiningEquations(Image(Wpiw,WJ,d)) : d in  [1..8]] >;
 ideal< CoordinateRing(PWw) | &cat[ DefiningEquations(Image(Wpiw,WJ,d)) : d in  [1..4]] >;
 MinimalBasis(Jw);
 X3w := Scheme(PWw,MinimalBasis(Jw));
 Wpiw := map<WJ->X3w | [k3, k1*k2, k1*k4, k3^2 + 2*v3, k3^2 + 2*v6, k2^3, k2^2*k4, k2*k4^2, k4^3, k1^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k1*k3 + k2*k4 + 2*v4), k1*(k2*k3 + 2*v2), k1*(k3*k4 + 2*v5), k4^2*(k3*k4 + 2*v5)]>;
 HSX3w  := HilbertSeries(Ideal(DefiningPolynomials(X3w)));
 [Wpiw(pts3WJ[7]), Wpiw(pts3WJ[8])];
 P4<qww, q1v, q2v, tv1, tv2> := ProjectiveSpace(Q0,4);
 proj2w := map<X3w->P4 |[cw^2, c1v, c2v, dv1, dv2]>;
 Proj2X3w := Scheme(P4,MinimalBasis(ideal< CoordinateRing(P4) | &cat[ DefiningEquations(Image(proj2w,X3w,d)) : d in  [1..4]] >));
 proj2w := map<X3w->Proj2X3w |[cw^2, c1v, c2v, dv1, dv2]>;
 IrreducibleComponents(BaseScheme(proj2w));
 [<Dimension(V), ArithmeticGenus(V),Degree(V)>: V in IrreducibleComponents(BaseScheme(proj2w))];
 P<qww, q1v, q2v, tv1, tv2> :=Ambient(P4);
 Q<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv>:= Ambient(X3);
 proj2walt := map<X3->Proj2X3w |[cwww, c1*cvw, c2*cvw, cwdv, -(4*c2*d1 + 4*c1*d2 + 2*c2*d2 + c1*cvw + c2*cvw + cwdv)]>;
 ICpbw:=IrreducibleComponents(BaseScheme(proj2walt));
 [[Dimension(Intersection(V,W)): V in ICpbw]: W in ICpbw];
 [IrreducibleComponents(Intersection(V,ICpbw[1])): V in ICpbw];
 [<Dimension(V), ArithmeticGenus(V)>: V in IrreducibleComponents(BaseScheme(proj2walt))];
 [[Dimension(Intersection(V,W)): W in setPtX3[1..5] cat setPtX3[7..10]]:  V in ICpbw];
 [[Dimension(Intersection(V,W)): W in setPtX3[1..5] cat setPtX3[7..10]]:  V in ICpb cat ICpbw];
 sPtX3w := [Pullback(proj2walt, pt): pt in SingularPoints(Proj2X3w)];
 setPtX3w:= SetToSequence(SequenceToSet(&cat [IrreducibleComponents(ReducedSubscheme(pt)): pt in sPtX3w]));
 [ADEtype(pt): pt in SingularPoints(Proj2X3w)];
 [[<Dimension(V), ArithmeticGenus(V)>: V in IrreducibleComponents(pt)]: pt in sPtX3w];
 [[<Dimension(V), ArithmeticGenus(V)>: V in IrreducibleComponents(ReducedSubscheme(pt))]: pt in sPtX3w];

 // <c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv>

 // Fibration
 Q0t<t> :=FunctionField(Q0);
 P2X3wt<qww, q1v, q2v, tv1, tv2> := BaseExtend(Proj2X3w,Q0t);
 P1t<x1,x2> := ProjectiveSpace(Q0t,1);
 Fib1t := map<P2X3wt->P1t |[q1v,q2v]>;
 coordt :=[1,t];
 G1t := Pullback(Fib1t,P1t!coordt);
 Cut := IrreducibleComponents(G1t)[1];
 linest:= IrreducibleComponents(Scheme(P2X3wt,[qww]));
 ptO := Cut!Points(IrreducibleComponents(Intersection(linest[2],Cut))[1])[1];
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
 MordellWeilGroup(E);
 AnalyticInformation(E);




 sPtX3w := [Pullback(proj2w, pt): pt in SingularPoints(Proj2X3w)];
 setPtX3w:= SetToSequence(SequenceToSet(&cat [IrreducibleComponents(ReducedSubscheme(pt)): pt in sPtX3w]));
 spptsw := [X3w![0, 0, 0, 0, 0, 4, 0, 0, 0, 0, -1, 0, 0, 0, 0], X3w![0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0], X3w![0, 0, 0, 0, 0, 4, 4, 4, 4, 0, 4, -1, 0, 0, 0], X3w![0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 4, -1, 0, 0, 0], X3w![0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, -1, 0, 0], X3w![0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, -1, 0], X3w![0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], X3w![1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]; //One missing cause defined over GF(125)
 ICw := IrreducibleComponents(Scheme(X3w,[cw, c1v, c2v, dv1, dv2]));
 IrreducibleComponents(Scheme(X3w,[cw, c1v, c2v, dv1, dv2, c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2]));
 rpt := IrreducibleComponents(ReducedSubscheme(Scheme(Proj2X3w,[qww-q1v, q2v+3*tv1])));
 ptsP2X3w:=Points(Proj2X3w); 
 IrreducibleComponents(Pullback(proj2w,ptsP2X3w[2]));
 [Degree(pt): pt in IrreducibleComponents(ReducedSubscheme(Pullback(proj2w,ptsP2X3w[2])))];
 P13<kwww, k1wv, k2wv, kwdv1, kwdv2, k111, k112, k122, k222, kvvv, k1d1, k1d2, kvdw1, kvdw2> := ProjectiveSpace(Q0,13);
 proj3w := map<X3w->P13 |[cw^3, cw*c1v, cw*c2v, cw*dv1, cw*dv2,c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2]>;
 Proj3X3w := Scheme(P13,MinimalBasis(ideal< CoordinateRing(P13) | &cat[ DefiningEquations(Image(proj3w,X3w,d)) : d in  [1..2]] >));
 proj3w := map<X3w->Proj3X3w |[cw^3, cw*c1v, cw*c2v, cw*dv1, cw*dv2,c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2]>;
 [Degree(pt): pt in IrreducibleComponents(ReducedSubscheme(Pullback(proj2w,ptsP2X3w[2])))];
 [IsSingular(proj3w(pt)): pt in spptsw];
 ptsP3X3w:=Points(Proj3X3w); 
 proj23w := map<Proj3X3w->Proj2X3w | [kwww, k1wv, k2wv, kwdv1, kwdv2]>;
 IrreducibleComponents(BaseScheme(proj23w));
[<Dimension(V), ArithmeticGenus(V),Degree(V)>: V in IrreducibleComponents(BaseScheme(proj23w))];
[<Dimension(V), ArithmeticGenus(V),Degree(V)>: V in IrreducibleComponents(ReducedSubscheme(BaseScheme(proj23w)))];
 blowuplines := IrreducibleComponents(Pullback(proj23w, ptsP2X3w[1]));
 singptsP3w := [Proj3X3w!Points(IrreducibleComponents(Pullback(proj23w, pt))[3])[1]:  pt in SingularPoints(Proj2X3w)];
 schP3w := [IrreducibleComponents(Pullback(proj23w, pt))[3]:  pt in SingularPoints(Proj2X3w)];
 [IsSingular(pt): pt in singptsP3w];
 [Degree(Intersection(blowuplines[1],sch)): sch in schP3w];
 [Degree(Intersection(blowuplines[2],sch)): sch in schP3w];
 [ADEtype(pt): pt in SingularPoints(Proj2X3w)];
 D1 := Divisor(Proj3X3w, blowuplines[1]);
 SelfIntersection(D1); //-2
 D2 := Divisor(Proj3X3w, blowuplines[2]);
 SelfIntersection(D2); //-2
 // <kwww, k1wv, k2wv, kwdv1, kwdv2, k111, k112, k122, k222, kvvv, k1d1, k1d2, kvdw1, kvdw2>
 // <c111, c112, c122, c222, c1d1, c1d2, c1vw, c2d1, c2d2, c2vw, kvvv, kwww, kvdw, kwdv>
 //cvdw2eq := -(4*c2*d1 + c1*d2 + cvvv + cvdw1);
 // cwdv2eq := -(4*c2*d1 + 4*c1*d2 + 2*c2*d2 + c2*cvw + cwww + cwdv1);
 -1*(4*c2*d1 + c1*d2 + cvvv + cvdw),-1*(4*c2*d1 + 4*c1*d2 + 2*c2*d2 + c2*cvw + cwww + cwdv);
 P<c111, c112, c122, c222, c1d1, c1d2, c1vw, c2d1, c2d2, c2vw, kvvv, kwww, kvdw, kwdv> := Ambient(Proj3X3);
 P3iso := map<Proj3X3->Proj3X3w |[kwww, c1vw, c2vw, kwdv, -(4*c2d1 + 4*c1d2 + 2*c2d2 + c2vw + kwww + kwdv), c111, c112, c122, c222, kvvv, c1d1,c1d2,kvdw, -(4*c2d1 + c1d2 + kvvv + kvdw)]>;


 PWv<cv, c1w, c2w, dw1, dw2, c111, c112, c122, c222, cwww, c1d1, c1d2, cwdv1, cwdv2, c22dv2>  := ProjectiveSpace(Q0, [1,2,2,2,2,3,3,3,3,3,3,3,3,3,4]);
 //Probably wrong Wpiv := map<WJ->PWv | [k1, k2*k3, k3*k4, k2*k3 + 2*v2, k3*k4 + 2*v5, k2^3, k2^2*k4, k2*k4^2, k4^3, k3^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k2*k4 + 2*v4), k3*(k3^2 + 2*v3), 2*k3*v6, 2*k4^2*v6]>;
 Jv := ideal< CoordinateRing(PWv) | &cat[ DefiningEquations(Image(Wpiv,WJ,d)) : d in  [1..8]] >;
 ideal< CoordinateRing(PWv) | &cat[ DefiningEquations(Image(Wpiv,WJ,d)) : d in  [1..4]] >;
 MinimalBasis(Jv);
 X3v := Scheme(PWv,MinimalBasis(Jv));
 Wpiv := map<WJ->PWv | [k1, k2*k3, k3*k4, k2*k3 + 2*v2, k3*k4 + 2*v5, k2^3, k2^2*k4, k2*k4^2, k4^3, k3^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k2*k4 + 2*v4), k3*(k3^2 + 2*v3), 2*k3*v6, 2*k4^2*v6]>;
 HSX3v  := HilbertSeries(Ideal(DefiningPolynomials(X3v)));
 [Wpiw(pts3WJ[6]), Wpiw(pts3WJ[9])];
 P4<qvv, q1w, q2w, tw1, tw2> := ProjectiveSpace(Q0,4);
 proj2v := map<X3v->P4 |[cv^2, c1w, c2w, dw1, dw2]>;
 Proj2X3v := Scheme(P4,MinimalBasis(ideal< CoordinateRing(P4) | &cat[ DefiningEquations(Image(proj2v,X3v,d)) : d in  [1..4]] >));
 proj2v := map<X3v->Proj2X3v |[cv^2, c1w, c2w, dw1, dw2]>;
 IrreducibleComponents(BaseScheme(proj2v));
 P<c1, c2, d1, d2, cvw, cvvv, cwww, cvdw, cwdv> := Ambient(X3);
 proj2valt := map<X3->Proj2X3v |[cvvv, c1*cvw, c2*cvw, cvdw, -(4*c2*d1 + c1*d2 + 2*cvvv + cvdw)]>;
 ICpbv:=IrreducibleComponents(BaseScheme(proj2valt));
 [[Dimension(Intersection(V,W)): V in ICpbv]: W in ICpbv];
 [IrreducibleComponents(Intersection(V,ICpbv[1])): V in ICpbv];
 [<Dimension(V), ArithmeticGenus(V)>: V in IrreducibleComponents(BaseScheme(proj2valt))];
 [[Dimension(Intersection(V,W)): W in setPtX3[1..5] cat setPtX3[7..10]]:  V in ICpbv];
 [[Dimension(Intersection(V,W)): W in setPtX3[1..5] cat setPtX3[7..10]]:  V in ICpb cat ICpbw cat ICpbv];
 [[Dimension(Intersection(V,W)): W in setPtX3[1..5] cat setPtX3[7..10]]:  V in ICpb cat ICpbw[1..2] cat ICpbv[1..2]];
 [[Dimension(Intersection(V,W)): W in setPtX3[1..5] cat setPtX3[7..10]]:  V in ICpb cat ICpbw[3..4] cat ICpbv[3..4]];


 HSP12 := 1/((1-t)^4*(1-t^2)^6);

 // The map into the desingularisation
 P13<c1dw1, c2dw1, c1dw2, c1vv, c2vv, c11w, c12w, c22w, cvww, cwd1, cwd2, cvdv1, cvdv2, c2dw2> := ProjectiveSpace(Q0,13);
 pi := map<WJ->P13 | [k2*(k2*k3 + 2*v2), k4*(k2*k3 + 2*v2),  k2*(k3*k4 + 2*v5), k1^2*k2, k1^2*k4, k2^2*k3, k2*k3*k4, k3*k4^2, k1*k3^2, k3*(k2^2 + k1*k3 + 2*v1), k3*(k1*k3 + k2*k4 + 2*v4), k1*(k3^2 + 2*v3), k1*(k3^2 + 2*v6), k4*(k3*k4 + 2*v5)]>;
 ideal< CoordinateRing(P13) | &cat[ DefiningEquations(Image(pi,WJ,d)) : d in  [1..1]] >;
 P10<c1vv, c2vv, c11w, c12w, c22w, cvww, cwd1, cwd2, cvdv1, cvdv2, c2dw2> := ProjectiveSpace(Q0,10);
 pi := map<WJ->P10 | [k1^2*k2, k1^2*k4, k2^2*k3, k2*k3*k4, k3*k4^2, k1*k3^2, k3*(k2^2 + k1*k3 + 2*v1), k3*(k1*k3 + k2*k4 + 2*v4), k1*(k3^2 + 2*v3), k1*(k3^2 + 2*v6), k4*(k3*k4 + 2*v5)]>;
 J10 := ideal< CoordinateRing(P10) | &cat[ DefiningEquations(Image(pi,WJ,d)) : d in  [1..2]] >;
 X3des := Scheme(P10,MinimalBasis(J10));
 // ptsX3des := Points(X3des);
 pides := map<WJ->X3des | [k1^2*k2, k1^2*k4, k2^2*k3, k2*k3*k4, k3*k4^2, k1*k3^2, k3*(k2^2 + k1*k3 + 2*v1), k3*(k1*k3 + k2*k4 + 2*v4), k1*(k3^2 + 2*v3), k1*(k3^2 + 2*v6), k4*(k3*k4 + 2*v5)]>;
 //IrreducibleComponents(Pullback(pides, ptsX3des[1]));
 //IrreducibleComponents(Pullback(pides, ptsX3des[7]));
 [<Dimension(V), Degree(V)>: V in IrreducibleComponents(ReducedSubscheme(Pullback(pides, ptsX3des[2])))];
 [<IsSingular(pt), #IrreducibleComponents(ReducedSubscheme(Pullback(pides, pt)))>: pt in ptsX3des];
 Proj2X3w;
 //<qww, q1v, q2v, tv1, tv2>;
 wmapw := map<X3des->Proj2X3w | [cvww, c1vv, c2vv, cvdv1, cvdv2]>;
 IC23w := IrreducibleComponents(BaseScheme(wmapw));
 //IsSingular(ptsP2X3w[1]);
 //[<ArithmeticGenus(V),Degree(V),Dimension(V)>: V in IrreducibleComponents(ReducedSubscheme(Pullback(wmapv,ptsP2X3w[1])))];
 //[<V, ArithmeticGenus(V),Degree(V),Dimension(V)>: V in IrreducibleComponents(ReducedSubscheme(Pullback(wmapv,ptsP2X3w[2])))];
 //[<pt, IsSingular(pt),[<Degree(V),Dimension(V)>: V in IrreducibleComponents(ReducedSubscheme(Pullback(wmapv,pt)))]>: pt in ptsP2X3w];
 //[<pt, ADEtype(pt),[<Degree(V),Dimension(V)>: V in IrreducibleComponents(ReducedSubscheme(Pullback(wmapv,pt)))]>: pt in SingularPoints(Proj2X3w)];
// redsch:=ReducedSubscheme(Pullback(wmapv,ptsP2X3w[2]));

 wmap:= map<X3des->Proj2X3 | [c11w, c12w, c22w, cwd1, cwd2, cvww]>;
 IC23 := IrreducibleComponents(BaseScheme(wmap));
 [<Dimension(V), Degree(V)>: V in IC23];
 Dimension(Intersection(IC23[1],IC23[2]));
 Degree(Intersection(IC23[1],IC23[2]));
 [<ADEtype(pt),[<Dimension(V),Degree(V)>: V in IrreducibleComponents(Pullback(wmap,pt))]>:pt in SingularPoints(Proj2X3)];
 [<ADEtype(pt),[<Dimension(V),Degree(V)>: V in IrreducibleComponents(ReducedSubscheme(Pullback(wmap,pt)))]>:pt in SingularPoints(Proj2X3)];
 pre5A2 := [IrreducibleComponents(Pullback(wmap,pt))[3]: pt in SingularPoints(Proj2X3)[6..10]];
 [[Dimension(Intersection(V,W)): V in pre5A2]: W in pre5A2];
 conicdes := IrreducibleComponents(Pullback(wmap,SingularPoints(Proj2X3)[5]))[3];
 [Dimension(Intersection(l,conicdes)): l in pre5A2];
 [Dimension(Intersection(l,IC23[1])): l in pre5A2];
 [Dimension(Intersection(l,IC23[2])): l in pre5A2];
 [Dimension(Intersection(l,Intersection(IC23[1],IC23[2]))): l in pre5A2];
 [Dimension(Intersection(l,conicdes)): l in IC23];
 point2A1 := [IrreducibleComponents(Pullback(wmap,pt))[3]: pt in SingularPoints(Proj2X3)[[1,4]]];
 [Dimension(Intersection(l,conicdes)): l in point2A1 ];
 [Dimension(Intersection(l,IC23[1])): l in point2A1 ];
 [Dimension(Intersection(l,IC23[2])): l in point2A1 ];
 line2A1 := [IrreducibleComponents(Pullback(wmap,pt))[3]: pt in SingularPoints(Proj2X3)[[2,3]]];
[Dimension(Intersection(l,conicdes)): l in line2A1 ];
 [Dimension(Intersection(l,IC23[1])): l in line2A1 ];
 [Dimension(Intersection(l,IC23[2])): l in line2A1 ];
 [[Dimension(Intersection(l,w)): l in line2A1 ]: w in pre5A2];
 [X3des!Points(pt)[1]: pt in point2A1];
 [X3des!Points(Intersection(l,conicdes))[1]: l in line2A1];
 [X3des!Points(Intersection(l,conicdes))[1]: l in IC23];
 

 // <qvv, q1w, q2w, tw1, tw2>
 P<c1vv, c2vv, c11w, c12w, c22w, cvww, cwd1, cwd2, cvdv1, cvdv2, c2dw2> := Ambient(X3des);
 wmapv2:= map<X3des->Proj2X3v | [c2vv, c12w, c22w, -(2*c1vv + 2*c2vv + 4*cwd2 + 2*cvdv1 + cvdv2), c2dw2]>;
 IC232v := IrreducibleComponents(BaseScheme(wmapv2));
 [<Dimension(V), Degree(V)>: V in IC232v];
 [[<Dimension(Intersection(l1,l2)),Degree(Intersection(l1,l2))>: l1 in IC232v]: l2 in IC232v];
 [[Points(Intersection(l1,l2)): l1 in IC232v[1..2]]: l2 in IC232v[3..4]];
 [<ADEtype(pt),[<Dimension(V),Degree(V)>: V in IrreducibleComponents(Pullback(wmapv2,pt))]>:pt in SingularPoints(Proj2X3v)];
 [<ADEtype(pt),[<Dimension(V),Degree(V)>: V in IrreducibleComponents(ReducedSubscheme(Pullback(wmapv2,pt)))]>:pt in SingularPoints(Proj2X3v)];
 pre2A2


 
 wmapv2:= map<X3des->Proj2X3v | [c2vv, c12w, c22w, -(2*c1vv + 2*c2vv + 4*cwd2 + 2*cvdv1 + cvdv2), c2dw2]>;
 IC232v := IrreducibleComponents(BaseScheme(wmapv2));
 
 //[<pt, IsSingular(pt),[<Degree(V),Dimension(V)>: V in IrreducibleComponents(ReducedSubscheme(Pullback(wmap,pt)))]>: pt in ptsP2X3];
 
 // Fibration desingularisation
 Q0t<t> :=FunctionField(Q0);
 X3dest<c1vv, c2vv, c11w, c12w, c22w, cvww, cwd1, cwd2, cvdv1, cvdv2, c2dw2> := BaseExtend(X3des,Q0t);
 P1t<x1,x2> := ProjectiveSpace(Q0t,1);
 Fib1t := map<X3dest->P1t |[c1vv,c2vv]>;
 Fib2t := map<X3dest->P1t |[c11w,c12w]>;
 Fib3t := map<X3dest->P1t |[c12w,c22w]>;
 coordt :=[1,t];
 G1t := Intersection(Intersection(Pullback(Fib1t,P1t!coordt),Pullback(Fib2t,P1t!coordt)),Pullback(Fib3t,P1t!coordt));
 Cut := IrreducibleComponents(G1t)[1];
 linest:= IrreducibleComponents(Scheme(X3dest,[cvww]));
 [<Dimension(l),ArithmeticGenus(l)>: l in linest];
 ptO := Cut!Points(IrreducibleComponents(Intersection(linest[4],Cut))[1])[1];
 D := Divisor(ptO);
 phi3t:= DivisorMap(3*D);
 E1, mapE1 := EllipticCurve(Image(phi3t), phi3t(ptO));
 E, mapE2 := MinimalModel(E1);
 mapE:=mapE1*mapE2;
 invmapE := Inverse(mapE);
 BadPlaces(E);
 LocalInformation(E);
 G, mapG := TorsionSubgroup(E);
 P2 := mapG(G.1);
 Inverse(mapE)(P2);



// This is how I found the invariant of degree 4 that I was missing
 //pols := [k1^2*k2^2, k1^2*k2*k4, k1^2*k4^2, k2^3*k3, k2^2*k3*k4, k2*k3*k4^2, k3*k4^3, k1^3*k3, k1*k2*k3^2, k1*k3^2*k4, k3^4, k1^2*(k2^2 + k1*k3 + 2*v1), k2*k3*(k2^2 + k1*k3 + 2*v1),  k3*k4*(k2^2 + k1*k3 + 2*v1), k1^2*(k2*k4 + 2*v4), k2*k3*(k2*k4 + 2*v4), k3*k4*(k2*k4 + 2*v4), k1*k2*(k3^2 + 2*v3), k1*k4*(k3^2 + 2*v3), k3^2*(k3^2 + 2*v3), (k3^2 + 2*v3)^2,  2*k1*k2*v6, 2*k1*k4*v6, 2*k3^2*v6, 2*(k3^2 + 2*v3)*v6, 4*v6^2, k2^2*(k2*k3 + 2*v2), k2*k4*(k2*k3 + 2*v2), k4^2*(k2*k3 + 2*v2), k1*k3*(k2*k3 + 2*v2),  (k2^2 + k1*k3 + 2*v1)*(k2*k3 + 2*v2), (k2*k3 + 2*v2)*(k2*k4 + 2*v4), k2^2*(k3*k4 + 2*v5), k2*k4*(k3*k4 + 2*v5), k4^2*(k3*k4 + 2*v5), k1*k3*(k3*k4 + 2*v5),  (k2^2 + k1*k3 + 2*v1)*(k3*k4 + 2*v5), (k2*k4 + 2*v4)*(k3*k4 + 2*v5)];
 //schpols := [Scheme(WJ,pol): pol in pols];
 //schpolsim := [Wpi(sc): sc in schpols];
 //PWw<cw, c1v, c2v, dv1, dv2, c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2, a> := ProjectiveSpace(Q0, [1,2,2,2,2,3,3,3,3,3,3,3,3,3,4]);
 //colmaps:=[map<WJ->PWw | [k3, k1*k2, k1*k4, k3^2 + 2*v3, 2*v6, k2^3, k2^2*k4, k2*k4^2, k4^3, k1^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k2*k4 + 2*v4), k1*(k2*k3 + 2*v2), k1*(k3*k4 + 2*v5),  pol]>: pol in pols];
// [<#Basis(ideal< CoordinateRing(PWw) | &cat[ DefiningEquations(Image(map,WJ,d)) : d in  [1..4]] >),DefiningPolynomials(map)[15]>: map in colmaps];
// Wpi := map<WJ->PWw | [k3, k1*k2, k1*k4, k3^2 + 2*v3, 2*v6, k2^3, k2^2*k4, k2*k4^2, k4^3, k1^3, k2*(k2^2 + k1*k3 + 2*v1), k2*(k2*k4 + 2*v4), k1*(k2*k3 + 2*v2), k1*(k3*k4 + 2*v5),  k1^2*(k2^2 + k1*k3 + 2*v1)]>;
 //ideal< CoordinateRing(PWw) | &cat[ DefiningEquations(Image(Wpi,WJ,d)) : d in  [1..4]] >;
 //Jw := ideal< CoordinateRing(PWw) | &cat[ DefiningEquations(Image(Wpi,WJ,d)) : d in  [1..4]] >;
 //MinimalBasis(Jw);
 //X3w := Scheme(PWw,MinimalBasis(Jw));
 //HSX3w  := HilbertSeries(Ideal(DefiningPolynomials(HSX3w)));

//HSKum  := HilbertSeries(Ideal(DefiningPolynomials(GeneralKummerSurface(J3))));
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


// Inverse birational transformation
 PWrev<cwdv, cvdw, cwww, cvvv, cvw, d2, d1, c2, c1> := ProjectiveSpace(Q0, [3, 3, 3, 3, 2, 2, 2, 1, 1]);
 Wpirev := map<WJ->PWrev | [k3*(k3^2 + 2*v3), k1*(k2*k3 + 2*v2), k3^3, k1^3, k1*k3, k1*k3 + k2*k4 + 2*v4, k2^2 + k1*k3 + 2*v1, k4, k2]>;
 Jrev := ideal< CoordinateRing(PWrev) | &cat[ DefiningEquations(Image(Wpirev,WJ,d)) : d in  [1..6]] >;
 MinimalBasis(Jrev);
 X3rev := Scheme(PWrev,MinimalBasis(Jrev));
 simpmap := map<X3rev->PWrev | [cwdv, cvdw, cwww, cvvv, cvw, d2, d1, c2, c1]>;
 Image(simpmap,X3rev,4);
 P5<c11, c12, c22, b1, b2, kvw>:= Ambient(Proj2X3);
 P<cvw, d2, d1, c2, c1> := ProjectiveSpace(Q0, [2,2,2,1,1]);
 PWrev<cwdv, cvdw, cwww, cvvv, cvw, d2, d1, c2, c1> := Ambient(X3rev);
 phi32 := map<X3rev->P | [cvw, d2, d1, c2, c1]>;
 X2 := Image(phi32,X3rev,4);
 P<cvw, d2, d1, c2, c1> :=Ambient(X2);
 phi23 := map<X2 ->X3rev | [-(c1+c2)^2*(3*c2^2*c1^2 + 4*c2*c1^3 + 3*c1^4 + 3*cvw*c2^2 + 2*cvw*c2*c1 + 4*d2*c2*c1 + 3*d1*c2*c1 + 4*cvw*c1^2 + 3*d2*c1^2 + 2*d2^2 + 4*cvw*d1 + d2*d1 + 2*d1^2),-(c1+c2)^2*(c2^2*c1^2 + 2*c2*c1^3 + cvw*c2^2 + 2*cvw*c2*c1 + 2*d1*c2*c1 + cvw*c1^2 + 3*d2*c1^2 + 4*d2^2 + 3*d2*d1),-(c1+c2)^2*(4*c2^2*c1^2 + 2*c2*c1^3 + 3*c1^4 + 4*cvw*c2^2 + cvw*c2*c1 + 2*cvw*c1^2 + 4*cvw*d1 + 4*d2*d1 + 2*d1^2), -(c1+c2)^2*(c2*c1^3 + c1^4 + 2*cvw*c2*c1 + 2*cvw*c1^2 + 4*cvw*d1 + 4*d2*d1 + 4*d1^2) ,(c1+c2)^2*cvw, (c1+c2)^2*d2, (c1+c2)^2*d1, (c1+c2)*c2, (c1+c2)*c1 ]>;
 phi23 := map<X2 ->X3 |[c1*(c1 + c2), c2*(c1 + c2), (c1 + c2)^2*d1, (c1 + c2)^2*d2,  (c1 + c2)^2*cvw, -((c1 + c2)^2*(c1^4 + c1^3*c2 + 2*c1^2*cvw +     2*c1*c2*cvw + 4*cvw*d1 + 4*d1^2 + 4*d1*d2)),  -((c1 + c2)^2*(3*c1^4 + 2*c1^3*c2 + 4*c1^2*c2^2 + 2*c1^2*cvw +     c1*c2*cvw + 4*c2^2*cvw + 4*cvw*d1 + 2*d1^2 + 4*d1*d2)),  -((c1 + c2)^2*(2*c1^3*c2 + c1^2*c2^2 + c1^2*cvw + 2*c1*c2*cvw +     c2^2*cvw + 2*c1*c2*d1 + 3*c1^2*d2 + 3*d1*d2 + 4*d2^2)),  -((c1 + c2)^2*(3*c1^4 + 4*c1^3*c2 + 3*c1^2*c2^2 + 4*c1^2*cvw +     2*c1*c2*cvw + 3*c2^2*cvw + 3*c1*c2*d1 + 4*cvw*d1 + 2*d1^2 +     3*c1^2*d2 + 4*c1*c2*d2 + d1*d2 + 2*d2^2))] >;
 
 [<V,Dimension(V)>: V in IrreducibleComponents(ReducedSubscheme(BaseScheme(phi23)))];
 IrreducibleComponents(Scheme(Proj2X3,[c11+c12,c12+c22]));
 P32 := map<Proj3X3->Proj2X3 | [c111,c112, c122,c1d1, c1d2, c1vw]>;
 P5<c11, c12, c22, b1, b2, kvw>:= Ambient(Proj2X3);
 P13<c111, c112, c122, c222, c1d1, c1d2, c1vw, c2d1, c2d2, c2vw, kvvv, kwww, kvdw, kwdv> := ProjectiveSpace(Q0, 13);
 P23 := map<Proj2X3->Proj3X3 | [c11*(c11+c12), c12*(c11+c12), c22*(c11+c12),c22*(c12+c22),b1*(c11+c12),b2*(c11+c12),kvw*(c11+c12),b1*(c12+c22),b2*(c12+c22),kvw*(c12+c22),-(c12*c11 + c11^2 + 2*kvw*c12 + 2*kvw*c11 + 4*b2*b1 + 4*b1^2),-(c12*c11 + c11^2 + 2*kvw*c12 + 2*kvw*c11 + 4*b2^2 + 4*b1^2), -(4*c12*c11 + 4*c11^2 + 3*kvw*c12 + 2*b1*c12 + 3*kvw*c11 + 3*b2*c11 + 2*kvw*b2 + 2*kvw*b1 + b2*b1 + b1^2) ,-(kvw*c12 + 4*b2*c12 + 3*b1*c12 + kvw*c11 + 3*b2*c11 + kvw*b2 + 3*kvw*b1)]>;
 BaseScheme(P23);
 P4<qww, q1v, q2v, tv1, tv2> := Ambient(Proj2X3w);
 P2w3 := map<Proj2X3w->Proj3X3 | [c11*(c11+c12), c12*(c11+c12), c22*(c11+c12),c22*(c12+c22),b1*(c11+c12),b2*(c11+c12),kvw*(c11+c12),b1*(c12+c22),b2*(c12+c22),kvw*(c12+c22),-(c12*c11 + c11^2 + 2*kvw*c12 + 2*kvw*c11 + 4*b2*b1 + 4*b1^2),-(c12*c11 + c11^2 + 2*kvw*c12 + 2*kvw*c11 + 4*b2^2 + 4*b1^2), -(4*c12*c11 + 4*c11^2 + 3*kvw*c12 + 2*b1*c12 + 3*kvw*c11 + 3*b2*c11 + 2*kvw*b2 + 2*kvw*b1 + b2*b1 + b1^2) ,-(kvw*c12 + 4*b2*c12 + 3*b1*c12 + kvw*c11 + 3*b2*c11 + kvw*b2 + 3*kvw*b1)]>;
 P<cw, c1v, c2v, dv1, dv2, c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2, c22dw2> := Ambient(X3w);
 PWwrev<c22dw2, cvdw2, cvdw1, c1d2, c1d1, cvvv, c222, c122, c112, c111, dv2, dv1, c2v, c1v, cw> := ProjectiveSpace(Q0, [4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1]);
 P<cw, c1v, c2v, dv1, dv2, c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2, c22dw2> := Ambient(X3w);
 simpmapw := map<X3w->PWwrev | [c22dw2, cvdw2, cvdw1, c1d2, c1d1, cvvv, c222, c122, c112, c111, dv2, dv1, c2v, c1v, cw]>;
 Image(simpmapw,X3w,4);
 spP< x, dv2, dv1, c2v, c1v, cw> := ProjectiveSpace(Q0, [3, 2, 2, 2, 2, 1]);
 P<cw, c1v, c2v, dv1, dv2, c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2, c22dw2> := Ambient(X3w);
 spimpmapw := map<X3w->spP | [cvvv, dv2, dv1, c2v, c1v, cw]>;
 Image(spimpmapw,X3w,6);
 Image(spimpmapw,X3w,3);
 [<P.n,#DefiningEquations(Image(map<X3w->spP | [P.n, dv2, dv1, c2v, c1v, cw]>,X3w,4))>:n in [6..14]];
 P[1];
 P5<x, tv2, tv1, q2v, q1v, qww> := ProjectiveSpace(Q0,[4,1,1,1,1,1]);
 DefiningEquations(Proj2X3w);
 aidsch := Scheme(P5,[qww^2 + 3*qww*q1v + 4*qww*q2v + 3*q1v*tv1 + 4*tv1^2, qww*q1v^2 + 3*q1v^3 + 3*qww*q1v*q2v + q1v^2*q2v + 2*qww*q2v^2 + 3*q1v*q2v^2 + 2*q1v^2*tv1 + q1v*q2v*tv1 + 4*q2v^2*tv1 + q1v*tv1^2 + 4*q2v*tv1^2 + 3*q1v^2*tv2 + 3*q1v*q2v*tv2 + 3*q1v*tv1*tv2 + 4*q2v*tv1*tv2 + 2*q1v*tv2^2]);
 repmap := map<X3w->P5 | [c112*c2v*cw + 4*c112*c1v*cw, dv2, dv1, c2v, c1v, cw^2]>;
 [<P.n,DefiningEquations(Image(map<X3w->P5 | [c2v*(c1v*cw-c2v*cw)*P.n, dv2, dv1, c2v, c1v, cw^2]>,X3w,4))[1]>:n in [6..14]];
 [<P.n,DefiningEquations(Image(map<X3w->P5 | [c2v*(c1v*cw+c2v*cw+3*cw^3)*P.n, dv2, dv1, c2v, c1v, cw^2]>,X3w,4))[1]>:n in [6..14]];
 //P5<x, tv2, tv1, q2v, q1v, qww> := ProjectiveSpace(Q0,[5,1,1,1,1,1]);
 //eqs := [<P.n,DefiningEquations(Image(map<X3w->P5 | [c2v*(c1v*cw-c2v*cw)*(c1v+c2v+3*cw^2)*P.n, dv2, dv1, c2v, c1v, cw^2]>,X3w,5))[1]>:n in [6..14]];
 
 P13<kwww, k1wv, k2wv, kwdv1, kwdv2, k111, k112, k122, k222, kvvv, k1d1, k1d2, kvdw1, kvdw2> := ProjectiveSpace(Q0,13);
 P<cw, c1v, c2v, dv1, dv2, c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2, c22dw2> := Ambient(X3w);
 proj3w := map<X3w->P13 |[cw^3, cw*c1v, cw*c2v, cw*dv1, cw*dv2,c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2]>;
 Proj3X3w := Scheme(P13,MinimalBasis(ideal< CoordinateRing(P13) | &cat[ DefiningEquations(Image(proj3w,X3w,d)) : d in  [1..2]] >));
 proj3w := map<X3w->Proj3X3w |[cw^3, cw*c1v, cw*c2v, cw*dv1, cw*dv2, c111, c112, c122, c222, cvvv, c1d1, c1d2, cvdw1, cvdw2]>;
 [x-eqs[n,2]: n in [1..9]];
 P4<qww, q1v, q2v, tv1, tv2> := Ambient(Proj2X3w);
 eqslist := [q2v*(q1v-q2v)*(q1v+q2v+3*qww)*qww^2,
    q2v*(q1v-q2v)*(q1v+q2v+3*qww)*qww*q1v,
    q2v*(q1v-q2v)*(q1v+q2v+3*qww)*qww*q2v,
    q2v*(q1v-q2v)*(q1v+q2v+3*qww)*qww*tv1,
    q2v*(q1v-q2v)*(q1v+q2v+3*qww)*qww*tv2,
    3*tv2*tv1*q2v*q1v^2 + 3*tv1*q2v^2*q1v^2 + 3*tv2*tv1*q1v^3 + 3*tv2*q2v*q1v^3 
        + 3*tv1*q2v*q1v^3 + q2v^2*q1v^3 + 3*tv2*q1v^4 + 2*q2v*q1v^4 + q1v^5 + 
        4*tv2*tv1*q1v^2*qww + 4*tv1*q2v*q1v^2*qww + q2v^2*q1v^2*qww + 
        4*tv2*q1v^3*qww + 4*q1v^4*qww + q2v*q1v^2*qww^2 + q1v^3*qww^2 + 
        4*q1v^2*qww^3,
    3*tv2*tv1*q2v^2*q1v + 3*tv1*q2v^3*q1v + 3*tv2*tv1*q2v*q1v^2 + 
        3*tv2*q2v^2*q1v^2 + 3*tv1*q2v^2*q1v^2 + q2v^3*q1v^2 + 3*tv2*q2v*q1v^3 + 
        2*q2v^2*q1v^3 + q2v*q1v^4 + 4*tv2*tv1*q2v*q1v*qww + 4*tv1*q2v^2*q1v*qww 
        + q2v^3*q1v*qww + 4*tv2*q2v*q1v^2*qww + 4*q2v*q1v^3*qww + 
        q2v^2*q1v*qww^2 + q2v*q1v^2*qww^2 + 4*q2v*q1v*qww^3,
    3*tv2*tv1*q2v^3 + 3*tv1*q2v^4 + 3*tv2*tv1*q2v^2*q1v + 3*tv2*q2v^3*q1v + 
        3*tv1*q2v^3*q1v + q2v^4*q1v + 3*tv2*q2v^2*q1v^2 + 2*q2v^3*q1v^2 + 
        q2v^2*q1v^3 + 4*tv2*tv1*q2v^2*qww + 4*tv1*q2v^3*qww + q2v^4*qww + 
        4*tv2*q2v^2*q1v*qww + 4*q2v^2*q1v^2*qww + q2v^3*qww^2 + q2v^2*q1v*qww^2 
        + 4*q2v^2*qww^3,
    tv2^2*q2v^3 + 2*tv2*tv1*q2v^3 + 2*tv2*q2v^4 + 2*tv1*q2v^4 + 2*tv2*q2v^3*q1v 
        + 3*tv2^2*q2v^2*qww + 2*tv2*tv1*q2v^2*qww + tv2*q2v^3*qww + 
        2*tv1*q2v^3*qww + 4*q2v^4*qww + 2*tv2*q2v^2*q1v*qww + 2*q2v^2*q1v^2*qww 
        + q2v^2*q1v*qww^2 + 4*q2v^2*qww^3,
    3*tv2*tv1*q2v^2*q1v + 3*tv1*q2v^3*q1v + 2*tv2*tv1*q2v*q1v^2 + 
        3*tv2*q2v^2*q1v^2 + 2*tv1*q2v^2*q1v^2 + q2v^3*q1v^2 + 2*tv2*q2v*q1v^3 + 
        4*q2v*q1v^4 + 4*q2v^4*qww + q2v^3*q1v*qww + q2v^3*qww^2 + 
        3*q2v^2*q1v*qww^2 + q2v*q1v^2*qww^2,
    tv1*q2v^3*q1v + tv2*q2v^2*q1v^2 + 2*q2v^3*q1v^2 + 3*tv2*q2v*q1v^3 + 
        3*tv1*q2v*q1v^3 + 3*q2v^2*q1v^3 + 2*tv2*q1v^4 + 4*tv1*q1v^4 + 4*q1v^5 + 
        3*tv1*q2v^2*q1v*qww + tv2*q2v*q1v^2*qww + 4*q2v^2*q1v^2*qww + 
        4*tv2*q1v^3*qww + 4*q2v*q1v^3*qww + 4*tv2*q1v^2*qww^2 + 
        4*tv1*q1v^2*qww^2 + 4*q2v*q1v^2*qww^2 + 4*q1v^3*qww^2,
    3*tv1*q2v^3*q1v + 3*tv2*q2v^2*q1v^2 + 4*tv1*q2v^2*q1v^2 + q2v^3*q1v^2 + 
        3*tv2*q2v*q1v^3 + tv1*q2v*q1v^3 + 2*q2v^2*q1v^3 + q2v*q1v^4 + 
        3*tv2*q2v^2*q1v*qww + 2*tv1*q2v^2*q1v*qww + 3*q2v^3*q1v*qww + 
        2*tv2*q2v*q1v^2*qww + tv1*q2v*q1v^2*qww + 4*q2v^2*q1v^2*qww + 
        q2v*q1v^3*qww + 4*tv2*q2v*q1v*qww^2 + 4*tv1*q2v*q1v*qww^2 + 
        4*q2v^2*q1v*qww^2 + 4*q2v*q1v^2*qww^2,
    2*tv2*tv1*q2v^2*q1v + 3*tv2*tv1*q2v*q1v^2 + tv1*q2v^2*q1v^2 + 
        4*tv1*q2v*q1v^3 + 2*tv1*q2v^3*qww + q2v^4*qww + tv2*q2v^2*q1v*qww + 
        3*tv1*q2v^2*q1v*qww + 2*q2v^3*q1v*qww + 4*tv2*q2v*q1v^2*qww + 
        2*q2v^2*q1v^2*qww + 4*q2v^3*qww^2 + 2*q2v^2*q1v*qww^2 + 
        4*q2v*q1v^2*qww^2,
    tv1*q2v^4 + 2*tv2*tv1*q2v^2*q1v + tv2*q2v^3*q1v + tv1*q2v^3*q1v + 
        2*q2v^4*q1v + 3*tv2*tv1*q2v*q1v^2 + 4*tv2*q2v^2*q1v^2 + 
        4*tv1*q2v^2*q1v^2 + 4*tv1*q2v*q1v^3 + 3*q2v^2*q1v^3 + tv1*q2v^3*qww + 
        q2v^4*qww + 2*tv2*q2v^2*q1v*qww + 2*q2v^3*q1v*qww + 3*tv2*q2v*q1v^2*qww 
        + 4*tv1*q2v*q1v^2*qww + 3*q2v^2*q1v^2*qww + 4*q2v*q1v^3*qww + 
        4*q2v^3*qww^2 + 2*q2v^2*q1v*qww^2 + 4*q2v*q1v^2*qww^2 ];
 P23w := iso<Proj2X3w->Proj3X3w | eqslist, [kwww, k1wv, k2wv, kwdv1, kwdv2]>;
 IP23w := Inverse(P23w);
 IP23w :=map<Proj3X3w->Proj2X3w | [kwww, k1wv, k2wv, kwdv1, kwdv2]>;
 BaseScheme(P23w);
 GCD(eqslist);
 ICBSw := IrreducibleComponents(BaseScheme(P23w));
 [<Dimension(V), Degree(V)>: V in ICBSw];
 [<Points(pt)[1],IsSingular(Proj2X3w!Points(pt)[1])>: pt in ICBSw[5..10]];
 [<Points(pt)[1],ADEtype(Proj2X3w!Points(pt)[1])>: pt in ICBSw[[5,6,7,8,10]]];
 [<Points(pt)[1],[<Dimension(V), Degree(V), ArithmeticGenus(V)>: V in IrreducibleComponents(ReducedSubscheme(Pullback(IP23w, Points(pt)[1])))]>: pt in ICBSw[5..10]];

 Image(repmap,X3w,3);
 PWrev<cwdv, cvdw, cwww, cvvv, cvw, d2, d1, c2, c1> := ProjectiveSpace(Q0, [3, 3, 3, 3, 2, 2, 2, 1, 1]);
 Wpirev := map<WJ->PWrev | [k3*(k3^2 + 2*v3), k1*(k2*k3 + 2*v2), k3^3, k1^3, k1*k3, 
 k2*k4 + 2*v4, k2^2 + k1*k3 + 2*v1, k4, k2]>;
 Jrev := ideal< CoordinateRing(PWrev) | &cat[ DefiningEquations(Image(Wpirev,WJ,d)) : d in  [1..6]] >;
 MinimalBasis(Jrev);
 X3rev := Scheme(PWrev,MinimalBasis(Jrev));

 P1122<cc1, cc2, dd1, dd2> := ProjectiveSpace(Q0,[1,1,2,2]);
 pi1122 := map<WJ->P1122 | [k2, k4, k2^2 + k1*k3 + 2*v1, k1*k3 + k2*k4 + 2*v4]>;
 ideal< CoordinateRing(P1122) | &cat[ DefiningEquations(Image(pi1122,WJ,d)) : d in  [1..3]] >; // This stops working for bigger d

 P11222<c1, c2, d1, d2, cvw> := ProjectiveSpace(Q0,[1,1,2,2,2]);
 pi11222 := map<WJ->P11222 | [k2, k4, k2^2 + k1*k3 + 2*v1, k1*k3 + k2*k4 + 2*v4, k1*k3]>;
 IrreducibleComponents(BaseScheme(pi11222)); //The 4 singular points that get blown-up
 X11222 := Scheme(P11222,MinimalBasis(ideal< CoordinateRing(P11222) | &cat[ DefiningEquations(Image(pi11222,WJ,d)) : d in  [1..5]] >)); 
 ICsingX12:=  IrreducibleComponents(ReducedSubscheme(JacobianSubrankScheme(X11222)));
 [<Dimension(V),ArithmeticGenus(V)>: V in ICsingX12];
 [<pt,ADEtype(pt)>: pt in SingularPoints(Proj2X3)];
 IrreducibleComponents(Scheme(X11222,[c1,c2]));
 proj1122 := map<X11222->P1122 | [c1,c2,d1,d2]>;
 X1122 := Image(proj1122);
 IrreducibleComponents(BaseScheme(proj1122));
 SX := JacobianSubrankScheme(X1122);
 SXred := ReducedSubscheme(SX);  
 ICSX := IrreducibleComponents(SXred);
 ICP:=[IrreducibleComponents(Pullback(proj1122,V)): V in ICSX];
 [IrreducibleComponents(Intersection(V[1],V[2])): V in ICP[1..3]];
 [[<Dimension(Intersection(V,W))>: V in &cat ICP]:W in &cat ICP];
 [[<Dimension(Intersection(V,W)),IrreducibleComponents(Intersection(V,W))>: V in &cat ICP]:W in &cat ICP];

 IrreducibleComponents(ReducedSubscheme(Scheme(X1122,[cc1,cc2])));



 P12w<cw, cv1, cv2, dv1, dv2> := ProjectiveSpace(Q0,[1,2,2,2,2]);
 pi12w := map<WJ->P12w | [k3, k1*k2, k1*k4, k3^2 + 2*v3, k3^2 + 2*v6]>;
 #IrreducibleComponents(BaseScheme(pi12w));
 X12w := Scheme(P12w,MinimalBasis(ideal< CoordinateRing(P12w) | &cat[ DefiningEquations(Image(pi12w,WJ,d)) : d in  [1..6]] >));
 ICsingX12w := IrreducibleComponents(Intersection(ReducedSubscheme(JacobianSubrankScheme(X12w)),X12w));
 [<Dimension(V),ArithmeticGenus(V)>: V in ICsingX12w];
 [IrreducibleComponents(ReducedSubscheme(Pullback(pi12w,Scheme(X12w, DefiningEquations(V))))): V in ICsingX12w];
 [<ReducedSubscheme(V),Dimension(V),ArithmeticGenus(V)>: V in IrreducibleComponents(Scheme(WJ,[k3,k1*k2 + k1*k4 + 3*(k3^2 + 2*v6),k3^2 + 2*v3]))];
 [<pt,ADEtype(pt)>: pt in SingularPoints(Proj2X3w)];
 IrreducibleComponents(ReducedSubscheme(Scheme(X12w,[cw]))) eq ICsingX12w;
 [[<Dimension(Intersection(W,V)), ReducedSubscheme(Intersection(W,V))> : V in ICsingX12w[1..3]]: W in ICsingX12w[1..3]];
 Ppw := <cw, cv1, cv2, dv1, dv2> := ProjectiveSpace(Q0,[1,2,2,2,2]);

p4w := map<WJ->Proj2X3w | [k3^2, k1*k2, k1*k4, k3^2 + 2*v3, k3^2 + 2*v6]>;
[<IrreducibleComponents(ReducedSubscheme(Pullback(p4w,pt))),ADEtype(pt)>: pt in SingularPoints(Proj2X3w)];


P54<c11, c12, c22, b1, b2, kvw, qww, q1v, q2v, tv1, tv2>:=ProductProjectiveSpace(Q0,[5,4]);
picomb := map<WJ->P54 | [k2^2, k2*k4, k4^2, k2^2 + k1*k3 + 2*v1, k1*k3 + k2*k4 + 2*v4, k1*k3, k3^2, k1*k2, k1*k4, k3^2 + 2*v3, k3^2 + 2*v6]>;
// 
I0 := ideal< CoordinateRing(P54) | MinimalBasis(&cat[ DefiningEquations(Image(picomb,WJ,[0,d])) : d in  [0..3]]) >;
// I1 := ideal< CoordinateRing(P54) | MinimalBasis(&cat[ DefiningEquations(Image(picomb,WJ,[1,d])) : d in  [0..3]] cat Basis(I0)) >;
deg2 :=[[i,2-i]: i in [0..2]];
I2 := ideal< CoordinateRing(P54) | MinimalBasis(&cat[ DefiningEquations(Image(picomb,WJ,d)) : d in deg2]) >;
deg3 :=[[i,3-i]: i in [0..3]];
I3 := ideal< CoordinateRing(P54) | MinimalBasis(&cat[ DefiningEquations(Image(picomb,WJ,d)) : d in deg3]) >;
deg4 :=[[i,4-i]: i in [0..4]];
I4 := ideal< CoordinateRing(P54) | MinimalBasis(&cat[ DefiningEquations(Image(picomb,WJ,d)) : d in deg4]) >;
Bt:=MinimalBasis(Basis(I2) cat Basis(I3) cat Basis(I4));
[#[V: V in Bt |  (Degree(V) eq i)]: i in [1..4]];
[[#[V: V in Bt |  (Degrees(P54,V) eq [i,j])]: i in [0..2]]:j in [0..3]];

X1w :=Scheme(P54,Bt);
proj1 := map<X1w->Proj2X3 | [c11, c12, c22, b1, b2, kvw]>;
[<Dimension(Pullback(proj1,pt)),ADEtype(pt)>: pt in SingularPoints(Proj2X3)];
pb1 := [Pullback(proj1,pt): pt in SingularPoints(Proj2X3)];
IrreducibleComponents(BaseScheme(proj1));
projw := map<X1w->Proj2X3w | [qww, q1v, q2v, tv1, tv2]>;
[<Dimension(Pullback(projw,pt)),ADEtype(pt)>: pt in SingularPoints(Proj2X3w)];
IrreducibleComponents(BaseScheme(projw));
pbw := [Pullback(projw,pt): pt in SingularPoints(Proj2X3w)];
[[Dimension(Intersection(V,W)): V in pb1]: W in pbw];
listpts:=&cat[&cat[IrreducibleComponents(Intersection(V,W)): V in pb1]: W in pbw];
ptscoord :=[[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1], [4, 4, 4, 4, -1, 0, 0, 4, 4, 0, -1], [1, 0, 0, -1, 0, 0, 0, 3, 0, 0, -1], [4, 0, 0, -1, 0, 0, 0, 2, 0, 1, -1], [0, 0, 0, 3, 1, -1, 0, 2, 3, 1, -1], [4, 1, 4, 0, 0, -1, 0, 2, 3, 1, -1], [1, 1, 1, 4, -1, 0, 0, 2, 2, 1, -1], [0, 0, 0, 4, 1, -1, 4, 0, 0, 1, -1], [4, 1, 4, 0, 0, -1, 4, 0, 0, 1, -1], [0, 0, 0, 1, 4, -1, 1, 0, 0, 1, -1], [4, 1, 4, 0, 0, -1, 1, 0, 0, 1, -1], [0, 0, 0, 2, 4, -1, 0, 1, -1, 0, 0], [4, 1, 4, 0, 0, -1, 0, 1, -1, 0, 0]];
pts := [X1w!pt: pt in ptscoord];
[<pt,IsSingular(pt)> : pt in pts];
singpts:=[pt: pt in pts | IsSingular(pt)];
[ADEtype(pt): pt in singpts];
IrreducibleComponents(BaseScheme(projw));


P544<c11, c12, c22, b1, b2, kvw, qww, q1v, q2v, tv1, tv2, qvv, q1w, q2w, tw1, tw2> := ProductProjectiveSpace(Q0,[5,4,4]);
pitot:= map<WJ->P544 | [k2^2, k2*k4, k4^2, k2^2 + k1*k3 + 2*v1, k1*k3 + k2*k4 + 2*v4, k1*k3, k3^2, k1*k2, k1*k4, k3^2 + 2*v3, k3^2 + 2*v6, k1^2, k2*k3, k3*k4, k2*k3 + 2*v2, k3*k4 + 2*v5]>;
deg2 :=&cat[[[2-i-j,i,j]: i in [0..2-j]]: j in [0..2]];
I2 := ideal< CoordinateRing(P544) | MinimalBasis(&cat[ DefiningEquations(Image(pitot,WJ,d)) : d in deg2]) >;
deg3 :=&cat[[[3-i-j,i,j]: i in [0..3-j]]: j in [0..3]];
I3 := ideal< CoordinateRing(P544) | MinimalBasis(&cat[ DefiningEquations(Image(pitot,WJ,d)) : d in deg3]) >;
Bt := MinimalBasis(Basis(I2) cat Basis(I3));
[[[#[V: V in Bt |  (Degrees(P544,V) eq [n-i-j,i,j])]: i in [0..n-j]]:j in [0..n]]: n in [1..3]];
deg4 :=&cat[[[4-i-j,i,j]: i in [0..4-j]]: j in [0..4]];
// I4 := ideal< CoordinateRing(P544) | MinimalBasis(&cat[ DefiningEquations(Image(pitot,WJ,d)) : d in deg4]) >; //We can see that basically all the relations are in degree less than 4
// Bt4 := MinimalBasis(Basis(I2) cat Basis(I3) cat Basis(I4));
// Bt eq Bt4;
Y3 := Scheme(P544,Bt);
proj1 := map<Y3->Proj2X3 | [c11, c12, c22, b1, b2, kvw]>;
[<Dimension(Pullback(proj1,pt)),#IrreducibleComponents(Pullback(proj1,pt)), ADEtype(pt)>: pt in SingularPoints(Proj2X3)];
pb1 := [Pullback(proj1,pt): pt in SingularPoints(Proj2X3)];
projw := map<Y3->Proj2X3w | [qww, q1v, q2v, tv1, tv2]>;
[<Dimension(Pullback(projw,pt)),#IrreducibleComponents(Pullback(projw,pt)),ADEtype(pt)>: pt in SingularPoints(Proj2X3w)];
pbw := [Pullback(projw,pt): pt in SingularPoints(Proj2X3w)];
[Dimension(BaseScheme(proj1)),Dimension(BaseScheme(projw)),Dimension(BaseScheme(projv))];
projv := map<Y3->Proj2X3v | [qvv, q1w, q2w, tw1, tw2]>;
[<Dimension(Pullback(projv,pt)),#IrreducibleComponents(Pullback(projv,pt)),ADEtype(pt)>: pt in SingularPoints(Proj2X3v)];
pbv := [Pullback(projv,pt): pt in SingularPoints(Proj2X3v)];
listpts1w:=&cat[&cat[IrreducibleComponents(Intersection(V,W)): V in pb1]: W in pbw];
listptstot:=&cat[&cat[&cat[IrreducibleComponents(Intersection(Intersection(V,W),U)): V in pbv]: W in pbw]: U in pb1];
listptstot:=&cat[&cat[IrreducibleComponents(Intersection(V,W)): V in pbv]: W in listpts1w];
ptscoord := [[0, 0, 0, 4, 1, -1, 4, 0, 0, 1, -1, 0, 4, 1, 1, -1], [0, 0, 0, 3, 1, -1, 0, 2, 3, 1, -1, 3, 0, 0, -1, 0], [0, 0, 0, 2, 4, -1, 0, 1, -1, 0, 0, 3, 0, 0, 0, -1], [0, 0, 0, 1, 4, -1, 1, 0, 0, 1, -1, 0, 1, 4, 1, -1], [4, 1, 4, 0, 0, -1, 0, 2, 3, 1, -1, 3, 0, 0, -1, 0], [4, 1, 4, 0, 0, -1, 4, 0, 0, 1, -1, 0, 4, 1, 1, -1], [4, 1, 4, 0, 0, -1, 1, 0, 0, 1, -1, 0, 1, 4, 1, -1], [4, 1, 4, 0, 0, -1, 0, 1, -1, 0, 0, 3, 0, 0, 0, -1], [4, 4, 4, 4, -1, 0, 0, 4, 4, 0, -1, 0, 4, 4, 4, -1], [1, 1, 1, 4, -1, 0, 0, 2, 2, 1, -1, 0, 1, 1, 4, -1], [4, 0, 0, -1, 0, 0, 0, 2, 0, 1, -1, 0, 4, 0, -1, 0], [1, 0, 0, -1, 0, 0, 0, 3, 0, 0, -1, 0, 1, 0, -1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]];
pts := [Y3!pt: pt in ptscoord];
[<pt,IsSingular(pt)> : pt in pts];
[IrreducibleComponents(proj1(W)): W in pbw];