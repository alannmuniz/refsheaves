--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
---- Setup for the examples in "Rank Two Reflexive Sheaves on P^3 with c_2 =4" 
---- by Marcos Jardim and Alan Muniz 
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- SETUP
-------------------------
restart;
loadPackage "SpaceCurves"
--kk = ZZ/1009;
kk = QQ;
R = kk[x_0..x_3]
O = OO_(Proj R);
-- computes the Rao function of v(I)
rao = I ->( 
    L :=  toList(-3..5);
    M := raoModule I;
    return matrix{L,apply(L, i-> hilbertFunction(i,M))};
    );
------ computes the cohomology table of F
cohTab = F ->( L:={}; 
for i from -3 to 5 do L = append(L, {i, rank  HH^0(F(i)), rank  HH^1(F(i)), rank HH^2(F(i)), rank HH^3(F(i)) });
return transpose matrix L);
------ computes the HH^1 dimensions of F
dimHH1 = F ->( L:={};  
for i from -3 to 3 do L = append(L, {i, rank  HH^1(F(i)) });
return transpose matrix L);
--- dualizing sheaf of C
omg = C-> prune sheaf Ext^1(C,R^{-4}); 
--- random complete intersenction containing Z
rand = (Z,a,b) -> (
    zz := gens Z;
    ma := matrix basis(a, Z); 
    mb := matrix basis(b, Z); 
    sa :=rank source ma;
    sb :=rank source mb;
    return trim ideal(zz*ma*random(kk^(sa), kk^1),zz*mb*random(kk^(sb), kk^1) );
   )
------------ make a random complete intersection (a,b)
cint = (a,b) -> ideal random(R^{a,b},R^1);
--- make disjoint lines
dsLns = k->  intersect apply( entries random(R^k,R^{2:-1}), ideal);  
--- make disjoint points
pts = k ->   intersect apply( entries random(R^k,R^{3:-1}), ideal);   
--- values of the Hilbert polynomials
chi0 = (c3,l) -> l^3//3 + 2*l^2 - l//3 + c3//2 - 6 ;
chi1 = (c3,l) -> l^3//3 + 3*l^2//2 -11*l//6 + c3//2 -5 ;
chiRef = (c1,c2,c3,l) ->2 + binomial(c1+ 2*l+3,3) - 2*(c2+ l^2 + c1*l)+(c3-(c1+l*2)*(c2+ l^2 + c1*l))//2 -1;
-----------  Chern classes of a given rank-2 sheaf
chern = F ->(
g = hilbertPolynomial(F,Projective=>false);
i := (vars ring g)_0_0;
X := coefficient(i^2,g);
Y := coefficient(i,g);
Z := coefficient(i^0,g);
return (2*X-4,  2*X^2-4*X-Y+11/3,  4/3*X^3-2*X*Y+2*Z);
);
--- Serre Correspondence Ext^1(C,R^{-k})_0 = H^0(omg_C(4-k))
Serre = (C,k) -> ( -- compute an extension of C by O(-k) and normalizes it
    if rank source basis(0,Ext^1(C,R^{-k})) == 0 then error "No extensions in this degree!";
    resC := res C; -- resolution of R/C
    I := ideal vars R;
    E := cokernel(resC.dd_3);
    H := Hom(E, I*(R^{-k}));
    r := rank source basis(0, H);
    A0 := mingens(H)*matrix basis(0,H)*random(kk^r, kk^1);
    A1 := transpose matrix adjoint'(A0,cover E, R^{-k}); --extension part
    A2 := A1 % transpose resC.dd_2;
    A :=  (resC.dd_2)||transpose A2;
    F0 := coker A;
    if codim annihilator Ext^1(F0,R^1) < 3 then error "F is not reflexive";
    F:= prune sheaf F0;
    return F(k//2);
    );
------------------------------------------------
----END OF SETUP
----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

