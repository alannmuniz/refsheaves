--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
---- Examples for "Rank Two Reflexive Sheaves on P^3 with c_2 =4" 
---- by Marcos Jardim and Alan Muniz 
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--
--EVEN CASES 
--
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- c_1 = 0, c_3 = 2, spec = {-1,0,0,0}, h^0(F(1)) = 0
-- deg(C) = 8 and  g(C) = 2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2";
Z = pts 1;
F = flatten entries ((gens Z)*(matrix basis(2,Z))*random(R^(10 -degree Z),R^4));
C = intersect(ideal(F_0,F_1),ideal(F_2,F_3));
F = Serre(C,4);
chern F
Ext^2(F,F) -- unobstucted


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 2, spec = {-1,0,0,0}, h^0(F(1)) > 0
-- deg(C) = 5 and  g(C) = -3
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
C = intersect(dsLns 3, ideal(x_0, random(2,R)));
degree C, genus C
F = Serre(C,2);
chern F
Ext^2(F,F) -- unobstructed


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 2, spec = {-1,-1,0,1}, h^0(F(1)) = 0
-- deg(C) = 8 and  g(C) = 2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--- Since h^1(I_C) = 1 either C is a disjoint union of two curves or it is nonreduced. 
--- But F corresponds to a section of H^0(omg_C) and C does not lie on a cubic, this 
--- eliminates every possibility for C reduced. So C is nonreduced, if it exists. 

--- double elliptic quartic gives the other spectrum
--- maybe some quadruple structure 

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 2, spec = {-1,-1,0,1}, h^0(F(1)) > 0
-- deg(C) = 5 and  g(C) = -3
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
q = x_0*x_3 + x_1^2+ x_2^2 + x_2*x_3;-- smooth quadric through (0:0:0:1)
D = topComponents((ideal(x_0,q))^2 + ideal(x_0*x_3^3-q*x_2^2)); --double conic
C = intersect(D, ideal(x_1,x_2));
degree C, genus C
F = Serre(C,2);
chern F
time Ext^2(F,F) --unobstructed


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 4, spec = {-1,-1,0,0}, h^0(F(1)) = 0
-- deg(C) = 8 and  g(C) = 3
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2";
Z = pts 2;
F = flatten entries ((gens Z)*(matrix basis(2,Z))*random(R^(10 -degree Z),R^4));
C = intersect(ideal(F_0,F_1),ideal(F_2,F_3));
F = Serre(C,4);
chern F
Ext^2(F,F) -- unobstucted

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 4, spec = {-1,-1,0,0}, h^0(F(1)) > 0
-- deg(C) = 5 and  g(C) = -2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
C' = ideal(x_0^2, x_0*x_3,x_3^2, x_0*x_2-x_1*x_3);
C'' = minors(2, matrix{{x_0,x_1,x_2},{x_1,x_2,x_3}});
C = intersect(C',C'');
degree C, genus C
F = Serre(C,2); 
chern F
Ext^2(F,F) --unobstructed

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 4, spec = {-2,-1,0,1}, h^0(F(1)) = 0
-- deg(C) = 8 and  g(C) = 3
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

--- 1 is in the spectrum, so h^1(I_C) = h^1(F(-2))  = 1. So either C is nonreduced or 
--- a disjoint union. The later cannot afford a reflexive extension. So C must be nonreduced.

------- double elliptic quartic gives the other spectrum.


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 4, spec = {-2,-1,0,1}, h^0(F(1)) > 0
-- deg(C) = 5 and  g(C) = -2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
q = x_0*x_3 + x_1^2+ x_2^2 + x_2*x_3;
D = ideal(x_0^2, x_0*q, q^2, x_0*random(3,R)-q*random(2,R));
C = intersect(D, ideal(x_0,x_1));
degree C, genus C
F = Serre(C,2);
Ext^2(F,F) --obstructed

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 6, spec = {-1,-1,-1,0}, h^0(F(1)) = 0
-- deg(C) = 8 and  g(C) = 4
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2";
Z = pts 3;
F = flatten entries ((gens Z)*(matrix basis(2,Z))*random(R^(10 -degree Z),R^4));
C = intersect(ideal(F_0,F_1),ideal(F_2,F_3));
F = Serre(C,4);
Ext^2(F,F) -- unobstucted

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 6, spec = {-1,-1,-1,0}, h^0(F(1)) > 0
-- deg(C) = 5 and  g(C) = -1
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
C' = minors(2, matrix{{x_0,x_1,x_2},{x_1,x_2,x_3}});
C'' = ideal(x_0, random(2,R));
C = intersect(C',C'');
degree C, genus C
F = Serre(C,2);
Ext^2(F,F) -- unobstructed


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 6, spec = {-2,-1,0,0}, h^0(F(1)) = 0
-- deg(C) = 8 and  g(C) = 4
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

--????????????



--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 6, spec = {-2,-1,0,0}, h^0(F(1)) > 0
-- deg(C) = 5 and  g(C) = -1
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
q = x_0*x_3 + x_1^2+ x_2^2 + x_2*x_3;
p1 = x_1*random(1,R)+x_2*random(1,R)+x_3*random(1,R);
p2 = x_1*random(2,R)+x_2*random(2,R)+x_3*random(2,R);
D = topComponents ideal(x_1^2, x_1*q, q^2, x_1*p2-q*p1);
degree D, genus D
C = intersect(D, ideal(x_0,x_1)); degree C, genus C
F = Serre(C,2);
chern F
Ext^2(F,F) --unobstructed

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 8, spec = {-1,-1,-1,-1}, h^0(F(1)) = 0
-- deg(C) = 8 and  g(C) = 5
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
Z = pts 4;
C = intersect(rand(Z,2,2), rand(Z,2,2));
degree C, genus C
minimalBetti  C
F = Serre(C,4);
chern F
Ext^2(F,F)


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 8, spec = {-1,-1,-1,-1}, h^0(F(1)) =1
-- deg(C) = 5 and  g(C) = 0
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
Y  = intersect( ideal curve(3,0,R), ideal(x_0,x_1));
C = quotient(rand(Y,3,3),Y);
degree C, genus C, isSmooth C
rao C
F = Serre(C,2);
Ext^2(F,F) --unobstructed



    
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 8, spec = {-1,-1,-1,-1}, h^0(F(1)) = 2
-- deg(C) = 5 and  g(C) = 0
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
Y = intersect(ideal(x_0,x_1), ideal(x_2,x_3), ideal(x_0+x_2,x_1-x_3))
C = quotient(rand(Y,2,4),Y);
degree C, genus C
rao C
F = Serre(C,2);
chern F
Ext^2(F,F)


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 8, spec = {-2,-1,-1,0}, h^0(F(1)) = 1
-- deg(C) = 5 and  g(C) = 0
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
Y = intersect(ideal random(R^{2,1}, R^1),ideal random(R^{2,1}, R^1));
C = quotient(rand(Y,3,3),Y);
degree C, genus C, isSmooth C
rao C
F = Serre(C,2);
Ext^2(F,F) --unobstructed
----
p = ideal(x_0,x_1,x_3);
C = intersect(rand(p,1,3), ideal(x_0,x_1), ideal(x_2,x_3));
degree C, genus C
rao C
F = Serre(C,2);
Ext^2(F,F)

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 8, spec = {-2,-1,-1,0}, h^0(F(1)) = 2
-- deg(C) = 5 and  g(C) = 0
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
C = intersect(ideal(x_0, random(2,R)), ideal(x_1,random(3,R)));
degree C, genus C
rao C
F = Serre(C,2);
Ext^2(F,F)

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 10, spec = {-2,-1,-1,-1}, h^0(F(1)) = 1
-- deg(C) = 5 and  g(C) = 1
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
Y = kernel map(kk[s,t],R, {s^4,s^3*t,s*t^3,t^4});
C = quotient(rand(Y,3,3),Y);
degree C, genus C, isSmooth C
minimalBetti C
F = Serre(C,2);
Ext^2(F,F)


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 10, spec = {-2,-1,-1,-1}, h^0(F(1)) = 2
-- deg(C) = 5 and  g(C) = 1
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
p = ideal(x_0,x_1,x_2);
C = intersect(rand(p,1,3), rand(p,1,2));
degree C, genus C
rao C
F = Serre(C,2);
Ext^2(F,F)

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 12, spec = {-2,-2,-1,-1}, h^0(F(1)) = 2
-- deg(C) = 5 and  g(C) = 2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
Y = ideal(x_0,x_1);
C = quotient(rand(Y,3,2),Y);
degree C, genus C
minimalBetti C
F = Serre(C,2);
chern F
Ext^2(F,F)

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = 0, c_3 = 12, spec = {−3, −2, −1, 0}, h^0(F(1)) = 3
-- deg(C) = 5 and  g(C) = 2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
Y = ideal(x_3, random(4,R));
C = intersect(Y, ideal(x_0,x_1));
degree C, genus C 
rao C
F = Serre(C,2);
Ext^2(F,F)



--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--ODD CASES 
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 2, spec = {-1,-1,-1,0}, h^0(F(1)) = 0
-- deg(C) = 6 and g(C) = -1
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
C1 = minors(2, random(R^3,R^{-1,-1}));
C2 = minors(2, random(R^3,R^{-1,-1}));
C = intersect(C1,C2);
degree C, genus C
F = Serre(C,3);
Ext^2(F,F)
chern F
---
C1' = kernel map(kk[s,t], R, {s^4,s^3*t, s*t^3,t^4});
C2' = ideal random(R^1, R^{-1,-2});
C' = intersect(C1',C2');
degree C', genus C'
F' = Serre(C',3);
Ext^2(F',F') 



--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 2, spec = {-1,-1,-1,0}, h^0(F(1)) > 0
-- deg(C) = 4 and g(C) = -4
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

--????????

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 2, spec = {-2,-1,0,0}, h^0(F(1)) = 0
-- deg(C) = 6 and g(C) = -1
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

--???????

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 2, spec = {-2,-1,0,0}, h^0(F(1)) > 0
-- deg(C) = 4 and g(C) = -4
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
f = x_3*x_2^2 - x_1*(x_1-x_3)*(x_1-7*x_3);
isSmooth ideal(x_0,f)
a = random(3,R)*x_0+random(3,R)*x_1+random(3,R)*x_2;
b = random(1,R)*x_0+random(1,R)*x_1+random(1,R)*x_2;
C = topComponents ideal( x_0^2, x_0*f, f^2, x_0*a-f*b);
degree C, genus C
rao C
F = Serre(C,3);
HH^0(F(1))
---
q = x_0^2-x_1*x_2;
a = x_3^3*x_0+x_0^3*x_1+x_2^3*x_3;
b = (x_2^2+x_1^2)*x_0+x_1^2*x_1+x_0^2*x_3;
C = topComponents ideal( x_3^2, x_3*q, q^2, x_3*a-q*b);
degree C, genus C
F = Serre(C,1);
time Ext^2(F,F) --unobstructed
chern F

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 4, spec = {-1,-1,-1,-1}, h^0(F(1)) = 0
-- deg(C) = 6 and g(C) = 0
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
S = kk[y_0..y_6];
Y = kernel map(kk[s,t], S, {t^6, t^5*s, t^4*s^2, t^3*s^3, t^2*s^4, t*s^5, s^6});
phi = map(S,R, random(S^{1}, S^4));
C = preimage(phi,Y);
degree C, genus C
F = Serre(C,3);
Ext^2(F,F)



--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 4, spec = {-1,-1,-1,-1}, h^0(F(1)) > 0
-- deg(C) = 4 and g(C) = -3
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
C = intersect(ideal(x_0,x_2), ideal(x_1,x_3), ideal(x_0-x_1,x_2-x_3), ideal(x_0+x_1,x_2+x_3));
degree C, genus C
minimalBetti C
F = Serre(C,1);
time Ext^2(F,F)


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 4, spec = {-2,-1,-1,0}, h^0(F(1)) = 0
-- deg(C) = 6 and g(C) = 0
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
C1 = minors(2,matrix{{x_0,x_1,x_2},{x_1,x_2,x_3}});
C2 = ideal(x_0, random(3,R));
C = intersect(C1,C2);
degree C, genus C
F = Serre(C,3);
chern F
time Ext^2(F,F)

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 4, spec = {-2,-1,-1,0}, h^0(F(1)) > 0
-- deg(C) = 4 and g(C) = -3
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

--?????

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 6, spec = {-2,-1,-1,-1}, h^0(F(1)) = 0
-- deg(C) = 6 and g(C) = 1
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
Y = dsLns 3; 
C = quotient(rand(Y,3,3),Y);
degree C, genus C
isSmooth C
F = Serre(C,3);
Ext^2(F,F)

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 6, spec = {-2,-1,-1,-1}, h^0(F(1)) > 0
-- deg(C) = 4 and g(C) = -2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

--?????

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 6, spec = {-2,-2,-1,0}, h^0(F(1)) = 0
-- deg(C) = 6 and g(C) = 1
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

--??????


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 6, spec = {-2,-2,-1,0}, h^0(F(1)) > 0
-- deg(C) = 4 and g(C) = -2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
C = intersect(ideal(x_0, random(3,R)), ideal(x_1, random(3,R)));
degree C, genus C
F = Serre(C,3);
Ext^2(F,F)
HH^0(F(1))
---
C1 = ideal(x_0^2,x_0*x_1,x_1^2,x_0*x_3^4-x_1*x_2^4);
C2 = ideal(x_0, x_2);
C3 = ideal(x_1, x_3);
C = intersect(C1,C2,C3);
degree C, genus C
F = Serre(C,1);
chern F
Ext^2(F,F)
rao C

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 8, spec = {-2,-2,-1,-1}, h^0(F(1)) = 0
-- deg(C) = 6 and g(C) = 2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
--h^0 IC for a linked curve to a 3,-1
X = sheaf module ideal random(R^{3,3},R^1);
f = n -> max(rank HH^0(X(n)) -( 3*(2-n) + 2),0)
g = n -> (
    j := 1;
    if n < 0 then j = 0;
    return f(n) -3*f(n-1)+3*f(n-2)-f(n-3) - j;
    );
rc = n -> g(n) - g(n-1);
for i from -3 to 8 do print (i, g i, rc i)
--------------
-- linked to a 3,-1 curve 
Y = intersect(ideal random(R^{2,1}, R^1), ideal random(R^{1,1},R^1));
C = quotient(rand(Y,3,3),Y);
degree C, genus C
isSmooth C
F = Serre(C,3);
Ext^2(F,F)
---------- twc + plane cubic 
Z = pts 2;
Y = intersect(ideal(x_0,x_1), Z);
C' = quotient(rand(Y,2,2),Y); --twc through Z
C = intersect(C', rand(Z,1,3));
degree C, genus C
F = Serre(C,3);
Ext^2(F,F)


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 8, spec = {-2,-2,-1,-1}, h^0(F(1)) > 0
-- deg(C) = 4 and g(C) = -1
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
C =intersect(ideal(x_0, random(2,R)), ideal(x_1, random(2,R)));
degree C, genus C
F = Serre(C,1);
Ext^2(F,F) --unobstructed
----
-- unobstructed example 
p = random(2,R);
q = random(2,R);
r = random(1,R);
s = random(2,R);
t = x_1*random(1,R);
A = map(R^{-4,-3,-3,-2,-1}, R^{-5,-5,-4,-4},
     matrix{{x_1 , -x_0, 0, 0},
    	   {0 , p, -x_1, 0 },
	   {p-q, 0 , -x_0, x_0},
	   {r*p, 0 , -r*x_0-q, p},
	   {s*p, -t*q , -s*x_0 , t*x_1}});
F = sheaf coker A;
chern F
codim ann sheafExt^1(F,O)
Ext^2(F,F)

-- another unobstructed example 
f = random(2,R);
g = random(2,R);
i = random(1,R);
j = random(1,R);
i' = random(2,R);
j' = 0
k' = x_0+x_1
A = map(R^{-4,-3,-3,-2,-1}, R^{-5,-5,-4,-4},
     matrix{{x_1 , -x_0, 0, 0},
    	   {f , 0, -x_0, 0 },
	   {-g, f , -x_1, x_0},
	   {0, i*f-j*g , -i*x_1-g, f+j*x_1},
	   {0, i'*f-j'*g , -i'*x_1-k'*g, k'*f+j'*x_1}
	   });
A*transpose matrix{{x_0,x_1,f,g}}
A%ideal(x_0,x_1,f)
F = sheaf coker A;
chern F
codim ann sheafExt^1(F,O)
Ext^2(F,F)
OZ = sheaf(R^1/ideal(x_0,x_1,f));
Ext^1(F,OZ)

----------
C = intersect(minors(2,matrix{{x_0,x_1,x_2},{x_1,x_2,x_3}}), ideal(x_0,x_3));
degree C, genus C
F = Serre(C,1);
Ext^2(F,F)
rao C

--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 8, spec = {-3,-2,-1,0}, h^0(F(1)) = 0
-- deg(C) = 6 and g(C) = 2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
p = ideal(x_0,x_1,x_2);
C' = intersect(ideal(x_0,x_1), ideal(x_2,x_3));
C = intersect(C',rand(p,1,4));
degree C, genus C
F = Serre(C,3) -- not reflexive
----
rC = res C;
m = random(rC_2,rC_3);
E0 = coker m;
H = Hom(E0, rC_1++R^{-3});
r = rank source basis(0, H)
A0 = mingens(H)*matrix basis(0,H)*random(kk^r, kk^1);
A1 = matrix adjoint'(A0,rC_2,rC_1++R^{-3});
F = sheaf prune coker A1;
codim ann sheafExt^1(F,O) --not reflexive




--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 8, spec = {-3,-2,-1,0}, h^0(F(1)) > 0
-- deg(C) = 4 and g(C) = -1
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
C = ideal(x_0^2,x_0*x_1,x_1^4, x_0*random(4,R)-x_1^3*x_2^2 );
degree C, genus C
rao C
F = Serre(C,1);
Ext^2(F,F)


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 10, spec = {-2,-2,-2,-1}, h^0(F(1)) = 0
-- deg(C) = 6 and g(C) = 2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
Y = minors(2, matrix{{x_0,x_1,x_2},{x_1,x_2,x_3}});
C = quotient(rand(Y,3,3),Y);
degree C, genus C
minimalBetti C
------
F = sheaf prune  cokernel random(R^{5:-2},R^{3:-3});
codim annihilator sheafExt^1(F,O)
cohTab F
Ext^2(F,F)


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 10, spec = {-2,-2,-2,-1}, h^0(F(1)) > 0
-- deg(C) = 4 and g(C) = 0
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"

C = ker map(kk[s,t],R,matrix{{t^4,s*t^3,s^3*t, s^4}}) --rational quartic
rao C
minimalBetti C
F = Serre(C,1);
chern F
cohTab F
Ext^2(F,F) --unobstructed


----
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 10, spec = {-3,-2,-1,-1}, h^0(F(1)) = 0
-- deg(C) = 6 and g(C) = 2
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
Y = intersect(ideal(x_0,x_1),ideal(x_2,x_3));
P = Y + ideal random(1,R);
C = intersect(Y,rand(P,4,1));
degree C, genus C
rao C
minimalBetti C
F = Serre(C,3);
HH^0(F(1))
cohTab F
Ext^2(F,F) --unobstructed


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 10, spec = {-3,-2,-1,-1}, h^0(F(1)) >  0
-- deg(C) = 4 and g(C) = 0
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
C = intersect( ideal(x_0,x_1),ideal(x_2, random(3,R)))
degree C, genus C
minimalBetti C
rao C
F = Serre(C,1);
cohTab F
Ext^2(F,F) --unobstructed

------ obstructed case
m = matrix{{x_0,x_1,x_2,x_3^3}}
l1 = random(1,R);
l2 = random(1,R);
f = random(3,R);
A = map(R^{2:-4,2:-2,-1}, R^{3:-5,-3},
     matrix{{-x_2,0,x_0,0}, 
	    {-x_1,x_0,0,0},
	    {0,-x_3^3,0,x_1},
	    {-x_3^3, -x_1^2*x_2,x_1^3,x_0},
	    {0, f*x_2, -f*x_1 -x_2*x_3^3, x_2^2}})
A*transpose m
F = sheaf coker A;
chern F 
cohTab F
codim ann sheafExt^1(F,O) >= 3 -- check it is reflexive
codim ann sheafExt^2(F,O) >= 4
codim ann sheafExt^3(F,O) >= 4
Ext^2(F,F)


-----special 4-line
C = ideal(x_0^2,x_0*x_1,x_1^4, x_1^3*x_2+x_0*x_3^3)
degree C, genus C
F = Serre(C,1);
F' = Serre2(C,1)
F
tabCoh F'
Ext^2(F,F)

----------
d = 4;
f = random(d-2,R);
g = random(d-2,R);
f'= random(d-1,R);
A = map(R^{-1-d,-d,2:-2,-1}, R^{2:-2-d,-1-d,-3},
matrix{{x_1 , -x_0 ,0 ,0},
    {x_2*x_3, 0 , -x_0, 0},
    {x_2^d+x_3^d, f*x_2*x_3, -f*x_1, -x_0},
    {0,x_2^d+x_3^d+g*x_2*x_3, -g*x_1,-x_1},
    {0, f'*x_2*x_3,(x_2^d+x_3^d) -f'*x_1, -x_2*x_3}}); 
F = sheaf coker A;
rank F, codim ann sheafExt^1(F,O), degree ann sheafExt^1(F,O)
Ext^2(F,F) 



--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Ex. c_1 = -1, c_3 = 12, spec = {-3,-2,-2,-1}, h^0(F(1)) > 0
-- deg(C) = 4 and g(C) = 1
--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart;
load "refshsetup.m2"
--- generic
F = sheaf cokernel random(R^{-3,-2,-2,-1}, R^{-4,-3});
codim annihilator sheafExt^1(F,O)
chern F
cohTab F
Ext^2(F,F)

---non generic 
I = ideal vars R
H = Hom(R^{-3,-4},I*R^{-3,-2,-2,-1});
r = rank source basis(0, H)
A0 = mingens(H)*matrix basis(0,H)*random(kk^r, kk^1);
A1 = matrix adjoint'(A0,R^{-3,-4},R^{-3,-2,-2,-1});
F = sheaf prune coker A1;
codim annihilator sheafExt^1(F,O)
chern F
Ext^2(F,F)

---- plane cubic + line
C = ideal(x_0*x_1,x_0*x_2,x_1*random(2,R)-x_2*random(2,R)) ;
degree C, genus C
minimalBetti C
F = Serre(C,1);
Ext^2(F,F)
---- conic+ double line
C = intersect(ideal(x_0^2,x_1), ideal(x_0,random(2,R)));
degree C, genus C
minimalBetti C
F = Serre(C,1);
Ext^2(F,F)
chern F