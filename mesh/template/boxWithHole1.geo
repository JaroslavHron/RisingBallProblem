// This is file 'boxWithHole1.geo' is template .geo file, used for automatic mesh generation
// used for automatic mesh generation
// You have to define following variables
// L1 = 40.0;
// L2 = 10.0;
// s  = 0.5;
// a  = 1.0;
// alpha = 1.0;
// d1 = 1.0;
// d2 = 2.0;
// theta = pi/2.0

lcar1 = 100;
lcar3 = 100;
/////////////////////////////////
////////////// BOX //////////////
cp1 = newp; Point(cp1) = {+L2,+L2,+L2,lcar1};
cp2 = newp; Point(cp2) = {-L2,+L2,+L2,lcar1};
cp3 = newp; Point(cp3) = {-L2,-L2,+L2,lcar1};
cp4 = newp; Point(cp4) = {+L2,-L2,+L2,lcar1};
cp5 = newp; Point(cp5) = {+L2,+L2,-L1,lcar1};
cp6 = newp; Point(cp6) = {-L2,+L2,-L1,lcar1};
cp7 = newp; Point(cp7) = {-L2,-L2,-L1,lcar1};
cp8 = newp; Point(cp8) = {+L2,-L2,-L1,lcar1};

cl1 = newreg; Line(cl1) = {cp1,cp2};
cl2 = newreg; Line(cl2) = {cp2,cp3};
cl3 = newreg; Line(cl3) = {cp3,cp4};
cl4 = newreg; Line(cl4) = {cp4,cp1};
cl5 = newreg; Line(cl5) = {cp5,cp1};
cl6 = newreg; Line(cl6) = {cp6,cp2};
cl7 = newreg; Line(cl7) = {cp7,cp3};
cl8 = newreg; Line(cl8) = {cp8,cp4};
cl9 = newreg; Line(cl9) = {cp5,cp6};
cl10= newreg; Line(cl10)= {cp6,cp7};
cl11= newreg; Line(cl11)= {cp7,cp8};
cl12= newreg; Line(cl12)= {cp8,cp5};

cll1 = newreg; Line Loop(cll1) = {+cl1,+cl2,+cl3,+cl4};
cll2 = newreg; Line Loop(cll2) = {+cl5,-cl4,-cl8,+cl12};
cll3 = newreg; Line Loop(cll3) = {+cl6,-cl1,-cl5,+cl9};
cll4 = newreg; Line Loop(cll4) = {+cl7,-cl2,-cl6,+cl10};
cll5 = newreg; Line Loop(cll5) = {+cl8,-cl3,-cl7,+cl11};
cll6 = newreg; Line Loop(cll6) = {-cl9,-cl10,-cl11,-cl12};

cps1 = newreg; Plane Surface(cps1) = {cll1};
cps2 = newreg; Plane Surface(cps2) = {cll2};
cps3 = newreg; Plane Surface(cps3) = {cll3};
cps4 = newreg; Plane Surface(cps4) = {cll4};
cps5 = newreg; Plane Surface(cps5) = {cll5};
cps6 = newreg; Plane Surface(cps6) = {cll6};

csl = newreg; Surface Loop(csl) = {cps2,cps3,cps4,cps5,cps6,cps1};

////////////////////////////////////
////////////// SPHERE //////////////

x = 0 ; y = 0.0 ; z = 0.0 ; r = 1.0 ;

sp1 = newp; Point(sp1) = {x,  y,  z} ;
sp2 = newp; Point(sp2) = {x+r,y,  z} ;
sp3 = newp; Point(sp3) = {x,  y+r,z} ;
sp4 = newp; Point(sp4) = {x,  y,  z+r} ;
sp5 = newp; Point(sp5) = {x-r,y,  z} ;
sp6 = newp; Point(sp6) = {x,  y-r,z} ;
sp7 = newp; Point(sp7) = {x,  y,  z-r} ;

sc1 = newreg; Circle(sc1) = {sp2,sp1,sp7}; sc2 = newreg; Circle(sc2) = {sp7,sp1,sp5};
sc3 = newreg; Circle(sc3) = {sp5,sp1,sp4}; sc4 = newreg; Circle(sc4) = {sp4,sp1,sp2};
sc5 = newreg; Circle(sc5) = {sp2,sp1,sp3}; sc6 = newreg; Circle(sc6) = {sp3,sp1,sp5};
sc7 = newreg; Circle(sc7) = {sp5,sp1,sp6}; sc8 = newreg; Circle(sc8) = {sp6,sp1,sp2};
sc9 = newreg; Circle(sc9) = {sp7,sp1,sp3}; sc10 = newreg; Circle(sc10) = {sp3,sp1,sp4};
sc11 = newreg; Circle(sc11) = {sp4,sp1,sp6}; sc12 = newreg; Circle(sc12) = {sp6,sp1,sp7};

sll1 = newreg; Line Loop(sll1) = {sc5,sc10,sc4};   
sll2 = newreg; Line Loop(sll2) = {sc9,-sc5,sc1};   
sll3 = newreg; Line Loop(sll3) = {sc12,-sc8,-sc1}; 
sll4 = newreg; Line Loop(sll4) = {sc8,-sc4,sc11};  
sll5 = newreg; Line Loop(sll5) = {-sc10,sc6,sc3};  
sll6 = newreg; Line Loop(sll6) = {-sc11,-sc3,sc7}; 
sll7 = newreg; Line Loop(sll7) = {-sc2,-sc7,-sc12};
sll8 = newreg; Line Loop(sll8) = {-sc6,-sc9,sc2};  

srs1 = newreg; Ruled Surface(srs1) = {sll1};
srs2 = newreg; Ruled Surface(srs2) = {sll2};
srs3 = newreg; Ruled Surface(srs3) = {sll3};
srs4 = newreg; Ruled Surface(srs4) = {sll4};
srs5 = newreg; Ruled Surface(srs5) = {sll5};
srs6 = newreg; Ruled Surface(srs6) = {sll6};
srs7 = newreg; Ruled Surface(srs7) = {sll7};
srs8 = newreg; Ruled Surface(srs8) = {sll8};

// We then store the surface loops identification numbers in a list for later
// reference (we will need these to define the final volume):
ssl = newreg;

Surface Loop(ssl) = {srs8,srs5,srs1,srs2,srs3,srs7,srs6,srs4};

///////////////////////////////////////
////////// FINAL VOLUME ///////////////
cv  = newreg; Volume(cv) = {csl,ssl};

noslip = newreg;
outflow= newreg;
ball   = newreg;
Physical Surface(1) = {cps1,cps2,cps3,cps4,cps5};
Physical Surface(2) = {cps6};
Physical Surface(3) = {srs8,srs5,srs1,srs2,srs3,srs7,srs6,srs4};

Physical Volume(1) = {cv};

Field[2] = MathEval;
Field[2].F = "sqrt(x^2+y^2+z^2)";

Field[3] = MathEval;
Field[3].F = "acos(z/F2)";

Field[4] = MathEval;
Field[4].F = Sprintf("%e* min( 1.0/cos( min( max(F3-%e, 0), pi/2.0 - 0.0001)  ), %e)", d1,theta,d2/d1);

Field[1] = MathEval;
Field[1].F = Sprintf("%e + %e * max( F2-F4,0 )^%e + %e * (F2-1.0)^%e",s,a,alpha,b,alpha);
//Field[1].F = "0.1 + F4/pi";
Background Field = 1;






