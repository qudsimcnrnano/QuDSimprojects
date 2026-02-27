cl1=8;
ox1=8;


hl=18.92; // hartree factor

cr=4; // Channel radius in nm
t1=4;  // oxide1 thickness in nm


R2=cr+t1;


Point(1)={0,0,0,cl1};
Point(2)={-cr*hl,0,0,cl1};  // Radius of channel
Point(3)={cr*hl,0,0,cl1};   // Radius of channel
Point(4)={-R2*hl,0,0,ox1};  // Radius of oxide 1
Point(5)={R2*hl,0,0,ox1};   // Radius of oxide 1



Circle(20) = {2,1,3};
Circle(21) = {3,1,2};
Circle(22) = {4,1,5};
Circle(23) = {5,1,4};



Line Loop(41) = {20,21};
Line Loop(42) = {22,23};



Plane Surface(100)={41,42};
Plane Surface(101)={41};


Physical Surface(200)={100}; //oxide1 IL
Physical Surface(201)={101}; //channel
