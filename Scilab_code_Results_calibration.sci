// Optimal proteome allocation and the temperature dependence of microbial growth laws
// Francis Mairet,Jean-Luc Gouzé, Hidde de Jong
// npj Systems Biology and Applications
//
// SI: Scilab code for Figure 2 (calibration)
xdel(winsid())
clear
//////////////////////       Parameters     //////////////
R=8.314;
Tref=310;
rq= 0.138;  
q=0.46;
K=0.003;  

// param=[ k_Rref  (k_f/k_u)ref   E_R    d     E_u   E_f  ]
param=[   8.42057  356.8579     41921  -0.0944764  575311  488858 ];

k_Mref= 1.54389974;
k_Mref_rich=  9.98920676;

//////////////////   Experimental data   ///////////////////////
///  thermal growth curve in minimal medium (Smirnova et al., 2016)
T_E1=273.15+[20 24 30 37 40 42 44 46];
mu_E1=[0.1384 0.2431 0.4120 0.631 0.7025 0.6980 0.5383 0.07928];

/////  thermal growth curve in rich medium (Barber, 1908)
T_E2=[283.15 285.54 287.04 289.55 292.65 294.75 296.75 298.65 300.75 302.75...
 304.75 306.75 308.15 310.85 311.65 313.95 315.95 317.75 318.35 319.95 320.84];
mu_E2=[0.053 0.0792 0.2296 0.34428 0.4181 0.6807 0.8362 1.0119 1.1382 1.2997...
 1.5620 1.8082 1.8070 2.0832 2.25306 2.34 2.04117 1.9236 1.6049 0.8546 0.016];

///// Proteomic data from Schmidt et al. (2016)
// Batch 37 °C
r37_E=0.2678;
c37_E=.0372;  
m37_E=0.374;

// Batch 42 °C 
r42_E=0.2535;
c42_E=0.0663;
m42_E=0.365;

// Batch 37 °C in rich medium (LB)
rLB_E=0.394
cLB_E=0.0628;
mLB_E=0.2591;

// Chemostat 37 °C
muChem_E=[0.5 0.35  0.2 0.12];
rChem_E=[0.1956 0.1729 0.1607 0.1529];
cChem_E=[0.0327 0.0318 0.0287 0.028];
mChem_E=[0.443 0.477 0.49 0.498];

////////  Function to compute the optimal allocation (Eq. (6))  ///////////////
function [p,r,m,c,mu]=allocation(TT,pa,kMref_g)
kR=pa(1)* exp(-pa(3)/R*(1/TT-1/Tref));;
kM_g=kMref_g* exp(-pa(3)/R*(1/TT-1/Tref));
kf_ku= pa(2) * exp(-pa(6)/R*(1/TT-1/Tref)) / (  (1-pa(4)*pa(5)/R*(1./TT-1/Tref))^(1/pa(4))  );    
p=(kM_g*K/kR)^0.5;
c=1/kf_ku * (-1 +  (  1 + (1-q)*kf_ku)^0.5 );
r=(1-q-c)* p*(K+p)/(K+2*K*p+p^2);
m=1-q-r-c;
mu=kR*p/(K+p)*r* c/ (1/kf_ku + c);
endfunction
   
/////// Computation of resource allocation in the different cases    //////////

T=280:0.2:322;

//// thermal growth curve in minimal medium
popt=[];ropt=[];mopt=[];muopt=[];copt=[];
for i=1:length(T)
[popt(i),ropt(i),mopt(i),copt(i),muopt(i)]=allocation(T(i),param,k_Mref);
end
mu_1=muopt;

// allocation in batch at 37 and 42 °C
r37=(rq+ropt(find(T>37+273.15,1)))
c37=(copt(find(T>37+273.15,1)))
m37=(mopt(find(T>37+273.15,1)))

r42=(rq+ropt(find(T>42+273.15,1)))
c42=(copt(find(T>42+273.15,1)))
m42=(mopt(find(T>42+273.15,1)))

// allocation in chemostat
g=0.01:0.01:1;
popt=[];ropt=[];mopt=[];muopt=[];copt=[];  
for i=1:length(g)
[popt(i),ropt(i),mopt(i),copt(i),muopt(i)]=allocation(Tref,param,g(i)*k_Mref);
end

rChem=rq+interp1(muopt,ropt,muChem_E);
cChem=interp1(muopt,copt,muChem_E);
mChem=interp1(muopt,mopt,muChem_E);

//// thermal growth curve in rich medium
popt=[];ropt=[];mopt=[];muopt=[];copt=[];
for i=1:length(T)
[popt(i),ropt(i),mopt(i),copt(i),muopt(i)]=allocation(T(i),param,k_Mref_rich);
end
mu_2=muopt;

//  allocation in batch in rich medium
rLB=(rq+ropt(find(T>37+273.15,1)))
cLB=(copt(find(T>37+273.15,1)))
mLB=(mopt(find(T>37+273.15,1)))

///////////////////////////       Figures           ///////////////////////////
//// Thermal growth curves
scf(1)
plot(T,mu_1,T,mu_2)
plot(T_E1,mu_E1,'o',T_E2,mu_E2,'x')
xlabel('Temperature (K)')
ylabel('Specific growth rate (/h)')

////  Proteome allocation
A=[];
A(1,:)= [rChem_E(4) mChem_E(4) cChem_E(4) 1-rChem_E(4)-mChem_E(4)-cChem_E(4)]; 
A(2,:)= [rChem_E(3) mChem_E(3) cChem_E(3) 1-rChem_E(3)-mChem_E(3)-cChem_E(3)]; 
A(3,:)= [rChem_E(2) mChem_E(2) cChem_E(2) 1-rChem_E(2)-mChem_E(2)-cChem_E(2)];
A(4,:)= [rChem_E(1) mChem_E(1) cChem_E(1) 1-rChem_E(1)-mChem_E(1)-cChem_E(1)];
A(5,:)= [r37_E m37_E c37_E 1-r37_E-m37_E-c37_E]; 
A(6,:)= [r42_E m42_E c42_E 1-r42_E-m42_E-c42_E]; 
A(7,:)= [rLB_E mLB_E cLB_E 1-rLB_E-mLB_E-cLB_E];
A=100*A;

B=[];
B(1,:)= [rChem(4) mChem(4) cChem(4) q-rq]; 
B(2,:)= [ rChem(3) mChem(3) cChem(3) q-rq]; 
B(3,:)= [ rChem(2) mChem(2) cChem(2) q-rq];
B(4,:)= [rChem(1) mChem(1) cChem(1) q-rq];
B(5,:)= [ r37 m37 c37 q-rq]; 
B(6,:)= [r42 m42 c42 q-rq]; 
B(7,:)= [rLB mLB cLB q-rq];
B=100*B;

scf(3)
subplot(311)
bar([A(:,1) B(:,1)])
legend(['Experiments','Model'],2)
ylabel('Ribosomal proteins (%)')
a=gca();
a.x_ticks.labels=["Chem 0.12/h";"Chem 0.2/h";"Chem 0.35/h";"Chem 0.5/h";...
"Batch 37";"Batch 42";"Batch 37 rich"]

subplot(312)
bar([A(:,2) B(:,2)])
ylabel('Metabolic proteins (%)')
a=gca();
a.x_ticks.labels=["Chem 0.12/h";"Chem 0.2/h";"Chem 0.35/h";"Chem 0.5/h";...
"Batch 37";"Batch 42";"Batch 37 rich"]

subplot(313)
bar([A(:,3) B(:,3)])
ylabel('Chaperones (%)')
a=gca();
a.x_ticks.labels=["Chem 0.12/h";"Chem 0.2/h";"Chem 0.35/h";"Chem 0.5/h";...
"Batch 37";"Batch 42";"Batch 37 rich"]


