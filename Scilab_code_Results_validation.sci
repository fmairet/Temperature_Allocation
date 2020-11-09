// Optimal proteome allocation and the temperature dependence of microbial growth laws
// Francis Mairet,Jean-Luc Gouzé, Hidde de Jong
// npj Systems Biology and Applications
//
// SI: Scilab code for Figure 3 (validation)
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
 
/////////   Experimental data (Herendeen et al., 1979; Pedersen et al., 1978)////////
// Effect of temperature
T_t=273.15+[13.5 15 23 30 37 42 46];
r_t=[0.73 0.78 0.97 1 1 0.97 0.46];
mu_t=[0.15 0.19 0.58 1.17 1.77 1.84 0.97];
c_t=[0.94 0.85 0.61 0.76 1 1.65 7.22];
 
// Effect of substrate
mu_s=[0.38 0.77 1.01 1.5 1.98];
r_s=[.58 .84 1 1.33 1.74]; // S1      
c_s=[0.94 1 1 1.18 1.22]
 
////////  Function to compute the optimal allocation (Eq. (6))  ///////////////

function [p,r,m,c,mu]=allocation(TT,pa,kMref_g)
kR=pa(1)* exp(-pa(3)/R*(1/TT-1/Tref));;
kM_g=kMref_g* exp(-pa(3)/R*(1/TT-1/Tref));
kf_ku= pa(2) * exp(-pa(6)/R*(1/TT-1/Tref)) /  (  (1-pa(4)*pa(5)/R*(1./TT-1/Tref))^(1/pa(4))  );      
p=(kM_g*K/kR)^0.5;
c=1/kf_ku * (-1 +  (  1 + (1-q)*kf_ku)^0.5 );
r=(1-q-c)* p*(K+p)/(K+2*K*p+p^2);
m=1-q-r-c;
mu=kR*p/(K+p)*r* c/ (1/kf_ku + c);
endfunction
   
/////// Computation of resource allocation in the different cases    //////////

// Effect of temperature
T=290:0.1:320;
popt=[];ropt=[];mopt=[];muopt=[];copt=[];
for i=1:length(T)
[popt(i),ropt(i),mopt(i),copt(i),muopt(i)]=allocation(T(i),param,k_Mref_rich);
end

// Effect of substrate limitation
g=0.05:0.01:1;
popt2=[];ropt2=[];mopt2=[];muopt2=[];copt2=[];
 for i=1:length(g)
[popt2(i),ropt2(i),mopt2(i),copt2(i),muopt2(i)]=allocation(Tref,param,g(i)*k_Mref_rich);
end

/////////////          Figures         /////////////////////////////////////
// Colorbar 
Tmin=min(min(T),min(T_t));
Tmax=max(max(T),max(T_t));
TN=(T-Tmin)/(Tmax-Tmin);
cmap=jetcolormap(64);

//  Ribosomal protein content
scf(1)
subplot(221)
for j=1:length(T)
 plot(muopt(j),(ropt(j)+rq)*100,'.','MarkerSize',7,'Color',cmap(1+floor(63*TN(j)),: ))    
end
a=gca();
a.zoom_box=[0,15,2.4,42];
ylabel('Optimal mass fraction r (%)')
title("Effect of temperature",'fontsize',2)

subplot(222)
for j=1:length(muopt2)
 plot(muopt2(j),(ropt2(j)+rq)*100,'.','MarkerSize',7,'Color',cmap(1+floor(63*TN(find(T>Tref,1))),:))   
end
a=gca();
a.zoom_box=[0,15,2.4,42];
title("Effect of substrate limitation",'fontsize',2)
 
subplot(223)
for j=1:length(mu_t)
  plot(mu_t(j),r_t(j),'.','MarkerSize',7,'Color',cmap(1+floor(63*(T_t(j)-Tmin)/(Tmax-Tmin)),:))  
end
a=gca();
a.zoom_box=[0,0.2,2.4,1.1];
xlabel('Specific growth rate (/h)')
ylabel('Ribosomal protein S1 (relative level)')

subplot(224)
for j=1:length(mu_s)
  plot(mu_s(j),r_s(j)/r_s(length(mu_s)),'.','MarkerSize',7,'Color',cmap(1+floor(63*(Tref-Tmin)/(Tmax-Tmin)),:))  
end
a=gca();
a.zoom_box=[0,0.2,2.4,1.1];
xlabel('Specific growth rate (/h)')

//  Chaperone content
scf(2)
subplot(221)
for j=1:length(T)
 plot(muopt(j),100*copt(j),'.','MarkerSize',7,'Color',cmap(1+floor(63*TN(j)),:))    
end
a=gca();
a.zoom_box=[0,0,2.4,25];
ylabel('Optimal mass fraction c (%)')
title("Effect of temperature",'fontsize',2)

subplot(222)
for j=1:length(muopt2)
 plot(muopt2(j),100*copt2(j),'.','MarkerSize',7,'Color',cmap(1+floor(63*TN(find(T>Tref,1))),:))   
end
a=gca();
a.zoom_box=[0,0,2.4,25];
title("Effect of substrate limitation",'fontsize',2)
    
subplot(223)
for j=1:length(mu_t)
 plot(mu_t(j),c_t(j),'.','MarkerSize',7,'Color',cmap(1+floor(63*(T_t(j)-Tmin)/(Tmax-Tmin)),:))  
end
a=gca();
a.zoom_box=[0,0,2.4,8];
xlabel('Specific growth rate (/h)')
ylabel('Chaperone groEL (relative level) ')

subplot(224)
for j=1:length(mu_s)
  plot(mu_s(j),c_s(j)/c_s(length(mu_s)),'.','MarkerSize',7,'Color',cmap(1+floor(63*(Tref-Tmin)/(Tmax-Tmin)),:))  
end
a=gca();
a.zoom_box=[0,0,2.4,8];
xlabel('Specific growth rate (/h)')








