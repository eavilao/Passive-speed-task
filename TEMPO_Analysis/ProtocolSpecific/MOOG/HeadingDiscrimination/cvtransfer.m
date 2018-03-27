function  cvstar=cvtransfer(cv,meanisi)

isi=[5.0  
5.5  
6.0 
6.5  
7.0  
7.5  
12.5  
17.5  
22.5  
27.5 
32.5  
37.5  
42.5  
47.5  
52.5  ];
a=[ 0.36  
 0.40  
 0.46  
 0.53  
 0.55  
 0.56  
  0.84 
  1.15  
 1.49 
 1.66  
 1.68  
  1.80  
 1.82  
 1.88  
 1.93 ];
b=[0.63 
0.66 
0.73 
0.79 
0.80  
0.81  
0.97 
1.02  
1.04  
1.01  
0.96  
0.93  
0.91 
0.90  
0.89  ];


a_inter= interp1(isi,a,meanisi);
b_inter= interp1(isi,b,meanisi);
 

if(meanisi>52.5)
    meanisi=52.5;
     a_inter=1.93;
     b_inter=0.89;,
 
    
end
 

cvstar=(cv/a_inter)^(1/b_inter); 
 