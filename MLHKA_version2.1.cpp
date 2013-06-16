
/******************************************************************************
 *  Copyright (C) 2004-2010 Stephen Wright <stephen.wright@utoronto.ca>                          *
 *  Portions of the code (C) 2010 Jeffrey Ross-Ibarra <rossibarra@ucdavis.edu>   *
 *                                                                              *
 *  This program is free software: you can redistribute it and/or modify                *
 *  it under the terms of the GNU General Public License as published by        *
 *  the Free Software Foundation, either version 3 of the License, or           *
 *  (at your option) any later version.                                         *
 *                                                                              *
 *  This program is distributed in the hope that it will be useful,             *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
 *  GNU General Public License for more details.                                *
 *                                                                              *
 *  You should have received a copy of the GNU General Public License           *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.       *
 *****************************************************************************/


#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string>
using namespace std;  //introduces namespace std
long double thetalikeli (long int, long int, long double);
long double divlikeli (long int, long double, long double);
long double factorial(long double);
long double HKAlikeli (long double[52]);

#define LOCI 100


double drand48(void);
void srand48(long int seedval);



int main( int argc, char *argv[]) 
{ 
     const char * ofn;
	ofn="likelihoods.txt";
        const char * ifn;
	ifn="infile.txt";
int opt=1, errors=0;
    long int SEED=int(rand()*1000000000), CHAIN=100000;
        // get variables
        while( opt<argc ){
                switch( argv[opt][1] ){
                        case 'o': //output file name
                                opt++;
                        	ofn = argv[opt++];
                                break;
                        case 'i': //input file name
                                opt++;
                        	ifn = argv[opt++];
                                break;
			  case 'c': // chain length
                                opt++;
                                CHAIN = atoi(argv[opt++]);
                                break;
			  case 's': // seed
                                opt++;
                                SEED = atoi(argv[opt++]);
                                break;
                        default :
                                ++errors;
                                break;
                }
        }

cout << SEED << endl;
ofstream outfile (ofn);   
ifstream infile (ifn);
if (!infile) {
       //print error message to standard error stream:
       cerr << "Can not open input file" << endl;
       exit(1);}
    long double L[LOCI];

    long double k[LOCI],maxmaxlik, maxmaxk[LOCI],d[LOCI], maxmaxtime,maxmaxth[LOCI],kmax[LOCI],scalar[LOCI],th[LOCI],Like, time, maxlik,timemax,thmax[LOCI];
    /*int name, run, start; -mpc */
    long int test, m, w, x, r, p, locus[LOCI], mod, n[LOCI], S[LOCI], y, b, h, zz;
    string locusname[LOCI], locussel[LOCI];

   infile >> r; 
   if (r>100)
   { cout << "Error: maximum 100 loci.";
     exit(1);
   }
          
   infile >> mod;
  
  
   if (mod>0)
   {
   for (zz=0;zz<mod;zz++)
    infile>>
    locussel[zz];
    }
  infile >> time; 

srand48(SEED);
   
   maxlik=-100000000;
  
   
   cout <<" number of loci: " <<r << endl;
  if (mod==0)
   {cout << " neutral model (all k=1)"<<endl;}
 
 
  else 
  
   { if (mod==r)
     { cout<< " error- number of non neutral loci=number of loci!" << endl;
       }
   else
   {
   outfile << " Testing for departure from neutrality at "; 
      
    cout << "Testing for departure from neutrality at ";
     for (m=0;m<mod;m++)
       { outfile<< locussel[m] << ", ";
         cout << locussel[m] << ", ";
 
      
       }
       
       outfile << endl;
       cout << endl;
  }
  }
  
  maxlik=-1000000000;
   maxmaxlik=-1000000000;
    
    long int kor;
    kor=0;      
 
   for (h=0; h<r; h++)
     {thmax[h]=0;
    maxmaxth[h]=0;
     k[h]=1;
     
     infile >>locusname[h]>> L[h]>>S[h] >>n[h]>>d[h]>>th[h] >> scalar[h];
  //  if (th[h]==0) {
//	  cerr << "Locus " << locusname[h] << " has a theta of 0 and will be reset to 0.0001." << endl; 
//	  th[h]=0.0001;
//	}
     
     for ( test = 0; test < mod; test++ )
     { if ( locusname[h]==locussel[test] )
     {locus[test]=h+1;
    
	
 
     }    
  }
}

long int counter,set;
set=1;
counter=0;

  timemax=0;
  for (p=0; p<CHAIN; p++)
  {counter=counter+1;
  if (counter==1000)
  { 
  long double num;
  
  num=counter*set;
  cout <<num<<endl;
    set=set+1;
    counter=0;}
    
     
  if (kor==11000)
  {kor=0;}

   Like=0; 
  
   b=0;
   
  for (w=0; w<r; w++)
 { 
 
 Like = Like+log(thetalikeli(S[w],n[w],(k[w]*L[w]*th[w]*scalar[w])))
 + (long double)(divlikeli(    (long int)(d[w]), (long double)(L[w] * th[w] * scalar[w]), (long double)(time / scalar[w])    ) );
 
 }


if (p==0)
{

if (Like<maxlik)
{outfile << "error: starting parameter values too far from true values!"<< endl <<" Try new values!";
cout << "error: starting parameter values too far from true values!"<< endl <<" Try new values!";

 exit(1);
 
 }
}

 if (Like>maxlik)
     {
    if (Like>maxmaxlik)
      {maxmaxlik=Like;
       
       maxmaxtime=time;
       for (y=0; y<r; y++)
       {
       maxmaxth[y]=th[y];
       maxmaxk[y]=k[y];
       }
      }   
    
 for (y=0; y<r; y++)
  {thmax[y]=th[y];
  kmax[y]=k[y];
  } 

      maxlik=Like;
      
  timemax=time;

  } 
     
  else 
    {
    if (kor>10000 && kor<11000)
    {
    if (drand48()<exp(0.5*(Like-maxlik)))
       {
       
       maxlik=Like;  
  
  timemax=time;
 
 for (y=0; y<r; y++)
  {thmax[y]=th[y];
  
  kmax[y]=k[y]; } 
  
  }
    
  else 
  
  {time=timemax;
   for (y=0; y<r; y++)
      {th[y]=thmax[y];}
      
  k[y]=kmax[y];
  
  }  

    }
    
    else
    {
    
    if (drand48()<exp(50*(Like-maxlik)))
       {

       maxlik=Like;  
  
  timemax=time;
 
 for (y=0; y<r; y++)
  {thmax[y]=th[y];
  
  kmax[y]=k[y]; } 
 
  }
    
else 
  
  {time=timemax;
   for (y=0; y<r; y++)
      {th[y]=thmax[y];}
      
  k[y]=kmax[y];
  
  } 
 
  }

} 
 
   
long int h;

h=(long int)(drand48()*(r+1+mod));

if (h<(1+mod))
   {
   if (h==0)
   
   {
   if (drand48()<0.5)
   {
   time=timemax-(drand48());}
   
   else
   {time=timemax+(drand48());}
   
   if(time<0)
     { time=0;}
   }
   
   
   else
   
   {
   if (drand48()<0.5)
   {
   k[locus[h-1]-1]=kmax[locus[h-1]-1]-(drand48());}
   
   else
   {k[locus[h-1]-1]=kmax[locus[h-1]-1]+(drand48());}
   
   if(k[locus[h-1]-1]<0)
     { k[locus[h-1]-1]=0;}
   }
   
   
   }
   
 else 
 {
 if (drand48()<0.5)
 {
 th[h-1-mod]=thmax[h-1-mod]-(0.01*drand48());
 
 if(th[h-1-mod]<0)
 {th[h-1-mod]=0;}
 
 }

else 
{th[h-1-mod]=thmax[h-1-mod]+(0.01*drand48());
if(th[h-1-mod]<0)
 {th[h-1-mod]=0;}

}

}


kor=kor+1;
}

long int work;
  outfile << "ML " << "T ";

for (x=0; x<r; x++)
 { outfile <<"theta(" << locusname[x] << ") k(" << locusname[x] << ") ";}
outfile << endl;
 outfile << maxmaxlik << " " << maxmaxtime<< " ";
      for (work=0; work<r; work++)
         { outfile << maxmaxth[work]<< " "<< maxmaxk[work] << " ";}
   outfile << endl;

	return 0;
}


long double thetalikeli (long int S, long int n, long double T) {  
    
   long int i, x, y, z;
   long double Q,P[500][500],X;
   
   for (i=0; i<(S+1); i++) {  
   
    P[2][i]=(pow ((T/(1+T)),i)) * (1/(1+T));

   }
   
   for (x=3; x<(n+1); x++)
    { for (y=0; y<(S+1); y++)
        { 
        
        X=0;
        for (z=0; z<(y+1); z++)
        
        {
         P[x][y]= X +( P[(x-1)][y-z]*(pow(T/(T-1+x),z))*((x-1)/(x+T-1)));
  		 X=P[x][y];  
  		
        }
      }
      
  }
      
      
	
	Q=P[n][S];
	
	return Q;
    
}

long double divlikeli (long int diverge, long double T, long double tim)



{    

     long double  K, J,C,F,S;

     int r,q;





     K=(T)*(tim);

    

    

    

    C=0;

    for (r=0; r<(diverge+1); r++)

     { 
	 F=0;
	 S=0;
	 for (q=1; q<(r+1); q++)
	 {F=F-log(q);}
	 
	 S=r*((log(K)+log(1+T)-log(T)));
	 
	 C=C+exp(S+F-K);
	 
	 
	

	
	 
	 
	 
	 }

	
J=log(1-(T/(1+T)))+(diverge*(log(T)-log(T+1)))+log(C);
    

    

    

    

 

 if (J==0)

 {cout << "Error, probability of 0!";

 

 exit(1);

 }

 

 else{

   

  return J;

    } 

}




long double factorial(long double numbe) {

	long double temp;

	if(numbe <= 1) return 1;
    
    temp=numbe;
    while (numbe>1)
    {
	temp = temp * (numbe-1);
	numbe=(numbe-1);
	}
	
	return temp;
}

