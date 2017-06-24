#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include "/usr/local/MATLAB/R2017a/extern/include/engine.h"

using namespace std;
#define BUFFSIZE 256

float function_z(float y, float t);
float y_analytical(float t);


//Declaring function z (dy/dt = z)
float diffeq(float y, float t){

    float z = (t+1)/(y);
    return z;
}

//Declare Analytical Solution

float ComputeAnalytical(float t){

    float y_real = sqrt(pow(t,2)+(2*t)+9);
    return y_real;
}






int main(){

		//START OF C PROGRAMMING
		//---------------------------------------------------------------------------------

		int t_final = 3;
		double y_zero = 3;
		int t_zero = 0;
		float delta_t = 0.1;
		int length_t = ((t_final - t_zero)/(delta_t))+2;
		double half_dt = 0.5*delta_t;
		double time[length_t];
		double y_euler[length_t];
		double y_analytical[length_t];
		double y_midpoint[length_t];
		double y_ab2[length_t];
		double vector_z[length_t]; //for midpoint

			

		//ICS

		time[0] = 0;
		y_euler[0] = y_zero; 
		y_midpoint[0] = y_zero;
		y_ab2[0] = y_zero;
		

		//Generate Time Vector and Y analytical
   		for (int i = 0; i<length_t; i++){

   				time[i+1] = time[i] + delta_t; //populate time vector;
   				
         		y_analytical[i] = ComputeAnalytical(time[i]);
         		//cout << time[i] << "\t" << y_analytical[i] << endl;

       	}

   		
   		
       	// Euler method
		
		for (int i = 0; i<length_t; i++){
			
				y_euler[i+1] = y_euler[i] + delta_t*diffeq(y_euler[i],time[i]);


   		}

   		//Midpoint

   		for (int i = 0; i< length_t; i++){
        		vector_z[i] = y_midpoint[i] + (0.5*delta_t) * diffeq(y_midpoint[i],time[i]); //the half dt times function evaluated at tn part
        		y_midpoint[i+1] = y_midpoint[i] + delta_t*diffeq(vector_z[i], time[i] + (0.5*delta_t));

         }



         //Multistep ADAMS BADSFORTH 2 

         for (int i = 0; i< length_t; i++){

		        float fk = diffeq(y_ab2[i],time[i]);
		        float fk2;
		        if (i == 0){
		            fk2 = diffeq(y_ab2[0], time[0]);
        		}
        		else {	
            	fk2 = diffeq(y_ab2[i-1], time[i-1]);	
        		}
       			y_ab2[i+1] = y_ab2[i] + (delta_t/2)*(3*fk - fk2);

   	 	}

	
   		//-----------------------------------------------------------------------------------
		
		//START OF MATLAB PROGRAMMING

		Engine *ep;

		//Create pointer arrays

		mxArray *time_pointer = NULL;
		mxArray *y_analytical_pointer = NULL;
		mxArray *y_midpoint_pointer = NULL;
		mxArray *y_ab2_pointer = NULL;
		mxArray *y_euler_pointer = NULL, *result = NULL;


		char buffer[BUFFSIZE+1];


		//check for error
		if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return EXIT_FAILURE;
		}


		//PUT C INTO MATLAB
		time_pointer = mxCreateDoubleMatrix(1,length_t, mxREAL);
		memcpy((void *)mxGetPr(time_pointer), (void *)time,sizeof(time));
		engPutVariable(ep,"time_matlab",time_pointer);

		

		y_analytical_pointer = mxCreateDoubleMatrix(1,length_t, mxREAL);	
		memcpy((void *)mxGetPr(y_analytical_pointer), (void *)y_analytical,sizeof(y_analytical));
		engPutVariable(ep,"y_analytical_matlab",y_analytical_pointer);
		

		y_euler_pointer = mxCreateDoubleMatrix(1,length_t, mxREAL);
		memcpy((void *)mxGetPr(y_euler_pointer), (void *)y_euler,sizeof(y_euler));
		engPutVariable(ep,"y_euler_matlab",y_euler_pointer);


		y_midpoint_pointer = mxCreateDoubleMatrix(1,length_t, mxREAL);
		memcpy((void *)mxGetPr(y_midpoint_pointer), (void *)y_midpoint,sizeof(y_midpoint));
		engPutVariable(ep,"y_midpoint_matlab",y_midpoint_pointer);

		y_ab2_pointer = mxCreateDoubleMatrix(1,length_t, mxREAL);
		memcpy((void *)mxGetPr(y_ab2_pointer), (void *)y_ab2,sizeof(y_ab2));
		engPutVariable(ep,"y_ab2_matlab",y_ab2_pointer);


		
		engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure
		

		engEvalString(ep,"f = plot(time_matlab,y_analytical_matlab,'-k','Linewidth',2);");
		engEvalString(ep,"title('All numerical methods vs Analytical solution');");
		engEvalString(ep,"xlabel('Time [seconds]');");
		engEvalString(ep,"ylabel('Function values');");
		engEvalString(ep,"hold on;");
		engEvalString(ep,"grid on;");
		engEvalString(ep,"box on;");
		engEvalString(ep,"axis([0 3 3 5])");
		engEvalString(ep,"f = plot(time_matlab,y_euler_matlab,'xg','Linewidth',0.5);");
		engEvalString(ep,"f = plot(time_matlab,y_midpoint_matlab,'or','Linewidth',0.5);");
		engEvalString(ep,"f = plot(time_matlab,y_ab2_matlab,'*b','Linewidth',0.5);");
		engEvalString(ep,"legend('Analytical Soln ','Euler Method','Midpoint Method','Adams-Badsforth Multistep Method','Location','northwest');");
		engEvalString(ep,"saveas(f,'allMethod.png')");
		engEvalString(ep,"close");
		

		return EXIT_SUCCESS;

}

   			
   		


   	


