#include <cmath>
using namespace std;
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
float function_z(float y, float t);
float y_analytical(float t);


//Declaring function z (dy/dt = z)
float function_z(float y, float t){

    float z = (t+1)/(y);
    return z;
}

//Declare Analytical Solution

float y_analytical(float t){

    float y_analytical = sqrt(pow(t,2)+(2*t)+9);
    return y_analytical;
}


class MethodClass{
	public:

		float t_final = 3;
		float y_zero = 3;
		float t_zero = 0;
		float delta_t = 0.1;
		float length_t = ((t_final - t_zero)/(delta_t))+2;
		float half_dt = 0.5*delta_t;
		vector<float> y_euler = vector<float>(length_t);
		vector<float> y_midpoint = vector<float>(length_t);
		vector<float> vector_z = vector<float>(length_t);
		vector<float> y_real = vector<float>(length_t);
		vector<float> y_ab2= vector<float>(length_t);
		vector<float> t = vector<float>(length_t);

		float AnalyticalSolution(){
   			for (int i = 0; i<length_t-1; i++){

   				t.at(i+1) = t.at(i) + delta_t; //populate time vector;
         		y_real.at(i) = y_analytical(t.at(i));

       		}

   		} 		

		
		float EulerMethod(){
			for (int i = 0; i<length_t-1; i++){

				y_euler.at(0) = y_zero; //initial condition 
				t.at(0) = t_zero;
				t.at(i+1) = t.at(i) + delta_t; //populate time vector;
				y_euler.at(i+1) = y_euler.at(i) + delta_t*function_z(y_euler.at(i),t.at(i));
								
   			}
   			
   		}


   		float MidPointMethod(){
   			 
   			 for (int i = 0; i< length_t-1; i++){
        		y_midpoint.at(0) = y_zero;
        		vector_z.at(i) = y_midpoint.at(i) + (0.5*delta_t) * function_z(y_midpoint.at(i),t.at(i)); //the half dt times function evaluated at tn part
        		y_midpoint.at(i+1) = y_midpoint.at(i) + delta_t*function_z(vector_z.at(i), t.at(i) + (0.5*delta_t));
          	}

   		}


   		float MultiStepMethod(){
   			for (int i = 0; i< length_t - 1; i++){

		        y_ab2.at(0) = y_zero;
		        float fk = function_z(y_ab2.at(i),t.at(i));
		        float fk2;
		        if (i == 0){
		            fk2 = function_z(y_ab2.at(0), t.at(0));
        		}
        		else {	
            	fk2 = function_z(y_ab2.at(i-1), t.at(i-1));	
        		}
       			y_ab2.at(i+1) = y_ab2.at(i) + (delta_t/2)*(3*fk - fk2);

   	 		}  
   		}

   		float MethodPlotting(){
   			plt::title("Numerical Methods vs Real Solution");
   			plt::plot(t,y_real,"ok");
   			plt::plot(t,y_euler,"-r");
   			plt::plot(t,y_ab2, "-g");
   			plt::plot(t,y_midpoint,"sb");
   			plt::grid(true);
   			plt::xlim(0,3);
   			plt::ylim(3,5);
   			plt::named_plot("Euler",t,y_euler,"-r");
		    plt::named_plot("Analytical", t, y_real,"ok");
		    plt::named_plot("Midpoint", t, y_midpoint,"-b");
		    plt::named_plot("AB2", t, y_ab2,"-g");
		    plt::legend();
		    plt::show();
   		}

};
