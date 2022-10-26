
//***************************************************************************************

// Code developed by Mojtaba Madadi Asl.
// Email: mojtabamadadi7@gmail.com

//***************************************************************************************

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "random_number.cpp"
#include "spike_statistics.cpp"

using namespace std;

//*****************************************************************

int main(){

double PI = 4. * atan(1);
double dt = 0.1;
int T = 250000;
long iseed = -64L;

bool STIMULATION;
int t_ON1, t_ON2, t_OFF, inter_step, intra_step, i, j, k, t, m1, m2, tf1, tf2, tff1, tff2, tbb1, tbb2, nnb1, nnb2, nnf1, nnf2, RANGE, c, counter1, counter2, f1, f2, win, RANGE1;
float V1[T], V2[T], S1[T], S2[T], g12, g21, tau, Vrest, Vth, tau_m, IEXT1, IEXT2, ISTIM1, ISTIM2, ta, td, noise1, noise2, gmin, gmax, A_plus, A_minus, tau_plus, tau_minus, delta_g, timelag, Tmax, Tmin, dx; 

vector<float> ISI1_ON;
vector<float> ISI2_ON;
vector<float> F1_ON;
vector<float> F2_ON;
vector<float> DT1_ON;
vector<float> DT2_ON;

vector<float> ISI1_OFF;
vector<float> ISI2_OFF;
vector<float> F1_OFF;
vector<float> F2_OFF;
vector<float> DT1_OFF;
vector<float> DT2_OFF;

ofstream outfile0("activity_weights.txt");
ofstream outfile1("firing_rates.txt");
ofstream outfile2("ISI_ON.txt");
ofstream outfile3("ISI_OFF.txt");

//*****************************************************************

V1[0] = ran2(&iseed);
V2[0] = ran2(&iseed);
Vrest = 0.0;
Vth = 1.0;
tau_m = 10.;
g12 = 0.5;
g21 = 0.5;
tf1 = 0;
tf2 = 0;
tff1 = 0;
tff2 = 0;
tbb1 = 0;
tbb2 = 0;
nnf1 = -2000;
nnf2 = -2000;
nnb1 = -2000;
nnb2 = -2000;
IEXT1 = 0.882;
IEXT2 = 0.882;
gmin = 0.05;
gmax = 1.0;
A_plus = 0.008;
A_minus = 0.005;
tau_plus = 10;
tau_minus = 20;
c = 0;
counter1 = 0;
counter2 = 0;
f1 = 0;
f2 = 0;
win = 5000;

inter_step = 10000;
intra_step = 300;
t_ON1 = 100000;
t_ON2 = 100150;
t_OFF = 155000;
m1 = 0;
m2 = 0;
STIMULATION = false;

td = 0.5;
ta = 10.5;
tau = int((td + ta) / dt);

//*******************************************************

Tmin = 0;
Tmax = 30;
dx = 0.5;

RANGE = ((Tmax - Tmin) / dx) + 1;
float CC[RANGE], XX[RANGE];
for (i = 0; i < RANGE; i++) XX[i] = i * dx + Tmin;

//****************************************************************************


for (t = 0; t < T - 1; t++){

	noise1 = 0.2 * gasdev(&iseed);

    if (t == t_ON1) STIMULATION = true;
    if (t == t_OFF) STIMULATION = false;

    if (STIMULATION == true){
        if (t == (t_ON1 + (m1 * intra_step))){
            ISTIM1 = 100.0;
            m1++;
            if (m1 == 5){
                m1 = 0;
                t_ON1 = t_ON1 + inter_step;
            }
        }
    }


	V1[t + 1] = V1[t] + dt * (1. / tau_m) * (Vrest - V1[t] + (noise1 / sqrt(dt)) + g12 * Delta_Function(t - tf2 - tau, dt) + IEXT1 + ISTIM1);

	ISTIM1 = 0.0;

	if (V1[t] > Vth){
		V1[t + 1] = 0.0;
	        S1[t] = 1.0;
        	tf1 = t;

		if (t < t_ON1){
			F1_ON.push_back(t * dt);
			DT1_ON.push_back(abs(tf1 - tf2) * dt);
		}

		if (t > t_OFF){
			F1_OFF.push_back(t * dt);
			DT1_OFF.push_back(abs(tf1 - tf2) * dt);
		}		
	}

	tff1 = tf1 + int(ta / dt);
	tbb1 = tf1 + int(td / dt);

    if (t == tff1){
        timelag = abs(tff1 - nnb2) * dt;
        nnf1 = tff1;
        if (timelag != 0) delta_g = A_minus * exp(- timelag / tau_minus);
        else delta_g = 0.0;
        g21 = g21 - delta_g;
    }

    if (t == tbb1){
        timelag = abs(tbb1 - nnf2) * dt;
        nnb1 = tbb1;
        if (timelag != 0) delta_g = A_plus * exp(- timelag / tau_plus);
        else delta_g = 0.0;
        g12 = g12 + delta_g;
    }
                    
    if (g12 > gmax) g12 = gmax;
    if (g21 > gmax) g21 = gmax;
    if (g12 < gmin) g12 = gmin;
    if (g21 < gmin) g21 = gmin;	

	noise2 = 0.2 * gasdev(&iseed);

    if (STIMULATION == true){
        if (t == (t_ON2 + (m2 * intra_step))){
            ISTIM2 = 100.0;
            m2++;
            if (m2 == 5){
                m2 = 0;
                t_ON2 = t_ON2 + inter_step;
            }
        }
    }


	V2[t + 1] = V2[t] + dt * (1. / tau_m) * (Vrest - V2[t] + (noise2 / sqrt(dt)) + g21 * Delta_Function(t - tf1 - tau, dt) + IEXT2 + ISTIM2);

	ISTIM2 = 0.0;

	if (V2[t] > Vth){
		V2[t + 1] = 0.0;
        	S2[t] = 1.0;
        	tf2 = t;

		if (t < t_ON1){
			F2_ON.push_back(t * dt);
			DT2_ON.push_back(abs(tf1 - tf2) * dt);
		}

		if (t > t_OFF){
			F2_OFF.push_back(t * dt);
			DT2_OFF.push_back(abs(tf1 - tf2) * dt);
		}		
	}

	tff2 = tf2 + int(ta / dt);
	tbb2 = tf2 + int(td / dt);

    if (t == tff2){
        timelag = abs(tff2 - nnb1) * dt;
        nnf2 = tff2;
        if (timelag != 0) delta_g = A_minus * exp(- timelag / tau_minus);
        else delta_g = 0.0;
        g12 = g12 - delta_g;
    }

    if (t == tbb2){
        timelag = abs(tbb2 - nnf1) * dt;
        nnb2 = tbb2;
        if (timelag != 0) delta_g = A_plus * exp(- timelag / tau_plus);
        else delta_g = 0.0;
        g21 = g21 + delta_g;
    }
                    
    if (g12 > gmax) g12 = gmax;
    if (g21 > gmax) g21 = gmax;
    if (g12 < gmin) g12 = gmin;
    if (g21 < gmin) g21 = gmin;	

    outfile0 << t * dt / 1000. << '\t' << S1[t] << '\t' << S2[t] << '\t' << g21 << '\t' << g12 << endl;
}


// Firing rate ***********************************************************************************************

RANGE1 = (T / win) + 1;
float CC1[RANGE1], CC2[RANGE1], XX1[RANGE1];

	for (i = 0; i < RANGE1; i++){
		XX1[i] = i * win;
		
		for (j = XX1[i]; j < XX1[i] + win; j++){
			if (S1[j] == 1) f1++;
			if (S2[j] == 1) f2++;
		}

		CC1[i] = f1 / (win * dt / 1000.);
		CC2[i] = f2 / (win * dt / 1000.);

		f1 = 0;
		f2 = 0;
	}

	for(i = 0; i < RANGE1; i++) outfile1 << XX1[i] * dt / 1000. << '\t' << CC1[i] << '\t' << CC2[i] << endl;


// ISI *******************************************************************************************************

        for (j = 0; j < F1_ON.size(); j++) ISI1_ON.push_back(F1_ON[j+1] - F1_ON[j]);
        for (j = 0; j < F2_ON.size(); j++) ISI2_ON.push_back(F2_ON[j+1] - F2_ON[j]);
        for (j = 0; j < F1_OFF.size(); j++) ISI1_OFF.push_back(F1_OFF[j+1] - F1_OFF[j]);
        for (j = 0; j < F2_OFF.size(); j++) ISI2_OFF.push_back(F2_OFF[j+1] - F2_OFF[j]);

	for (j = 0; j < F1_ON.size(); j++) outfile2 << DT1_ON[j] << '\t' << ISI1_ON[j] << '\t' << DT2_ON[j] << '\t' << ISI2_ON[j] << endl;
	for (j = 0; j < F1_OFF.size(); j++) outfile3 << DT1_OFF[j] << '\t' << ISI1_OFF[j] << '\t' << DT2_OFF[j] << '\t' << ISI2_OFF[j] << endl;


//************************************************************************************************************

outfile0.close();
outfile1.close();
outfile2.close();
outfile3.close();

return 0;

}

// The End *****************************************************

