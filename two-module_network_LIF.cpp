
//***************************************************************************************

// Code developed by Mojtaba Madadi Asl.
// Email: mojtabamadadi7@gmail.com

//***************************************************************************************

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "random_number.cpp"
#include "spike_statistics.cpp"

using namespace std;

//**************************************************************************

const int NODE = 200;
const int T = 250000;

bool STIMULATION;
double I1[NODE], I2[NODE], V1[NODE][T], V2[NODE][T], g_intra11[NODE][NODE], g_intra22[NODE][NODE], g_inter21[NODE][NODE], g_inter12[NODE][NODE], r1[NODE][T], r2[NODE][T], td_intra[NODE][NODE], ta_intra[NODE][NODE], tau_intra[NODE][NODE], td_inter[NODE][NODE], ta_inter[NODE][NODE], tau_inter[NODE][NODE], g_mean21[T], g_mean12[T], LOOP[NODE][NODE], kint[NODE][NODE], NET1_bin[T], NET2_bin[T], R1[T], R2[T], Phi1[T], Phi2[T];
int t_ON1, t_ON2, t_OFF, inter_step, intra_step, m1, m2, t, i, j, k, m, Lnum, r33, c, c1, c2, RANGE, RANGE1, adj_intra11[NODE][NODE], adj_intra22[NODE][NODE], adj_inter12[NODE][NODE], adj_inter21[NODE][NODE], tf1[NODE], tf2[NODE], inh_num1, inh_num2, t_ref1[NODE], t_ref2[NODE], tff1[NODE], tbb1[NODE], nnf1[NODE], nnb1[NODE], tff2[NODE], tbb2[NODE], nnf2[NODE], nnb2[NODE], counter1, counter2, counter3, counter4, counter_syn1, counter_syn2, counter_syn3, counter_syn4, total, A1[NODE][T], A2[NODE][T];
double dx, dx1, dummy1, dummy11, dummy2, dummy22, r11, r22, TT, tau_plus, tau_minus, A_plus, A_minus, noise, ISYN, IEXT, ISTIM1, ISTIM2, xrand, rate, gmax, gmin, Vrest, Vth, NET1, NET2, tau_syn, tau_m, tau_r, MEANex, MEANin, SD, delta_g, timelag, dummy3, dummy4, irr_min, irr_max, sum, h, loops, dG, rx1, ry1, rx2, ry2, sum_phi1, sum_phi2;

vector<double> signal1;
vector<double> signal2;

vector<double> ISI1;
vector<double> ISI2;
vector<vector<double> > SPIKE_TRAIN1;
vector<vector<double> > SPIKE_TRAIN2;

//**************************************************************************

int main(){

bool FLAG;
long iseed = -64L;
long iseed1 = -70L;
long iseed2 = -84L;
const double PI = 4. * atan(1.);
const double dt = 0.1;
const double rate = 1.0;
const double PP = rate * 0.001;
const double P_INTRA = 0.15;
const double P_INTER = 0.15;

ofstream outfile0("initial_strengths_M1.txt");
ofstream outfile1("initial_strengths_M1.txt");
ofstream outfile2("phase_distribution.txt");
ofstream outfile3("mean_coupling.txt");
ofstream outfile4("raster_M1.txt");
ofstream outfile5("activity_M1.txt");
ofstream outfile6("raster_M2.txt");
ofstream outfile7("activity_M2.txt");
ofstream outfile8("final_strengths.txt");
ofstream outfile9("irregularity.txt");
ofstream outfile10("bin_population_activity.txt");
ofstream outfile11("fft.txt");
ofstream outfile12("correlation.txt");

tau_plus = 10;
tau_minus = 20;
A_plus = 0.008;
A_minus = 0.005;
gmax = 1.0;
gmin = 0.05;
Vrest = 0.0;
Vth = 1.0;
tau_m = 10;
tau_r = 10.0;
tau_syn = 5.0;
ISYN = 0.0;
IEXT = 0.0;
ISTIM1 = 0.0;
ISTIM2 = 0.0;

inh_num1 = 0;
inh_num2 = 0;
SD = 0.05;
MEANex = 0.2;
MEANin = 0.8;
NET1 = 0;
NET2 = 0;
counter1 = 0;
counter2 = 0;
counter3 = 0;
counter4 = 0;
counter_syn1 = 0;
counter_syn2 = 0;
counter_syn3 = 0;
counter_syn4 = 0;
dummy1 = 0.0;
dummy11 = 0.0;
dummy2 = 0.0;
dummy22 = 0.0;
dummy3 = 0.0;
dummy4 = 0.0;
rx1 = 0.0;
ry1 = 0.0;
rx2 = 0.0;
ry2 = 0.0;
sum = 0.;
sum_phi1 = 0;
sum_phi2 = 0;
h = 0.3;
loops = 0.;
Lnum = 0;
total = 0;
dx = 0.05;
c = 0;

inter_step = 10000;
intra_step = 300;
t_ON1 = 100000;
t_ON2 = 100150;
t_OFF = 155000;
m1 = 0;
m2 = 0;

STIMULATION = false;

RANGE = ((gmax - gmin) / dx) + 2;
double CC[RANGE], CC0[RANGE], XX[RANGE];

// Zero initialization *****************************************************************************

    for (i = 0; i < NODE; i++){
        V1[i][0] = ran2(&iseed) * Vth;
        V2[i][0] = ran2(&iseed) * Vth;
        r1[i][0] = 5.0;
        r2[i][0] = 5.0;
        I1[i] = 0.0;
        I2[i] = 0.0;
        t_ref1[i] = 0;
        t_ref2[i] = 0;
        tf1[i] = -2000;
        tf2[i] = -2000;
	tff1[i] = 0;
        tbb1[i] = 0;
        nnf1[i] = -2000;
        nnb1[i] = -2000;
	tff2[i] = 0;
        tbb2[i] = 0;
        nnf2[i] = -2000;
        nnb2[i] = -2000;

        for (j = 0; j < NODE; j++){
            g_intra11[i][j] = 0.0;
            g_intra22[i][j] = 0.0;
            g_inter12[i][j] = 0.0;
            g_inter21[i][j] = 0.0;
            adj_intra11[i][j] = 0;
            adj_intra22[i][j] = 0;
            adj_inter12[i][j] = 0;
            adj_inter21[i][j] = 0;
	    ta_intra[i][j] = 0.0;
            td_intra[i][j] = 0.0;
            tau_intra[i][j] = ta_intra[i][j] + td_intra[i][j];
	    ta_inter[i][j] = 10.5;
            td_inter[i][j] = 0.5;
            tau_inter[i][j] = ta_inter[i][j] + td_inter[i][j];
            }
        }

    for (i = 0; i < NODE; i++){
        for (j = 0; j < T; j++){
            A1[i][j] = 0;
            A2[i][j] = 0;
        }
    }



// Initializing random network ****************************************************************************

    
        for (i = 0; i < NODE; i++){
            for (j = 0; j < i; j++){

                r11 = ran2(&iseed);

                if (r11 <= P_INTRA){
                     adj_intra11[i][j] = 1;
                     if (r11 <= P_INTRA * P_INTRA) adj_intra11[j][i] = 1;
                }

                r11 = ran2(&iseed);

                if (r11 <= P_INTRA){
                     adj_intra22[i][j] = 1;
                     if (r11 <= P_INTRA * P_INTRA) adj_intra22[j][i] = 1;
                }
            }
        }

        for (i = 0; i < NODE; i++){
            for (j = 0; j < NODE; j++){

                r11 = ran2(&iseed);

                if (r11 <= P_INTER){
                     adj_inter21[i][j] = 1;
                     if (r11 <= P_INTER * P_INTER) adj_inter12[j][i] = 1;
                }

                r11 = ran2(&iseed);

                if (r11 <= P_INTER){
                     adj_inter12[i][j] = 1;
                     if (r11 <= P_INTER * P_INTER) adj_inter21[j][i] = 1;
                }
            }
        }

            while(inh_num1 < (0.2 * NODE)){
                r33 = ran2(&iseed) * NODE;
                FLAG = false;

            for (j = 0; j < NODE; j++){
                if (adj_intra11[j][r33] == 1){
                    adj_intra11[j][r33] = -1;
                    adj_inter21[j][r33] = 0;
                    FLAG = true;
                }
            }
                    if (FLAG = true) inh_num1++;
            }


            while(inh_num2 < (0.2 * NODE)){
                r33 = ran2(&iseed) * NODE;
                FLAG = false;

            for (j = 0; j < NODE; j++){
                if (adj_intra22[j][r33] == 1){
                    adj_intra22[j][r33] = -1;
                    adj_inter12[j][r33] = 0;
                    FLAG = true;
                }
            }
                    if (FLAG = true) inh_num2++;
            }


// Initializing synaptic strength: Guassian distribution ************************************************************


         for (i = 0; i < NODE; i++){
             for (j = 0; j < NODE; j++){

                 if (adj_intra11[i][j] == 1){
                    g_intra11[i][j] = (gasdev(&iseed) * SD) + MEANex;
                    while (g_intra11[i][j] < gmin || g_intra11[i][j] > gmax) g_intra11[i][j] = (gasdev(&iseed) * SD) + MEANex;
                 }

                 if (adj_intra11[i][j] == -1){
                    g_intra11[i][j] = (gasdev(&iseed) * SD) + MEANin;
                    while (g_intra11[i][j] < gmin || g_intra11[i][j] > gmax) g_intra11[i][j] = (gasdev(&iseed) * SD) + MEANin;
                 }

                 if (adj_inter21[i][j] == 1){
                    g_inter21[i][j] = (gasdev(&iseed) * SD) + MEANex;
                    while (g_inter21[i][j] < gmin || g_inter21[i][j] > gmax) g_inter21[i][j] = (gasdev(&iseed) * SD) + MEANex;
                 }

                 if (adj_intra22[i][j] == 1){
                    g_intra22[i][j] = (gasdev(&iseed) * SD) + MEANex;
                    while (g_intra22[i][j] < gmin || g_intra22[i][j] > gmax) g_intra22[i][j] = (gasdev(&iseed) * SD) + MEANex;
                 }

                 if (adj_intra22[i][j] == -1){
                    g_intra22[i][j] = (gasdev(&iseed) * SD) + MEANin;
                    while (g_intra22[i][j] < gmin || g_intra22[i][j] > gmax) g_intra22[i][j] = (gasdev(&iseed) * SD) + MEANin;
                 }

                 if (adj_inter12[i][j] == 1){
                    g_inter12[i][j] = (gasdev(&iseed) * SD) + MEANex;
                    while (g_inter12[i][j] < gmin || g_inter12[i][j] > gmax) g_inter12[i][j] = (gasdev(&iseed) * SD) + MEANex;
                 }
            }
         }

// Initial Strengths Distribution: Bin diagrams ***********************************************************************

    for (i = 0; i < RANGE; i++) XX[i] = i * dx + gmin;

    for (i = 0; i < NODE; i++){
        for (j = 0; j < NODE; j++){
		    if (adj_inter12[i][j] != 0){
                c = int(g_inter12[i][j] / dx);
                CC0[c - 1]++;
			    counter4++;
			}
			if (adj_inter21[i][j] != 0){
			    c = int(g_inter21[i][j] / dx);
                CC[c - 1]++;
			    counter3++;
			}
        }
    }

    for(i = 0; i < RANGE; i++) outfile0 << XX[i] << '\t' << CC [i] / double(counter3) << endl;
    for(i = 0; i < RANGE; i++) outfile1 << XX[i] << '\t' << CC0[i] / double(counter4) << endl;

//**************************************************************************************************
//**************************************************************************************************

    for (t = 0; t < T; t++){


// Calculating mean coupling and phase diagrams *************************************


		for (i = 0; i < NODE; i++){
                    rx1 = rx1 + (1. / double(NODE)) * cos(V1[i][t] * 2.0 * PI);
                    ry1 = ry1 + (1. / double(NODE)) * sin(V1[i][t] * 2.0 * PI);
                    rx2 = rx2 + (1. / double(NODE)) * cos(V2[i][t] * 2.0 * PI);
                    ry2 = ry2 + (1. / double(NODE)) * sin(V2[i][t] * 2.0 * PI);

		    sum_phi1 = sum_phi1 + ((1.0 / double(NODE)) * V1[i][t] * 2.0 * PI);
		    sum_phi2 = sum_phi2 + ((1.0 / double(NODE)) * V2[i][t] * 2.0 * PI);

                }

                    R1[t] = sqrt(rx1 * rx1 + ry1 * ry1);
                    R2[t] = sqrt(rx2 * rx2 + ry2 * ry2);
		    Phi1[t] = sum_phi1;
		    Phi2[t] = sum_phi2;

                    if (t == 10000) outfile2 << R1[t] << '\t' << Phi1[t] << '\t' << R2[t] << '\t' << Phi2[t] << endl;

                    rx1 = 0.;
                    ry1 = 0.;
                    rx2 = 0.;
                    ry2 = 0.;
		    sum_phi1 = 0.0;
		    sum_phi2 = 0.0;


                for (i = 0; i < NODE; i++){
                    for (j = 0; j < NODE; j++){
                        if (adj_inter21[i][j] == 1){
                            dummy3 = dummy3 + g_inter21[i][j];
                            counter1++;
                        }

                        if (adj_inter12[i][j] == 1){
                            dummy4 = dummy4 + g_inter12[i][j];
                            counter2++;
                        }
                    }
                }

                    g_mean21[t] = (1. / double(counter1)) * dummy3; 
                    g_mean12[t] = (1. / double(counter2)) * dummy4; 
                    dG = abs(g_mean21[t] - g_mean12[t]);

                    outfile3 << t * dt / 1000. << '\t' << g_mean21[t] << '\t' << g_mean12[t] << '\t' << dG << endl;

                    dummy3 = 0.;
                    dummy4 = 0.;
                    counter1 = 0;
                    counter2 = 0;

// Module 1 *****************************************************************


        for (i = 0; i < NODE; i++){

            xrand = ran2(&iseed1);

            if (xrand <= PP){
                I1[i] = r1[i][0];
                t_ref1[i] = t;
            }

            else I1[i] = r1[i][0] * exp(- ((t - t_ref1[i]) * dt) / tau_r);
            IEXT = I1[i];

            for (j = 0; j < NODE; j++){
                if (t == (tf1[j] + int(tau_intra[i][j] / dt))) dummy1 = g_intra11[i][j] * adj_intra11[i][j] * r1[j][0];
                else dummy1 = 0;

		if (t == (tf2[j] + int(tau_inter[i][j] / dt))) dummy2 = g_inter12[i][j] * adj_inter12[i][j] * r2[j][0];
                else dummy2 = 0;

                dummy11 = dummy11 + dummy1;
                dummy22 = dummy22 + dummy2;
            }

                ISYN = dummy11 + dummy22;

            dummy11 = 0.0;
            dummy22 = 0.0;

            ISTIM1 = 0.0;

            if (t == t_ON1) STIMULATION = true;
            if (t == t_OFF) STIMULATION = false;

            if (STIMULATION == true){
                if (t == (t_ON1 + (m1 * intra_step))){
                    ISTIM1 = 40.0;
                    if (i == NODE - 1) m1++;
                    if (m1 == 5){
                        m1 = 0;
                        t_ON1 = t_ON1 + inter_step;
                    }
                }
            }


            V1[i][t + 1] = V1[i][t] + dt * (1. / tau_m) * (Vrest - V1[i][t] + ISYN + IEXT + ISTIM1);

            ISYN = 0.;

            		if (V1[i][t] > Vth){
                	V1[i][t + 1] = 0.0;
                        A1[i][t] = 1;
                	tf1[i] = t;
                        SPIKE_TRAIN1.push_back(vector<double>());
                        if (0.001 * T < t) SPIKE_TRAIN1[i].push_back(t);

                	NET1++;
                	outfile4 << t * dt / 1000. << '\t' << i + 1 << endl;
                    }

		for (j = 0; j < NODE; j++){

	    tff1[i] = tf1[i] + int(ta_inter[j][i] / dt);
            tbb1[i] = tf1[i] + int(td_inter[j][i] / dt);

            if (t == tff1[i]){
		        if (adj_inter21[j][i] == 1){
                        timelag = abs(tff1[i] - nnb2[j]) * dt;
                        nnf1[i] = tff1[i];
                        if (timelag =! 0) delta_g = A_minus * exp(- timelag / tau_minus);
                        else delta_g = 0.0;
                        g_inter21[j][i] = g_inter21[j][i] - delta_g;
                }
            }

            if (t == tbb1[i]){
		        if (adj_inter12[i][j] == 1){
                        timelag = abs(tbb1[i] - nnf2[j]) * dt;
                        nnb1[i] = tbb1[i];
                        if (timelag =! 0) delta_g = A_plus * exp(- timelag / tau_plus);
                        else delta_g = 0.0;
                        g_inter12[i][j] = g_inter12[i][j] + delta_g;
                }
            }
                    
                    if (g_inter12[i][j] > gmax) g_inter12[i][j] = gmax;
                    if (g_inter21[j][i] > gmax) g_inter21[j][i] = gmax;
                    if (g_inter12[i][j] < gmin) g_inter12[i][j] = gmin;
                    if (g_inter21[j][i] < gmin) g_inter21[j][i] = gmin;
		}
        }

        outfile5 << t * dt / 1000. << '\t' << NET1 / double(NODE) << endl;
        NET1_bin[t] = NET1;
        NET1 = 0;

// Module 2 **************************************************************************************************
  

        for (i = 0; i < NODE; i++){

            xrand = ran2(&iseed2);

            if (xrand <= PP){
                I2[i] = r2[i][0];
                t_ref2[i] = t;
            }

            else I2[i] = r2[i][0] * exp(- ((t - t_ref2[i]) * dt) / tau_r);
            IEXT = I2[i];

            for (j = 0; j < NODE; j++){
                if (t == (tf2[j] + int(tau_intra[i][j] / dt))) dummy1 = g_intra22[i][j] * adj_intra22[i][j] * r2[j][0];
                else dummy1 = 0;

        	if (t == (tf1[j] + int(tau_inter[i][j] / dt))) dummy2 = g_inter21[i][j] * adj_inter21[i][j] * r1[j][0];
                else dummy2 = 0;

                dummy11 = dummy11 + dummy1;
                dummy22 = dummy22 + dummy2;
            }

                ISYN = dummy11 + dummy22;

            dummy11 = 0.0;
            dummy22 = 0.0;

            ISTIM2 = 0.0;

            if (STIMULATION == true){
                    if (t == (t_ON2 + (m2 * intra_step))){
                        ISTIM2 = 40.0;
                        if (i == NODE - 1) m2++;
                        if (m2 == 5){
                            m2 = 0;
                            t_ON2 = t_ON2 + inter_step;
                    }
                }
            }

            V2[i][t + 1] = V2[i][t] + dt * (1. / tau_m) * (Vrest - V2[i][t] + ISYN + IEXT + ISTIM2);

            ISYN = 0.;

            		if (V2[i][t] > Vth){
                	V2[i][t + 1] = 0.0;
                        A2[i][t] = 1;
                	tf2[i] = t;

                        SPIKE_TRAIN2.push_back(vector<double>());
                        if (0.001 * T < t) SPIKE_TRAIN2[i].push_back(t);

                	NET2++;
                	outfile6 << t * dt / 1000. << '\t' << i + 1 << endl;
                    }


		for (j = 0; j < NODE; j++){

	    tff2[i] = tf2[i] + int(ta_inter[j][i] / dt);
            tbb2[i] = tf2[i] + int(td_inter[j][i] / dt);

            if (t == tff2[i]){
		        if (adj_inter12[j][i] == 1){
                        timelag = abs(tff2[i] - nnb1[j]) * dt;
                        nnf2[i] = tff2[i];
                        if (timelag =! 0) delta_g = A_minus * exp(- timelag / tau_minus);
                        else delta_g = 0.0;
                        g_inter12[j][i] = g_inter12[j][i] - delta_g;
                }
            }

            if (t == tbb2[i]){
		        if (adj_inter21[i][j] == 1){
                        timelag = abs(tbb2[i] - nnf1[j]) * dt;
                        nnb2[i] = tbb2[i];
                        if (timelag =! 0) delta_g = A_plus * exp(- timelag / tau_plus);
                        else delta_g = 0.0;
                        g_inter21[i][j] = g_inter21[i][j] + delta_g;
                }
            }

                    if (g_inter21[i][j] > gmax) g_inter21[i][j] = gmax;
                    if (g_inter12[j][i] > gmax) g_inter12[j][i] = gmax;
                    if (g_inter21[i][j] < gmin) g_inter21[i][j] = gmin;
                    if (g_inter12[j][i] < gmin) g_inter12[j][i] = gmin;
		}
        }

        outfile7 << t * dt / 1000. << '\t' << NET2 / double(NODE) << endl;
        NET2_bin[t] = NET2;
        NET2 = 0;
    } // End of time loop


// Final Strengths Distribution: Bin diagrams ***********************************************************************

    for (i = 0; i < RANGE; i++) XX[i] = i * dx + gmin;

    for (i = 0; i < NODE; i++){
        for (j = 0; j < NODE; j++){
		    if (adj_inter12[i][j] != 0){
                c = int(g_inter12[i][j] / dx);
                CC0[c - 1]++;
			    counter4++;
			}
			if (adj_inter21[i][j] != 0){
			    c = int(g_inter21[i][j] / dx);
                CC[c - 1]++;
			    counter3++;
			}
        }
    }

    for(i = 0; i < RANGE; i++) outfile8 << XX[i] << '\t' << CC [i] / double(counter3) << '\t' << CC0[i] / double(counter4) << endl;

// CV *********************************************************************************************

dx1 = 0.05;
c1 = 0;
c2 = 0;
irr_min = 0.0;
irr_max = 1.0;

RANGE1 = ((irr_max - irr_min) / dx1) + 1;
double CC1[RANGE1], CC2[RANGE1], XX1[RANGE1];

 	for (i = 0; i < RANGE1; i++) XX1[i] = i * dx1 + irr_min;

    for (i = 0; i < NODE; i++){ 
        for (j = 0; j < (SPIKE_TRAIN1[i].size()) - 1; j++) ISI1.push_back(SPIKE_TRAIN1[i][j+1] - SPIKE_TRAIN1[i][j]);
        for (j = 0; j < (SPIKE_TRAIN2[i].size()) - 1; j++) ISI2.push_back(SPIKE_TRAIN2[i][j+1] - SPIKE_TRAIN2[i][j]);

		c1 = int(Coefficient_Variation(ISI1, ISI1.size()) / dx1);
        CC1[c1]++;

		c2 = int(Coefficient_Variation(ISI2, ISI2.size()) / dx1);
        CC2[c2]++;

        ISI1.clear();            
        ISI2.clear();            
    }

    for (i = 0; i < RANGE1; i++) outfile9 << XX1[i] << '\t' << CC1[i] / double(NODE) << '\t' << CC2[i] / double(NODE) << endl;


// Population Activity *************************************************************************

int dx2 = 100;
int RANGE2 = (T / dx2) + 1;
int XX2[RANGE2];

double sum_M1 = 0, sum_M2 = 0, ACT1[RANGE2], ACT2[RANGE2];

 	for (i = 0; i < RANGE2; i++){
        XX2[i] = i * dx2;

        for (j = XX2[i]; j < (XX2[i] + dx2); j++){
            sum_M1 = sum_M1 + NET1_bin[j];
            sum_M2 = sum_M2 + NET2_bin[j];
        }

        ACT1[i] = sum_M1;
        ACT2[i] = sum_M2;

        sum_M1 = 0;
        sum_M2 = 0;
    }

    for (i = 0; i < RANGE2; i++) outfile10 << XX2[i] * dt / 1000. << '\t' << ACT1[i] / (double(NODE) * dx2) << '\t' << ACT2[i] / (double(NODE) * dx2) << endl;


// FFT **************************************************************************************************

        for (i = 0; i < T; i++){
            signal1.push_back(NET1_bin[i]);
            signal2.push_back(NET2_bin[i]);
        }

    int NUM_POINTS = signal1.size();

    vector<double> freq1(NUM_POINTS);
    vector<double> freq2(NUM_POINTS);
    vector<double> freqAmplitudes1(NUM_POINTS);
    vector<double> freqAmplitudes2(NUM_POINTS);

    fft(signal1, freq1, freqAmplitudes1, 0.01);
    fft(signal2, freq2, freqAmplitudes2, 0.01);

    for (i = 0; i < NUM_POINTS / 2.0; i++) outfile11 << freq1[i] << '\t' << freqAmplitudes1[i] << '\t' << freq2[i] << '\t' << freqAmplitudes2[i] << endl;


// Pearson Correlation Coefficient *************************************************************************

double dx3 = 0.1, bin = 100;
int n11, n12, n21, n22, c5 = 0, c6 = 0;
int RANGE3 = (2.0 / double(dx3)) + 1, RANGE4 = (T / bin) + 1; 
int CC5[RANGE3], CC6[RANGE3], XX4[RANGE4], A1_bin[NODE][RANGE4], A2_bin[NODE][RANGE4], arr11[RANGE4], arr12[RANGE4], arr21[RANGE4], arr22[RANGE4];
double XX3[RANGE3];

    for (i = 0; i < NODE; i++){
        for (j = 0; j < RANGE4; j++){
            A1_bin[i][j] = 0;
            A2_bin[i][j] = 0;
        }
    }

 	for (i = 0; i < RANGE3; i++) XX3[i] = i * dx3;
 

    for (i = 0; i < RANGE4; i++){
        XX4[i] = i * bin;

        for (k = 0; k < NODE; k++){
            for (j = XX4[i]; j < (XX4[i] + bin); j++){
                if (A1[k][j] == 1) A1_bin[k][i] = 1;
                if (A2[k][j] == 1) A2_bin[k][i] = 1;
            }
        }
    }

    for (i = 0; i < NODE; i++){
        for (k = 0; k < NODE; k++){
            for (j = 0; j < RANGE4; j++){
                if (i != k){
                    arr11[j] = A1_bin[i][j];
                    arr12[j] = A1_bin[k][j];
                    arr21[j] = A2_bin[i][j];
                    arr22[j] = A2_bin[k][j];
                }
            }

        n11 = sizeof(arr11) / sizeof(arr11[0]);
        n11 = sizeof(arr12) / sizeof(arr12[0]);
        n21 = sizeof(arr21) / sizeof(arr21[0]);
        n22 = sizeof(arr22) / sizeof(arr22[0]);

        c5 = Pearson_Correlation_Array(arr11, arr12, n11, n12) / dx3;
        c6 = Pearson_Correlation_Array(arr21, arr22, n21, n22) / dx3;

        CC5[c5 - 1]++;
        CC6[c6 - 1]++;

        }
    }


    for (i = 0; i < RANGE3; i++) outfile12 << XX3[i] - 1.0 << '\t' << CC5[i] / double((NODE)*(NODE-1)) << '\t' << CC6[i] / double((NODE)*(NODE-1)) << endl;


// *********************************************************************************************
// *********************************************************************************************

    outfile0.close();
    outfile1.close();
    outfile2.close();
    outfile3.close();
    outfile4.close();
    outfile5.close();
    outfile6.close();
    outfile7.close();
    outfile8.close();
    outfile9.close();
    outfile10.close();
    outfile11.close();
    outfile12.close();

return 0;

}

// The End *****************************************************

