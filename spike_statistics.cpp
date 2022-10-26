
//**********************************************************************************************
//**********************************************************************************************

// Code developed by Mojtaba Madadi Asl.
// Email: mojtabamadadi7@gmail.com
// Necessary functions for spike count statistics.

//**********************************************************************************************
//**********************************************************************************************

#include <bits/stdc++.h>
#include <iostream>
#include <vector> 
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <array>
#include <assert.h>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <fftw3.h>
#include <cstdio>
#include <assert.h>
#include <cmath>

#define REAL 0
#define IMAG 1

using namespace std; 


//**********************************************************************************************
//**********************************************************************************************

// Delta Function 

double Delta_Function(int x, double dt){
    
	double flag;

	if (x == 0) flag = 1. / dt;
	else flag = 0;               
	return(flag);
}


//**********************************************************************************************
//**********************************************************************************************

// Mean (Vector)

double Mean(vector<double> &vec, int n){ 
	
	double sum = 0.0; 

	for (int i = 0; i < n; i++) sum = sum + vec[i]; 
	return sum / n; 
} 

// Mean (Array) 

double Mean_Array(int arr[], int n){ 

	double sum = 0.0; 

	for (int i = 0; i < n; i++) sum = sum + arr[i]; 
	return sum / n; 
} 


//**********************************************************************************************
//**********************************************************************************************

// Standard Deviation (Vector)

double Standard_Deviation(vector<double> &vec, int n){ 
	
	double sum = 0.0;

	for (int i = 0; i < n; i++) sum = sum + (vec[i] - Mean(vec, n)) * (vec[i] - Mean(vec, n)); 
	return sqrt(sum / (n - 1)); 
} 

// Standard Deviation (Array)

double Standard_Deviation_Array(int arr[], int n){ 

	double sum = 0.0;

	for (int i = 0; i < n; i++) sum = sum + (arr[i] - Mean_Array(arr, n)) * (arr[i] - Mean_Array(arr, n)); 
	return sqrt(sum / (n - 1)); 
} 


//**********************************************************************************************
//**********************************************************************************************

// Coefficient of Variation (Vector)

double Coefficient_Variation(vector<double> &vec, int n){ 

	return Standard_Deviation(vec, n) / Mean(vec, n); 
} 

// Coefficient of Variation (Array)

double Coefficient_Variation_Array(int arr[], int n){ 

	return Standard_Deviation_Array(arr, n) / Mean_Array(arr, n); 

} 

//**********************************************************************************************
//**********************************************************************************************

// Pearson Correlation Coefficient (Array)

double Pearson_Correlation_Array(int arr1[], int arr2[], int n1, double n2){ 

	double sum = 0.0;
    int n;

    if (n1 != n2){
        n = 0;
        cout << "Error: Arrays dimentions do not match!" << endl;
    }

    else n = n1;

	for (int i = 0; i < n; i++) sum = sum + ((arr1[i] - Mean_Array(arr1, n1)) * (arr2[i] - Mean_Array(arr2, n2))); 

	return sum / (Standard_Deviation_Array(arr1, n1) * Standard_Deviation_Array(arr2, n2)); 
} 

//**********************************************************************************************
//**********************************************************************************************

// Fast Fourier Transform (by using FFTW3)

void fft(vector<double> &signal, vector<double> &freq, vector<double> &freqAmplitudes, double dt)
{
    const int NUM_POINTS = signal.size();
    const double fs = 1.0 / (dt + 0.0);
    double *signalArray = new double[NUM_POINTS];
    unsigned flags{0};

    fftw_complex result[NUM_POINTS / 2 + 1];
    fftw_plan plan = fftw_plan_dft_r2c_1d(NUM_POINTS,
                                          signalArray,
                                          result,
                                          flags);
    for (int i = 0; i < NUM_POINTS; i++)
        signalArray[i] = signal[i];

    fftw_execute(plan);
    for (int i = 0; i < (0.5 * NUM_POINTS); i++)
    {
        freqAmplitudes[i] = 2.0 / (double)(NUM_POINTS) * sqrt(result[i][REAL] * result[i][REAL] + result[i][IMAG] * result[i][IMAG]);
        freq[i] = i / double(NUM_POINTS) * fs;
    }

    fftw_destroy_plan(plan);
    delete[] signalArray;
}

//**********************************************************************************************
//**********************************************************************************************

