#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

int binomial_simple(
    double S,                   // Stock Price
    double K,                   // Strike Price
    double r,                   // Interest Rate
    double q,                   // Dividend Yield
    double sigma,               // Stock Volatility
    double T,                   // Time of Expiration 
    double t0,                  // Initial Time
    bool call,                  // True - Call, False - Put
    bool American,              // True - American, False - European
    int n,                      // Number Timesteps
    double &V                   // Option Fair Value
    ) {

	// Validation
	if (n < 1) return 1;
	if (S <= 0) return 1;	
	if (T <= t0) return 1;		
	if (sigma <= 0) return 1;		
	
	// Parameter Calculation
	double dt = (T - t0) / double(n);
	double df = exp(-r * dt);
	double growth = exp((r - q) * dt);
	double u = exp(sigma * sqrt(dt));
	double d = 1.0 / u;
	double p_prob = (growth - d) / (u - d);
	double q_prob = 1.0 - p_prob;
 
	// Additional Validation   
	if (p_prob < 0.0) return 1;		
	if (p_prob > 1.0) return 1;	
	
    // Memory Allocation
	double** stock_nodes = new double*[n + 1];
	double** option_nodes = new double*[n + 1];
	double* S_tmp = NULL;
	double* V_tmp = NULL;
	for (int i = 0 ; i <= n ; i++) {
		stock_nodes[i] = new double[n + 1];
		option_nodes[i] = new double[n + 1];
		S_tmp = stock_nodes[i];
		V_tmp = option_nodes[i];
		for (int j = 0; j <= n ; j++) {
			S_tmp[j] = 0;
			V_tmp[j] = 0;
		}
	}
	S_tmp = stock_nodes[0];
	S_tmp[0] = S;

	// Indexing the Nodes (For Stock Price)
	for (int i = 1 ; i <= n ; i++) {
		double* prev = stock_nodes[i - 1];
		S_tmp = stock_nodes[i];
		S_tmp[0] = prev[0] * d;
		for (int j = 1; j <= n ; j++) {
			S_tmp[j] = S_tmp[j - 1] * u * u;
		}
	}

    // Setting the Terminal Payoffs
    for (int i = 0 ; i <= n ; i++) {
        S_tmp = stock_nodes[i];
        V_tmp = option_nodes[i];
        for (int j = 0 ; j <= n ; j++) {
            double intrinsic = 0;
            if (call && S_tmp[j] > K) intrinsic = S_tmp[j] - K;
            if (!call && K > S_tmp[j]) intrinsic = K - S_tmp[j];
            V_tmp[j] = intrinsic;
        }
    }

    // Main Valuation
    for (int i = n-1 ; i >= 0 ; i--) {
        S_tmp = stock_nodes[i];
        V_tmp = option_nodes[i];
        double* V_next = option_nodes[i+1];
        for (int j = 0 ; j <= i ; j++) {
            double tmp = df * (p_prob * V_next[j + 1] + q_prob * V_next[j]);
            if (American) {
                V_tmp[j] = max(V_tmp[j] , tmp);
            }
            else V_tmp[j] = tmp;
        }
    }

    // Set V - Fair Value
    int i = 0;
    V_tmp = option_nodes[i];
    V = V_tmp[0];

    // Memory Deallocation
    for (int k = 0; k <= n; k++) {
    	delete [] stock_nodes[k];
    	delete [] option_nodes[k];
    }
    delete [] stock_nodes;
    delete [] option_nodes;

    return 0;
}

int main() {
    // Test Cases
    double S = 100.0;
    double K = 100.0;
    double r = 0.1;
    double q = 0.0;
    double sigma = 0.5;
    double T = 0.3;
    double t0 = 0.0;
    double V;
    int n = 3;
    bool call;
    bool American;

    // Output File
    ofstream myfile("result.txt");

    // European Call
    call = true;
    American = false;
	binomial_simple(S, K, r, q, sigma, T, t0, call, American, n, V);
    myfile << " European Call " << V << endl;

    // European Put
    call = false;
    American = false;
	binomial_simple(S, K, r, q, sigma, T, t0, call, American, n, V);
    myfile << " European Put " << V << endl;

    // American Call
    call = true;
    American = true;
	binomial_simple(S, K, r, q, sigma, T, t0, call, American, n, V);
    myfile << " American Call " << V << endl;

    // American Put
    call = false;
    American = true;
	binomial_simple(S, K, r, q, sigma, T, t0, call, American, n, V);
    myfile << " American Put " << V << endl;

    myfile.close();
	return 0;
}

/* COMMENTED OUT FUNCTIONS RELATING TO VISUALIZATION & DEBUGGING */

// void print_tree(double** nodes, int steps, ofstream &myfile) {
//     // Fillers for the parts without the nodes
//     string tabs = " x    ";

//     // Parameters
//     int rows = steps;
//     int cols = steps;

//     // Set Precision
//     myfile << setprecision(2) << fixed;

//     // Print Upper Half & Mid
//     for (int i = steps ; i >= 0 ; i--) {
//         int temp = i;
//         for (int j = 0 ; j < temp ; j++) {
//             myfile << tabs;
//         }
//         if (nodes[rows][cols] < 99.99) myfile << "0";
//         if (nodes[rows][cols] < 9.99) myfile << "0";
//         myfile << nodes[rows][cols];
//         int temp2 = i;
//         while (steps - temp2 >= 2) {
//             int rows_temp = rows + 2;
//             int cols_temp = cols + 1;
//             myfile << tabs;
//             if (nodes[rows_temp][cols_temp] < 99.99) myfile << "0";
//             if (nodes[rows_temp][cols_temp] < 9.99) myfile << "0";
//             myfile << nodes[rows_temp][cols_temp];
//             temp2 = temp2 + 2;
//         }
//         rows--;
//         cols--;
//         myfile << endl;
//     }

//     // Reset Parameters
//     rows = 1;
//     cols = 0;

//     // Print Lower Half
//     for (int i = 1 ; i <= steps ; i++) {
//         int temp = i;
//         for (int j = 0 ; j < temp ; j++) {
//             myfile << tabs;
//         }
//         if (nodes[rows][cols] < 99.99) myfile << "0";
//         if (nodes[rows][cols] < 9.99) myfile << "0";
//         myfile << nodes[rows][cols];
//         int temp2 = steps;
//         while (temp2 - i >= 2) {
//             int rows_temp = rows + 2;
//             int cols_temp = cols + 1;
//             myfile << tabs;
//             if (nodes[rows_temp][cols_temp] < 99.99) myfile << "0";
//             if (nodes[rows_temp][cols_temp] < 9.99) myfile << "0";
//             myfile << nodes[rows_temp][cols_temp];
//             temp2 = temp2 - 2;
//         }
//         rows++;
//         myfile << endl;
//     }   
// }

// // Print Stock Nodes
// myfile << " --- STOCK PRICE --- " << endl;
// print_tree(stock_nodes, n, myfile);
// myfile << endl;

// // Print Terminal Payoff
// myfile << " --- TERMINAL PAYOFF --- " << endl;
// print_tree(option_nodes, n, myfile);
// myfile << endl;

// // Print Main Valuation
// myfile << " --- MAIN VALUATION --- " << endl;
// print_tree(option_nodes, n, myfile);
// myfile << endl;