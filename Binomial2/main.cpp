#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

class Derivative {
    public:
        virtual ~Derivative() {}

        virtual double TerminalPayoff(double S) const {return 0;}
        virtual int ValuationTests(double S , double &V) const {return 0;}

        double r;
        double q;
        double sigma;
        double T;

    protected:
        Derivative() {r = 0; q = 0; sigma = 0; T = 0;}
};

class Option : public Derivative {
    public:
        Option() {K = 0; isCall = false; isAmerican = false;}
        virtual ~Option() {}

        virtual double TerminalPayoff(double S) const;
        virtual int ValuationTests(double S , double &V) const;

        double K;
        bool isCall;
        bool isAmerican;
};

double Option::TerminalPayoff(double S) const {
    double intrinsic = 0;
    if (isCall && S > K) intrinsic = S - K;
    if (!isCall && K > S) intrinsic = K - S;
    return intrinsic;
}

int Option::ValuationTests(double S, double &V) const {
    if (isAmerican) {
        double temp = 0;
        if (isCall && S > K) temp = S - K;
        if (!isCall && K > S) temp = K - S;
        V = max(V, temp);
    } 
    return 0;
}

class BinomialModel {
    public:
        BinomialModel(int n);
        ~BinomialModel();

        int FairValue(int n , Derivative* p_derivative, double S, double t0, double &V);

    private:
        void Clear();
        int Allocate(int n);

        int n_tree;
        double** stock_nodes;
        double** value_nodes;
};

BinomialModel::BinomialModel(int n) {
    n_tree = 0;
    stock_nodes = NULL;
    value_nodes = NULL;
    Allocate(n);
}

BinomialModel::~BinomialModel() {
    Clear();
}

void BinomialModel::Clear() {
    if (n_tree <= 0) return;
    for (int i = 0 ; i <= n_tree ; i++) {
		delete[] stock_nodes[i];
		delete[] value_nodes[i];
    }
    delete [] stock_nodes;
    delete [] value_nodes;
    stock_nodes = NULL;
    value_nodes = NULL;
}

int BinomialModel::Allocate(int n) {
    if (n <= n_tree) return 0;

    Clear();

    n_tree = n;
    
    stock_nodes = new double*[n + 1];
	value_nodes = new double*[n + 1];
	double* S_tmp = NULL;
	double* V_tmp = NULL;
	for (int i = 0 ; i <= n ; i++) {
		stock_nodes[i] = new double[n + 1];
		value_nodes[i] = new double[n + 1];
		S_tmp = stock_nodes[i];
		V_tmp = value_nodes[i];
		for (int j = 0; j <= n ; j++) {
			S_tmp[j] = 0;
			V_tmp[j] = 0;
		}
	}
    return 0;
}

int BinomialModel::FairValue(int n, Derivative* p_derivative, double S, double t0, double &V) {
    V = 0.0;

    // Validation
    if (n < 1) return 1;
    if (S <= 0) return 1;
    if (p_derivative == NULL) return 1;
    if (p_derivative -> T <= t0) return 1;
    if (p_derivative -> sigma <= 0) return 1;

    // Local Variable Declaration
	double* S_tmp = NULL;
	double* V_tmp = NULL;

    // Calculate Parameters
    double dt = (p_derivative -> T - t0) / double(n);
    double df = exp(-1 * p_derivative -> r * dt);
    double growth = exp((p_derivative -> r - p_derivative -> q) * dt);
    double u = exp(p_derivative -> sigma * sqrt(dt));
    double d = 1.0 / u;
    double p_prob = (growth - d) / (u - d);
    double q_prob = 1.0 - p_prob;

    // Additional Validation
    if (p_prob < 0.0) return 1;
    if (p_prob > 1.0) return 1;

    // Allocate Memory
    Allocate(n);
    S_tmp = stock_nodes[0];
    S_tmp[0] = S;

    // Set Stock Prices
    for (int i = 1 ; i <= n ; i++) {
		double* prev = stock_nodes[i - 1];
		S_tmp = stock_nodes[i];
		S_tmp[0] = prev[0] * d;
		for (int j = 1; j <= n ; j++) {
			S_tmp[j] = S_tmp[j - 1] * u * u;
		}
	}

    // Set Terminal Payoff
    int i = n;
    S_tmp = stock_nodes[i];
    V_tmp = value_nodes[i];
    for (int j = 0 ; j <= n ; j++) {
        V_tmp[j] = p_derivative -> TerminalPayoff(S_tmp[j]);
    }

    // Main Valuation
    for (int i = n - 1 ; i >= 0 ; i--) {
        S_tmp = stock_nodes[i];
        V_tmp = value_nodes[i];
        double* V_next = value_nodes[i + 1];
        for (int j = 0 ; j <= i ; j++) {
            V_tmp[j] = df * (p_prob * V_next[j + 1] + q_prob * V_next[j]);
            p_derivative -> ValuationTests(S_tmp[j], V_tmp[j]);
        }
    }

    // Set Fair Value
    V_tmp = value_nodes[0];
    V = V_tmp[0];

    // Deallocate
    // Clear();

    return 0;
}

int main() {

    // Test Cases
    int rc = 0;

    // Output File
    ofstream myfile("result.txt");

    // Test Cases
    double S = 100.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.01;
    double sigma = 0.5;
    double T = 1.0;
    double t0 = 0.0;

    // Class - European Put
    Option Eur_put;
    Eur_put.r = r;
    Eur_put.q = q;
    Eur_put.sigma = sigma;
    Eur_put.T = T;
    Eur_put.K = K;
    Eur_put.isCall = false;
    Eur_put.isAmerican = false;

    // Class - American Put
    Option Am_put;
    Am_put.r = r;
    Am_put.q = q;
    Am_put.sigma = sigma;
    Am_put.T = T;
    Am_put.K = K;
    Am_put.isCall = false;
    Am_put.isAmerican = true;

    // Class - European Call
    Option Eur_call;
    Eur_call.r = r;
    Eur_call.q = q;
    Eur_call.sigma = sigma;
    Eur_call.T = T;
    Eur_call.K = K;
    Eur_call.isCall = true;
    Eur_call.isAmerican = false;

    // Class - American Call
    Option Am_call;
    Am_call.r = r;
    Am_call.q = q;
    Am_call.sigma = sigma;
    Am_call.T = T;
    Am_call.K = K;
    Am_call.isCall = true;
    Am_call.isAmerican = true;

    // Fair Values
    double FV_Am_put = 0;
    double FV_Eur_put = 0;
    double FV_Am_call = 0;
    double FV_Eur_call = 0;

    // Steps
    int n = 100;

    // Class for Binomial Model
    BinomialModel binom(n);

    // Test Variables
    double dS = 0.1;
    int imax = 2000;
    int i;
    for (i = 1 ; i <= imax ; i++) {
        S = i * dS;
        rc = binom.FairValue(n, &Am_put, S, t0, FV_Am_put);
        rc = binom.FairValue(n, &Eur_put, S, t0, FV_Eur_put);
        rc = binom.FairValue(n, &Am_call, S, t0, FV_Am_call);
        rc = binom.FairValue(n, &Eur_call, S, t0, FV_Eur_call);

        myfile << setw(16) << S << " ";
        myfile << setw(16) << FV_Am_put << " ";
        myfile << setw(16) << FV_Eur_put << " ";
        myfile << setw(16) << FV_Am_call << " ";
        myfile << setw(16) << FV_Eur_call << " ";
        myfile << endl;
    }

    myfile.close();
	return 0;
}
