#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

class Derivative {
    public:
        // Destructor
        virtual ~Derivative() {}

        // Functions - Override via Option Class
        virtual double TerminalPayoff(double S) const {return 0;}
        virtual int ValuationTests(double S, double &V) const {return 0;}

        // Class Variables
        double r;
        double q;
        double sigma;
        double T;
    
    protected:
        // Default Constructor
        Derivative() {r = 0; q = 0; sigma = 0; T = 0;}
};

class BinaryOption : public Derivative {
    public:
        // Default Constructor & Desctructor
        BinaryOption() {K = 0; isCall = false;}
        virtual ~BinaryOption() {}

        // Functions - Definition to Follow
        virtual double TerminalPayoff(double S) const;

        // Class member variables
        double K;
        bool isCall;
};

double BinaryOption::TerminalPayoff(double S) const {
    // Default payout value
    double payout = 0;

    // For either call or put, if conditions are met, pay $1
    if (isCall && S >= K) payout = 1;
    if (!isCall && K > S) payout = 1;

    // Return the calculated value (if condition is met)
    return payout;
}

class Option : public Derivative {
    public:
        // Default Constructor & Destructor
        Option() {K = 0; isCall = false; isAmerican = false;}
        virtual ~Option() {}

        // Functions - Definition to Follow
        virtual double TerminalPayoff(double S) const;
        virtual int ValuationTests(double S, double &V) const;

        // Class member variables
        double K;
        bool isCall;
        bool isAmerican;
};

double Option::TerminalPayoff(double S) const {
    // Default intrinsic value
    double intrinsic = 0;

    // For either call or put, if conditions are met, set intrinsic
    if (isCall && S > K) intrinsic = S - K;
    if (!isCall && K > S) intrinsic = K - S;

    // Return the calculated value (if condition is met)
    return intrinsic;
}

int Option::ValuationTests(double S, double &V) const {
    // Default intrinsic value
    double temp = 0;

    // Perform the valuation test only for American options
    if (isAmerican) {

        // For either call or put, if conditions are met, set intrinsic
        if (isCall && S > K) temp = S - K;
        if (!isCall && K > S) temp = K - S;

        // Set the value to max of the existing node or the intrinsic value
        V = max(V, temp);
    }

    // Return 0 if successfully done
    return 0;
}

class Straddle : public Derivative {
    public:
        // Default Constructor & Destructor
        Straddle() {K = 0; isAmerican = false;}
        virtual ~Straddle() {}

        // Functions - Definition to Follow
        virtual double TerminalPayoff(double S) const;
        virtual int ValuationTests(double S, double &V) const;

        // Class member variables
        double K;
        bool isAmerican;
};

double Straddle::TerminalPayoff(double S) const {
    // Initialize
    double intrinsic = 0;

    // Intrinsic value will always be absolute value of S - K
    intrinsic = fabs(S - K);

    // Return the calculated value
    return intrinsic;
}

int Straddle::ValuationTests(double S, double &V) const {
    // Initialize
    double temp = 0;

    // Perform the valuation test for American straddle
    if (isAmerican) {

        // Intrinsic value will always be absolute value of S - K
        temp = fabs(S - K);

        // Set the higher one of the exisiting value or the intrinsic value to V
        V = max(V, temp);
    }

    // Return 0 if successfully done
    return 0;
}

class UpOutBarrierCallOption : public Derivative {
    public:
        // Default Constructor & Destructor
        UpOutBarrierCallOption() {K = 0; isAmerican = false;}

        // Functions - Definition to Follow
        virtual double TerminalPayoff(double S) const;
        virtual int ValuationTests(double S, double &V) const;

        // Class member variables
        double K;
        double B;
        bool isAmerican;
};

double UpOutBarrierCallOption::TerminalPayoff(double S) const {
    // Initialize
    double rebate = 0;

    // Rebate calculation
    if (S < K) rebate = 0;
    if (K <= S && S < B) rebate = S - K;
    if (S >= B) rebate = B - K;

    // Return the calculated rebate
    return rebate;
}

int UpOutBarrierCallOption::ValuationTests(double S, double &V) const {
    // Valuation testing
    if (S > B) V = B - K;

    // Additional case for American
    if (isAmerican) {
        if (B > S && S > K && S - K > V) {
            V = S - K;
        }
    }

    // Return 0 on success
    return 0;
}

class BinomialModel {
    public:
        // Constructor & Destructor
        BinomialModel(int n);
        ~BinomialModel(); 

        // Functions - Definition to Follow
        int FairValue(int n, Derivative* p_derivative, double S, double t0, double &V);
        int ImpliedVolatility(int n, Derivative* p_derivative, double S, double t0, double target, double &implied_vol, int &num_iter);

    private:
        // Functions - Definition to Follow
        void Clear();
        int Allocate(int n);
        int ImpliedVolatilityPrivate(int n, Derivative* p_derivative, double S, double t0, double target, double &implied_vol, int &num_iter);

        // Class Variables
        int n_tree;
        double** stock_nodes;
        double** value_nodes;
};

int BinomialModel::Allocate(int n) {
    // Validation
    if (n <= n_tree) return 0;

    // Clear existing memory, if any
    Clear();

    // Assign number steps to the BM class variable
    n_tree = n;

    // Allocate the memory
    stock_nodes = new double*[n + 1];
    value_nodes = new double*[n + 1];
    double* S_tmp = NULL;
    double* V_tmp = NULL;
    for (int i = 0 ; i <= n ; i++) {
        stock_nodes[i] = new double[n + 1];
        value_nodes[i] = new double[n + 1];
        S_tmp = stock_nodes[i];
        V_tmp = value_nodes[i];
        for (int j = 0 ; j <= n ; j++) {
            S_tmp[j] = 0;
            V_tmp[j] = 0;
        }
    }

    // Successful allocation - Return 0
    return 0;
}

BinomialModel::BinomialModel(int n) {
    // Constructor, sets n_tree to 0 and pointers to NULL
    // Invokes Allocate() - Definition to Follow
    n_tree = 0;
    stock_nodes = NULL;
    value_nodes = NULL;
    Allocate(n);
}

BinomialModel::~BinomialModel() {
    // Destructor
    // Invokes Clear() - Definition to Follow
    Clear();
}

void BinomialModel::Clear() {
    // Validation
    if (n_tree <= 0) return;

    // Delete inner array pointers
    for (int i = 0 ; i <= n_tree ; i++) {
        delete[] stock_nodes[i];
        delete[] value_nodes[i];
    }

    // Delete the array pointers
    delete[] stock_nodes;
    delete[] value_nodes;

    // Clean up dangling pointers
    stock_nodes = NULL;
    value_nodes = NULL;
}

int BinomialModel::FairValue(int n, Derivative* p_derivative, double S, double t0, double &V) {
    // Initialize fair value variable
    V = 0.0;

    // Validation
    if (n < 1 || S <= 0) return 1;
    if (p_derivative == NULL || p_derivative -> T <= t0 || p_derivative -> sigma <= 0) return 1;

    // Declare the local pointers
    double* S_tmp = NULL;
    double* V_tmp = NULL;

    // Calculate the parameters
    double dt = (p_derivative -> T - t0) / double(n);
    double df = exp(-1 * p_derivative -> r * dt);
    double growth = exp((p_derivative -> r - p_derivative -> q) * dt);
    double u = exp(p_derivative -> sigma * sqrt(dt));
    double d = 1.0 / u;
    double p_prob = (growth - d) / (u - d);
    double q_prob = 1.0 - p_prob;

    // Validation
    if (p_prob < 0.0 || p_prob > 1.0) return 1;

    // Allocate memory
    Allocate(n);
    S_tmp = stock_nodes[0];
    S_tmp[0] = S;

    // Set stock prices
    for (int i = 1 ; i <= n ; i++) {
        double* prev = stock_nodes[i - 1];
        S_tmp = stock_nodes[i];
        S_tmp[0] = prev[0] * d;
        for (int j = 1 ; j <= n ; j++) {
            S_tmp[j] = S_tmp[j - 1] * u * u;
        }
    }

    // Set terminal payoff - on the outermost nodes
    int i = n;
    S_tmp = stock_nodes[i];
    V_tmp = value_nodes[i];
    for (int j = 0 ; j <= n ; j++) {
        V_tmp[j] = p_derivative -> TerminalPayoff(S_tmp[j]);
    }

    // Main valuation
    for (int i = n - 1 ; i >= 0 ; i--) {
        S_tmp = stock_nodes[i];
        V_tmp = value_nodes[i];
        double* V_next = value_nodes[i + 1];
        for (int j = 0 ; j <= i ; j++) {
            V_tmp[j] = df * (p_prob * V_next[j + 1] + q_prob * V_next[j]);
            p_derivative -> ValuationTests(S_tmp[j], V_tmp[j]);
        }
    }

    // Set the fair value
    V_tmp = value_nodes[0];
    V = V_tmp[0];

    // Explicit de-allocation of memory, if needed
    // Clear();

    // If successful operation until here, return 0
    return 0;
}

int BinomialModel::ImpliedVolatility(int n, Derivative* p_derivative, double S, double t0, double target, double &implied_vol, int &num_iter) {
    // Test Case
    int rc = 0;

    // Container for the original volatility
    const double saved_vol = p_derivative -> sigma;

    // Run through the private function
    rc = ImpliedVolatilityPrivate(n, p_derivative, S, t0, target, implied_vol, num_iter);
    p_derivative -> sigma = saved_vol;

    // Return the result
    return rc;
}

int BinomialModel::ImpliedVolatilityPrivate(int n, Derivative* p_derivative, double S, double t0, double target, double &implied_vol, int &num_iter) {
    // Constants
    const double tol = 1.0e-4;
    const int max_iter = 100;

    // Initialize Referenced Variables
    implied_vol = 0;
    num_iter = 0;

    // Local Variables
    double sigma_low = 0.01;
    double sigma_high = 3.0;
    double FV_low = 0;
    double FV_high = 0;
    double FV = 0;

    // Validation for Sigma Low (Success)
    p_derivative -> sigma = sigma_low;
    FairValue(n, p_derivative, S, t0, FV_low);
    double diff_FV_low = FV_low - target;
    if (fabs(diff_FV_low) <= tol) {
        implied_vol = p_derivative -> sigma;
        return 0;
    }

    // Validation for Sigma High (Success)
    p_derivative -> sigma = sigma_high;
    FairValue(n, p_derivative, S, t0, FV_high);
    double diff_FV_high = FV_high - target;
    if (fabs(diff_FV_high) <= tol) {
        implied_vol = p_derivative -> sigma;
        return 0;
    }

    // Validation for Out of Reach (Fail)
    if (diff_FV_low * diff_FV_high > 0) return 1;

    // Iteration (Success or Loop Over)
    for (int i = 0 ; i < max_iter ; i++) {
        p_derivative -> sigma = 0.5 * (sigma_low + sigma_high);
        FairValue(n, p_derivative, S, t0, FV);
        double diff_FV = FV - target;
        if (fabs(diff_FV) <= tol) {
            implied_vol = p_derivative -> sigma;
            num_iter = i;
            return 0;
        }
        else if ((diff_FV_low * diff_FV) > 0)
            sigma_low = p_derivative -> sigma;
        else
            sigma_high = p_derivative -> sigma;
    }

    // If Iteration Fails (Fail)
    num_iter = max_iter;
    implied_vol = 0;
    return 1;
}

void price_from_yield (double F, int ID, double y, double &B) {

    /*
     *  Coupon rate is replaced as it will depend on the student ID
     *  Number of steps also replaced as it will always be 8
     *  Arguments:
     *  F - Face Value
     *  ID - Student ID
     *  y - Yield
     *  B - Bond Fair Value
     */
    
    // Convert ID into an integer array
    int coupons[8];
    int id = ID;
    for (int i = 7; i >= 0; i--) {
        coupons[i] = id % 10;
        id /= 10;
    }

    // Bond FV Calculation
    y = y * 0.01;
    B = 0;
    for (int i = 0 ; i < 7 ; i++) {
        B += ((double)coupons[i] / 2) / pow((1 + (y / 2)), i + 1);
    }
    B += (F + ((double)coupons[7] / 2)) / pow((1 + (y / 2)), 8);

}

int yield_from_price (double F, int ID, double B, double tol, int max_iter, double &y) {

    /*
     *  Coupon rate is replced as it will depend on the student ID
     *  Number of steps also replaced as will will be 8
     *  Arguments:
     *  F - Face Value
     *  ID - Student ID
     *  B - Bond Market Price
     *  tol - Tolerance
     *  max_iter - Max Iterations
     *  y - Yield
     */

    // Convert ID into an integer array
    int coupons[8];
    int id = ID;
    for (int i = 7; i >= 0; i--) {
        coupons[i] = id % 10;
        id /= 10;
    }

    // Local Parameters
    double y_low = 0; double y_high = 100;
    double b_low = 0; double b_high = 0; double b_mid = 0;

    // Obtain Bond Fair Value (Lower Bound)
    price_from_yield(F, ID, y_low, b_low);
    if (fabs(b_low - B) <= tol) {
        y = y_low;
        return 0;
    }
    if (b_low < B) {
        y = 0;
        return 1;
    }

    // Obtain Bond Fair Value (Upper Bound)
    price_from_yield(F, ID, y_high, b_high);
    if (fabs(b_high - B) <= tol) {
        y = y_high;
        return 0;
    }
    if (b_high > B) {
        y = 0;
        return 1;
    }

    // Iteration Using Bisection
    for (int i = 0 ; i < max_iter ; i++) {
        double y_mid = (y_low + y_high) / 2;
        price_from_yield(F, ID, y_mid, b_mid);
        if (fabs(b_mid - B) <= tol) {
            y = y_mid;
            return 0;
        }
        else if (b_mid > B) {
            y_low = y_mid;
        }
        else {
            y_high = y_mid;
        }
        if ((y_high - y_low) <= tol) {
            y = y_mid;
            return 0;
        }
    }

    // If Iteration Exits (Failed)
    return 1;
}

int main() {
    // Test Result
    int rc = 0;

    // Set precision to 5 decimal places
    cout << setprecision(5) << fixed;

    // Fair value container
    double U = 0.0;
    double U_left = 0.0;
    double U_right = 0.0;

    // Params
    const double K = 100.0;           // Same for all 3 days
    const double B = 130.0;           // Same for all 3 days
    const double r = 0.0;             // Same for all 3 days
    const double q = 0.0;             // Same for all 3 days
    const double T = 0.5;             // Same for all 3 days
    const bool isAmerican = true;     // Same for all 3 days

    // Fluctuations
    double ds1 = 1224 / 1e4;    // First 4 digits of student id
    double ds2 = -1635 / 1e4;   // Last 4 digits of student id

    // Binomial Model
    int n = 1000;               // 1000 steps
    BinomialModel binom(n);

    // UpOutBarrierCallOption Object & set the constants
    UpOutBarrierCallOption Am_call;
    Am_call.K = K;
    Am_call.B = B;
    Am_call.r = r;
    Am_call.q = q;
    Am_call.T = T;
    Am_call.isAmerican = isAmerican;

    // Day 0 Params
    double S0 = 90.0;
    double t0_0 = 0.0;
    double M0 = 5.0;

    // Day 0 Results
    double sig_0 = 0.0;
    int sig_0_iter = 0;
    double Delta_0 = 0.0;
    double money_0 = 0.0;

    // Day 0 Calculations for implied volatility
    rc = binom.ImpliedVolatility(n, &Am_call, S0, t0_0, M0, sig_0, sig_0_iter);
    Am_call.sigma = sig_0;
    cout << "Day 0 Volatility: " << sig_0 << endl;

    // Day 0 Calculation for fair value using above volatility
    rc = binom.FairValue(n, &Am_call, S0, t0_0, U);
    cout << "Day 0 Option Fair Value: " << U << endl;

    // Day 0 Calculation for Delta_0
    rc = binom.FairValue(n, &Am_call, S0 + 1, t0_0, U_left);
    rc = binom.FairValue(n, &Am_call, S0 - 1, t0_0, U_right);
    Delta_0 = (U_left - U_right) / 2;
    cout << "Day 0 Delta: " << Delta_0 << endl;

    // Day 0 Calculation for money_0
    money_0 = (Delta_0 * S0) - M0;
    cout << "Day 0 Money: " << money_0 << endl;
    cout << endl;

    // Day 1 Params
    double S1 = S0 + ds1;
    double t0_1 = 0.01;
    double M1 = 5.2;

    // Day 1 Results
    double sig_1 = 0.0;
    int sig_1_iter = 0;
    double Delta_1 = 0.0;
    double money_1 = 0.0;

    // Day 1 Calculations for implied volatility
    rc = binom.ImpliedVolatility(n, &Am_call, S1, t0_1, M1, sig_1, sig_1_iter);
    Am_call.sigma = sig_1;
    cout << "Day 1 Volatility: " << sig_1 << endl;

    // Day 1 Calculation for fair value using above volatility
    rc = binom.FairValue(n, &Am_call, S1, t0_1, U);
    cout << "Day 1 Option Fair Value: " << U << endl;

    // Day 1 Calculation for Delta_1
    rc = binom.FairValue(n, &Am_call, S1 + 1, t0_1, U_left);
    rc = binom.FairValue(n, &Am_call, S1 - 1, t0_1, U_right);
    Delta_1 = (U_left - U_right) / 2;
    cout << "Day 1 Delta: " << Delta_1 << endl;

    // Day 1 Calculation for money_1
    money_1 = money_0 + ((Delta_1 - Delta_0) * S1);
    cout << "Day 1 Money: " << money_1 << endl;
    cout << endl;

    // Day 2 Params
    double S2 = S1 + ds2;
    double t0_2 = 0.02;
    double M2 = 5.15;

    // Day 2 Results
    double sig_2 = 0.0;
    int sig_2_iter = 0;
    double Delta_2 = 0.0;
    double money_2 = 0.0;

    // Day 2 Calculations for implied volatility
    rc = binom.ImpliedVolatility(n, &Am_call, S2, t0_2, M2, sig_2, sig_2_iter);
    Am_call.sigma = sig_2;
    cout << "Day 2 Volatility: " << sig_2 << endl;

    // Day 2 Calculation for fair value using above volatility
    rc = binom.FairValue(n, &Am_call, S2, t0_2, U);
    cout << "Day 2 Option Fair Value: " << U << endl;

    // Day 2 Calculation for Delta_2
    rc = binom.FairValue(n, &Am_call, S2 + 1, t0_2, U_left);
    rc = binom.FairValue(n, &Am_call, S2 - 1, t0_2, U_right);
    Delta_2 = (U_left - U_right) / 2;
    cout << "Day 2 Delta: " << Delta_2 << endl;

    // Day 2 Calculation for money_2
    money_2 = money_1 + ((Delta_2 - Delta_1) * S2);
    cout << "Day 2 Money: " << money_2 << endl;
    cout << endl;

    // Closing Profit Params
    double profit = 0.0;

    // Closing Profit Calculation
    profit = money_2 + M2 - (Delta_2 * S2);
    cout << "Closing Profit: " << profit;

    return 0;
}