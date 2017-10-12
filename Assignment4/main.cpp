#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;

// Name : Eric Kwon
// Assignment 2
// Due : October 13, 2017

// 4.1 Forward Price
void forward_price (double s0, double t0, double r, double q, double tf, double &fv) {
    double r_dec = r * 0.01;
    double q_dec = q * 0.01;
    double param = (r_dec - q_dec) * (tf - t0);
    fv = s0 * exp(param);
}

// 4.2 Forwards: Arbitrage
void forward_profit (double fv, double ffv, double &profit, string &strategy) {
    if (fv > ffv) {
        strategy = "(Short Sell) ";
        profit = fv - ffv;
    }
    else {
        strategy = "(Buy Low Sell High) ";
        profit = ffv - fv;
    }
}

// 4.3 Futures: Mark to Market
void random_walk (double si[], double fi[], double (&money)[7], int size, ofstream& myfile) {
    for (int i=1 ; i <= size-2 ; i++) {
        if (fi[i] > fi[i-1]) {
            money[i] = fi[i] - fi[i-1];
        }
        else if (fi[i] < fi[i-1]) {
            money[i] = -1 * (fi[i-1] - fi[i]);
        }
        else {
            money[i] = 0;
        }
        money[6] += money[i];
    }
    money[5] = -1 * si[5];
    money[6] += money[5];

    myfile << left << setw(3) << setfill(' ') << "i";
    myfile << left << setw(7) << setfill(' ') << "si";
    myfile << left << setw(7) << setfill(' ') << "fi";
    myfile << left << setw(7) << setfill(' ') << "recvd";
    myfile << left << setw(7) << setfill(' ') << "paid";
    myfile << endl;
    for (int i=1 ; i<6 ; i++) {
        myfile << left << setw(3) << setfill(' ') << i;
        myfile << left << setw(7) << setfill(' ') << si[i];
        myfile << left << setw(7) << setfill(' ') << fi[i];
        if (money[i] > 0) {
            myfile << left << setw(7) << setfill(' ') << money[i];
            myfile << left << setw(7) << setfill(' ') << 0;
        }
        else if (money[i] < 0) {
            myfile << left << setw(7) << setfill(' ') << 0;
            myfile << left << setw(7) << setfill(' ') << (-1 * money[i]);
        }
        else {
            myfile << left << setw(7) << setfill(' ') << 0;
            myfile << left << setw(7) << setfill(' ') << 0;
        }
        myfile << endl;
    }
    myfile << "Total Money Paid: " << (-1 * money[6]) << endl;
}
int main() {
    
    // Output file
    ofstream myfile;
    myfile.open("solutions.txt");

    // 4.1 Forward Price Variables
    double fv = 0;
    double s0 = 0;
    double t0 = 0;
    double r = 0;
    double q = 0;
    double tf = 0;

    // 4.1.1 Fair Value for Below Parameters
    s0 = 100.5;
    t0 = 0;
    r = 5.5;
    q = 0;
    tf = 0.75;
    forward_price (s0, t0, r, q, tf, fv);
    myfile << "4.1.1 Solution: " << fv << endl << endl;

    // 4.1.2 Fair Value for Below Parameters
    s0 = 95.5;
    t0 = 0;
    r = 5.1;
    q = 1.1;
    tf = 0.65;
    forward_price (s0, t0, r, q, tf, fv);
    myfile << "4.1.2 Solution: " << fv << endl << endl;

    // 4.2 Forwards: Arbitrage Variables
    double ffv = 0;
    double profit = 0;
    string strategy = "";

    // 4.2 Paramaters
    s0 = 100.0;
    t0 = 0;
    r = 5.0;
    q = 0.0;
    tf = 1.0;
    forward_price (s0, t0, r, q, tf, fv);

    //4.2.1 Strategy & Profit for F=105.0
    ffv = 105.0;
    forward_profit(fv, ffv, profit, strategy);
    myfile << "4.2.1 Solution: " << strategy << profit << endl << endl;

    //4.2.2 Strategy & Profit for F=106.0
    ffv = 106.0;
    forward_profit(fv, ffv, profit, strategy);
    myfile << "4.2.2 Solution: " << strategy << profit << endl << endl;
    
    //4.3.1 Random Walk #1 Variables
    double rw1_s[6] = {100.0, 99.5, 101.3, 101.3, 100.2, 99.3};
    double rw1_f[6] = {105.5, 103.3, 104.1, 102.1, 101.3, 99.3};
    double rw1_m[7];
    myfile << "4.3.1 Solution: " << endl;
    random_walk(rw1_s, rw1_f, rw1_m, 6, myfile);
    myfile << endl;

    //4.3.2 Random Walk #2 Variables
    double rw2_s[6] = {100.0, 100.9, 103.8, 106.1, 107.5, 108.3};
    double rw2_f[6] = {105.5, 106.3, 108.7, 109.2, 108.3, 108.3};
    double rw2_m[7];
    myfile << "4.3.2 Solution: " << endl;
    random_walk(rw2_s, rw2_f, rw2_m, 6, myfile);
    myfile << endl;

    // 4.3.3 Summary
    myfile << "4.3.3 Solution: ";
    myfile << "Yes the random walk affect the total amount paid by the investor" << endl;

    myfile.close();
    return 0;
}