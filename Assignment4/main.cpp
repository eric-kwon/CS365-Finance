// Version Update (10-12-17)
// - Edited random_walk function to calculate the mark to market value for day 5
// - Edited random_walk to allow partial runs
// - Added function that compares initial future value to the value after random walks
// - Added question solutions from 4.3.4 Addendum

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
void forward_profit (double fv, double ffv, ofstream &file) {
    if (fv > ffv) {
        file << "Strategy: Short Selling for profit of: " << fv - ffv << endl;
        file << "1. Enter into forward contract to BUY" << endl;
        file << "2. Short sell the promised stock to gain 100.0 cashflow" << endl;
        file << "3. Using this cashflow, accrue interest until maturity (value): " << fv << endl;
        file << "4. Upon maturity use the " << fv << " to cover the forward price of " << ffv << endl;
        file << "5. Then, the remaining value is net profit: " << fv - ffv << endl << endl;
    }
    else {
        file << "Strategy: Buy Low and Sell High for profit of: " << ffv - fv << endl;
        file << "1. Enter into forward contract to SELL" << endl;
        file << "2. Buy a stock by borrowing 100.0" << endl;
        file << "3. Interest will accrue on the money borrowed to become " << fv << " at maturity" << endl;
        file << "4. Upon maturity sell the stock that was bought and receive " << ffv << " and pay " << fv << " to the lender" << endl;
        file << "5. Then, the remaining value is net profit: " << ffv - fv << endl << endl;
    }
}

// 4.3 Futures: Mark to Market
void random_walk (double si[], double fi[], double (&money)[7], int size, ofstream &myfile) {
    for (int i=1 ; i < size ; i++) {
        if (fi[i] > fi[i-1]) {
            money[i] = fi[i] - fi[i-1];
        }
        else if (fi[i] < fi[i-1]) {
            money[i] = -1 * (fi[i-1] - fi[i]);
        }
        else {
            money[i] = 0;
        }
        money[size] += money[i];
    }
    money[size] += (-1 * fi[size-1]);

    myfile << left << setw(3) << setfill(' ') << "i";
    myfile << left << setw(7) << setfill(' ') << "si";
    myfile << left << setw(7) << setfill(' ') << "fi";
    myfile << left << setw(7) << setfill(' ') << "recvd";
    myfile << left << setw(7) << setfill(' ') << "paid";
    myfile << endl;
    for (int i=1 ; i<size ; i++) {
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
    myfile << "Total Money Paid: " << (-1 * money[size]) << endl;
}

// 4.3 Futures: Premature Sale
void premature_sale (double si[], double fi[], double (&money)[7], int size, ofstream &myfile) {
    for (int i=1 ; i < size ; i++) {
        if (fi[i] > fi[i-1]) {
            money[i] = fi[i] - fi[i-1];
        }
        else if (fi[i] < fi[i-1]) {
            money[i] = -1 * (fi[i-1] - fi[i]);
        }
        else {
            money[i] = 0;
        }
        money[size] += money[i];
    }

    myfile << left << setw(3) << setfill(' ') << "i";
    myfile << left << setw(7) << setfill(' ') << "si";
    myfile << left << setw(7) << setfill(' ') << "fi";
    myfile << left << setw(7) << setfill(' ') << "recvd";
    myfile << left << setw(7) << setfill(' ') << "paid";
    myfile << endl;
    for (int i=1 ; i<size ; i++) {
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
    if (money[size] > 0)
        myfile << "Account Gains: " << money[size] << endl;
    else if (money[size] < 0)
        myfile << "Account Losses: " << (-1 * money[size]) << endl;
    else
        myfile << "No Gain/Loss" << endl;
}

// 4.3.3 Discrepancy Check - Random Walk
string random_walk_check (double f0, double ff) {
    ff *= -1;
    string result = "";
    if (f0 == ff)
        return result = "Random Walk DOES NOT Affect Amount Paid";
    else
        return result = "Random Walk DOES Affect Amount Paid";
}

// 4.3.7 Discrepancy Check - Premature Sale
string premature_sale_check (double pf1, double pf2) {
    string result = "";
    if (pf1 == pf2)
        return result = "Profit/Loss IS NOT Affected";
    else   
        return result = "Profit/Loss IS Affected";
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
    forward_price(s0, t0, r, q, tf, fv);
    myfile << "4.1.1 Solution: " << fv << endl << endl;

    // 4.1.2 Fair Value for Below Parameters
    s0 = 95.5;
    t0 = 0;
    r = 5.1;
    q = 1.1;
    tf = 0.65;
    forward_price(s0, t0, r, q, tf, fv);
    myfile << "4.1.2 Solution: " << fv << endl << endl;

    // 4.2 Forwards: Arbitrage Variables
    double ffv = 0;
    string strategy = "";

    // 4.2 Paramaters
    s0 = 100.0;
    t0 = 0;
    r = 5.0;
    q = 0.0;
    tf = 1.0;
    forward_price(s0, t0, r, q, tf, fv);

    //4.2.1 Strategy & Profit for F=105.0
    ffv = 105.0;
    myfile << "4.2.1 Solution: " << endl;
    forward_profit(fv, ffv, myfile);

    //4.2.2 Strategy & Profit for F=106.0
    ffv = 106.0;
    myfile << "4.2.2 Solution: " << endl;
    forward_profit(fv, ffv, myfile);
    
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
    myfile << "4.3.3 Solution: " << endl;
    myfile << "For 4.3.1: " << random_walk_check(rw1_f[0], rw1_m[6]) << endl;
    myfile << "For 4.3.2: " << random_walk_check(rw2_f[0], rw2_m[6]) << endl << endl;

    //4.3.5 Random Walk #3 Variables
    double rw3_s[6] = {100.0, 99.5, 101.3, 101.3, 100.1, 102.3};
    double rw3_f[6] = {105.5, 106.3, 105.1, 105.8, 104.2, 102.3};
    double rw3_m[7];
    myfile << "4.3.5 Solution: Full Run" << endl;
    random_walk(rw3_s, rw3_f, rw3_m, 6, myfile);
    myfile << endl;

    //4.3.5.1 Random Walk #3 Until Day 2
    double rw3_m_2[7];
    myfile << "4.3.5 Solution: Sell on Day 2" << endl;
    premature_sale(rw3_s, rw3_f, rw3_m_2, 3, myfile);
    myfile << endl;

    //4.3.5.2 Random Walk #3 Until Day 3
    double rw3_m_3[7];
    myfile << "4.3.5 Solution: Sell on Day 2" << endl;
    premature_sale(rw3_s, rw3_f, rw3_m_3, 4, myfile);
    myfile << endl;

    //4.3.6 Random Walk #4 Variables
    double rw4_s[6] = {100.0, 99.7, 101.2, 101.5, 102.8, 102.3};
    double rw4_f[6] = {105.5, 106.3, 105.1, 105.8, 104.2, 102.3};
    double rw4_m[7];
    myfile << "4.3.6 Solution: Full Run" << endl;
    random_walk(rw4_s, rw4_f, rw4_m, 6, myfile);
    myfile << endl;

    //4.3.6.1 Random Walk #3 Until Day 2
    double rw4_m_2[7];
    myfile << "4.3.6 Solution: Sell on Day 2" << endl;
    premature_sale(rw4_s, rw4_f, rw4_m_2, 3, myfile);
    myfile << endl;

    //4.3.6.2 Random Walk #3 Until Day 3
    double rw4_m_3[7];
    myfile << "4.3.6 Solution: Sell on Day 2" << endl;
    premature_sale(rw4_s, rw4_f, rw4_m_3, 4, myfile);
    myfile << endl;

    //4.3.7 Summary
    myfile << "4.3.7 Solution:" << endl;
    myfile << "For Sale on Day 2: " << premature_sale_check(rw3_m_2[3], rw4_m_2[3]) << endl;
    myfile << "For Sale on Day 3: " << premature_sale_check(rw3_m_3[4], rw4_m_3[4]);

    myfile.close();
    return 0;
}