
// CSCI365 Homework Assignment 1
// Name: Eric Kwon

#include <iostream>
#include <math.h>
using namespace std;

// Problem #1 - Future Value
double future_value (double f0, double t0, double t1, double r) {
    double r_decimal = 0.01 * r;
    double f1 = f0 * exp(r_decimal * (t1-t0));
    return f1;
}

void load_future_value () {
    double f0, t0, t1, r;
    cout << "Enter principal amount (USD): ";
    cin >> f0;
    cout << "Enter today's year (yyyy): ";
    cin >> t0;
    cout << "Enter future year (yyyy): ";
    cin >> t1;
    cout << "Enter the interest rate (%): ";
    cin >> r;
    cout << "Future value (USD): " << future_value(f0, t0, t1, r) << endl;
}

// Problem #2 - Discount Factor and Rate
int df_and_r (double f0, double f1, double t0, double t1, double &df, double &r) {
    if (t1 - t0 <= 0.0) {
        df = 0;
        r = 0;
        return -1;
    }
    if ((f0 < 0.0) || (f1 < 0.0)) {
        df = 0;
        r = 0;
        return -2;
    }
    r = log(f1 / f0) / (t1 - t0);
    double param = -1 * r * (t1 - t0);
    df = exp(param);
    return 0;
}

void load_df_and_r () {
    double df, r, f0, f1, t0, t1;
    int result = 0;
    cout << "Enter today's cashflow: ";
    cin >> f0;
    cout << "Enter future's cashflow: ";
    cin >> f1;
    cout << "Enter today's year: ";
    cin >> t0;
    cout << "Enter future's year: ";
    cin >> t1;
    result = df_and_r(f0, f1, t0, t1, df, r);
    if (result == -1)
        cout << "Enter a valid date";
    if (result == -2)
        cout << "Enter a valid cashflow";
    if (result == 0) {
        df = df * 100;
        r = r * 100;
        cout << "Discount Factor: " << df << "%\n";
        cout << "Rate: " << r << "%\n";
    }
}

// Problem #3 - Price from Yield
void price_from_yield (double f, double c, double y, int n, double &b) {
    y = y * 0.01;
    b = 0;
    for (int i=1 ; i < n ; i++) {
        b += (c / 2) / pow((1 + (y / 2)), i);
    }
    b += (f + (c / 2)) / pow((1 + (y / 2)), n);
}

void load_price_from_yield () {
    double f, c, y, n, b;
    cout << "Enter the face value (USD): ";
    cin >> f;
    cout << "Enter the coupon value (USD): ";
    cin >> c;
    cout << "Enter the yield rate (%): ";
    cin >> y;
    cout << "Number of payments per year: ";
    cin >> n;
    price_from_yield(f, c, y, n, b);
    cout << "Bond Price (USD): " << b << endl;
}

// Problem #4 - Yield from Price
int yield_from_price (double f, double c, int n, double b, double tol, int max_iter, double &y) {
    double y_low = 0; double y_high = 100;
    double b_low = 0; double b_high = 0;
    price_from_yield(f, c, y_low, n, b_low);
    if ((abs(b_low - b) <= tol)) {
        y = y_low;
        return 0;
    }
    if (b_low < b) {
        y = 0;
        return 1;
    }
    price_from_yield(f, c, y_high, n, b_high);
    if ((abs(b_high - b) <= tol)) {
        y = y_high;
        return 0;
    }
    if (b_high > b) {
        y = 0;
        return 1;
    }
    double b_mid = 0;
    for (int i = 0 ; i < max_iter ; i++) {
        double y_mid = (y_low + y_high) / 2;
        price_from_yield (f, c, y_mid, n, b_mid);
        if (abs(b_mid - b) <= tol) {
            y = y_mid;
            return 0;
        }
        else if (b_mid > b) {
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
    return -1;
}

void load_yield_from_price () {
    double f, c, b, tol, y;
    int n, max_iter;
    cout << "Enter the face value (USD): ";
    cin >> f;
    cout << "Enter the coupon value (USD): ";
    cin >> c;
    cout << "Enter bond market price (USD):";
    cin >> b;
    cout << "Number of payments per year: ";
    cin >> n;
    cout << "Enter the tolerance (%): ";
    cin >> tol;
    tol = tol * 0.01;
    cout << "Enter maximum number of iterations: ";
    cin >> max_iter;
    int result = yield_from_price(f, c, n, b, tol, max_iter, y);
    if (result == 0)
        cout << "Yield (%): " << y << endl;
    else
        cout << "Yield is out of normal bounds" << endl;
}

int main() {
//    Uncomment function that needs to be called
//    load_future_value();
//    load_df_and_r();
//    load_price_from_yield();
//    load_yield_from_price();
    return 0;
}