#include <cmath>
#include <iostream>
#include <nlohmann/json.hpp>
#include <random>

class Bond {
public:
  double price, coupon_rate, discount_rate, face_value, periods_to_maturity;
  std::string name, currency;

  Bond(double cr, double dr, double fv, double ptm) {
    coupon_rate = cr;
    discount_rate = dr;
    face_value = fv;
    periods_to_maturity = ptm;
  }

  Bond(std::string n, std::string c, double cr, double dr, double fv,
       double ptm) {
    name = n;
    currency = c;
    coupon_rate = cr;
    discount_rate = dr;
    face_value = fv;
    periods_to_maturity = ptm;
  }

  void update_price() {
    // P = ((FC)((1 - (1 + D) ^ -N) / D)) + ((F)(1 + D) ^ -N)
    price = ((this->face_value * this->coupon_rate) * ((1 - std::pow(1 + this->discount_rate, (-1 * this->periods_to_maturity))) / this->discount_rate)) + ((this->face_value) * (std::pow(1 + this->discount_rate, (-1 * this->periods_to_maturity))));
  }
};

double bond_price(Bond bond) {
  // calculate the discounted value of the face value
  double price = bond.face_value /
                 std::pow(1 + bond.discount_rate, bond.periods_to_maturity);

  // add the value of the periodic coupon payments
  for (int i = 1; i <= (int)bond.periods_to_maturity; i++) {
    price += (bond.face_value * bond.coupon_rate) /
             std::pow(1 + bond.discount_rate, i);
  }

  return price;
}

class RandDouble {
  private:
  std::random_device rd;
    std::default_random_engine re;
    std::uniform_real_distribution<double> rng;

  public:
  RandDouble(double min, double max) {
    re = std::default_random_engine (rd());
    rng = std::uniform_real_distribution<double> (min, max);
  }

  double next() { // I miss you C# ♡ ♡ ♡
    return rng(re);
  }

};

int main() {
  std::cout << "Bond Pricer - Joshua O'Riordan" << std::endl << std::endl;

  if (true) { // comments
    // Bond Attributes:

    // Face Value - The price that the bond issuer pays upon bond maturation.
    // Coupon rate - A fraction of the face value that dictates the value of coupon payments (if applicable). 
    // Discount rate - A risk based value that takes into account: inflation rate, risk of default, interest rate volatility etc.

    // Definitions of simple bond valuation formula - https://en.wikipedia.org/wiki/Bond_valuation#Present_value_approach

    // P = bond price
    // F = face value
    // D = discount rate
    // C = coupon rate
    // N = number of periods
    // t = current period

    // P = ∑((FC) / (1 + D) ^ t) + (F / (1 + D) ^ N)

    // remove ∑ for easier rearranging of P and F
    // P = ((FC)((1 - (1 + D) ^ -N) / D)) + ((F)(1 + D) ^ -N)
  }

  if (true) { // simple price calculations
    Bond example (0.05, 0.04, 1000, 10);
    double example_price = bond_price(example);

    std::cout << "Bond Price: " << example_price << std::endl << std::endl;

    // Alternative price calculation implementations
    double tfv =
        ((example.face_value * example.coupon_rate) *
         ((1 - std::pow(1 + example.discount_rate,
                        (-1 * example.periods_to_maturity))) /
          example.discount_rate)) +
        ((example.face_value) * (std::pow(1 + example.discount_rate,
                                          (-1 * example.periods_to_maturity))));
    std::cout << "TFV: " << tfv << std::endl;

    tfv = example.face_value *
          ((example.coupon_rate *
            ((1 - std::pow(1 + example.discount_rate,
                           (-1 * example.periods_to_maturity))) /
             example.discount_rate)) +
           std::pow(1 + example.discount_rate,
                    (-1 * example.periods_to_maturity)));
    std::cout << "TFV II: " << tfv << std::endl;

    double face_val =
        tfv / ((example.coupon_rate *
                ((1 - std::pow(1 + example.discount_rate,
                               (-1 * example.periods_to_maturity))) /
                 example.discount_rate)) +
               std::pow(1 + example.discount_rate,
                        (-1 * example.periods_to_maturity)));
    std::cout << "FACV: " << face_val << std::endl << std::endl;

    // Coupon yield as a decimal, coupon / recalculated price
    for (int i = 0; i < 10; i++) {
      double cur_bond_price = bond_price(example);
      example.periods_to_maturity -= 1;

      std::cout << "Coupon Yield (" << i << "): "
                << (example.face_value * example.coupon_rate) / cur_bond_price
                << std::endl;
    }
    std::cout << std::endl;
  }

  if (true) { // testing discount rate variations effects on pricing
    Bond dr_var(0.05, 0.00, 1000, 10);
    for (int i = 0; i < 5; i++) {
      std::cout << "DR = " << dr_var.discount_rate
                << " -> Price = " << bond_price(dr_var) << std::endl;
      dr_var.discount_rate += 0.01; // vary by 1%
    }
    std::cout << std::endl;
  }

  if (true) { // comparing a bond investment against cash left to inflation
    double inflation_rate = 0.034; // 3.4% per year, pulled from a google search (USA)

    Bond inflation_test(0.05, 0.04, 1000, 10);
    double price = bond_price(inflation_test);
    double compare = price;
    for (int i = 0; i < 10; i++) {
      compare *= (1 + inflation_rate);
    }
    std::cout << "Example Bond Price: " << price << std::endl;
    std::cout << "Bond Price adjusted to maturity date: " << compare << std::endl;
    std::cout << "Bond Payout: " << inflation_test.face_value + ((inflation_test.coupon_rate * inflation_test.face_value) * inflation_test.periods_to_maturity) << std::endl;
    std::cout << std::endl;

  }

  if (true) { // calculate discount rate based on expected default rate
    // Df = discount factor
    // Dr = discount rate
    // p = period

    // Df = 1 / ((1 + Dr) ^ p)
    // Df = (1 + Dr) ^ -p
    
    // Dr = (Df ^ (1 / -p)) - 1
    Bond defaulted(0.05, 0.04, 1000, 10);

    double unadjusted_price = bond_price(defaulted);
    unadjusted_price = 1000;

    for (int i = 1; i <= 5; i++) {
      double def_risk = i * 0.01; // i% risk of default
      
      std::cout << def_risk << "% Default risk: " << unadjusted_price * (1 / (1 + def_risk)) << std::endl; // treats def_risk as a lifetime value

      double Df = 1 / std::pow(1 + def_risk, defaulted.periods_to_maturity);
      // double Df = std::pow(1 + def_risk, -1 * defaulted.periods_to_maturity); // would also work
      std::cout << unadjusted_price * Df << std::endl; // treats def_risk as a yearly

      double Dr = std::pow(Df, 1 / (-1 * defaulted.periods_to_maturity)) - 1; // already have def risk, this is just a test of my equation rearrangement
      // std::cout << Dr << std::endl;
    }
    /* interesting note here, the two formulae implemented return different results,
    one formula assumes def_risk is the lifetime default risk of the bond, the other 
    treats it as a yearly (periodic) risk of default */
    std::cout << std::endl;
  }

  if (true) { // calculate desired discount factor using forecasted inflation and default rates (rng)
    Bond var_rates (0.05, 0.04, 1000, 10);
    std::vector<double> inf_rates (var_rates.periods_to_maturity);
    std::vector<double> def_rates (var_rates.periods_to_maturity);

    RandDouble rng(0, 0.05);

    for (int i = 0; i < var_rates.periods_to_maturity; i++) {
      // all values for inflation rate and default probability are randomly generated numbers between 0.0% and 5.0%
      inf_rates[i] = rng.next();
      def_rates[i] = rng.next();

      // std::cout << inf_rates[i] << "   " << def_rates[i] << std::endl;
    }
    // std::cout << std::endl;

    double face_discount_factor = 1;
    double inf_only_fdf = 1; // inflation only face discount factor

    for (int i = 0; i < var_rates.periods_to_maturity; i++) {
      face_discount_factor *= 1 / (1 + inf_rates[i] + def_rates[i]);

      inf_only_fdf *= 1 / (1 + inf_rates[i]);
    }

    std::cout << "Face Value Discount Factor: " << face_discount_factor << std::endl;
    std::cout << "Inflation Only FVDF: " << inf_only_fdf << std::endl;

    // average inr_rate + def_rate for each period is 0.05, so use that to form a typical discount rate for the bond duration
    double tmp = 1 / std::pow(1 + 0.05, var_rates.periods_to_maturity);
    std::cout << "Avg Disc: " << tmp << std::endl;
  }

  return 0;
}
