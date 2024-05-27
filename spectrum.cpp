#include<cstdio>
#include<cmath>
#include<thread>
#include<future>

double pi = 3.1415926535897931e+00;

double epsilon_c = 3.1, epsilon_e = .1, H_con_m = -1e-4, H_exp_m = 1e-4, Delta_eta_B= 1e-4, c_s=1.0;
//These parameters are wanted. They should be global parameters, but can be modified in functions.
//Except these, of cause, we still need the wavenumber *k*

inline double Bessel_1(double nu, double x){
	if(x>3e3){
		return sqrt(2 / pi / x) * cos(x - nu * pi / 2. - pi / 4.);
	}	
    double result = .0, ak = 1 / (pow(2., nu) * tgamma(nu + 1.));
    int kn;
    double y0, y1, y2, dx, x0;
    if(x > 1.){
        for(kn = 0; ak == .0 || ak+result!=result;){
            result += ak;
            kn += 2;
            if(nu != - 2.){
                ak *= - 1. / (kn * (kn + 2 * nu));
            }
            else{
                ak = pow(-1 , kn / 2.) / tgamma(kn / 2.) / tgamma(nu + kn / 2. + 1.);
            }
        }
        y0 = result * pow(1., nu);
        result = .0;
        ak = 1 / (pow(2., nu) * tgamma(nu + 1.));
        for(kn = 0; ak == .0 || ak+result!=result;){
            result += ak * (kn + nu);
            kn += 2;
            if(nu != - 2.){
                ak *= - 1. / (kn * (kn + 2 * nu));
            }
            else{
                ak = pow(-1 , kn / 2.) / tgamma(kn / 2.) / tgamma(nu + kn / 2. + 1.);
            }
        }
        y1 = result * pow(1., nu - 1);
        dx = 1e-3;
        for(x0 = 1.; x0 < x; x0 += dx){
            y2 = - (x0 * y1 + (x0 * x0 - nu * nu) * y0) / pow(x0, 2.);
            y1 += y2 * dx;
            y0 += y1 * dx;
        }
        x0 -=dx;
        y2 = - (x0 * y1 + (x0 * x0 - nu * nu) * y0) / pow(x0, 2.);
        y1 += y2 * (x - x0);
        y0 += y1 * (x - x0);
        return y0;
    }
    else{
        for(kn = 0; ak == .0 || fabs(ak)>fabs(result*1e-5)||result==.0;){
            result += ak;
            kn += 2;
            if(nu != - 1.){
                ak *= - 1. * x * x / (kn * (kn + 2 * nu));
            }
            else{
                ak = pow(-1 , kn / 2.) * x * x / tgamma(kn / 2.) / tgamma(nu + kn / 2. + 1.);
            }
        }
        return result * pow(x, nu);
    }
}

//double bessel_value[12];

inline double bessel_2_value(int n, double nu,double*bessel_value){
    double a0 = sin(nu * pi), a1 = cos(nu * pi);
    return bessel_value[n - 6] * a1 / a0 - bessel_value[n] / a0;
}

#define Hankel_1_re(n) (bessel_value[n])
#define Hankel_1_im(n) (bessel_value[n + 6])
#define Hankel_2_re(n) (bessel_value[n])
#define Hankel_2_im(n) (bessel_value[n + 6] * -1)

inline double c5_re(double k,double*bessel_value){
	double l = c_s * k;
    return - (( - sin((epsilon_c - 2) * pi / (2 * (epsilon_c - 1))) * (epsilon_e - 1) *
                    /*R*/
    sqrt(1 / (H_con_m - epsilon_c * H_con_m)) * H_exp_m * sqrt(1 / (H_exp_m - epsilon_e * H_exp_m)) *
    pow(pi, 1.5) * (k * (Hankel_1_re(1) - Hankel_1_re(2)) * (((epsilon_e - 1) * H_exp_m - 2 * k /
                            /*RR*/
    tan(Delta_eta_B * k)) * Hankel_2_re(3) + k * (Hankel_2_re(4) - Hankel_2_re(5))) - k *
                            /*RRR*/
    (Hankel_1_im(1) - Hankel_1_im(2)) * (((epsilon_e - 1) * H_exp_m - 2 * l / tan(Delta_eta_B * l)) *
        /*RI*/
    Hankel_2_im(3) + k * (Hankel_2_im(4) - Hankel_2_im(5))) + Hankel_1_re(0) * (((epsilon_c - 1) *
    /*RII*/                                                     /*RR*/
    (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 * ((epsilon_c - 1) * H_con_m + H_exp_m -
    epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_2_re(3) + k * ((epsilon_c - 1) *
                                                        /*RRR*/
    H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_2_re(4) - Hankel_2_re(5))) - Hankel_1_im(0) *
                                                                                    /*RI*/
    (((epsilon_c - 1) * (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 * ((epsilon_c - 1) *
    H_con_m + H_exp_m - epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_2_im(3) + k *
                                                                            /*RII*/
    ((epsilon_c - 1) * H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_2_im(4) - Hankel_2_im(5)))) -
    cos((epsilon_c - 2) * pi / (2 * (epsilon_c - 1))) * (epsilon_e - 1) *
    /*I*/
    sqrt(1 / (H_con_m - epsilon_c * H_con_m)) * H_exp_m * sqrt(1 / (H_exp_m - epsilon_e * H_exp_m)) *
    pow(pi, 1.5) * (k * (Hankel_1_im(1) - Hankel_1_im(2)) * (((epsilon_e - 1) * H_exp_m - 2 * k /
                            /*II*/
    tan(Delta_eta_B * k)) * Hankel_2_re(3) + k * (Hankel_2_re(4) - Hankel_2_re(5))) + k *
                            /*IIR*/
    (Hankel_1_re(1) - Hankel_1_re(2)) * (((epsilon_e - 1) * H_exp_m - 2 * l / tan(Delta_eta_B * l)) *
        /*IR*/
    Hankel_2_im(3) + k * (Hankel_2_im(4) - Hankel_2_im(5))) + Hankel_1_im(0) * (((epsilon_c - 1) *
    /*IRI*/                                                     /*II*/
    (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 * ((epsilon_c - 1) * H_con_m + H_exp_m -
    epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_2_re(3) + k * ((epsilon_c - 1) *
                                                        /*IIR*/
    H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_2_re(4) - Hankel_2_re(5))) + Hankel_1_re(0) *
                                                                                    /*IR*/
    (((epsilon_c - 1) * (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 * ((epsilon_c - 1) *
    H_con_m + H_exp_m - epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_2_im(3) + k *
                                                                            /*IRI*/
    ((epsilon_c - 1) * H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_2_im(4) - Hankel_2_im(5))))) *
    sin(Delta_eta_B * l) / (8 * l * (2 * (epsilon_e - 1) *
    H_exp_m + k * pi * bessel_value[5] * bessel_value[9] - k * pi * bessel_value[3] *
    bessel_value[11])));
}

inline double c5_im(double k,double*bessel_value){
	double l = c_s * k;
    return  - (( - sin((epsilon_c - 2) * pi / (2 * (epsilon_c - 1))) * (epsilon_e - 1) *
                    /*R*/
    sqrt(1 / (H_con_m - epsilon_c * H_con_m)) * H_exp_m * sqrt(1 / (H_exp_m - epsilon_e * H_exp_m)) *
    pow(pi, 1.5) * (k * (Hankel_1_re(1) - Hankel_1_re(2)) * (((epsilon_e - 1) * H_exp_m - 2 * k /
                            /*RR*/
    tan(Delta_eta_B * k)) * Hankel_2_im(3) + k * (Hankel_2_im(4) - Hankel_2_im(5))) + k *
                            /*RRI*/
    (Hankel_1_im(1) - Hankel_1_im(2)) * (((epsilon_e - 1) * H_exp_m - 2 * l / tan(Delta_eta_B * l)) *
        /*RI*/
    Hankel_2_re(3) + k * (Hankel_2_re(4) - Hankel_2_re(5))) + Hankel_1_re(0) * (((epsilon_c - 1) *
    /*RIR*/                                                     /*RR*/
    (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 * ((epsilon_c - 1) * H_con_m + H_exp_m -
    epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_2_im(3) + k * ((epsilon_c - 1) *
                                                        /*RRI*/
    H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_2_im(4) - Hankel_2_im(5))) + Hankel_1_im(0) *
                                                                                    /*RI*/
    (((epsilon_c - 1) * (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 * ((epsilon_c - 1) *
    H_con_m + H_exp_m - epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_2_re(3) + k *
                                                                            /*RIR*/
    ((epsilon_c - 1) * H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_2_re(4) - Hankel_2_re(5)))) +
    cos((epsilon_c - 2) * pi / (2 * (epsilon_c - 1))) * (epsilon_e - 1) *
    /*I*/
    sqrt(1 / (H_con_m - epsilon_c * H_con_m)) * H_exp_m * sqrt(1 / (H_exp_m - epsilon_e * H_exp_m)) *
    pow(pi, 1.5) * (k * (Hankel_1_re(1) - Hankel_1_re(2)) * (((epsilon_e - 1) * H_exp_m - 2 * k /
                            /*IR*/
    tan(Delta_eta_B * k)) * Hankel_2_re(3) + k * (Hankel_2_re(4) - Hankel_2_re(5))) - k *
                            /*IRR*/
    (Hankel_1_im(1) - Hankel_1_im(2)) * (((epsilon_e - 1) * H_exp_m - 2 * l / tan(Delta_eta_B * l)) *
        /*II*/
    Hankel_2_im(3) + k * (Hankel_2_im(4) - Hankel_2_im(5))) + Hankel_1_re(0) * (((epsilon_c - 1) *
    /*III*/                                                     /*IR*/
    (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 * ((epsilon_c - 1) * H_con_m + H_exp_m -
    epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_2_re(3) + k * ((epsilon_c - 1) *
                                                        /*IRR*/
    H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_2_re(4) - Hankel_2_re(5))) - Hankel_1_im(0) *
                                                                                    /*II*/
    (((epsilon_c - 1) * (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 * ((epsilon_c - 1) *
    H_con_m + H_exp_m - epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_2_im(3) + k *
                                                                            /*III*/
    ((epsilon_c - 1) * H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_2_im(4) - Hankel_2_im(5))))) *
    sin(Delta_eta_B * l) / (8 * l * (2 * (epsilon_e - 1) *
    H_exp_m + k * pi * bessel_value[5] * bessel_value[9] - k * pi * bessel_value[3] *
    bessel_value[11])));
}

inline double c6_re(double k,double*bessel_value){
	double l = c_s * k;
    return ( - sin((epsilon_c - 2) * pi / (2 * (epsilon_c - 1))) * (epsilon_e - 1) * sqrt(1 / (H_con_m -
                /*R*/
    epsilon_c * H_con_m)) * H_exp_m * sqrt(1 / (H_exp_m - epsilon_e * H_exp_m)) * pow(pi, 1.5) * (k *
    (Hankel_1_re(1) - Hankel_1_re(2)) * (((epsilon_e - 1) * H_exp_m - 2 * l / tan(Delta_eta_B * l)) *
        /*RR*/
    Hankel_1_re(3) + k * (Hankel_1_re(4) - Hankel_1_re(5))) - k *(Hankel_1_im(1) - Hankel_1_im(2)) *
    /*RRR*/                                                         /*RI*/
    (((epsilon_e - 1) * H_exp_m - 2 * l / tan(Delta_eta_B * l)) * Hankel_1_im(3) + k *
                                                                    /*RII*/
    (Hankel_1_im(4) - Hankel_1_im(5))) + Hankel_1_re(0) * (((epsilon_c - 1) * (epsilon_e - 1) *
                                            /*RR*/
    H_exp_m * H_con_m + 4 * l * l - 2 * ((epsilon_c - 1) * H_con_m + H_exp_m - epsilon_e * H_exp_m) *
    l / tan(Delta_eta_B * l)) * Hankel_1_re(3) + k * ((epsilon_c - 1) * H_con_m + 2 * k /
                                /*RRR*/
    tan(Delta_eta_B * k)) * (Hankel_1_re(4) - Hankel_1_re(5))) - Hankel_1_im(0) * (((epsilon_c - 1) *
                                                                    /*RI*/
    (epsilon_e - 1) * H_exp_m * H_con_m + 4 * l * l - 2 * ((epsilon_c - 1) * H_con_m + H_exp_m -
    epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_1_im(3) + k * ((epsilon_c - 1) *
                                                        /*RII*/
    H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_1_im(4) - Hankel_1_im(5)))) -
    cos((epsilon_c - 2) * pi / (2 * (epsilon_c - 1))) *
    /*I*/
    (epsilon_e - 1) * sqrt(1 / (H_con_m - epsilon_c * H_con_m)) * H_exp_m *
    sqrt(1 / (H_exp_m - epsilon_e * H_exp_m)) * pow(pi, 1.5) * (k * (Hankel_1_im(1) - Hankel_1_im(2)) *
                                                                        /*II*/
    (((epsilon_e - 1) * H_exp_m - 2 * l / tan(Delta_eta_B * l)) * Hankel_1_re(3) + k * (Hankel_1_re(4) -
                                                                    /*IIR*/
    Hankel_1_re(5))) + k * (Hankel_1_re(1) - Hankel_1_re(2)) * (((epsilon_e - 1) * H_exp_m - 2 * k /
                            /*IR*/
    tan(Delta_eta_B * k)) * Hankel_1_im(3) + k * (Hankel_1_im(4) - Hankel_1_im(5))) + Hankel_1_im(0) *
                            /*IRI*/                                                     /*II*/
    (((epsilon_c - 1) * (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 *((epsilon_c - 1) *
    H_con_m + H_exp_m - epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_1_re(3) + k *
                                                                            /*IIR*/
    ((epsilon_c - 1) * H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_1_re(4) - Hankel_1_re(5))) +
    Hankel_1_re(0) * (((epsilon_c - 1) * (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 *
    /*IR*/
    ((epsilon_c - 1) * H_con_m + H_exp_m - epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) *
    Hankel_1_im(3) + k * ((epsilon_c - 1) * H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_1_im(4) -
    /*IRI*/
    Hankel_1_im(5))))) *
    sin(Delta_eta_B * l) / (8 * l * (2 * (epsilon_e - 1) * H_exp_m + k * pi *
    bessel_value[5] * bessel_value[9] - k * pi * bessel_value[3] * bessel_value[11]));
}

inline double c6_im(double k,double*bessel_value){
	double l = c_s * k;
    return ( - sin((epsilon_c - 2) * pi / (2 * (epsilon_c - 1))) * (epsilon_e - 1) * sqrt(1 / (H_con_m -
                /*R*/
    epsilon_c * H_con_m)) * H_exp_m * sqrt(1 / (H_exp_m - epsilon_e * H_exp_m)) * pow(pi, 1.5) * (k *
    (Hankel_1_re(1) - Hankel_1_re(2)) * (((epsilon_e - 1) * H_exp_m - 2 * l / tan(Delta_eta_B * l)) *
        /*RR*/
    Hankel_1_im(3) + k * (Hankel_1_im(4) - Hankel_1_im(5))) + k *(Hankel_1_im(1) - Hankel_1_im(2)) *
    /*RRI*/                                                         /*RI*/
    (((epsilon_e - 1) * H_exp_m - 2 * l / tan(Delta_eta_B * l)) * Hankel_1_re(3) + k *
                                                                    /*RIR*/
    (Hankel_1_re(4) - Hankel_1_re(5))) + Hankel_1_re(0) * (((epsilon_c - 1) * (epsilon_e - 1) *
                                            /*RR*/
    H_exp_m * H_con_m + 4 * l * l - 2 * ((epsilon_c - 1) * H_con_m + H_exp_m - epsilon_e * H_exp_m) *
    l / tan(Delta_eta_B * l)) * Hankel_1_im(3) + k * ((epsilon_c - 1) * H_con_m + 2 * k /
                                /*RRI*/
    tan(Delta_eta_B * k)) * (Hankel_1_im(4) - Hankel_1_im(5))) + Hankel_1_im(0) * (((epsilon_c - 1) *
                                                                    /*RI*/
    (epsilon_e - 1) * H_exp_m * H_con_m + 4 * l * l - 2 * ((epsilon_c - 1) * H_con_m + H_exp_m -
    epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_1_re(3) + k * ((epsilon_c - 1) *
                                                        /*RIR*/
    H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_1_re(4) - Hankel_1_re(5)))) +
    cos((epsilon_c - 2) * pi / (2 * (epsilon_c - 1))) *
    /*I*/
    (epsilon_e - 1) * sqrt(1 / (H_con_m - epsilon_c * H_con_m)) * H_exp_m *
    sqrt(1 / (H_exp_m - epsilon_e * H_exp_m)) * pow(pi, 1.5) * (k * (Hankel_1_re(1) - Hankel_1_re(2)) *
                                                                        /*IR*/
    (((epsilon_e - 1) * H_exp_m - 2 * l / tan(Delta_eta_B * l)) * Hankel_1_re(3) + k * (Hankel_1_re(4) -
                                                                    /*IRR*/
    Hankel_1_re(5))) - k * (Hankel_1_im(1) - Hankel_1_im(2)) * (((epsilon_e - 1) * H_exp_m - 2 * k /
                            /*II*/
    tan(Delta_eta_B * k)) * Hankel_1_im(3) + k * (Hankel_1_im(4) - Hankel_1_im(5))) + Hankel_1_re(0) *
                            /*III*/                                                     /*IR*/
    (((epsilon_c - 1) * (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 *((epsilon_c - 1) *
    H_con_m + H_exp_m - epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) * Hankel_1_re(3) + k *
                                                                            /*IRR*/
    ((epsilon_c - 1) * H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_1_re(4) - Hankel_1_re(5))) -
    Hankel_1_im(0) * (((epsilon_c - 1) * (epsilon_e - 1) * H_con_m * H_exp_m + 4 * l * l - 2 *
    /*II*/
    ((epsilon_c - 1) * H_con_m + H_exp_m - epsilon_e * H_exp_m) * l / tan(Delta_eta_B * l)) *
    Hankel_1_im(3) + k * ((epsilon_c - 1) * H_con_m + 2 * l / tan(Delta_eta_B * l)) * (Hankel_1_im(4) -
    /*IRI*/
    Hankel_1_im(5))))) *
    sin(Delta_eta_B * l) / (8 * l * (2 * (epsilon_e - 1) * H_exp_m + k * pi *
    bessel_value[5] * bessel_value[9] - k * pi * bessel_value[3] * bessel_value[11]));
}

inline double power_spectrum_bounce_inflation(double k){
	double bessel_value[12];
	if(epsilon_c==1.)epsilon_c+=1e-5;
	if(epsilon_e==1.)epsilon_e+=1e-5;

	double nu1 = .5+1/(1-epsilon_c),nu2=1.5+1/(1-epsilon_c),nu3=(1+epsilon_c)/(2-2*epsilon_c), 
    nu4 = .5+1/(1-epsilon_e),nu5=1.5+1/(1-epsilon_e),nu6=(1+epsilon_e)/(2-2*epsilon_e);
    if(sin(nu1*pi)==.0||sin(nu2*pi)==.0||sin(nu3*pi)==.0){
        epsilon_c+=epsilon_c*1e-5;
    }

    bessel_value[0] = Bessel_1(.5 + 1 / (1 - epsilon_c), k / (H_con_m - epsilon_c * H_con_m));
    bessel_value[1] = Bessel_1(1.5 + 1 / (1 - epsilon_c), k / (H_con_m - epsilon_c * H_con_m));
    bessel_value[2] = Bessel_1((1 + epsilon_c) / (2 - 2 * epsilon_c), k / (H_con_m - epsilon_c * H_con_m));

    bessel_value[3] = Bessel_1(.5 + 1 / (1 - epsilon_e), k / (H_exp_m - epsilon_e * H_exp_m));
    bessel_value[4] = Bessel_1(1.5 + 1 / (1 - epsilon_e), k / (H_exp_m - epsilon_e * H_exp_m));
    bessel_value[5] = Bessel_1((1 + epsilon_e) / (2 - 2 * epsilon_e), k / (H_exp_m - epsilon_e * H_exp_m));
    
    bessel_value[6] = Bessel_1(-1 * (.5 + 1 / (1 - epsilon_c)), k / (H_con_m - epsilon_c * H_con_m));
    bessel_value[7] = Bessel_1(-1 * (1.5 + 1 / (1 - epsilon_c)), k / (H_con_m - epsilon_c * H_con_m));
    bessel_value[8] = Bessel_1(-1 * (1 + epsilon_c) / (2 - 2 * epsilon_c), k / (H_con_m - epsilon_c * H_con_m));
    
    bessel_value[9] = Bessel_1(-1 * (.5 + 1 / (1 - epsilon_e)), k / (H_exp_m - epsilon_e * H_exp_m));
    bessel_value[10] = Bessel_1(-1 * (1.5 + 1 / (1 - epsilon_e)), k / (H_exp_m - epsilon_e * H_exp_m));
    bessel_value[11] = Bessel_1(-1 * (1 + epsilon_e) / (2 - 2 * epsilon_e), k / (H_exp_m - epsilon_e * H_exp_m));

    bessel_value[6] = bessel_2_value(6, .5 + 1 / (1 - epsilon_c),bessel_value);
    bessel_value[7] = bessel_2_value(7, 1.5 + 1 / (1 - epsilon_c),bessel_value);
    bessel_value[8] = bessel_2_value(8, (1 + epsilon_c) / (2 - 2 * epsilon_c),bessel_value);

    bessel_value[9] = bessel_2_value(9, .5 + 1 / (1 - epsilon_e),bessel_value);
    bessel_value[10] = bessel_2_value(10, 1.5 + 1 / (1 - epsilon_e),bessel_value);
    bessel_value[11] = bessel_2_value(11, (1 + epsilon_e) / (2 - 2 * epsilon_e),bessel_value);

    /*for(int a=0;a<12;a++){
	    printf("%e ",bessel_value[a]);
    }*/

    return pow(c5_re(k,bessel_value) - c6_re(k,bessel_value), 2.) + pow(c5_im(k,bessel_value) - c6_im(k,bessel_value), 2.);
}

int main(int argc,char* argv[]){
    //int cpu_counts=std::thread::hardware_concurrency();
    int cpu_counts=20;
    double k[cpu_counts], ps, pst, delta_R;
    std::future<double>threads[cpu_counts];
    epsilon_c=atof(argv[1]);
    epsilon_e=atof(argv[2]);
    H_con_m=atof(argv[3]);
    H_exp_m=atof(argv[4]);
    Delta_eta_B=atof(argv[5]);
    //A_s=atof(argv[6]);
    //c_s=atof(argv[7]);
    if(Delta_eta_B==.0){
	    Delta_eta_B=1e-10;
    }
    //printf("%e\n",power_spectrum_bounce_inflation(1e-5));
    //printf("%e\n",Bessel_1(-2,1.1));
    for(int a=0;a<1000;a++){
        if(a>=cpu_counts){
		    delta_R=H_exp_m*H_exp_m/(8*pi*pi*epsilon_e);
            ps=delta_R*threads[a%cpu_counts].get();
            pst=ps*16*epsilon_e;
            printf("%e %e %e\n",k[a%cpu_counts],ps,pst);
        }
        k[a%cpu_counts] = pow(10.,a*7e-3+-6);
        threads[a%cpu_counts]=std::async(power_spectrum_bounce_inflation,k[a%cpu_counts]);
    }
    for(int a=1000;a<1000+cpu_counts;a++){
	delta_R=H_exp_m*H_exp_m/(8*pi*pi*epsilon_e);
        ps=delta_R*threads[a%cpu_counts].get();
        pst=ps*16*epsilon_e;
        printf("%e %e %e\n",k[a%cpu_counts],ps,pst);
    }
    return 0;
}

//g++ spectrum.cpp -o power_spectrum_bounce_inflation -lpthread

