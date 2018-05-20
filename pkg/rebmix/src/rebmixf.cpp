#include <math.h>
#include <float.h>

#include <stdio.h>
#include <ctype.h>
#include <time.h>

#include "base.h"
#include "rebmixf.h"

#if (_MEMORY_LEAK_SWITCH)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

// CompnentDistribution constructor.

CompnentDistribution::CompnentDistribution(Base *owner)
{
    owner_ = owner;
    pdf_ = NULL;
    Theta_ = NULL;
} // CompnentDistribution

// CompnentDistribution destructor.

CompnentDistribution::~CompnentDistribution()
{
    int i;

    if (Theta_) {
        for (i = 0; i < length_Theta_; i++) {
            if (Theta_[i]) free(Theta_[i]);
        }

        free(Theta_);
    }

    if (pdf_) free(pdf_);
} // ~CompnentDistribution

int CompnentDistribution::Realloc(int length_pdf, int length_Theta, int *length_theta)
{
    int i;
    int Error = 0;

    length_pdf_ = length_pdf;

    pdf_ = (ParametricFamilyType_e*)realloc(pdf_, length_pdf_ * sizeof(ParametricFamilyType_e));

    Error = NULL == pdf_; if (Error) goto E0;

    length_Theta_ = length_Theta;

    length_theta_ = (int*)realloc(length_theta_, length_Theta_ * sizeof(int));

    Error = NULL == length_theta_; if (Error) goto E0;

    Theta_ = (FLOAT**)calloc((size_t)length_Theta_, sizeof(FLOAT*));

    Error = NULL == Theta_; if (Error) goto E0;

    for (i = 0; i < length_Theta_; i++) {
        length_theta_[i] = (int)labs(length_theta[i]);

        if (length_theta[i] > 0) {
            Theta_[i] = (FLOAT*)realloc(Theta_[i], length_theta_[i] * sizeof(FLOAT));

            Error = NULL == Theta_[i]; if (Error) goto E0;

            memset(Theta_[i], 0, length_theta_[i] * sizeof(FLOAT));
        }
    }

    for (i = 0; i < owner_->length_Theta_; i++) {
        owner_->length_theta_[i] = length_theta_[i];
    }

E0: return Error;
} // Realloc

int CompnentDistribution::Memmove(CompnentDistribution *CmpTheta)
{
    int i;
    int Error = 0;

    memmove(pdf_, CmpTheta->pdf_, length_pdf_ * sizeof(ParametricFamilyType_e));

    for (i = 0; i < length_Theta_; i++) if (CmpTheta->Theta_[i]) {
        memmove(Theta_[i], CmpTheta->Theta_[i], length_theta_[i] * sizeof(FLOAT));
    }

    return Error;
} // Memmove

// Rebmix constructor.

Rebmix::Rebmix()
{
    p_value_ = (FLOAT)0.0;
    min_dist_mul_ = (FLOAT)0.0;
    var_mul_ = (FLOAT)0.0;
    kmax_ = (FLOAT)0.0;
    ChiSqr_ = (FLOAT)0.0;
    curr_ = NULL;
    o_ = 0;
    open_ = NULL;
    save_ = NULL;
    Preprocessing_ = poHistogram;
    cmax_ = 15;
    Criterion_ = icAIC;
    Variables_ = NULL;
    IniTheta_ = NULL;
    length_K_ = 0;
    K_ = NULL;
    y0_ = NULL;
    ymin_ = NULL;
    ymax_ = NULL;
    ar_ = (FLOAT)0.1;
    Restraints_ = rtLoose;
    n_ = 0;
    Y_ = NULL;
    X_ = NULL;
    W_ = NULL;
    MixTheta_ = NULL;
    memset(&summary_, 0, sizeof(SummaryParameterType));
    opt_length_ = 0;
    opt_c_ = NULL;
    opt_IC_ = NULL;
    opt_logL_ = NULL;
    opt_D_ = NULL;
    all_length_ = 0;
    all_I_ = NULL;
    all_K_ = NULL;
    all_IC_ = NULL;
    memset(&additional_, 0, sizeof(AdditionalParameterType));
} // Rebmix

// Rebmix destructor.

Rebmix::~Rebmix()
{
    int i;

    if (all_IC_) free(all_IC_);

    if (all_K_) free(all_K_);

    if (all_I_) free(all_I_);

    if (opt_D_) free(opt_D_);

    if (opt_logL_) free(opt_logL_);

    if (opt_IC_) free(opt_IC_);

    if (opt_c_) free(opt_c_);

    if (summary_.h) free(summary_.h);

    if (summary_.ymax) free(summary_.ymax);

    if (summary_.ymin) free(summary_.ymin);

    if (summary_.y0) free(summary_.y0);

    if (MixTheta_) {
        for (i = 0; i < cmax_; i++) {
            if (MixTheta_[i]) delete MixTheta_[i];
        }

        delete[] MixTheta_;
    }

    if (W_) free(W_);

    if (X_) {
        for (i = 0; i < n_; i++) {
            if (X_[i]) free(X_[i]);
        }

        free(X_);
    }

    if (Y_) {
        for (i = 0; i < n_; i++) {
            if (Y_[i]) free(Y_[i]);
        }

        free(Y_);
    }

    if (ymax_) free(ymax_);

    if (ymin_) free(ymin_);

    if (y0_) free(y0_);

    if (K_) free(K_);

    if (IniTheta_) delete IniTheta_;

    if (Variables_) free(Variables_);

    if (save_) free(save_);

    if (open_) {
        for (i = 0; i < o_; i++) {
            if (open_[i]) free(open_[i]);
        }

        free(open_);
    }
} // ~Rebmix

// Adds number of classes or k-nearest neighbours to be processed. Returns 0 if none is added. Otherwise 1 is returned.

int Rebmix::Golden()
{
    FLOAT ICopt;
    int   i, iopt, Stop = 0;

    if (additional_.Bracket) {
        ICopt = FLOAT_MAX; iopt = 0;

        for (i = 0; i < all_length_; i++) if (all_K_[i]) {
            if (all_IC_[i] < ICopt) {
                ICopt = all_IC_[i]; iopt = i;
            }
        }

        additional_.a = 0; additional_.b = all_length_ - 1;

        for (i = 0; i < all_length_; i++) if (all_K_[i]) {
            if (i < iopt) {
                additional_.a = i;
            }
            else
            if (i > iopt) {
                additional_.b = i; break;
            }
        }

        additional_.c = additional_.b - (int)ceil((additional_.b - additional_.a) / Phi);
        additional_.d = additional_.a + (int)ceil((additional_.b - additional_.a) / Phi);

        all_K_[additional_.c] = additional_.c + all_K_[0];
        all_K_[additional_.d] = additional_.d + all_K_[0];

        additional_.Bracket = 0;
    }
    else {
        if (all_IC_[additional_.c] <= all_IC_[additional_.d]) {
            additional_.b = additional_.d;
        }
        else {
            additional_.a = additional_.c;
        }

        additional_.c = additional_.b - (int)ceil((additional_.b - additional_.a) / Phi);
        additional_.d = additional_.a + (int)ceil((additional_.b - additional_.a) / Phi);

        Stop = additional_.b - additional_.a < 3;

        all_K_[additional_.a] = additional_.a + all_K_[0];
        all_K_[additional_.b] = additional_.b + all_K_[0];
        all_K_[additional_.c] = additional_.c + all_K_[0];
        all_K_[additional_.d] = additional_.d + all_K_[0];
    }

    return (Stop);
} // Golden

// Perform necessary initializations.

int Rebmix::Initialize()
{
    int Error = 0;

    p_value_ = (FLOAT)0.0001;

    min_dist_mul_ = (FLOAT)2.5;

    var_mul_ = (FLOAT)0.0625;

    kmax_ = ((FLOAT)1.0 + (FLOAT)1.0 / length_pdf_) * (FLOAT)pow((FLOAT)n_, (FLOAT)1.0 / ((FLOAT)1.0 + (FLOAT)1.0 / length_pdf_));

    Error = GammaInv((FLOAT)1.0 - (FLOAT)2.0 * p_value_, (FLOAT)2.0, length_pdf_ / (FLOAT)2.0, &ChiSqr_);
    
    return (Error);
} // Initialize

// Preprocessing of observations for k-nearest neighbour.

int Rebmix::PreprocessingKNN(int   k,    // k-nearest neighbours.
                             FLOAT *h,   // Normalizing vector.
                             FLOAT **Y)  // Pointer to the input array [y0,...,yd-1,kl,logV,R].
{
    FLOAT *Dk = NULL;
    FLOAT Dc, R, logV, logVn;
    int   i, j, l, m, q;
    int   Error = n_ < 1;

    if (Error) goto E0;

    if (k > 1) k -= 1; else k = 1;

    Dk = (FLOAT*)malloc(k * sizeof(FLOAT));

    Error = NULL == Dk; if (Error) goto E0;

    logVn = length_pdf_ * LogPi / (FLOAT)2.0 - Gammaln((FLOAT)1.0 + length_pdf_ / (FLOAT)2.0);

    for (i = 0; i < n_; i++) {
        Dk[0] = FLOAT_MAX; q = 0;

        for (j = 0; j < n_; j++) if (i != j) {
            Dc = (FLOAT)0.0;

            for (l = 0; l < length_pdf_; l++) {
                R = (Y[i][l] - Y[j][l]) / h[l]; Dc += R * R;
            }

            q += Dc <= FLOAT_MIN;

            for (l = 0; l < k; l++) {
                if (Dc < Dk[l]) {
                    for (m = k - 1; m > l; m--) Dk[m] = Dk[m - 1];

                    if ((Dc > FLOAT_MIN) || (l != k - 1)) Dk[l] = Dc;

                    break;
                }
            }
        }

        R = (FLOAT)sqrt(Dk[k - 1]);

        if (q >= k) R *= (FLOAT)exp((FLOAT)log((k + (FLOAT)1.0) / (q + (FLOAT)2.0)) / length_pdf_);

        logV = logVn + length_pdf_ * (FLOAT)log(R);

        for (j = 0; j < length_pdf_; j++) logV += (FLOAT)log(h[j]);

        Y[i][length_pdf_] = (FLOAT)1.0; Y[i][length_pdf_ + 1] = logV; Y[i][length_pdf_ + 2] = R;
    }

E0: if (Dk) free(Dk);

    return Error;
} // PreprocessingKNN

// Preprocessing of observations for Parzen window.

int Rebmix::PreprocessingPW(FLOAT *h,   // Sides of the hypersquare.
                            FLOAT **Y)  // Pointer to the input array [y0,...,yd-1,kl,k].
{
    int i, j, l;
    int Error = n_ < 1;

    if (Error) goto E0;

    for (i = 0; i < n_; i++) {
        Y[i][length_pdf_] = (FLOAT)1.0; Y[i][length_pdf_ + 1] = (FLOAT)0.0;
    }

    for (i = 0; i < n_; i++) {
        for (j = i; j < n_; j++) {
            for (l = 0; l < length_pdf_; l++) if ((FLOAT)fabs(Y[i][l] - Y[j][l]) > (FLOAT)0.5 * h[l]) goto S0;

            Y[i][length_pdf_ + 1] += (FLOAT)1.0; if (i != j) Y[j][length_pdf_ + 1] += (FLOAT)1.0;
S0:;
        }
    }

E0: return Error;
} // PreprocessingPW 

// Preprocessing of observations for histogram.

int Rebmix::PreprocessingH(FLOAT *h,  // Sides of the hypersquare.
                           FLOAT *y0, // Origins.
                           int   *k,  // Total number of bins.
                           FLOAT **Y) // Pointer to the input array [y0,...,yd-1,kl].
{
    int i, j, l;
    int Error = n_ < 1;

    if (Error) goto E0;

    *k = 0;

    for (i = 0; i < n_; i++) {
        for (j = 0; j < length_pdf_; j++) {
            l = (int)floor((Y_[i][j] - y0[j]) / h[j] + (FLOAT)0.5);

            Y[*k][j] = y0[j] + l * h[j];
        }

        for (j = 0; j < *k; j++) {
            for (l = 0; l < length_pdf_; l++) if ((FLOAT)fabs(Y[j][l] - Y[*k][l]) > (FLOAT)0.5 * h[l]) goto S0;

            Y[j][length_pdf_] += (FLOAT)1.0; goto S1;
S0:;
        }

        Y[*k][length_pdf_] = (FLOAT)1.0; (*k)++;
S1:;
    }
       
E0: return Error;
} // PreprocessingH

// Preprocessing of observations for histogram.

int Rebmix::PreprocessingH(FLOAT *h,     // Sides of the hypersquare.
                           FLOAT *y0,    // Origins.
                           int   *k,     // Total number of bins.
                           FLOAT **Y,    // Pointer to the input array [y0,...,yd-1,kl].
                           int   *State) // State variable.
{
    int i, j, l;
    int Error = n_ < 1;

    if (Error) goto E0;

    *k = 0;

    for (i = 0; i < n_; i++) {
        for (j = 0; j < length_pdf_; j++) {
            l = (int)floor((Y_[i][j] - y0[j]) / h[j] + (FLOAT)0.5);

            Y[*k][j] = y0[j] + l * h[j];
        }

        for (j = 0; j < *k; j++) {
            for (l = 0; l < length_pdf_; l++) if ((FLOAT)fabs(Y[j][l] - Y[*k][l]) >(FLOAT)0.5 * h[l]) goto S0;

            Y[j][length_pdf_] += (FLOAT)1.0; goto S1;
S0:;
        }

        Y[*k][length_pdf_] = (FLOAT)1.0; (*k)++;

        if ((*State) && (*k > kmax_)) {
            *State = 2; goto E0;
        }
S1:;
    }

E0: return Error;
} // PreprocessingH

// Global mode detection for k-nearest neighbour.

int Rebmix::GlobalModeKNN(int   *m,  // Global mode.
                          FLOAT **Y, // Pointer to the input array [y0,...,yd-1,kl].
                          FLOAT *h,  // Normalizing vector.
                          int   *I)  // Pointer to the inlier observations. 
{
    FLOAT Dc, R, cur, imax, omax, in, on;
    int   i, im, om, j, l;
    int   Error = 0, Stop = 0;

    im = om = -1; imax = omax = in = on = (FLOAT)0.0;

    for (i = 0; i < n_; i++) {
        cur = Y[i][length_pdf_] / (FLOAT)exp(Y[i][length_pdf_ + 1]);

        if (I[i]) {
            in += cur;

            if (cur > imax) {
                im = i; imax = cur;
            }
        }
        else {
            on += cur;

            if (cur > omax) {
                om = i; omax = cur;
            }
        }
    }

    if (im < 0) {
        *m = om;
    }
    else
    if (om < 0) {
        *m = im;
    }
    else
    if (omax > imax) {
        *m = om;
    }
    else {
        *m = om;       

        for (i = 0; i < n_; i++) if (I[i]) {
            for (j = 0; j < n_; j++) if (!I[j]) {
                Dc = (FLOAT)0.0;

                for (l = 0; l < length_pdf_; l++) {
                    R = (Y[i][l] - Y[j][l]) / h[l]; Dc += R * R;
                }

                R = (FLOAT)sqrt(Dc);

                Stop = R < min_dist_mul_ * (FLOAT)0.5 * (Y[i][length_pdf_ + 2] + Y[j][length_pdf_ + 2]);

                if (Stop) {
                    *m = im; break;
                }
            }

            if (Stop) break;
        }
    }

    return Error;
} // GlobalModeKNN 

// Global mode detection for Parzen window.

int Rebmix::GlobalModePW(int   *m,  // Global mode.
                         FLOAT **Y, // Pointer to the input array [y0,...,yd-1,kl].
                         FLOAT *h,  // Sides of the hypersquare.
                         int   *I)  // Pointer to the inlier observations.   
{
    FLOAT cur, imax, omax, in, on;
    int   i, im, om, j, l;
    int   Error = 0, Stop = 0;

    im = om = -1; imax = omax = in = on = (FLOAT)0.0;

    for (i = 0; i < n_; i++) {
        cur = Y[i][length_pdf_] * Y[i][length_pdf_ + 1];

        if (I[i]) {
            in += cur;

            if (cur > imax) {
                im = i; imax = cur;
            }
        }
        else {
            on += cur;

            if (cur > omax) {
                om = i; omax = cur;
            }
        }
    }

    if (im < 0) {
        *m = om;
    }
    else
    if (om < 0) {
        *m = im;
    }
    else
    if (omax > imax) {
        *m = om;
    }
    else {
        *m = om;       

        for (i = 0; i < n_; i++) if (I[i]) {
            for (j = 0; j < n_; j++) if (!I[j]) {
                Stop = 1;
        
                for (l = 0; l < length_pdf_; l++) {
                    Stop &= (FLOAT)fabs(Y[i][l] - Y[j][l]) < min_dist_mul_ * h[l];
                }

                if (Stop) {
                    *m = im; break;
                }
            }

            if (Stop) break;
        }
    }

    return Error;
} // GlobalModePW

// Global mode detection for histogram.

int Rebmix::GlobalModeH(int   *m,  // Global mode.
                        int   k,   // Total number of bins.
                        FLOAT **Y, // Pointer to the input array [y0,...,yd-1,kl].
                        FLOAT *h,  // Sides of the hypersquare.
                        int   *I)  // Pointer to the inlier observations. 
{
    FLOAT cur, imax, omax, in, on;
    int   i, im, om, j, l;
    int   Error = 0, Stop = 0;

    im = om = -1; imax = omax = in = on = (FLOAT)0.0;

    for (i = 0; i < k; i++) {
        cur = Y[i][length_pdf_];

        if (I[i]) {
            in += cur;

            if (cur > imax) {
                im = i; imax = cur;
            }
        }
        else {
            on += cur;

            if (cur > omax) {
                om = i; omax = cur;
            }
        }
    }

    if (im < 0) {
        *m = om;
    }
    else
    if (om < 0) {
        *m = im;
    }
    else
    if (omax > imax) {
        *m = om;
    }
    else {
        *m = om;       

        for (i = 0; i < k; i++) if (I[i]) {
            for (j = 0; j < k; j++) if (!I[j]) {
                Stop = 1;
        
                for (l = 0; l < length_pdf_; l++) {
                    Stop &= (FLOAT)fabs(Y[i][l] - Y[j][l]) < min_dist_mul_ * h[l];
                }

                if (Stop) {
                    *m = im; break;
                }
            }

            if (Stop) break;
        }
    }

    return Error;
} // GlobalModeH

// Returns rough normal parameters.

int RoughNormalParameters(FLOAT ym,
                          FLOAT fm,
                          FLOAT *Mean,
                          FLOAT *Stdev)
{
    int Error = 0;

    *Mean = ym; *Stdev = (FLOAT)1.0 / (SqrtPi2 * fm);

    return Error;
} // RoughNormalParameters

// Returns rough lognormal parameters.

int RoughLognormalParameters(FLOAT ym,
                             FLOAT fm,
                             FLOAT *Mean,
                             FLOAT *Stdev)
{
    FLOAT Lambda, dLambda;
    FLOAT A[3];
    int   i;
    int   Error = 0;

    Error = ym <= FLOAT_MIN; if (Error) goto E0;

    A[0] = (FLOAT)2.0 * (FLOAT)log(SqrtPi2 * ym * fm);

    if (Eps / ((FLOAT)1.0 + Eps) + (FLOAT)log(Eps * ((FLOAT)1.0 + Eps)) + A[0] >= (FLOAT)0.0) {
        Lambda = (FLOAT)1.0 + Eps;
    }
    else {
        i = 1; Lambda = (FLOAT)1.0 + Eps; Error = 1;
        while ((i <= ItMax) && Error) {
            A[1] = (FLOAT)1.0 / Lambda; A[2] = Lambda - (FLOAT)1.0;

            dLambda = ((FLOAT)1.0 - A[1] + (FLOAT)log(Lambda * A[2]) + A[0]) / (A[1] * ((FLOAT)1.0 + A[1]) + (FLOAT)1.0 / A[2]);

            Lambda -= dLambda;

            if (IsNan(dLambda) || IsInf(dLambda)) {
                Error = 1; goto E0;
            }

            if ((FLOAT)fabs(dLambda) < Max(Eps * (FLOAT)fabs(Lambda), Eps)) Error = 0;

            i++;
        }
    }

    *Mean = Lambda - (FLOAT)1.0 + (FLOAT)log(ym);

    *Stdev = (FLOAT)pow(Lambda * (Lambda - (FLOAT)1.0), (FLOAT)0.5);

E0: return Error;
} // RoughLognormalParameters

// Returns rough Weibull parameters.

int RoughWeibullParameters(FLOAT ym,
                           FLOAT fm,
                           FLOAT *Theta,
                           FLOAT *Beta)
{
    FLOAT Alpha, dAlpha;
    FLOAT A[4];
    int   i;
    int   Error = 0;

    Error = ym <= FLOAT_MIN; if (Error) goto E0;

    A[0] = (FLOAT)exp((FLOAT)1.0) * ym * fm;

    if ((FLOAT)0.064024 - A[0] >= (FLOAT)0.0) {
        Alpha = (FLOAT)1.234332;
    }
    else {
        i = 1; Alpha = (FLOAT)1.3349695; Error = 1;
        while ((i <= ItMax) && Error) {
            A[1] = Alpha - (FLOAT)1.0;

            A[2] = (FLOAT)1.0 + (Euler + (FLOAT)log(A[1] / Alpha)) / Alpha;

            A[3] = (FLOAT)exp((FLOAT)1.0 / Alpha);

            dAlpha = (A[2] * A[1] * A[3] - A[0]) / (A[3] * ((FLOAT)1.0 - (A[1] - A[2]) / Alpha / Alpha));

            Alpha -= dAlpha;

            if (IsNan(dAlpha) || IsInf(dAlpha)) {
                Error = 1; goto E0;
            }

            if ((FLOAT)fabs(dAlpha) < Max(Eps * (FLOAT)fabs(Alpha), Eps)) Error = 0;

            i++;
        }
    }

    if (Error) goto E0;

    *Beta = Alpha + Euler + (FLOAT)log((FLOAT)1.0 - (FLOAT)1.0 / Alpha);

    *Theta = ym * (FLOAT)pow(Alpha / (Alpha - (FLOAT)1.0), (FLOAT)1.0 / (*Beta));

E0: return Error;
} // RoughWeibullParameters

// Returns rough Gamma parameters.

int RoughGammaParameters(FLOAT ym,
                         FLOAT fm,
                         FLOAT *Theta,
                         FLOAT *Beta)
{

    FLOAT Alpha, dAlpha;
    FLOAT A[5];
    int   i;
    int   Error = 0;

    Error = ym <= FLOAT_MIN; if (Error) goto E0;

    A[0] = (FLOAT)log(ym * fm * SqrtPi2);

    if ((FLOAT)2.016083 + A[0] <= (FLOAT)0.0) {
        Alpha = (FLOAT)1.000299;
    }
    else {
        i = 1; Alpha = (FLOAT)1.000299; Error = 1;
        while ((i <= ItMax) && Error) {
            A[1] = (FLOAT)log((FLOAT)1.0 - (FLOAT)1.0 / Alpha);

            A[2] = A[1] + (FLOAT)1.0 / Alpha;

            A[3] = Euler * ((FLOAT)1.0 + Alpha) / (Euler - (FLOAT)1.0 - Alpha * A[1]);

            A[4] = A[3] * ((FLOAT)1.0 + A[3] * (A[1] + (FLOAT)1.0 / (Alpha - (FLOAT)1.0)) / Euler) / ((FLOAT)1.0 + Alpha);

            dAlpha = (A[3] * A[2] + (FLOAT)0.5 * (FLOAT)log(A[3]) - A[0]) / (A[4] * (A[2] + (FLOAT)0.5 / A[3]) + A[3] / (Alpha - (FLOAT)1.0) / Alpha / Alpha);

            Alpha -= dAlpha;

            if (IsNan(dAlpha) || IsInf(dAlpha)) {
                Error = 1; goto E0;
            }

            if ((FLOAT)fabs(dAlpha) < Max(Eps * (FLOAT)fabs(Alpha), Eps)) Error = 0;

            i++;
        }
    }

    if (Error) goto E0;

    *Beta = Euler * ((FLOAT)1.0 + Alpha) / (Euler - (FLOAT)1.0 - Alpha * (FLOAT)log((FLOAT)1.0 - (FLOAT)1.0 / Alpha));

    *Theta = ym * Alpha / (Alpha - (FLOAT)1.0) / (*Beta);

E0: return Error;
} // RoughGammaParameters

// Returns rough von Mises parameters.

int RoughvonMisesParameters(FLOAT h,
                            FLOAT ym,
                            FLOAT fm,
                            FLOAT *Mean,
                            FLOAT *Kappa)
{
    FLOAT dKappa;
    FLOAT A[3];
    int   i;
    int   Error = 0;

    A[0] = (FLOAT)0.5 * h; A[1] = Pi2 - A[0];

    if (ym < A[0]) {
        fm += (A[0] - ym) / (A[0] + ym) * fm;
    }
    else
    if (ym > A[1]) {
        fm += (A[0] + ym - Pi2) / (A[0] - ym + Pi2) * fm;
    }

    *Mean = ym;

    A[0] = (FLOAT)log(Pi2 * fm);

    if (A[0] <= (FLOAT)0.0) {
        *Kappa = (FLOAT)0.0;
    }
    else
    if (A[0] >= (FLOAT)3.75) {
        *Kappa = (FLOAT)288.0;
    }
    else {
        *Kappa = (FLOAT)0.0;

        i = 1; Error = 1;
        while ((i <= ItMax) && Error) {
            A[1] = BesselI0(*Kappa); A[2] = BesselI1(*Kappa);

            dKappa = (*Kappa - (FLOAT)log(A[1]) - A[0]) / ((FLOAT)1.0 - A[2] / A[1]);

            if (IsNan(dKappa) || IsInf(dKappa)) {
                Error = 1; goto E0;
            }

            *Kappa -= dKappa;

            if ((FLOAT)fabs(dKappa) < Max(Eps * (FLOAT)fabs(*Kappa), Eps)) Error = 0;

            i++;
        }
    }

E0: return Error;
} // RoughvonMisesParameters

// Returns rough binomial parameters.

int RoughBinomialParameters(FLOAT ym,
                            FLOAT fm,
                            FLOAT n,
                            FLOAT *p)
{
    int Error = 0;

    if ((int)ym == 0) {
        *p = (fm < (FLOAT)1.0) ? (FLOAT)1.0 - (FLOAT)pow(fm, (FLOAT)1.0 / n) : (FLOAT)0.0;
    }
    else
    if ((int)ym == (int)n) {
        *p = (fm < (FLOAT)1.0) ? (FLOAT)pow(fm, (FLOAT)1.0 / n) : (FLOAT)1.0;
    }
    else {
        *p = ym / n;
    }

    return Error;
} // RoughBinomialParameters

// Returns rough Poisson parameters.

int RoughPoissonParameters(FLOAT ym,
                           FLOAT fm,
                           FLOAT *Theta)
{
    int Error = 0;

    if ((int)ym == 0) {
        *Theta = (fm < (FLOAT)1.0) ? -(FLOAT)log(fm) : (FLOAT)0.0;
    }
    else {
        *Theta = ym;
    }

    return Error;
} // RoughPoissonParameters

// Returns component marginal p.d.f..

int ComponentMarginalDist(int                  i,           // Index of variable y.
                          FLOAT                *Y,          // Pointer to the input point [y0,...,yd-1].
                          CompnentDistribution *CmpTheta,   // Component distribution type.
                          FLOAT                *CmpMrgDist) // Component marginal distribution.
{
    FLOAT y, ypb, p, Theta;
    int   j, n;
    int   Error = 0;

    switch (CmpTheta->pdf_[i]) {
    case pfNormal:
        y = (Y[i] - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]);

        *CmpMrgDist = (FLOAT)exp(-(y * y)) / (SqrtPi2 * CmpTheta->Theta_[1][i]);

        break;
    case pfLognormal:
        if (Y[i] > FLOAT_MIN) {
            y = ((FLOAT)log(Y[i]) - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]);

            *CmpMrgDist = (FLOAT)exp(-(y * y)) / (SqrtPi2 * CmpTheta->Theta_[1][i]) / Y[i];
        }
        else {
            *CmpMrgDist = (FLOAT)0.0;
        }

        break;
    case pfWeibull:
        if (Y[i] > FLOAT_MIN) {
            ypb = (FLOAT)exp(CmpTheta->Theta_[1][i] * (FLOAT)log(Y[i] / CmpTheta->Theta_[0][i]));

            *CmpMrgDist = CmpTheta->Theta_[1][i] * ypb * (FLOAT)exp(-ypb) / Y[i];
        }
        else {
            *CmpMrgDist = (FLOAT)0.0;
        }

        break;
    case pfGamma:
        if (Y[i] > FLOAT_MIN) {
            ypb = Y[i] / CmpTheta->Theta_[0][i];

            *CmpMrgDist = (FLOAT)exp(CmpTheta->Theta_[1][i] * (FLOAT)log(ypb) - ypb - Gammaln(CmpTheta->Theta_[1][i])) / Y[i];
        }
        else {
            *CmpMrgDist = (FLOAT)0.0;
        }

        break;
    case pfvonMises:
        if ((Y[i] < (FLOAT)0.0) || (Y[i] > Pi2)) {
            *CmpMrgDist = (FLOAT)0.0;
        }
        else {
            *CmpMrgDist = (FLOAT)exp(CmpTheta->Theta_[1][i] * (FLOAT)cos(Y[i] - CmpTheta->Theta_[0][i])) / Pi2 / BesselI0(CmpTheta->Theta_[1][i]);
        }

        break;
    case pfBinomial:
        j = (int)Y[i]; n = (int)CmpTheta->Theta_[0][i]; p = CmpTheta->Theta_[1][i];

        if (j < 0)
            *CmpMrgDist = (FLOAT)0.0;
        else
        if (j == 0)
            *CmpMrgDist = (FLOAT)pow((FLOAT)1.0 - p, n);
        else
        if (j == n)
            *CmpMrgDist = (FLOAT)pow(p, n);
        else
        if (j > n)
            *CmpMrgDist = (FLOAT)0.0;
        else
            *CmpMrgDist = (FLOAT)exp(Gammaln(n + (FLOAT)1.0) - Gammaln(j + (FLOAT)1.0) - Gammaln(n - j + (FLOAT)1.0)) *
                          (FLOAT)pow(p, j) * (FLOAT)pow((FLOAT)1.0 - p, n - j);

        break;
    case pfPoisson:
        j = (int)Y[i]; Theta = CmpTheta->Theta_[0][i];

        *CmpMrgDist = (FLOAT)exp(j * (FLOAT)log(Theta) - Theta - Gammaln(j + (FLOAT)1.0));

        break;
    case pfDirac:
        if ((FLOAT)fabs(Y[i] - CmpTheta->Theta_[0][i]) > FLOAT_MIN) {
            *CmpMrgDist = (FLOAT)0.0;
        }
        else {
            *CmpMrgDist = (FLOAT)1.0;
        }
        
        break;
    case pfUniform:
        if ((Y[i] > CmpTheta->Theta_[1][i]) || (Y[i] < CmpTheta->Theta_[0][i])) {
            *CmpMrgDist = (FLOAT)0.0;
        }
        else {
            *CmpMrgDist = (FLOAT)1.0 / (CmpTheta->Theta_[1][i] - CmpTheta->Theta_[0][i]);
        }
    }

    return Error;
} // ComponentMarginalDist 

// Rough component parameter estimation for k-nearest neighbours.

int Rebmix::RoughEstimationKNN(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,logV,R].
                               int                  k,           // k-nearest neighbours.
                               FLOAT                *h,          // Normalizing vector.
                               FLOAT                nl,          // Total number of observations in class l.
                               int                  m,           // Mode index.
                               CompnentDistribution *RigidTheta, // Rigid parameters.
                               CompnentDistribution *LooseTheta) // Loose parameters.
{
    int                i, ii, j, l, *N = NULL;
    RoughParameterType *Mode = NULL;
    FLOAT              CmpMrgDist, epsilon, logflm, flm, flmin, flmax, Dlm, Dlmin, Dc, R, *D = NULL;
    int                Error = 0, Stop;

    Mode = (RoughParameterType*)malloc(length_pdf_ * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    N = (int*)malloc(length_pdf_ * sizeof(int));

    Error = NULL == N; if (Error) goto E0;

    D = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == D; if (Error) goto E0;

    // Rigid restraints.

    logflm = (FLOAT)0.0;

    for (i = 0; i < length_pdf_; i++) {
        N[i] = 0; D[i] = (FLOAT)2.0 * Y[m][length_pdf_ + 2] * h[i];

        if (length_pdf_ > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                Dc = (FLOAT)0.0;

                for (l = 0; l < length_pdf_; l++) if (i != l) {
                    R = (Y[j][l] - Y[m][l]) / h[l]; Dc += R * R;
                }

                R = (FLOAT)sqrt(Dc);

                if (R > Y[m][length_pdf_ + 2]) goto S0;

                Mode[i].klm += Y[j][length_pdf_];

                X_[N[i]][i] = Y[m][i] + (int)floor((Y[j][i] - Y[m][i]) / D[i] + (FLOAT)0.5) * D[i];

                for (ii = 0; ii < N[i]; ii++) {
                    if ((FLOAT)fabs(X_[N[i]][i] - X_[ii][i]) < (FLOAT)0.5 * D[i]) goto S0;
                }

                N[i] += 1;
S0:;
            } 
        }
        else {
            Mode[i].klm = nl;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                X_[N[i]][i] = Y[m][i] + (int)floor((Y[j][i] - Y[m][i]) / D[i] + (FLOAT)0.5) * D[i];

                for (ii = 0; ii < N[i]; ii++) {
                    if ((FLOAT)fabs(X_[N[i]][i] - X_[ii][i]) < (FLOAT)0.5 * D[i]) goto S1;
                }

                N[i] += 1;
S1:;
            }
        }

        Mode[i].h = (FLOAT)2.0 * Y[m][length_pdf_ + 2];  Mode[i].ym = Y[m][i]; Mode[i].flm = Y[m][length_pdf_] * k / (Mode[i].klm * (FLOAT)2.0 * Y[m][length_pdf_ + 2] * h[i]); logflm += (FLOAT)log(Mode[i].flm);
    }

    epsilon = (FLOAT)exp(((FLOAT)log(Y[m][length_pdf_] * k / nl) - Y[m][length_pdf_ + 1] - logflm) / length_pdf_);

    for (i = 0; i < length_pdf_; i++) {
        if (epsilon < (FLOAT)1.0) Mode[i].flm *= epsilon;

        switch (RigidTheta->pdf_[i]) {
        case pfNormal:
            Error = RoughNormalParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfLognormal:
            Error = RoughLognormalParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfWeibull:
            Error = RoughWeibullParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfGamma:
            Error = RoughGammaParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfvonMises:
            Error = RoughvonMisesParameters(Mode[i].h, Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfBinomial:
            Error = RoughBinomialParameters(Mode[i].ym, Mode[i].flm, RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfPoisson:
            Error = RoughPoissonParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i]);

            if (Error) goto E0;

            break;
        case pfDirac:
            RigidTheta->Theta_[0][i] = Mode[i].ym;

            break;
        case pfUniform:
            Error = 1; goto E0; 
        }
    }

    Error = LooseTheta->Memmove(RigidTheta);

    if (Error) goto E0;

    if (Restraints_ == rtRigid) goto E0;

    // Loose restraints.

    for (i = 0; i < length_pdf_; i++) if (N[i] > 1) {
        if (LooseTheta->pdf_[i] == pfDirac) goto E1;

        // Bracketing.

        Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

        for (j = 0; j < N[i]; j++) {
            Error = ComponentMarginalDist(i, X_[j], LooseTheta, &CmpMrgDist);

            if (Error) goto E0;

            Dlm -= CmpMrgDist * D[i];
        }

        if (Dlm > (FLOAT)0.0) goto E1;

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)1.0 - (FLOAT)2.0 * p_value_; flmax = Mode[i].flm;

        // Bisection.

        Stop = 0;

        while (!Stop) {
            flm = (flmax + flmin) / (FLOAT)2.0;

            switch (LooseTheta->pdf_[i]) {
            case pfNormal:
                Error = RoughNormalParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfLognormal:
                Error = RoughLognormalParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfWeibull:
                Error = RoughWeibullParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfGamma:
                Error = RoughGammaParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfvonMises:
                Error = RoughvonMisesParameters(Mode[i].h, Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfBinomial:
                Error = RoughBinomialParameters(Mode[i].ym, flm, LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfPoisson:
                Error = RoughPoissonParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i]);

                if (Error) goto E0;

                break;
            case pfDirac:
                break;
            case pfUniform:
                Error = 1; goto E0; 
            }

            Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

            for (j = 0; j < N[i]; j++) {
                Error = ComponentMarginalDist(i, X_[j], LooseTheta, &CmpMrgDist);

                if (Error) goto E0;

                Dlm -= CmpMrgDist * D[i];
            }

            if (((FLOAT)fabs(Dlm) < Eps) || (flmax - flmin < Eps)) {
                Stop = 1;
            }
            else {
                if (Dlm * Dlmin > (FLOAT)0.0) {
                    flmin = flm; Dlmin = Dlm;
                }
                else {
                    flmax = flm;
                }
            }
        }
E1:;
    }

E0: if (D) free(D);

    if (N) free(N);

    if (Mode) free(Mode);

    return Error;
} // RoughEstimationKNN 

// Rough component parameter estimation for Parzen window.

int Rebmix::RoughEstimationPW(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,k].
                              FLOAT                *h,          // Sides of the hypersquare.
                              FLOAT                nl,          // Total number of observations in class l.
                              int                  m,           // Mode index.
                              CompnentDistribution *RigidTheta, // Rigid parameters.
                              CompnentDistribution *LooseTheta) // Loose parameters.
{
    int                i, ii, j, l, *N = NULL;
    RoughParameterType *Mode = NULL;
    FLOAT              CmpMrgDist, epsilon, logflm, flm, flmin, flmax, logV, Dlm, Dlmin;
    int                Error = 0, Stop;

    Mode = (RoughParameterType*)malloc(length_pdf_ * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    N = (int*)malloc(length_pdf_ * sizeof(int));

    Error = NULL == N; if (Error) goto E0;

    // Rigid restraints.

    logflm = (FLOAT)0.0; logV = (FLOAT)0.0;

    for (i = 0; i < length_pdf_; i++) {
        logV += (FLOAT)log(h[i]); N[i] = 0;

        if (length_pdf_ > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                for (l = 0; l < length_pdf_; l++) if ((i != l) && ((FLOAT)fabs(Y[j][l] - Y[m][l]) > (FLOAT)0.5 * h[l])) goto S0;

                Mode[i].klm += Y[j][length_pdf_];

                X_[N[i]][i] = Y[m][i] + (int)floor((Y[j][i] - Y[m][i]) / h[i] + (FLOAT)0.5) * h[i];

                for (ii = 0; ii < N[i]; ii++) {
                    if ((FLOAT)fabs(X_[N[i]][i] - X_[ii][i]) < (FLOAT)0.5 * h[i]) goto S0;
                }

                N[i] += 1;
S0:;
            }
        }
        else {
            Mode[i].klm = nl;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                X_[N[i]][i] = Y[m][i] + (int)floor((Y[j][i] - Y[m][i]) / h[i] + (FLOAT)0.5) * h[i];

                for (ii = 0; ii < N[i]; ii++) {
                    if ((FLOAT)fabs(X_[N[i]][i] - X_[ii][i]) < (FLOAT)0.5 * h[i]) goto S1;
                }

                N[i] += 1;
S1:;
            }
        }

        Mode[i].h = h[i]; Mode[i].ym = Y[m][i]; Mode[i].flm = Y[m][length_pdf_] * Y[m][length_pdf_ + 1] / (Mode[i].klm * h[i]); logflm += (FLOAT)log(Mode[i].flm);
    }

    epsilon = (FLOAT)exp(((FLOAT)log(Y[m][length_pdf_] * Y[m][length_pdf_ + 1] / nl) - logV - logflm) / length_pdf_);

    for (i = 0; i < length_pdf_; i++) {
        if (epsilon < (FLOAT)1.0) Mode[i].flm *= epsilon;

        switch (RigidTheta->pdf_[i]) {
        case pfNormal:
            Error = RoughNormalParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfLognormal:
            Error = RoughLognormalParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfWeibull:
            Error = RoughWeibullParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfGamma:
            Error = RoughGammaParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfvonMises:
            Error = RoughvonMisesParameters(Mode[i].h, Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfBinomial:
            Error = RoughBinomialParameters(Mode[i].ym, Mode[i].flm, RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfPoisson:
            Error = RoughPoissonParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i]);

            if (Error) goto E0;

            break;
        case pfDirac:
            RigidTheta->Theta_[0][i] = Mode[i].ym;
            
            break;
        case pfUniform:
            Error = 1; goto E0; 
        }
    }

    Error = LooseTheta->Memmove(RigidTheta);

    if (Error) goto E0;

    if (Restraints_ == rtRigid) goto E0;

    // Loose restraints.

    for (i = 0; i < length_pdf_; i++) if (N[i] > 1) {
        if (LooseTheta->pdf_[i] == pfDirac) goto E1;

        // Bracketing.

        Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

        for (j = 0; j < N[i]; j++) {
            Error = ComponentMarginalDist(i, X_[j], LooseTheta, &CmpMrgDist);

            if (Error) goto E0;

            Dlm -= CmpMrgDist * h[i];
        }

        if (Dlm > (FLOAT)0.0) goto E1;

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)1.0 - (FLOAT)2.0 * p_value_; flmax = Mode[i].flm;

        // Bisection.

        Stop = 0;

        while (!Stop) {
            flm = (flmax + flmin) / (FLOAT)2.0;

            switch (LooseTheta->pdf_[i]) {
            case pfNormal:
                Error = RoughNormalParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfLognormal:
                Error = RoughLognormalParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfWeibull:
                Error = RoughWeibullParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfGamma:
                Error = RoughGammaParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfvonMises:
                Error = RoughvonMisesParameters(Mode[i].h, Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfBinomial:
                Error = RoughBinomialParameters(Mode[i].ym, flm, LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfPoisson:
                Error = RoughPoissonParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i]);

                if (Error) goto E0;

                break;
            case pfDirac:
                break;
            case pfUniform:
                Error = 1; goto E0; 
            }

            Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

            for (j = 0; j < N[i]; j++) {
                Error = ComponentMarginalDist(i, X_[j], LooseTheta, &CmpMrgDist);

                if (Error) goto E0;

                Dlm -= CmpMrgDist * h[i];
            }

            if (((FLOAT)fabs(Dlm) < Eps) || (flmax - flmin < Eps)) {
                Stop = 1;
            }
            else {
                if (Dlm * Dlmin > (FLOAT)0.0) {
                    flmin = flm; Dlmin = Dlm;
                }
                else {
                    flmax = flm;
                }
            }
        }
E1:;
    }

E0: if (N) free(N);

    if (Mode) free(Mode);

    return Error;
} // RoughEstimationPW 

// Rough component parameter estimation for histogram.

int Rebmix::RoughEstimationH(int                  k,           // Total number of bins.
                             FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl].
                             FLOAT                *h,          // Sides of the hypersquare.
                             FLOAT                nl,          // Total number of observations in class l.
                             int                  m,           // Mode index.
                             CompnentDistribution *RigidTheta, // Rigid parameters.
                             CompnentDistribution *LooseTheta) // Loose parameters.
{
    int                i, j, l, *N = NULL;
    RoughParameterType *Mode = NULL;
    FLOAT              CmpMrgDist, epsilon, logflm, flm, flmin, flmax, logV, Dlm, Dlmin;
    int                Error = 0, Stop;

    Mode = (RoughParameterType*)malloc(length_pdf_ * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    N = (int*)malloc(length_pdf_ * sizeof(int));

    Error = NULL == N; if (Error) goto E0;

    // Rigid restraints.

    logflm = (FLOAT)0.0; logV = (FLOAT)0.0;

    for (i = 0; i < length_pdf_; i++) {
        logV += (FLOAT)log(h[i]); N[i] = 0;

        if (length_pdf_ > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                for (l = 0; l < length_pdf_; l++) if ((i != l) && (Y[j][l] != Y[m][l])) goto S0;

                Mode[i].klm += Y[j][length_pdf_]; X_[N[i]][i] = Y[j][i]; N[i] += 1;
S0:;
            }
        }
        else {
            Mode[i].klm = nl;

            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                X_[N[i]][i] = Y[j][i]; N[i] += 1;
            }
        }

        Mode[i].h = h[i]; Mode[i].ym = Y[m][i]; Mode[i].flm = Y[m][length_pdf_] / (Mode[i].klm * h[i]); logflm += (FLOAT)log(Mode[i].flm);
    }

    epsilon = (FLOAT)exp(((FLOAT)log(Y[m][length_pdf_] / nl) - logV - logflm) / length_pdf_);

    for (i = 0; i < length_pdf_; i++) {
        if (epsilon < (FLOAT)1.0) Mode[i].flm *= epsilon;

        switch (RigidTheta->pdf_[i]) {
        case pfNormal:
            Error = RoughNormalParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfLognormal:
            Error = RoughLognormalParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfWeibull:
            Error = RoughWeibullParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfGamma:
            Error = RoughGammaParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfvonMises:
            Error = RoughvonMisesParameters(Mode[i].h, Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfBinomial:
            Error = RoughBinomialParameters(Mode[i].ym, Mode[i].flm, RigidTheta->Theta_[0][i], &RigidTheta->Theta_[1][i]);

            if (Error) goto E0;

            break;
        case pfPoisson:
            Error = RoughPoissonParameters(Mode[i].ym, Mode[i].flm, &RigidTheta->Theta_[0][i]);

            if (Error) goto E0;

            break;
        case pfDirac:
            RigidTheta->Theta_[0][i] = Mode[i].ym;

            break;
        case pfUniform:
            Error = 1; goto E0; 
        }
    }

    Error = LooseTheta->Memmove(RigidTheta);

    if (Error) goto E0;

    if (Restraints_ == rtRigid) goto E0;

    // Loose restraints.

    for (i = 0; i < length_pdf_; i++) if (N[i] > 1) {
        if (LooseTheta->pdf_[i] == pfDirac) goto E1;

        // Bracketing.

        Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

        for (j = 0; j < N[i]; j++) {
            Error = ComponentMarginalDist(i, X_[j], LooseTheta, &CmpMrgDist);

            if (Error) goto E0;

            Dlm -= CmpMrgDist * h[i];
        }

        if (Dlm > (FLOAT)0.0) goto E1;

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)1.0 - (FLOAT)2.0 * p_value_; flmax = Mode[i].flm;

        // Bisection.

        Stop = 0;

        while (!Stop) {
            flm = (flmax + flmin) / (FLOAT)2.0;

            switch (LooseTheta->pdf_[i]) {
            case pfNormal:
                Error = RoughNormalParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfLognormal:
                Error = RoughLognormalParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfWeibull:
                Error = RoughWeibullParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfGamma:
                Error = RoughGammaParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfvonMises:
                Error = RoughvonMisesParameters(Mode[i].h, Mode[i].ym, flm, &LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfBinomial:
                Error = RoughBinomialParameters(Mode[i].ym, flm, LooseTheta->Theta_[0][i], &LooseTheta->Theta_[1][i]);

                if (Error) goto E0;

                break;
            case pfPoisson:
                Error = RoughPoissonParameters(Mode[i].ym, flm, &LooseTheta->Theta_[0][i]);

                if (Error) goto E0;

                break;
            case pfDirac:
                break;
            case pfUniform:
                Error = 1; goto E0; 
            }

            Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

            for (j = 0; j < N[i]; j++) {
                Error = ComponentMarginalDist(i, X_[j], LooseTheta, &CmpMrgDist);

                if (Error) goto E0;

                Dlm -= CmpMrgDist * h[i];
            }

            if (((FLOAT)fabs(Dlm) < Eps) || (flmax - flmin < Eps)) {
                Stop = 1;
            }
            else {
                if (Dlm * Dlmin > (FLOAT)0.0) {
                    flmin = flm; Dlmin = Dlm;
                }
                else {
                    flmax = flm;
                }
            }
        }
E1:;
    }

E0: if (N) free(N);

    if (Mode) free(Mode);

    return Error;
} // RoughEstimationH 

// Returns component p.d.f..

int Rebmix::ComponentDist(FLOAT                *Y,        // Pointer to the input point [y0,...,yd-1].
                          CompnentDistribution *CmpTheta, // Component parameters.
                          FLOAT                *CmpDist,  // Component distribution.
                          int                  *Outlier)  // 1 if outlier otherwise 0.
{
    FLOAT y, ypb, p, Theta;
    int   i, j, n;
    int   Error = 0;

    *CmpDist = (FLOAT)1.0; if (Outlier) *Outlier = 0;   

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        switch (CmpTheta->pdf_[i]) {
        case pfNormal:
            y = (Y[i] - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]); y *= y;

            if (Outlier) {
                *Outlier |= (FLOAT)2.0 * y > ChiSqr_;
            }

            *CmpDist *= (FLOAT)exp(-y) / (SqrtPi2 * CmpTheta->Theta_[1][i]);

            break;
        case pfLognormal:
            if (Y[i] > FLOAT_MIN) {
                y = ((FLOAT)log(Y[i]) - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]); y *= y;

                if (Outlier) {
                    *Outlier |= (FLOAT)2.0 * y > ChiSqr_;
                }

                *CmpDist *= (FLOAT)exp(-y) / (SqrtPi2 * CmpTheta->Theta_[1][i]) / Y[i];
            }
            else {
                *CmpDist = (FLOAT)0.0; goto E0;
            }

            break;
        case pfWeibull:
            if (Y[i] > FLOAT_MIN) {
                if (Outlier) {
                    y = WeibullInv((FLOAT)1.0 - p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                    *Outlier |= Y[i] > y;

                    y = WeibullInv(p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                    *Outlier |= Y[i] < y;
                }

                ypb = (FLOAT)exp(CmpTheta->Theta_[1][i] * (FLOAT)log(Y[i] / CmpTheta->Theta_[0][i]));

                *CmpDist *= CmpTheta->Theta_[1][i] * ypb * (FLOAT)exp(-ypb) / Y[i];
            }
            else {
                *CmpDist = (FLOAT)0.0; goto E0;
            }

            break;
        case pfGamma:
            if (Y[i] > FLOAT_MIN) {
                if (Outlier) {
                    Error = GammaInv((FLOAT)1.0 - p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i], &y);

                    if (Error) goto E0;

                    *Outlier |= Y[i] > y;

                    Error = GammaInv(p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i], &y);

                    if (Error) goto E0;

                    *Outlier |= Y[i] < y;
                }

                ypb = Y[i] / CmpTheta->Theta_[0][i];

                *CmpDist *= (FLOAT)exp(CmpTheta->Theta_[1][i] * (FLOAT)log(ypb) - ypb - Gammaln(CmpTheta->Theta_[1][i])) / Y[i];
            }
            else {
                *CmpDist = (FLOAT)0.0; goto E0;
            }

            break;
        case pfvonMises:
            if (Outlier) {
                y = vonMisesInv((FLOAT)1.0 - p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                *Outlier |= Y[i] > y;

                y = vonMisesInv(p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                *Outlier |= Y[i] < y;
            }

            if ((Y[i] < (FLOAT)0.0) || (Y[i] > Pi2)) {
                *CmpDist = (FLOAT)0.0; goto E0;
            }
            else {
                *CmpDist *= (FLOAT)exp(CmpTheta->Theta_[1][i] * (FLOAT)cos(Y[i] - CmpTheta->Theta_[0][i])) / Pi2 / BesselI0(CmpTheta->Theta_[1][i]);
            }

            break;
        case pfBinomial:
            if (Outlier) {
                y = BinomialInv((FLOAT)1.0 - p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                *Outlier |= Y[i] > y;

                y = BinomialInv(p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                *Outlier |= Y[i] < y;
            }

            j = (int)Y[i]; n = (int)CmpTheta->Theta_[0][i]; p = CmpTheta->Theta_[1][i];

            if (j < 0) {
                *CmpDist = (FLOAT)0.0; goto E0;
            }
            else
            if (j == 0)
                *CmpDist *= (FLOAT)pow((FLOAT)1.0 - p, n);
            else
            if (j == n)
                *CmpDist *= (FLOAT)pow(p, n);
            else
            if (j > n) {
                *CmpDist = (FLOAT)0.0; goto E0;
            }
            else
                *CmpDist *= (FLOAT)exp(Gammaln(n + (FLOAT)1.0) - Gammaln(j + (FLOAT)1.0) - Gammaln(n - j + (FLOAT)1.0)) *
                            (FLOAT)pow(p, j) * (FLOAT)pow((FLOAT)1.0 - p, n - j);

            break;
        case pfPoisson:
            if (Outlier) {
                y = PoissonInv((FLOAT)1.0 - p_value_, CmpTheta->Theta_[0][i]);

                *Outlier |= Y[i] > y;

                y = PoissonInv(p_value_, CmpTheta->Theta_[0][i]);

                *Outlier |= Y[i] < y;
            }

            j = (int)Y[i]; Theta = CmpTheta->Theta_[0][i];

            *CmpDist *= (FLOAT)exp(j * (FLOAT)log(Theta) - Theta - Gammaln(j + (FLOAT)1.0));

            break;
        case pfDirac:
            if ((FLOAT)fabs(Y[i] - CmpTheta->Theta_[0][i]) > FLOAT_MIN) {
                *CmpDist = (FLOAT)0.0; goto E0;
            }
            else {
                *CmpDist *= (FLOAT)1.0;
            }

            break;
        case pfUniform: 
            if ((Y[i] > CmpTheta->Theta_[1][i]) || (Y[i] < CmpTheta->Theta_[0][i])) {
                *CmpDist = (FLOAT)0.0; goto E0;
            }
            else {
                *CmpDist *= (FLOAT)1.0 / (CmpTheta->Theta_[1][i] - CmpTheta->Theta_[0][i]);
            }
        }
    }

E0: return Error;
} // ComponentDist

// Returns logarithm of component p.d.f..

int Rebmix::LogComponentDist(FLOAT                *Y,        // Pointer to the input point [y0,...,yd-1].
                             CompnentDistribution *CmpTheta, // Component parameters.
                             FLOAT                *CmpDist,  // Component distribution.
                             int                  *Outlier)  // 1 if outlier otherwise 0.
{
    FLOAT y, ypb, p, Theta;
    int   i, j, n;
    int   Error = 0;

    *CmpDist = (FLOAT)0.0; if (Outlier) *Outlier = 0;

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        switch (CmpTheta->pdf_[i]) {
        case pfNormal:
            y = (Y[i] - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]); y *= y;

            if (Outlier) {
                *Outlier |= (FLOAT)2.0 * y > ChiSqr_;
            }

            *CmpDist += -y - LogSqrtPi2  - (FLOAT)log(CmpTheta->Theta_[1][i]);

            break;
        case pfLognormal:
            if (Y[i] > FLOAT_MIN) {
                y = ((FLOAT)log(Y[i]) - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]); y *= y;

                if (Outlier) {
                    *Outlier |= (FLOAT)2.0 * y > ChiSqr_;
                }

                *CmpDist += -y - LogSqrtPi2 - (FLOAT)log(CmpTheta->Theta_[1][i]) - (FLOAT)log(Y[i]);
            }
            else {
                *CmpDist = -FLOAT_MAX; goto E0;
            }

            break;
        case pfWeibull:
            if (Y[i] > FLOAT_MIN) {
                if (Outlier) {
                    y = WeibullInv((FLOAT)1.0 - p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                    *Outlier |= Y[i] > y;

                    y = WeibullInv(p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                    *Outlier |= Y[i] < y;
                }

                ypb = (FLOAT)exp(CmpTheta->Theta_[1][i] * (FLOAT)log(Y[i] / CmpTheta->Theta_[0][i]));

                *CmpDist += (FLOAT)log(CmpTheta->Theta_[1][i]) + (FLOAT)log(ypb) - ypb - (FLOAT)log(Y[i]);
            }
            else {
                *CmpDist = -FLOAT_MAX; goto E0;
            }

            break;
        case pfGamma:
            if (Y[i] > FLOAT_MIN) {
                if (Outlier) {
                    Error = GammaInv((FLOAT)1.0 - p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i], &y);

                    if (Error) goto E0;

                    *Outlier |= Y[i] > y;

                    Error = GammaInv(p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i], &y);

                    if (Error) goto E0;

                    *Outlier |= Y[i] < y;
                }

                ypb = Y[i] / CmpTheta->Theta_[0][i];

                *CmpDist += CmpTheta->Theta_[1][i] * (FLOAT)log(ypb) - ypb - Gammaln(CmpTheta->Theta_[1][i]) - (FLOAT)log(Y[i]);
            }
            else {
                *CmpDist = -FLOAT_MAX; goto E0;
            }

            break;
        case pfvonMises:
            if (Outlier) {
                y = vonMisesInv((FLOAT)1.0 - p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                *Outlier |= Y[i] > y;

                y = vonMisesInv(p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                *Outlier |= Y[i] < y;
            }

            if ((Y[i] < (FLOAT)0.0) || (Y[i] > Pi2)) {
                *CmpDist = -FLOAT_MAX; goto E0;
            }
            else {
                *CmpDist += CmpTheta->Theta_[1][i] * (FLOAT)cos(Y[i] - CmpTheta->Theta_[0][i]) - LogPi2 - (FLOAT)log(BesselI0(CmpTheta->Theta_[1][i]));
            }

            break;
        case pfBinomial:
            if (Outlier) {
                y = BinomialInv((FLOAT)1.0 - p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                *Outlier |= Y[i] > y;

                y = BinomialInv(p_value_, CmpTheta->Theta_[0][i], CmpTheta->Theta_[1][i]);

                *Outlier |= Y[i] < y;
            }

            j = (int)Y[i]; n = (int)CmpTheta->Theta_[0][i]; p = CmpTheta->Theta_[1][i];

            if (j < 0) {
                *CmpDist = -FLOAT_MAX; goto E0;
            }
            else
            if (j == 0)
                *CmpDist += n * (FLOAT)log((FLOAT)1.0 - p);
            else
            if (j == n)
                *CmpDist += n * (FLOAT)log(p);
            else
            if (j > n) {
                *CmpDist = -FLOAT_MAX; goto E0;
            }
            else
                *CmpDist += Gammaln(n + (FLOAT)1.0) - Gammaln(j + (FLOAT)1.0) - Gammaln(n - j + (FLOAT)1.0) +
                            j * (FLOAT)log(p) + (n - j) * (FLOAT)log((FLOAT)1.0 - p);

            break;
        case pfPoisson:
            if (Outlier) {
                y = PoissonInv((FLOAT)1.0 - p_value_, CmpTheta->Theta_[0][i]);

                *Outlier |= Y[i] > y;

                y = PoissonInv(p_value_, CmpTheta->Theta_[0][i]);

                *Outlier |= Y[i] < y;
            }

            j = (int)Y[i]; Theta = CmpTheta->Theta_[0][i];

            *CmpDist += j * (FLOAT)log(Theta) - Theta - Gammaln(j + (FLOAT)1.0);

            break;
        case pfDirac:
            if ((FLOAT)fabs(Y[i] - CmpTheta->Theta_[0][i]) > FLOAT_MIN) {
                *CmpDist = -FLOAT_MAX; goto E0;
            }
            else {
                *CmpDist += (FLOAT)0.0;
            }

            break;
        case pfUniform:
            if ((Y[i] > CmpTheta->Theta_[1][i]) || (Y[i] < CmpTheta->Theta_[0][i])) {
                *CmpDist = -FLOAT_MAX; goto E0;
            }
            else {
                *CmpDist -= (FLOAT)log(CmpTheta->Theta_[1][i] - CmpTheta->Theta_[0][i]);
            }
        }
    }

E0: return Error;
} // LogComponentDist

// Enhanced component parameter estimation for k-nearest neighbours.

int Rebmix::EnhancedEstimationKNN(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,logV,R].
                                  FLOAT                nl,          // Total number of observations in class l.
                                  CompnentDistribution *RigidTheta, // Rigid parameters.
                                  CompnentDistribution *LooseTheta) // Loose parameters.
{
    CompnentDistribution *EnhanTheta = NULL;
    FLOAT                A[4], T[2];
    int                  i, j, l;
    FLOAT                dP, MrgVar, TmpVar;
    int                  Error = 0;

    EnhanTheta = new CompnentDistribution(this);

    Error = NULL == EnhanTheta; if (Error) goto E0;

    Error = EnhanTheta->Realloc(length_pdf_, length_Theta_, length_theta_);

    if (Error) goto E0;

    for (i = 0; i < length_pdf_; i++) {
        switch (RigidTheta->pdf_[i]) {
        case pfNormal:
            EnhanTheta->pdf_[i] = pfNormal;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] /= nl;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                T[0] = Y[j][i] - EnhanTheta->Theta_[0][i];

                EnhanTheta->Theta_[1][i] += Y[j][length_pdf_] * T[0] * T[0];
            }

            EnhanTheta->Theta_[1][i] /= nl;

            if (EnhanTheta->Theta_[1][i] <= FLOAT_MIN) {
                Error = 1; goto E0;
            }

            EnhanTheta->Theta_[1][i] = (FLOAT)sqrt(EnhanTheta->Theta_[1][i]);

            TmpVar = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[1][i];
            MrgVar = RigidTheta->Theta_[1][i] * RigidTheta->Theta_[1][i];

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfLognormal:
            EnhanTheta->pdf_[i] = pfLognormal;

            for (j = 0; j < n_; j++) {
                if ((Y[j][length_pdf_] > FLOAT_MIN) && (Y[j][i] > FLOAT_MIN)) {
                    T[0] = Y[j][length_pdf_] * (FLOAT)log(Y[j][i]);

                    EnhanTheta->Theta_[0][i] += T[0];
                    EnhanTheta->Theta_[1][i] += T[0] * (FLOAT)log(Y[j][i]);
                }
            }

            EnhanTheta->Theta_[0][i] /= nl;
            EnhanTheta->Theta_[1][i] = EnhanTheta->Theta_[1][i] / nl - EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[0][i];

            if (EnhanTheta->Theta_[1][i] <= FLOAT_MIN) {
                Error = 1; goto E0;
            }

            EnhanTheta->Theta_[1][i] = (FLOAT)sqrt(EnhanTheta->Theta_[1][i]);

            TmpVar = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[1][i];
            MrgVar = RigidTheta->Theta_[1][i] * RigidTheta->Theta_[1][i];

            TmpVar = ((FLOAT)exp(TmpVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * EnhanTheta->Theta_[0][i] + TmpVar);
            MrgVar = ((FLOAT)exp(MrgVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * RigidTheta->Theta_[0][i] + MrgVar);

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfWeibull:
            EnhanTheta->pdf_[i] = pfWeibull;

            EnhanTheta->Theta_[1][i] = RigidTheta->Theta_[1][i];

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                memset(&A, 0, 4 * sizeof(FLOAT));

                for (l = 0; l < n_; l++) {
                    if ((Y[l][length_pdf_] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                        T[0] = (FLOAT)log(Y[l][i]);
                        T[1] = (FLOAT)exp(T[0] * EnhanTheta->Theta_[1][i]);

                        A[0] += Y[l][length_pdf_] * T[0];
                        A[1] += Y[l][length_pdf_] * T[1] * T[0];
                        A[2] += Y[l][length_pdf_] * T[1];
                        A[3] += Y[l][length_pdf_] * T[1] * T[0] * T[0];
                    }
                }

                A[0] /= nl; T[0] = A[1] / A[2]; T[0] *= T[0]; T[1] = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[1][i];

                dP = ((FLOAT)1.0 / EnhanTheta->Theta_[1][i] + A[0] - A[1] / A[2]) / (T[0] - A[3] / A[2] - (FLOAT)1.0 / T[1]);

                EnhanTheta->Theta_[1][i] -= dP;

                if (IsNan(dP) || IsInf(dP)) {
                    Error = 1; goto E0;
                }

                if ((FLOAT)fabs(dP) < Max(Eps * (FLOAT)fabs(EnhanTheta->Theta_[1][i]), Eps)) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            EnhanTheta->Theta_[0][i] = (FLOAT)exp((FLOAT)log(A[2]) / EnhanTheta->Theta_[1][i]);

            TmpVar = EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[0][i];
            MrgVar = RigidTheta->Theta_[0][i] * RigidTheta->Theta_[0][i];

            TmpVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / EnhanTheta->Theta_[1][i])) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / EnhanTheta->Theta_[1][i]));
            MrgVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / RigidTheta->Theta_[1][i])) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / RigidTheta->Theta_[1][i]));

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfGamma:
            EnhanTheta->pdf_[i] = pfGamma;

            EnhanTheta->Theta_[1][i] = RigidTheta->Theta_[1][i];

            memset(&A, 0, 2 * sizeof(FLOAT));

            for (l = 0; l < n_; l++) {
                if ((Y[l][length_pdf_] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                    A[0] += Y[l][length_pdf_] * Y[l][i];
                    A[1] += Y[l][length_pdf_] * (FLOAT)log(Y[l][i]);
                }
            }

            A[0] /= nl; A[1] /= nl;

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                if (Digamma(EnhanTheta->Theta_[1][i], &T[0]) || Digamma(EnhanTheta->Theta_[1][i] + Eps, &T[1])) goto E0;

                dP = ((FLOAT)log(EnhanTheta->Theta_[1][i]) - T[0] - (FLOAT)log(A[0]) + A[1]) / ((FLOAT)1.0 / EnhanTheta->Theta_[1][i] - (T[1] - T[0]) / Eps);

                EnhanTheta->Theta_[1][i] -= dP;

                if (IsNan(dP) || IsInf(dP)) {
                    Error = 1; goto E0;
                }

                if ((FLOAT)fabs(dP) < Max(Eps * (FLOAT)fabs(EnhanTheta->Theta_[1][i]), Eps)) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            EnhanTheta->Theta_[0][i] = A[0] / EnhanTheta->Theta_[1][i];

            TmpVar = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[0][i];
            MrgVar = RigidTheta->Theta_[1][i] * RigidTheta->Theta_[0][i] * RigidTheta->Theta_[0][i];

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfvonMises:
            EnhanTheta->pdf_[i] = pfvonMises;

            memset(&A, 0, 3 * sizeof(FLOAT));

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                A[0] += Y[j][length_pdf_] * (FLOAT)cos(Y[j][i]);
                A[1] += Y[j][length_pdf_] * (FLOAT)sin(Y[j][i]);
            }

            A[0] /= nl; A[1] /= nl; A[2] = (FLOAT)sqrt((FLOAT)pow(A[0], (FLOAT)2.0) + (FLOAT)pow(A[1], (FLOAT)2.0));

            if (A[1] > FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] = (FLOAT)2.0 * (FLOAT)atan((A[2] - A[0]) / A[1]);
            }
            else
            if (A[1] < -FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] = (FLOAT)2.0 * (FLOAT)atan((A[2] - A[0]) / A[1]) + Pi2;
            }
            else
            if (A[0] > FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] = (FLOAT)0.0;
            }
            else
            if (A[0] < -FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] = Pi;
            }
            else {
                Error = 1; goto E0;
            }

            EnhanTheta->Theta_[1][i] = RigidTheta->Theta_[1][i];

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                A[0] = BesselI0(EnhanTheta->Theta_[1][i]); A[1] = BesselI1(EnhanTheta->Theta_[1][i]);

                dP = (A[1] - A[2] * A[0]) / (A[0] - (A[2] + (FLOAT)1.0 / EnhanTheta->Theta_[1][i]) * A[1]);

                if (IsNan(dP) || IsInf(dP)) {
                    Error = 1; goto E0;
                }

                EnhanTheta->Theta_[1][i] -= dP;

                if ((FLOAT)fabs(dP) < Max(Eps * (FLOAT)fabs(EnhanTheta->Theta_[1][i]), Eps)) Error = 0;

                j++;
            }

            if (Error) goto E0;

            TmpVar = (FLOAT)1.0 - BesselI1(EnhanTheta->Theta_[1][i]) / BesselI0(EnhanTheta->Theta_[1][i]);
            MrgVar = (FLOAT)1.0 - BesselI1(RigidTheta->Theta_[1][i]) / BesselI0(RigidTheta->Theta_[1][i]);

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfBinomial:
            EnhanTheta->pdf_[i] = pfBinomial;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                T[0] += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] = RigidTheta->Theta_[0][i];

            EnhanTheta->Theta_[1][i] = T[0] / EnhanTheta->Theta_[0][i] / nl;

            if ((EnhanTheta->Theta_[0][i] < (FLOAT)0.0) || (EnhanTheta->Theta_[1][i] < (FLOAT)0.0) || (EnhanTheta->Theta_[1][i] >(FLOAT)1.0)) {
                Error = 1; goto E0;
            }

            TmpVar = EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[1][i] * ((FLOAT)1.0 - EnhanTheta->Theta_[1][i]);
            MrgVar = RigidTheta->Theta_[0][i] * RigidTheta->Theta_[1][i] * ((FLOAT)1.0 - RigidTheta->Theta_[1][i]);

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfPoisson:
            EnhanTheta->pdf_[i] = pfPoisson;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                T[0] += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] = T[0] / nl;

            EnhanTheta->Theta_[1][i] = (FLOAT)0.0;

            if (EnhanTheta->Theta_[0][i] < (FLOAT)0.0) {
                Error = 1; goto E0;
            }

            TmpVar = EnhanTheta->Theta_[0][i];
            MrgVar = RigidTheta->Theta_[0][i];

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfDirac:
            EnhanTheta->pdf_[i] = pfDirac;

            EnhanTheta->Theta_[0][i] = RigidTheta->Theta_[0][i];

            break;
        case pfUniform:
            Error = 1; goto E0;
        }
    }

    Error = LooseTheta->Memmove(EnhanTheta);

    if (Error) goto E0;

E0: if (EnhanTheta) delete EnhanTheta;

    return Error;
} // EnhancedEstimationKNN

// Enhanced component parameter estimation for Parzen window.

int Rebmix::EnhancedEstimationPW(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,k].
                                 FLOAT                nl,          // Total number of observations in class l.
                                 CompnentDistribution *RigidTheta, // Rigid parameters.
                                 CompnentDistribution *LooseTheta) // Loose parameters.
{
    CompnentDistribution *EnhanTheta = NULL;
    FLOAT                A[4], T[2];
    int                  i, j, l;
    FLOAT                dP, MrgVar, TmpVar;
    int                  Error = 0;

    EnhanTheta = new CompnentDistribution(this);

    Error = NULL == EnhanTheta; if (Error) goto E0;

    Error = EnhanTheta->Realloc(length_pdf_, length_Theta_, length_theta_);

    if (Error) goto E0;

    for (i = 0; i < length_pdf_; i++) {
        switch (RigidTheta->pdf_[i]) {
        case pfNormal:
            EnhanTheta->pdf_[i] = pfNormal;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] /= nl;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                T[0] = Y[j][i] - EnhanTheta->Theta_[0][i];

                EnhanTheta->Theta_[1][i] += Y[j][length_pdf_] * T[0] * T[0];
            }

            EnhanTheta->Theta_[1][i] /= nl;

            if (EnhanTheta->Theta_[1][i] <= FLOAT_MIN) {
                Error = 1; goto E0;
            }

            EnhanTheta->Theta_[1][i] = (FLOAT)sqrt(EnhanTheta->Theta_[1][i]);

            TmpVar = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[1][i];
            MrgVar = RigidTheta->Theta_[1][i] * RigidTheta->Theta_[1][i];

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfLognormal:
            EnhanTheta->pdf_[i] = pfLognormal;

            for (j = 0; j < n_; j++) {
                if ((Y[j][length_pdf_] > FLOAT_MIN) && (Y[j][i] > FLOAT_MIN)) {
                    T[0] = Y[j][length_pdf_] * (FLOAT)log(Y[j][i]);

                    EnhanTheta->Theta_[0][i] += T[0];
                    EnhanTheta->Theta_[1][i] += T[0] * (FLOAT)log(Y[j][i]);
                }
            }

            EnhanTheta->Theta_[0][i] /= nl;
            EnhanTheta->Theta_[1][i] = EnhanTheta->Theta_[1][i] / nl - EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[0][i];

            if (EnhanTheta->Theta_[1][i] <= FLOAT_MIN) {
                Error = 1; goto E0;
            }

            EnhanTheta->Theta_[1][i] = (FLOAT)sqrt(EnhanTheta->Theta_[1][i]);

            TmpVar = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[1][i];
            MrgVar = RigidTheta->Theta_[1][i] * RigidTheta->Theta_[1][i];

            TmpVar = ((FLOAT)exp(TmpVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * EnhanTheta->Theta_[0][i] + TmpVar);
            MrgVar = ((FLOAT)exp(MrgVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * RigidTheta->Theta_[0][i] + MrgVar);

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfWeibull:
            EnhanTheta->pdf_[i] = pfWeibull;

            EnhanTheta->Theta_[1][i] = RigidTheta->Theta_[1][i];

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                memset(&A, 0, 4 * sizeof(FLOAT));

                for (l = 0; l < n_; l++) {
                    if ((Y[l][length_pdf_] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                        T[0] = (FLOAT)log(Y[l][i]);
                        T[1] = (FLOAT)exp(T[0] * EnhanTheta->Theta_[1][i]);

                        A[0] += Y[l][length_pdf_] * T[0];
                        A[1] += Y[l][length_pdf_] * T[1] * T[0];
                        A[2] += Y[l][length_pdf_] * T[1];
                        A[3] += Y[l][length_pdf_] * T[1] * T[0] * T[0];
                    }
                }

                A[0] /= nl; T[0] = A[1] / A[2]; T[0] *= T[0]; T[1] = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[1][i];

                dP = ((FLOAT)1.0 / EnhanTheta->Theta_[1][i] + A[0] - A[1] / A[2]) / (T[0] - A[3] / A[2] - (FLOAT)1.0 / T[1]);

                EnhanTheta->Theta_[1][i] -= dP;

                if (IsNan(dP) || IsInf(dP)) {
                    Error = 1; goto E0;
                }

                if ((FLOAT)fabs(dP) < Max(Eps * (FLOAT)fabs(EnhanTheta->Theta_[1][i]), Eps)) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            EnhanTheta->Theta_[0][i] = (FLOAT)exp((FLOAT)log(A[2]) / EnhanTheta->Theta_[1][i]);

            TmpVar = EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[0][i];
            MrgVar = RigidTheta->Theta_[0][i] * RigidTheta->Theta_[0][i];

            TmpVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / EnhanTheta->Theta_[1][i])) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / EnhanTheta->Theta_[1][i]));
            MrgVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / RigidTheta->Theta_[1][i])) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / RigidTheta->Theta_[1][i]));

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfGamma:
            EnhanTheta->pdf_[i] = pfGamma;

            EnhanTheta->Theta_[1][i] = RigidTheta->Theta_[1][i];

            memset(&A, 0, 2 * sizeof(FLOAT));

            for (l = 0; l < n_; l++) {
                if ((Y[l][length_pdf_] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                    A[0] += Y[l][length_pdf_] * Y[l][i];
                    A[1] += Y[l][length_pdf_] * (FLOAT)log(Y[l][i]);
                }
            }

            A[0] /= nl; A[1] /= nl;

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                if (Digamma(EnhanTheta->Theta_[1][i], &T[0]) || Digamma(EnhanTheta->Theta_[1][i] + Eps, &T[1])) goto E0;

                dP = ((FLOAT)log(EnhanTheta->Theta_[1][i]) - T[0] - (FLOAT)log(A[0]) + A[1]) / ((FLOAT)1.0 / EnhanTheta->Theta_[1][i] - (T[1] - T[0]) / Eps);

                EnhanTheta->Theta_[1][i] -= dP;

                if (IsNan(dP) || IsInf(dP)) {
                    Error = 1; goto E0;
                }

                if ((FLOAT)fabs(dP) < Max(Eps * (FLOAT)fabs(EnhanTheta->Theta_[1][i]), Eps)) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            EnhanTheta->Theta_[0][i] = A[0] / EnhanTheta->Theta_[1][i];

            TmpVar = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[0][i];
            MrgVar = RigidTheta->Theta_[1][i] * RigidTheta->Theta_[0][i] * RigidTheta->Theta_[0][i];

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfvonMises:
            EnhanTheta->pdf_[i] = pfvonMises;

            memset(&A, 0, 3 * sizeof(FLOAT));

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                A[0] += Y[j][length_pdf_] * (FLOAT)cos(Y[j][i]);
                A[1] += Y[j][length_pdf_] * (FLOAT)sin(Y[j][i]);
            }

            A[0] /= nl; A[1] /= nl; A[2] = (FLOAT)sqrt((FLOAT)pow(A[0], (FLOAT)2.0) + (FLOAT)pow(A[1], (FLOAT)2.0));

            if (A[1] > FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] = (FLOAT)2.0 * (FLOAT)atan((A[2] - A[0]) / A[1]);
            }
            else
            if (A[1] < -FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] = (FLOAT)2.0 * (FLOAT)atan((A[2] - A[0]) / A[1]) + Pi2;
            }
            else
            if (A[0] > FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] = (FLOAT)0.0;
            }
            else
            if (A[0] < -FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] = Pi;
            }
            else {
                Error = 1; goto E0;
            }

            EnhanTheta->Theta_[1][i] = RigidTheta->Theta_[1][i];

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                A[0] = BesselI0(EnhanTheta->Theta_[1][i]); A[1] = BesselI1(EnhanTheta->Theta_[1][i]);

                dP = (A[1] - A[2] * A[0]) / (A[0] - (A[2] + (FLOAT)1.0 / EnhanTheta->Theta_[1][i]) * A[1]);

                if (IsNan(dP) || IsInf(dP)) {
                    Error = 1; goto E0;
                }

                EnhanTheta->Theta_[1][i] -= dP;

                if ((FLOAT)fabs(dP) < Max(Eps * (FLOAT)fabs(EnhanTheta->Theta_[1][i]), Eps)) Error = 0;

                j++;
            }

            if (Error) goto E0;

            TmpVar = (FLOAT)1.0 - BesselI1(EnhanTheta->Theta_[1][i]) / BesselI0(EnhanTheta->Theta_[1][i]);
            MrgVar = (FLOAT)1.0 - BesselI1(RigidTheta->Theta_[1][i]) / BesselI0(RigidTheta->Theta_[1][i]);

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfBinomial:
            EnhanTheta->pdf_[i] = pfBinomial;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                T[0] += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] = RigidTheta->Theta_[0][i];

            EnhanTheta->Theta_[1][i] = T[0] / EnhanTheta->Theta_[0][i] / nl;

            if ((EnhanTheta->Theta_[0][i] < (FLOAT)0.0) || (EnhanTheta->Theta_[1][i] < (FLOAT)0.0) || (EnhanTheta->Theta_[1][i] >(FLOAT)1.0)) {
                Error = 1; goto E0;
            }

            TmpVar = EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[1][i] * ((FLOAT)1.0 - EnhanTheta->Theta_[1][i]);
            MrgVar = RigidTheta->Theta_[0][i] * RigidTheta->Theta_[1][i] * ((FLOAT)1.0 - RigidTheta->Theta_[1][i]);

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfPoisson:
            EnhanTheta->pdf_[i] = pfPoisson;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                T[0] += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] = T[0] / nl;

            EnhanTheta->Theta_[1][i] = (FLOAT)0.0;

            if (EnhanTheta->Theta_[0][i] < (FLOAT)0.0) {
                Error = 1; goto E0;
            }

            TmpVar = EnhanTheta->Theta_[0][i];
            MrgVar = RigidTheta->Theta_[0][i];

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfDirac:
            EnhanTheta->pdf_[i] = pfDirac;

            EnhanTheta->Theta_[0][i] = RigidTheta->Theta_[0][i];

            break;
        case pfUniform:
            Error = 1; goto E0;
        }
    }

    Error = LooseTheta->Memmove(EnhanTheta);

    if (Error) goto E0;

E0: if (EnhanTheta) delete EnhanTheta;

    return Error;
} // EnhancedEstimationPW

// Enhanced component parameter estimation for histogram.

int Rebmix::EnhancedEstimationH(int                  k,           // Total number of bins.
                                FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,k].
                                FLOAT                nl,          // Total number of observations in class l.
                                CompnentDistribution *RigidTheta, // Rigid parameters.
                                CompnentDistribution *LooseTheta) // Loose parameters.
{
    CompnentDistribution *EnhanTheta = NULL;
    FLOAT                A[4], T[2];
    int                  i, j, l;
    FLOAT                dP, MrgVar, TmpVar;
    int                  Error = 0;

    EnhanTheta = new CompnentDistribution(this);

    Error = NULL == EnhanTheta; if (Error) goto E0;

    Error = EnhanTheta->Realloc(length_pdf_, length_Theta_, length_theta_);

    if (Error) goto E0;

    for (i = 0; i < length_pdf_; i++) {
        switch (RigidTheta->pdf_[i]) {
        case pfNormal:
            EnhanTheta->pdf_[i] = pfNormal;

            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] /= nl;

            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                T[0] = Y[j][i] - EnhanTheta->Theta_[0][i];

                EnhanTheta->Theta_[1][i] += Y[j][length_pdf_] * T[0] * T[0];
            }

            EnhanTheta->Theta_[1][i] /= nl;

            if (EnhanTheta->Theta_[1][i] <= FLOAT_MIN) {
                Error = 1; goto E0;
            }

            EnhanTheta->Theta_[1][i] = (FLOAT)sqrt(EnhanTheta->Theta_[1][i]);

            TmpVar = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[1][i];
            MrgVar = RigidTheta->Theta_[1][i] * RigidTheta->Theta_[1][i];

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfLognormal:
            EnhanTheta->pdf_[i] = pfLognormal;

            for (j = 0; j < k; j++) {
                if ((Y[j][length_pdf_] > FLOAT_MIN) && (Y[j][i] > FLOAT_MIN)) {
                    T[0] = Y[j][length_pdf_] * (FLOAT)log(Y[j][i]);

                    EnhanTheta->Theta_[0][i] += T[0];
                    EnhanTheta->Theta_[1][i] += T[0] * (FLOAT)log(Y[j][i]);
                }
            }

            EnhanTheta->Theta_[0][i] /= nl;
            EnhanTheta->Theta_[1][i] = EnhanTheta->Theta_[1][i] / nl - EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[0][i];

            if (EnhanTheta->Theta_[1][i] <= FLOAT_MIN) {
                Error = 1; goto E0;
            }

            EnhanTheta->Theta_[1][i] = (FLOAT)sqrt(EnhanTheta->Theta_[1][i]);

            TmpVar = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[1][i];
            MrgVar = RigidTheta->Theta_[1][i] * RigidTheta->Theta_[1][i];

            TmpVar = ((FLOAT)exp(TmpVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * EnhanTheta->Theta_[0][i] + TmpVar);
            MrgVar = ((FLOAT)exp(MrgVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * RigidTheta->Theta_[0][i] + MrgVar);

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfWeibull:
            EnhanTheta->pdf_[i] = pfWeibull;

            EnhanTheta->Theta_[1][i] = RigidTheta->Theta_[1][i];

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                memset(&A, 0, 4 * sizeof(FLOAT));

                for (l = 0; l < k; l++) {
                    if ((Y[l][length_pdf_] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                        T[0] = (FLOAT)log(Y[l][i]);
                        T[1] = (FLOAT)exp(T[0] * EnhanTheta->Theta_[1][i]);

                        A[0] += Y[l][length_pdf_] * T[0];
                        A[1] += Y[l][length_pdf_] * T[1] * T[0];
                        A[2] += Y[l][length_pdf_] * T[1];
                        A[3] += Y[l][length_pdf_] * T[1] * T[0] * T[0];
                    }
                }

                A[0] /= nl; T[0] = A[1] / A[2]; T[0] *= T[0]; T[1] = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[1][i];

                dP = ((FLOAT)1.0 / EnhanTheta->Theta_[1][i] + A[0] - A[1] / A[2]) / (T[0] - A[3] / A[2] - (FLOAT)1.0 / T[1]);

                EnhanTheta->Theta_[1][i] -= dP;

                if (IsNan(dP) || IsInf(dP)) {
                    Error = 1; goto E0;
                }

                if ((FLOAT)fabs(dP) < Max(Eps * (FLOAT)fabs(EnhanTheta->Theta_[1][i]), Eps)) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            EnhanTheta->Theta_[0][i] = (FLOAT)exp((FLOAT)log(A[2]) / EnhanTheta->Theta_[1][i]);

            TmpVar = EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[0][i];
            MrgVar = RigidTheta->Theta_[0][i] * RigidTheta->Theta_[0][i];

            TmpVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / EnhanTheta->Theta_[1][i])) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / EnhanTheta->Theta_[1][i]));
            MrgVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / RigidTheta->Theta_[1][i])) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / RigidTheta->Theta_[1][i]));

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfGamma:
            EnhanTheta->pdf_[i] = pfGamma;

            EnhanTheta->Theta_[1][i] = RigidTheta->Theta_[1][i];

            memset(&A, 0, 2 * sizeof(FLOAT));

            for (l = 0; l < k; l++) {
                if ((Y[l][length_pdf_] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                    A[0] += Y[l][length_pdf_] * Y[l][i];
                    A[1] += Y[l][length_pdf_] * (FLOAT)log(Y[l][i]);
                }
            }

            A[0] /= nl; A[1] /= nl;

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                if (Digamma(EnhanTheta->Theta_[1][i], &T[0]) || Digamma(EnhanTheta->Theta_[1][i] + Eps, &T[1])) goto E0;

                dP = ((FLOAT)log(EnhanTheta->Theta_[1][i]) - T[0] - (FLOAT)log(A[0]) + A[1]) / ((FLOAT)1.0 / EnhanTheta->Theta_[1][i] - (T[1] - T[0]) / Eps);

                EnhanTheta->Theta_[1][i] -= dP;

                if (IsNan(dP) || IsInf(dP)) {
                    Error = 1; goto E0;
                }

                if ((FLOAT)fabs(dP) < Max(Eps * (FLOAT)fabs(EnhanTheta->Theta_[1][i]), Eps)) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            EnhanTheta->Theta_[0][i] = A[0] / EnhanTheta->Theta_[1][i];

            TmpVar = EnhanTheta->Theta_[1][i] * EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[0][i];
            MrgVar = RigidTheta->Theta_[1][i] * RigidTheta->Theta_[0][i] * RigidTheta->Theta_[0][i];

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfvonMises:
            EnhanTheta->pdf_[i] = pfvonMises;

            memset(&A, 0, 3 * sizeof(FLOAT));

            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                A[0] += Y[j][length_pdf_] * (FLOAT)cos(Y[j][i]);
                A[1] += Y[j][length_pdf_] * (FLOAT)sin(Y[j][i]);
            }

            A[0] /= nl; A[1] /= nl; A[2] = (FLOAT)sqrt((FLOAT)pow(A[0], (FLOAT)2.0) + (FLOAT)pow(A[1], (FLOAT)2.0));

            if (A[1] > FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] = (FLOAT)2.0 * (FLOAT)atan((A[2] - A[0]) / A[1]);
            }
            else
            if (A[1] < -FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] = (FLOAT)2.0 * (FLOAT)atan((A[2] - A[0]) / A[1]) + Pi2;
            }
            else
            if (A[0] > FLOAT_MIN) { 
                EnhanTheta->Theta_[0][i] = (FLOAT)0.0;
            }
            else
            if (A[0] < -FLOAT_MIN) {
                EnhanTheta->Theta_[0][i] = Pi;
            }
            else {
                Error = 1; goto E0;
            }
            
            EnhanTheta->Theta_[1][i] = RigidTheta->Theta_[1][i];

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                A[0] = BesselI0(EnhanTheta->Theta_[1][i]); A[1] = BesselI1(EnhanTheta->Theta_[1][i]);

                dP = (A[1] - A[2] * A[0]) / (A[0] - (A[2] + (FLOAT)1.0 / EnhanTheta->Theta_[1][i]) * A[1]);

                if (IsNan(dP) || IsInf(dP)) {
                    Error = 1; goto E0;
                }

                EnhanTheta->Theta_[1][i] -= dP;

                if ((FLOAT)fabs(dP) < Max(Eps * (FLOAT)fabs(EnhanTheta->Theta_[1][i]), Eps)) Error = 0;

                j++;
            }

            if (Error) goto E0;

            TmpVar = (FLOAT)1.0 - BesselI1(EnhanTheta->Theta_[1][i]) / BesselI0(EnhanTheta->Theta_[1][i]);
            MrgVar = (FLOAT)1.0 - BesselI1(RigidTheta->Theta_[1][i]) / BesselI0(RigidTheta->Theta_[1][i]);

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfBinomial:
            EnhanTheta->pdf_[i] = pfBinomial;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                T[0] += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] = RigidTheta->Theta_[0][i];

            EnhanTheta->Theta_[1][i] = T[0] / EnhanTheta->Theta_[0][i] / nl;

            if ((EnhanTheta->Theta_[0][i] < (FLOAT)0.0) || (EnhanTheta->Theta_[1][i] < (FLOAT)0.0) || (EnhanTheta->Theta_[1][i] >(FLOAT)1.0)) {
                Error = 1; goto E0;
            }

            TmpVar = EnhanTheta->Theta_[0][i] * EnhanTheta->Theta_[1][i] * ((FLOAT)1.0 - EnhanTheta->Theta_[1][i]);
            MrgVar = RigidTheta->Theta_[0][i] * RigidTheta->Theta_[1][i] * ((FLOAT)1.0 - RigidTheta->Theta_[1][i]);

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfPoisson:
            EnhanTheta->pdf_[i] = pfPoisson;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                T[0] += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] = T[0] / nl;

            EnhanTheta->Theta_[1][i] = (FLOAT)0.0;

            if (EnhanTheta->Theta_[0][i] < (FLOAT)0.0) {
                Error = 1; goto E0;
            }

            TmpVar = EnhanTheta->Theta_[0][i];
            MrgVar = RigidTheta->Theta_[0][i];

            if (TmpVar < MrgVar * var_mul_) {
                Error = 1; goto E0;
            }

            break;
        case pfDirac:
            EnhanTheta->pdf_[i] = pfDirac;

            EnhanTheta->Theta_[0][i] = RigidTheta->Theta_[0][i];
  
            break;
        case pfUniform:
            Error = 1; goto E0;
        }
    }

    Error = LooseTheta->Memmove(EnhanTheta);

    if (Error) goto E0;

E0: if (EnhanTheta) delete EnhanTheta;

    return Error;
} // EnhancedEstimationH

// Moments calculation.

int Rebmix::MomentsCalculation(CompnentDistribution *CmpTheta, // Component parameters.
                               FLOAT                *FirstM,   // First moment.
                               FLOAT                *SecondM)  // Second moment.
{
    FLOAT R;
    int   i;
    int   Error = 0;

    for (i = 0; i < length_pdf_; i++) {
        switch (CmpTheta->pdf_[i]) {
        case pfNormal:
            FirstM[i] = CmpTheta->Theta_[0][i];

            SecondM[i] = CmpTheta->Theta_[1][i] * CmpTheta->Theta_[1][i] + CmpTheta->Theta_[0][i] * CmpTheta->Theta_[0][i];

            break;
        case pfLognormal:
            FirstM[i] = (FLOAT)exp(CmpTheta->Theta_[0][i] + (FLOAT)0.5 * CmpTheta->Theta_[1][i] * CmpTheta->Theta_[1][i]);

            SecondM[i] = (FLOAT)exp((FLOAT)2.0 * (CmpTheta->Theta_[0][i] + CmpTheta->Theta_[1][i] * CmpTheta->Theta_[1][i]));

            break;
        case pfWeibull:
            FirstM[i] = CmpTheta->Theta_[0][i] * (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)1.0 / CmpTheta->Theta_[1][i]));

            SecondM[i] = CmpTheta->Theta_[0][i] * CmpTheta->Theta_[0][i] * (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / CmpTheta->Theta_[1][i]));

            break;
        case pfGamma:
            FirstM[i] = CmpTheta->Theta_[1][i] * CmpTheta->Theta_[0][i];

            SecondM[i] = CmpTheta->Theta_[1][i] * CmpTheta->Theta_[0][i] * CmpTheta->Theta_[0][i] * ((FLOAT)1.0 + CmpTheta->Theta_[1][i]);

            break;
        case pfvonMises:
            R = BesselI1(CmpTheta->Theta_[1][i]) / BesselI0(CmpTheta->Theta_[1][i]);

            FirstM[i] = R * (FLOAT)cos(CmpTheta->Theta_[0][i]);
            
            SecondM[i] = R * (FLOAT)sin(CmpTheta->Theta_[0][i]);

            break;
        case pfBinomial:
            FirstM[i] = CmpTheta->Theta_[0][i] * CmpTheta->Theta_[1][i];

            SecondM[i] = (FLOAT)0.0;

            break;
        case pfPoisson:
            FirstM[i] = CmpTheta->Theta_[0][i];

            SecondM[i] = (FLOAT)0.0;

            break;
        case pfDirac:
            FirstM[i] = CmpTheta->Theta_[0][i];

            SecondM[i] = (FLOAT)0.0;
         
            break;
        case pfUniform:
            Error = 1; goto E0;
        }
    }

E0: return Error;
} // MomentsCalculation

// Returns Bayes Weibull parameters.

void BayesWeibullParameters(FLOAT FirstM,  // First moment.
                            FLOAT SecondM, // Second moment.
                            FLOAT *Theta1, // Component parameter.
                            FLOAT *Theta2) // Component parameter.
{
    FLOAT A;
    FLOAT xl, xm, xh, fl, fm, fh;
    FLOAT i;
    int   Error = 0, Stop;

    A = (FLOAT)log(SecondM / FirstM / FirstM); xl = 0.001; xh = (FLOAT)10.0;

    fl = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xl) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xl);
    fh = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xh) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xh);

    i = 1; Error = 1;
    while ((i <= ItMax) && Error) {
        if (fl * fh < (FLOAT)0.0)
            Error = 0;
        else
        if ((FLOAT)fabs(fl) < (FLOAT)fabs(fh)) {
            xl += (FLOAT)1.6 * (xl - xh);
            fl = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xl) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xl);
        }
        else {
            xh += (FLOAT)1.6 * (xh - xl);
            fh = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xh) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xh);
        }

        i++;
    }

    if (Error) goto E0;

    // Bisection.

    Stop = 0;

    while (!Stop) {
        xm = (xh + xl) / (FLOAT)2.0;

        fm = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xm) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xm);

        if (((FLOAT)fabs(fm) < Eps) || (xh - xl < Eps)) {
            Stop = 1;
        }
        else {
            if (fm * fl > (FLOAT)0.0) {
                xl = xm; fl = fm;
            }
            else {
                xh = xm;
            }
        }
    }

    *Theta2 = xm;

    *Theta1 = FirstM / (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xm));

E0: return;
} // BayesWeibullParameters

// Returns Bayes von Mises parameters.

void BayesvonMisesParameters(FLOAT FirstM,  // First moment.
                             FLOAT SecondM, // Second moment.
                             FLOAT *Theta1, // Component parameter.
                             FLOAT *Theta2) // Component parameter.
{
    FLOAT A[3];
    FLOAT dP, theta1, theta2;
    FLOAT i;
    int   Error = 0;

    A[2] = (FLOAT)sqrt((FLOAT)pow(FirstM, (FLOAT)2.0) + (FLOAT)pow(SecondM, (FLOAT)2.0));

    if (SecondM > FLOAT_MIN) {
        theta1 = (FLOAT)2.0 * (FLOAT)atan((A[2] - FirstM) / SecondM);
    }
    else
    if (SecondM < -FLOAT_MIN) {
        theta1 = (FLOAT)2.0 * (FLOAT)atan((A[2] - FirstM) / SecondM) + Pi2;
    }
    else
    if (FirstM > FLOAT_MIN) {
        theta1 = (FLOAT)0.0;
    }
    else
    if (FirstM < -FLOAT_MIN) {
        theta1 = Pi;
    }
    else {
        Error = 1; goto E0;
    }

    theta2 = *Theta2;

    i = 1; Error = 1;
    while ((i <= ItMax) && Error) {
        A[0] = BesselI0(theta2); A[1] = BesselI1(theta2);

        dP = (A[1] - A[2] * A[0]) / (A[0] - (A[2] + (FLOAT)1.0 / theta2) * A[1]);

        if (IsNan(dP) || IsInf(dP)) {
            Error = 1; goto E0;
        }

        theta2 -= dP;

        if ((FLOAT)fabs(dP) < Max(Eps * (FLOAT)fabs(theta2), Eps)) Error = 0;

        i++;
    }

    if (Error) goto E0;
    
    *Theta1 = theta1; *Theta2 = theta2;

E0: return;
} // BayesvonMisesParameters

// Bayes classification of the remaining observations for k-nearest neighbour.

int Rebmix::BayesClassificationKNN(FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1].
                                   int                  c,          // Number of components.
                                   FLOAT                *W,         // Component weights.
                                   CompnentDistribution **MixTheta, // Mixture parameters.
                                   FLOAT                **FirstM,   // First moments.
                                   FLOAT                **SecondM)  // Second moments.
{
    int   i, j, l, outlier, Outlier = 0;
    FLOAT CmpDist, Max, Tmp, dW, N = (FLOAT)0.0;
    int   Error = 0;

    for (i = 0; i < n_; i++) {
        if (Y[i][length_pdf_] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(Y[i], MixTheta[l], &CmpDist, &outlier);

            if (Error) goto E0;

            Max = W[l] * CmpDist; Outlier = outlier;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(Y[i], MixTheta[j], &CmpDist, &outlier);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) {
                    l = j; Max = Tmp; Outlier = outlier;
                }
            }

            if (Outlier) {
                N += Y[i][length_pdf_];
            }
            else { 
                dW = Y[i][length_pdf_] / n_; W[l] += dW;

                 for (j = 0; j < length_pdf_; j++) {
                    switch (MixTheta[l]->pdf_[j]) {
                    case pfNormal: case pfLognormal: case pfWeibull: case pfGamma: case pfBinomial: case pfPoisson: case pfDirac: case pfUniform:
                        FirstM[l][j] += dW * (Y[i][j] - FirstM[l][j]) / W[l];

                        SecondM[l][j] += dW * (Y[i][j] * Y[i][j] - SecondM[l][j]) / W[l];

                        break;
                    case pfvonMises:
                        FirstM[l][j] += dW * ((FLOAT)cos(Y[i][j]) - FirstM[l][j]) / W[l];

                        SecondM[l][j] += dW * ((FLOAT)sin(Y[i][j]) - SecondM[l][j]) / W[l];
                    }
                }
            }
        }
    }

    for (i = 0; i < c; i++) {
        W[i] *= n_ / (n_ - N);

        for (j = 0; j < length_pdf_; j++) {
            switch (MixTheta[i]->pdf_[j]) {
            case pfNormal:
                MixTheta[i]->Theta_[0][j] = FirstM[i][j];

                MixTheta[i]->Theta_[1][j] = (FLOAT)sqrt(SecondM[i][j] - MixTheta[i]->Theta_[0][j] * MixTheta[i]->Theta_[0][j]);

                break;
            case pfLognormal:
                MixTheta[i]->Theta_[0][j] = (FLOAT)2.0 * (FLOAT)log(FirstM[i][j]) - (FLOAT)0.5 * (FLOAT)log(SecondM[i][j]);


                MixTheta[i]->Theta_[1][j] = (FLOAT)sqrt((FLOAT)log(SecondM[i][j]) - (FLOAT)2.0 * (FLOAT)log(FirstM[i][j]));

                break;
            case pfWeibull:
                BayesWeibullParameters(FirstM[i][j], SecondM[i][j], &MixTheta[i]->Theta_[0][j], &MixTheta[i]->Theta_[1][j]);

                break;
            case pfGamma:
                MixTheta[i]->Theta_[1][j] = (FLOAT)1.0 / (SecondM[i][j] / FirstM[i][j] / FirstM[i][j] - (FLOAT)1.0);

                MixTheta[i]->Theta_[0][j] = FirstM[i][j] / MixTheta[i]->Theta_[1][j];

                break;
            case pfvonMises:
                BayesvonMisesParameters(FirstM[i][j], SecondM[i][j], &MixTheta[i]->Theta_[0][j], &MixTheta[i]->Theta_[1][j]);

                break;
            case pfBinomial:
                MixTheta[i]->Theta_[1][j] = FirstM[i][j] / MixTheta[i]->Theta_[0][j];

                break;
            case pfPoisson:
                MixTheta[i]->Theta_[0][j] = FirstM[i][j];

                break;
            case pfDirac:
                break;
            case pfUniform:
                Error = 1; goto E0;
            }
        }
    }

E0: return Error;
} // BayesClassificationKNN 

// Bayes classification of the remaining observations for Parzen window.

int Rebmix::BayesClassificationPW(FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1].
                                  int                  c,          // Number of components.
                                  FLOAT                *W,         // Component weights.
                                  CompnentDistribution **MixTheta, // Mixture parameters.
                                  FLOAT                **FirstM,   // First moments.
                                  FLOAT                **SecondM)  // Second moments.
{
    int   i, j, l, outlier, Outlier = 0;
    FLOAT CmpDist, Max, Tmp, dW, N = (FLOAT)0.0;
    int   Error = 0;

    for (i = 0; i < n_; i++) {
        if (Y[i][length_pdf_] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(Y[i], MixTheta[l], &CmpDist, &outlier);

            if (Error) goto E0;

            Max = W[l] * CmpDist; Outlier = outlier;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(Y[i], MixTheta[j], &CmpDist, &outlier);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) {
                    l = j; Max = Tmp; Outlier = outlier;
                }
            }

            if (Outlier) {
                N += Y[i][length_pdf_];
            }
            else { 
                dW = Y[i][length_pdf_] / n_; W[l] += dW;

                for (j = 0; j < length_pdf_; j++) {
                    switch (MixTheta[l]->pdf_[j]) {
                    case pfNormal: case pfLognormal: case pfWeibull: case pfGamma: case pfBinomial: case pfPoisson: case pfDirac: case pfUniform:
                        FirstM[l][j] += dW * (Y[i][j] - FirstM[l][j]) / W[l];

                        SecondM[l][j] += dW * (Y[i][j] * Y[i][j] - SecondM[l][j]) / W[l];

                        break;
                    case pfvonMises:
                        FirstM[l][j] += dW * ((FLOAT)cos(Y[i][j]) - FirstM[l][j]) / W[l];

                        SecondM[l][j] += dW * ((FLOAT)sin(Y[i][j]) - SecondM[l][j]) / W[l];
                    }
                }
            }
        }
    }

    for (i = 0; i < c; i++) {
        W[i] *= n_ / (n_ - N);

        for (j = 0; j < length_pdf_; j++) {
            switch (MixTheta[i]->pdf_[j]) {
            case pfNormal:
                MixTheta[i]->Theta_[0][j] = FirstM[i][j];

                MixTheta[i]->Theta_[1][j] = (FLOAT)sqrt(SecondM[i][j] - MixTheta[i]->Theta_[0][j] * MixTheta[i]->Theta_[0][j]);

                break;
            case pfLognormal:
                MixTheta[i]->Theta_[0][j] = (FLOAT)2.0 * (FLOAT)log(FirstM[i][j]) - (FLOAT)0.5 * (FLOAT)log(SecondM[i][j]);


                MixTheta[i]->Theta_[1][j] = (FLOAT)sqrt((FLOAT)log(SecondM[i][j]) - (FLOAT)2.0 * (FLOAT)log(FirstM[i][j]));

                break;
            case pfWeibull:
                BayesWeibullParameters(FirstM[i][j], SecondM[i][j], &MixTheta[i]->Theta_[0][j], &MixTheta[i]->Theta_[1][j]);

                break;
            case pfGamma:
                MixTheta[i]->Theta_[1][j] = (FLOAT)1.0 / (SecondM[i][j] / FirstM[i][j] / FirstM[i][j] - (FLOAT)1.0);

                MixTheta[i]->Theta_[0][j] = FirstM[i][j] / MixTheta[i]->Theta_[1][j];

                break;
            case pfvonMises:
                BayesvonMisesParameters(FirstM[i][j], SecondM[i][j], &MixTheta[i]->Theta_[0][j], &MixTheta[i]->Theta_[1][j]);

                break;
            case pfBinomial:
                MixTheta[i]->Theta_[1][j] = FirstM[i][j] / MixTheta[i]->Theta_[0][j];

                break;
            case pfPoisson:
                MixTheta[i]->Theta_[0][j] = FirstM[i][j];

                break;
            case pfDirac:
                break;
            case pfUniform:
                Error = 1; goto E0; 
            }
        }
    }

E0: return Error;
} // BayesClassificationPW 

// Bayes classification of the remaining observations for histogram.

int Rebmix::BayesClassificationH(int                  k,          // Total number of bins.
                                 FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1].
                                 int                  c,          // Number of components.
                                 FLOAT                *W,         // Component weights.
                                 CompnentDistribution **MixTheta, // Mixture parameters.
                                 FLOAT                **FirstM,   // First moments.
                                 FLOAT                **SecondM)  // Second moments.
{
    int   i, j, l, outlier, Outlier = 0;
    FLOAT CmpDist, Max, Tmp, dW, N = (FLOAT)0.0;
    int   Error = 0;

    for (i = 0; i < k; i++) {
        if (Y[i][length_pdf_] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(Y[i], MixTheta[l], &CmpDist, &outlier);

            if (Error) goto E0;

            Max = W[l] * CmpDist; Outlier = outlier;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(Y[i], MixTheta[j], &CmpDist, &outlier);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) {
                    l = j; Max = Tmp; Outlier = outlier;
                }
            }

            if (Outlier) {
                N += Y[i][length_pdf_];
            }
            else { 
                dW = Y[i][length_pdf_] / n_; W[l] += dW;

                for (j = 0; j < length_pdf_; j++) {
                    switch (MixTheta[l]->pdf_[j]) {
                    case pfNormal: case pfLognormal: case pfWeibull: case pfGamma: case pfBinomial: case pfPoisson: case pfDirac: case pfUniform:
                        FirstM[l][j] += dW * (Y[i][j] - FirstM[l][j]) / W[l];

                        SecondM[l][j] += dW * (Y[i][j] * Y[i][j] - SecondM[l][j]) / W[l];

                        break;
                    case pfvonMises:
                        FirstM[l][j] += dW * ((FLOAT)cos(Y[i][j]) - FirstM[l][j]) / W[l];

                        SecondM[l][j] += dW * ((FLOAT)sin(Y[i][j]) - SecondM[l][j]) / W[l];
                    }
                }
            }
        }
    }

    for (i = 0; i < c; i++) {
        W[i] *= n_ / (n_ - N);

        for (j = 0; j < length_pdf_; j++) {
            switch (MixTheta[i]->pdf_[j]) {
            case pfNormal:
                MixTheta[i]->Theta_[0][j] = FirstM[i][j];

                MixTheta[i]->Theta_[1][j] = (FLOAT)sqrt(SecondM[i][j] - MixTheta[i]->Theta_[0][j] * MixTheta[i]->Theta_[0][j]);

                break;
            case pfLognormal:
                MixTheta[i]->Theta_[0][j] = (FLOAT)2.0 * (FLOAT)log(FirstM[i][j]) - (FLOAT)0.5 * (FLOAT)log(SecondM[i][j]);


                MixTheta[i]->Theta_[1][j] = (FLOAT)sqrt((FLOAT)log(SecondM[i][j]) - (FLOAT)2.0 * (FLOAT)log(FirstM[i][j]));

                break;
            case pfWeibull:
                BayesWeibullParameters(FirstM[i][j], SecondM[i][j], &MixTheta[i]->Theta_[0][j], &MixTheta[i]->Theta_[1][j]);

                break;
            case pfGamma:
                MixTheta[i]->Theta_[1][j] = (FLOAT)1.0 / (SecondM[i][j] / FirstM[i][j] / FirstM[i][j] - (FLOAT)1.0);

                MixTheta[i]->Theta_[0][j] = FirstM[i][j] / MixTheta[i]->Theta_[1][j];

                break;
            case pfvonMises:
                BayesvonMisesParameters(FirstM[i][j], SecondM[i][j], &MixTheta[i]->Theta_[0][j], &MixTheta[i]->Theta_[1][j]);

                break;
            case pfBinomial:
                MixTheta[i]->Theta_[1][j] = FirstM[i][j] / MixTheta[i]->Theta_[0][j];

                break;
            case pfPoisson:
                MixTheta[i]->Theta_[0][j] = FirstM[i][j];

                break;
            case pfDirac:
                break;
            case pfUniform:
                Error = 1; goto E0; 
            }
        }
    }

E0: return Error;
} // BayesClassificationH

int Rebmix::DegreesOffreedom(int                  c,          // Number of components.
                             CompnentDistribution **MixTheta, // Mixture parameters.
                             int                  *M)         // Degrees of freedom.
{
    int i, j;
    int Error = 0;

    *M = c - 1;

    for (i = 0; i < c; i++) {
        for (j = 0; j < length_pdf_; j++) {
            switch (MixTheta[i]->pdf_[j]) {
            case pfNormal:
                *M += 2;

            break;
            case pfLognormal:
                *M += 2;

                break;
            case pfWeibull:
                *M += 2;

                break;
            case pfGamma:
                *M += 2;

                break;
            case pfvonMises:
                *M += 2;

                break;
            case pfBinomial:
                *M += 1;

            break;
            case pfPoisson:
                *M += 1;

                break;
            case pfDirac:
                *M += 1;

                break;
            case pfUniform:
                Error = 1; goto E0; 
            }
        }
    }

E0: return Error;
} // DegreesOffreedom

// Returns mixture p.d.f..

int Rebmix::MixtureDist(FLOAT                *Y,         // Pointer to the input point [y0,...,yd-1].
                        int                  c,          // Number of components.
                        FLOAT                *W,         // Component weights.
                        CompnentDistribution **MixTheta, // Mixture parameters.
                        FLOAT                *MixDist)   // Mixture distribution.
{
    FLOAT CmpDist;
    int   i;
    int   Error = 0;

    *MixDist = (FLOAT)0.0;

    for (i = 0; i < c; i++) {
        Error = ComponentDist(Y, MixTheta[i], &CmpDist, NULL);

        if (Error) goto E0;

        *MixDist += W[i] * CmpDist;
    }

E0: return Error;
} // MixtureDist

// Returns mixture probability.

int Rebmix::MixtureDist(FLOAT                logV,       // Logarithm of volume of the hypersquare.
                        FLOAT                *Y,         // Pointer to the input point [y0,...,yd-1].
                        int                  c,          // Number of components.
                        FLOAT                *W,         // Component weights.
                        CompnentDistribution **MixTheta, // Mixture parameters.
                        FLOAT                *MixDist)   // Mixture distribution.
{
    FLOAT logCmpDist;
    int   i;
    int   Error = 0;

    *MixDist = (FLOAT)0.0;

    for (i = 0; i < c; i++) {
        Error = LogComponentDist(Y, MixTheta[i], &logCmpDist, NULL);

        if (Error) goto E0;

        *MixDist += W[i] * (FLOAT)exp(logCmpDist + logV);
    }

E0: return Error;
} // MixtureDist 

// Returns information criterion for k-nearest neighbour.

int Rebmix::InformationCriterionKNN(int                  k,          // k-nearest neighbours.
                                    FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1,kl,logV,R].
                                    int                  c,          // Number of components.
                                    FLOAT                *W,         // Component weights.
                                    CompnentDistribution **MixTheta, // Mixture parameters.
                                    FLOAT                *IC,        // Information criterion.
                                    FLOAT                *logL,      // log-likelihood.
                                    int                  *M,         // Degrees of freedom.
                                    FLOAT                *D)         // Total of positive relative deviations.
{
    int   i, j;
    FLOAT E, SSE, EN, PW, K, PC, CmpDist, MixDist, tau;
    int   Error = 0;

    Error = DegreesOffreedom(c, MixTheta, M);

    if (Error) goto E0;

    *IC = *logL = EN = *D = SSE = PW = K = PC = (FLOAT)0.0;

    for (i = 0; i < n_; i++) {
        Error = MixtureDist(Y[i][length_pdf_ + 1], Y[i], c, W, MixTheta, &MixDist);

        if (Error) goto E0;

        E = Y[i][length_pdf_] / n_ - MixDist / k;

        if (E > (FLOAT)0.0) {
            *D += E;
        }

        Error = MixtureDist(Y[i], c, W, MixTheta, &MixDist);

        if (Error) goto E0;

        if (MixDist > FLOAT_MIN) {
            *logL += (FLOAT)log(MixDist);
        }
        else {
            *logL += (FLOAT)log(FLOAT_MIN);
        }

        switch (Criterion_) {
        case icAWE: case icCLC: case icICL: case icPC: case icICLBIC:
            for (j = 0; j < c; j++) {
                Error = ComponentDist(Y[i], MixTheta[j], &CmpDist, NULL);

                if (Error) goto E0;

                if (MixDist > FLOAT_MIN) {
                    tau = W[j] * CmpDist / MixDist;

                    if (tau < FLOAT_MIN) tau = FLOAT_MIN;
                }
                else {
                    tau = (FLOAT)1.0;
                }

                EN -= tau * (FLOAT)log(tau); PC += tau * tau;
            }

            break;
        case icSSE:
            SSE += E * E;

            break;
        case icAIC: case icAIC3: case icAIC4: case icAICc: case icBIC: case icCAIC: case icHQC: case icMDL2: case icMDL5: case icD: default:
            break;
        }
    }

    switch (Criterion_) {
    case icAIC: // AIC - Akaike information criterion Akaike (1973).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M);

        break;
    case icAIC3: // AIC3 - Modified Akaike information criterion Smith & Spiegelhalter (1980).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)3.0 * (*M);

        break;
    case icAIC4: // AIC4 - Modified Akaike information criterion Smith & Spiegelhalter (1980).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)4.0 * (*M);

        break;
    case icAICc: // AICc - Akaike second-order corrected information criterion for small sample sizes Hurvich & Tsai (1989).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * ((FLOAT)1.0 + ((*M) + 1) / (n_ - (*M) - (FLOAT)1.0));

        break;
    case icBIC: // BIC - Bayesian information criterion Schwarz (1978).
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icCAIC: // CAIC - Consistent Akaike information criterion Bozdogan (1987).
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * ((FLOAT)log((FLOAT)n_) + (FLOAT)1.0);

        break;
    case icHQC: // HQC - Hannan-Quinn information criterion Hannan & Quinn (1979).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log((FLOAT)log((FLOAT)n_));

        break;
    case icMDL2: // MDL2 - Minimum description length Liang et al.(1992).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icMDL5: // MDL5 - Minimum description length Liang et al.(1992).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)5.0 * (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icAWE: // AWE - Approximate weight of evidence criterion Banfield & Raftery (1993).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * (*M) * ((FLOAT)3.0 / (FLOAT)2.0 + (FLOAT)log((FLOAT)n_));

        break;
    case icCLC: // CLC - Classification likelihood criterion Biernacki & Govaert (1997).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN;

        break;
    case icICL: // ICL - Integrated classification likelihood Biernacki et al.(1998).
        for (j = 0; j < c; j++) {
            PW += W[j] * (FLOAT)log(W[j]); K += (FLOAT)Gammaln(W[j] * n_ + (FLOAT)0.5);
        }

        K += Gammaln(c * (FLOAT)0.5) - c * Gammaln((FLOAT)0.5) - Gammaln(n_ + c * (FLOAT)0.5);

        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * n_ * PW - (FLOAT)2.0 * K + ((*M) - c + (FLOAT)1.0) * (FLOAT)log((FLOAT)n_);

        break;
    case icPC: // PC - Partition coeficient Bezdek (1981).
        *IC = PC;

        break;
    case icICLBIC: // ICL-BIC - Integrated classification likelihood criterion Biernacki et al.(1998).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icD: // D - Total of positive relative deviations Nagode & Fajdiga (2011).
        *IC = *D;

        break;
    case icSSE: // SSE - Sum of squares error Bishop (1998).
        *IC = (FLOAT)0.5 * SSE;
    }

E0: return Error;
} // InformationCriterionKNN 

// Returns information criterion for Parzen window.

int Rebmix::InformationCriterionPW(FLOAT                logV,       // Logarithm of volume of the hypersquare.
                                   FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1,kl,k].
                                   int                  c,          // Number of components.
                                   FLOAT                *W,         // Component weights.
                                   CompnentDistribution **MixTheta, // Mixture parameters.
                                   FLOAT                *IC,        // Information criterion.
                                   FLOAT                *logL,      // log-likelihood.
                                   int                  *M,         // Degrees of freedom.
                                   FLOAT                *D)         // Total of positive relative deviations.
{
    int   i, j;
    FLOAT E, SSE, EN, PW, K, PC, CmpDist, MixDist, tau;
    int   Error = 0;

    Error = DegreesOffreedom(c, MixTheta, M);

    if (Error) goto E0;

    *IC = *logL = EN = *D = SSE = PW = K = PC = (FLOAT)0.0;

    for (i = 0; i < n_; i++) {
        Error = MixtureDist(logV, Y[i], c, W, MixTheta, &MixDist);

        if (Error) goto E0;

        E = Y[i][length_pdf_] / n_ - MixDist / Y[i][length_pdf_ + 1];

        if (E > (FLOAT)0.0) {
            *D += E;
        }

        Error = MixtureDist(Y[i], c, W, MixTheta, &MixDist);

        if (Error) goto E0;

        if (MixDist > FLOAT_MIN) {
            *logL += (FLOAT)log(MixDist);
        }
        else {
            *logL += (FLOAT)log(FLOAT_MIN);
        }

        switch (Criterion_) {
        case icAWE: case icCLC: case icICL: case icPC: case icICLBIC:
            for (j = 0; j < c; j++) {
                Error = ComponentDist(Y[i], MixTheta[j], &CmpDist, NULL);

                if (Error) goto E0;

                if (MixDist > FLOAT_MIN) {
                    tau = W[j] * CmpDist / MixDist;

                    if (tau < FLOAT_MIN) tau = FLOAT_MIN;
                }
                else {
                    tau = (FLOAT)1.0;
                }

                EN -= tau * (FLOAT)log(tau); PC += tau * tau;
            }

            break;
        case icSSE:
            SSE += E * E;

            break;
        case icAIC: case icAIC3: case icAIC4: case icAICc: case icBIC: case icCAIC: case icHQC: case icMDL2: case icMDL5: case icD: default:
            break;
        }
    }

    switch (Criterion_) {
    case icAIC: // AIC - Akaike information criterion Akaike (1973).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M);

        break;
    case icAIC3: // AIC3 - Modified Akaike information criterion Smith & Spiegelhalter (1980).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)3.0 * (*M);

        break;
    case icAIC4: // AIC4 - Modified Akaike information criterion Smith & Spiegelhalter (1980).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)4.0 * (*M);

        break;
    case icAICc: // AICc - Akaike second-order corrected information criterion for small sample sizes Hurvich & Tsai (1989).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * ((FLOAT)1.0 + ((*M) + 1) / (n_ - (*M) - (FLOAT)1.0));

        break;
    case icBIC: // BIC - Bayesian information criterion Schwarz (1978).
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icCAIC: // CAIC - Consistent Akaike information criterion Bozdogan (1987).
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * ((FLOAT)log((FLOAT)n_) + (FLOAT)1.0);

        break;
    case icHQC: // HQC - Hannan-Quinn information criterion Hannan & Quinn (1979).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log((FLOAT)log((FLOAT)n_));

        break;
    case icMDL2: // MDL2 - Minimum description length Liang et al.(1992).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icMDL5: // MDL5 - Minimum description length Liang et al.(1992).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)5.0 * (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icAWE: // AWE - Approximate weight of evidence criterion Banfield & Raftery (1993).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * (*M) * ((FLOAT)3.0 / (FLOAT)2.0 + (FLOAT)log((FLOAT)n_));

        break;
    case icCLC: // CLC - Classification likelihood criterion Biernacki & Govaert (1997).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN;

        break;
    case icICL: // ICL - Integrated classification likelihood Biernacki et al.(1998).
        for (j = 0; j < c; j++) {
            PW += W[j] * (FLOAT)log(W[j]); K += (FLOAT)Gammaln(W[j] * n_ + (FLOAT)0.5);
        }

        K += Gammaln(c * (FLOAT)0.5) - c * Gammaln((FLOAT)0.5) - Gammaln(n_ + c * (FLOAT)0.5);

        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * n_ * PW - (FLOAT)2.0 * K + ((*M) - c + (FLOAT)1.0) * (FLOAT)log((FLOAT)n_);

        break;
    case icPC: // PC - Partition coeficient Bezdek (1981).
        *IC = PC;

        break;
    case icICLBIC: // ICL-BIC - Integrated classification likelihood criterion Biernacki et al.(1998).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icD: // D - Total of positive relative deviations Nagode & Fajdiga (2011).
        *IC = *D;

        break;
    case icSSE: // SSE - Sum of squares error Bishop (1998).
        *IC = (FLOAT)0.5 * SSE;
    }

E0: return Error;
} // InformationCriterionPW 

// Returns information criterion for histogram.

int Rebmix::InformationCriterionH(FLOAT                logV,       // Logarithm of volume of the hypersquare.
                                  int                  k,          // Total number of bins.
                                  FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1,kl].
                                  int                  c,          // Number of components.
                                  FLOAT                *W,         // Component weights.
                                  CompnentDistribution **MixTheta, // Mixture parameters.
                                  FLOAT                *IC,        // Information criterion.
                                  FLOAT                *logL,      // log-likelihood.
                                  int                  *M,         // Degrees of freedom.
                                  FLOAT                *D)         // Total of positive relative deviations.
{
    int   i, j;
    FLOAT E, SSE, EN, PW, K, PC, CmpDist, MixDist, tau;
    int   Error = 0;

    Error = DegreesOffreedom(c, MixTheta, M);

    if (Error) goto E0;

    *IC = *logL = EN = *D = SSE = PW = K = PC = (FLOAT)0.0;

    for (i = 0; i < k; i++) {
        Error = MixtureDist(logV, Y[i], c, W, MixTheta, &MixDist);

        if (Error) goto E0;

        E = Y[i][length_pdf_] / n_ - MixDist;

        if (E > (FLOAT)0.0) {
            *D += E;
        }

        SSE += E * E;
    }

    for (i = 0; i < n_; i++) {
        Error = MixtureDist(Y_[i], c, W, MixTheta, &MixDist);

        if (Error) goto E0;

        if (MixDist > FLOAT_MIN) {
            *logL += (FLOAT)log(MixDist);
        }
        else {
            *logL += (FLOAT)log(FLOAT_MIN);
        }

        switch (Criterion_) {
        case icAWE: case icCLC: case icICL: case icPC: case icICLBIC:
            for (j = 0; j < c; j++) {
                Error = ComponentDist(Y_[i], MixTheta[j], &CmpDist, NULL);

                if (Error) goto E0;

                if (MixDist > FLOAT_MIN) {
                    tau = W[j] * CmpDist / MixDist;

                    if (tau < FLOAT_MIN) tau = FLOAT_MIN;
                }
                else {
                    tau = (FLOAT)1.0;
                }

                EN -= tau * (FLOAT)log(tau); PC += tau * tau;
            }

            break;
        case icSSE:
            break;
        case icAIC: case icAIC3: case icAIC4: case icAICc: case icBIC: case icCAIC: case icHQC: case icMDL2: case icMDL5: case icD: default:
            break;
        }
    }

    switch (Criterion_) {
    case icAIC: // AIC - Akaike information criterion Akaike (1973).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M);

        break;
    case icAIC3: // AIC3 - Modified Akaike information criterion Smith & Spiegelhalter (1980).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)3.0 * (*M);

        break;
    case icAIC4: // AIC4 - Modified Akaike information criterion Smith & Spiegelhalter (1980).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)4.0 * (*M);

        break;
    case icAICc: // AICc - Akaike second-order corrected information criterion for small sample sizes Hurvich & Tsai (1989).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * ((FLOAT)1.0 + ((*M) + 1) / (n_ - (*M) - (FLOAT)1.0));

        break;
    case icBIC: // BIC - Bayesian information criterion Schwarz (1978).
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icCAIC: // CAIC - Consistent Akaike information criterion Bozdogan (1987).
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * ((FLOAT)log((FLOAT)n_) + (FLOAT)1.0);

        break;
    case icHQC: // HQC - Hannan-Quinn information criterion Hannan & Quinn (1979).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log((FLOAT)log((FLOAT)n_));

        break;
    case icMDL2: // MDL2 - Minimum description length Liang et al.(1992).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icMDL5: // MDL5 - Minimum description length Liang et al.(1992).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)5.0 * (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icAWE: // AWE - Approximate weight of evidence criterion Banfield & Raftery (1993).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * (*M) * ((FLOAT)3.0 / (FLOAT)2.0 + (FLOAT)log((FLOAT)n_));

        break;
    case icCLC: // CLC - Classification likelihood criterion Biernacki & Govaert (1997).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN;

        break;
    case icICL: // ICL - Integrated classification likelihood Biernacki et al.(1998).
        for (j = 0; j < c; j++) {
            PW += W[j] * (FLOAT)log(W[j]); K += (FLOAT)Gammaln(W[j] * n_ + (FLOAT)0.5);
        }

        K += Gammaln(c * (FLOAT)0.5) - c * Gammaln((FLOAT)0.5) - Gammaln(n_ + c * (FLOAT)0.5);

        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * n_ * PW - (FLOAT)2.0 * K + ((*M) - c + (FLOAT)1.0) * (FLOAT)log((FLOAT)n_);

        break;
    case icPC: // PC - Partition coeficient Bezdek (1981).
        *IC = PC;

        break;
    case icICLBIC: // ICL-BIC - Integrated classification likelihood criterion Biernacki et al.(1998).
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (*M) * (FLOAT)log((FLOAT)n_);

        break;
    case icD: // D - Total of positive relative deviations Nagode & Fajdiga (2011).
        *IC = *D;

        break;
    case icSSE: // SSE - Sum of squares error Bishop (1998).
        *IC = (FLOAT)0.5 * SSE;
    }

E0: return Error;
} // InformationCriterionH

int Rebmix::CombineComponentsKNN(FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1,kl,k].
                                 int                  c,          // Number of components.
                                 FLOAT                *W,         // Component weights.
                                 CompnentDistribution **MixTheta, // Mixture parameters.
                                 int                  *F,         // From components.
                                 int                  *T,         // To components.
                                 FLOAT                *EN,        // Entropy.
                                 FLOAT                *ED)        // Entropy decrease.
{
    int   *C = NULL, i, ii, j, jj, II, J, Tmpi;
    FLOAT CmpDist, ed, en, MixDist, **Tau = NULL, *Tmpf = NULL;
    int   Error = 0;

    Tau = (FLOAT**)malloc(c * sizeof(FLOAT*));

    Error = NULL == Tau; if (Error) goto E0;

    for (i = 0; i < c; i++) {
        Tau[i] = (FLOAT*)malloc(n_ * sizeof(FLOAT));

        Error = NULL == Tau[i]; if (Error) goto E0;
    }

    en = (FLOAT)0.0;

    for (i = 0; i < n_; i++) {
        Error = MixtureDist(Y[i], c, W, MixTheta, &MixDist);

        if (Error) goto E0;

        for (j = 0; j < c; j++) {
            Error = ComponentDist(Y[i], MixTheta[j], &CmpDist, NULL);

            if (Error) goto E0;

            if (MixDist > FLOAT_MIN) {
                Tau[j][i] = W[j] * CmpDist / MixDist;

                if (Tau[j][i] < FLOAT_MIN) Tau[j][i] = FLOAT_MIN;

                en -= Tau[j][i] * (FLOAT)log(Tau[j][i]);
            }
            else {
                Tau[j][i] = (FLOAT)1.0;
            }
        }
    }

    C = (int*)malloc(c * sizeof(int));

    Error = NULL == C; if (Error) goto E0;

    for (i = 0; i < c; i++) {
        C[i] = i; F[i] = T[i] = 0; ED[i] = (FLOAT)0.0; EN[i] = en;       
    }

    i = c;

    while (i > 1) {
        II = J = 0; ED[i - 2] = (FLOAT)0.0;
 
        for (ii = 0; ii < i - 1; ii++) {
            for (j = ii + 1; j < i; j++) {
                ed = (FLOAT)0.0;
 
                for (jj = 0; jj < n_; jj++) {
                    ed -= Tau[ii][jj] * (FLOAT)log(Tau[ii][jj]) + Tau[j][jj] * (FLOAT)log(Tau[j][jj]);
                
                    ed += (Tau[ii][jj] + Tau[j][jj]) * (FLOAT)log(Tau[ii][jj] + Tau[j][jj]);
                }

                if (ed > ED[i - 2]) {
                    ED[i - 2] = ed; II = ii; J = j;
                }
            }
        }

        F[i - 2] = C[J] + 1, T[i - 2] = C[II] + 1;

        for (j = 0; j < n_; j++) {
            Tau[II][j] += Tau[J][j];
        }

        Tmpf = Tau[J]; Tmpi = C[J];

        for (j = J; j < c - 1; j++) {
            Tau[j] = Tau[j + 1]; C[j] = C[j + 1];
        }

        Tau[c - 1] = Tmpf; C[c - 1] = Tmpi;

        EN[i - 2] = (FLOAT)0.0;

        for (ii = 0; ii < i - 1; ii++) {
            for (jj = 0; jj < n_; jj++) {
                EN[i - 2] -= Tau[ii][jj] * (FLOAT)log(Tau[ii][jj]);
            }
        }

        i--;
    }

E0: if (C) free(C);

    if (Tau) {
        for (i = 0; i < c; i++) {
            if (Tau[i]) free(Tau[i]);
        }

        free(Tau);
    }

    return Error;
} // CombineComponentsKNN

// Returns combined components for Parzen window.

int Rebmix::CombineComponentsPW(FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1,kl,k].
                                int                  c,          // Number of components.
                                FLOAT                *W,         // Component weights.
                                CompnentDistribution **MixTheta, // Mixture parameters.
                                int                  *F,         // From components.
                                int                  *T,         // To components.
                                FLOAT                *EN,        // Entropy.
                                FLOAT                *ED)        // Entropy decrease.
{
    int   *C = NULL, i, ii, j, jj, II, J, Tmpi;
    FLOAT CmpDist, ed, en, MixDist, **Tau = NULL, *Tmpf = NULL;
    int   Error = 0;

    Tau = (FLOAT**)malloc(c * sizeof(FLOAT*));

    Error = NULL == Tau; if (Error) goto E0;

    for (i = 0; i < c; i++) {
        Tau[i] = (FLOAT*)malloc(n_ * sizeof(FLOAT));

        Error = NULL == Tau[i]; if (Error) goto E0;
    }

    en = (FLOAT)0.0;

    for (i = 0; i < n_; i++) {
        Error = MixtureDist(Y[i], c, W, MixTheta, &MixDist);

        if (Error) goto E0;

        for (j = 0; j < c; j++) {
            Error = ComponentDist(Y[i], MixTheta[j], &CmpDist, NULL);

            if (Error) goto E0;

            if (MixDist > FLOAT_MIN) {
                Tau[j][i] = W[j] * CmpDist / MixDist;

                if (Tau[j][i] < FLOAT_MIN) Tau[j][i] = FLOAT_MIN;

                en -= Tau[j][i] * (FLOAT)log(Tau[j][i]);
            }
            else {
                Tau[j][i] = (FLOAT)1.0;
            }
        }
    }

    C = (int*)malloc(c * sizeof(int));

    Error = NULL == C; if (Error) goto E0;

    for (i = 0; i < c; i++) {
        C[i] = i; F[i] = T[i] = 0; ED[i] = (FLOAT)0.0; EN[i] = en;       
    }

    i = c;

    while (i > 1) {
        II = J = 0; ED[i - 2] = (FLOAT)0.0;
 
        for (ii = 0; ii < i - 1; ii++) {
            for (j = ii + 1; j < i; j++) {
                ed = (FLOAT)0.0;
 
                for (jj = 0; jj < n_; jj++) {
                    ed -= Tau[ii][jj] * (FLOAT)log(Tau[ii][jj]) + Tau[j][jj] * (FLOAT)log(Tau[j][jj]);
                
                    ed += (Tau[ii][jj] + Tau[j][jj]) * (FLOAT)log(Tau[ii][jj] + Tau[j][jj]);
                }

                if (ed > ED[i - 2]) {
                    ED[i - 2] = ed; II = ii; J = j;
                }
            }
        }

        F[i - 2] = C[J] + 1, T[i - 2] = C[II] + 1;

        for (j = 0; j < n_; j++) {
            Tau[II][j] += Tau[J][j];
        }

        Tmpf = Tau[J]; Tmpi = C[J];

        for (j = J; j < c - 1; j++) {
            Tau[j] = Tau[j + 1]; C[j] = C[j + 1];
        }

        Tau[c - 1] = Tmpf; C[c - 1] = Tmpi;

        EN[i - 2] = (FLOAT)0.0;

        for (ii = 0; ii < i - 1; ii++) {
            for (jj = 0; jj < n_; jj++) {
                EN[i - 2] -= Tau[ii][jj] * (FLOAT)log(Tau[ii][jj]);
            }
        }

        i--;
    }

E0: if (C) free(C);

    if (Tau) {
        for (i = 0; i < c; i++) {
            if (Tau[i]) free(Tau[i]);
        }

        free(Tau);
    }

    return Error;
} // CombineComponentsPW

// Returns combined components for histogram.

int Rebmix::CombineComponentsH(int                  k,          // Total number of bins.
                               FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1,kl].
                               int                  c,          // Number of components.
                               FLOAT                *W,         // Component weights.
                               CompnentDistribution **MixTheta, // Mixture parameters.
                               int                  *F,         // From components.
                               int                  *T,         // To components.
                               FLOAT                *EN,        // Entropy.
                               FLOAT                *ED)        // Entropy decrease.
{
    int   *C = NULL, i, ii, j, jj, II, J, Tmpi;
    FLOAT CmpDist, ed, en, MixDist, **Tau = NULL, *Tmpf = NULL;
    int   Error = 0;

    Tau = (FLOAT**)malloc(c * sizeof(FLOAT*));

    Error = NULL == Tau; if (Error) goto E0;

    for (i = 0; i < c; i++) {
        Tau[i] = (FLOAT*)malloc(k * sizeof(FLOAT));

        Error = NULL == Tau[i]; if (Error) goto E0;
    }

    en = (FLOAT)0.0;

    for (i = 0; i < k; i++) {
        Error = MixtureDist(Y[i], c, W, MixTheta, &MixDist);

        if (Error) goto E0;

        for (j = 0; j < c; j++) {
            Error = ComponentDist(Y[i], MixTheta[j], &CmpDist, NULL);

            if (Error) goto E0;

            if (MixDist > FLOAT_MIN) {
                Tau[j][i] = W[j] * CmpDist / MixDist;

                if (Tau[j][i] < FLOAT_MIN) Tau[j][i] = FLOAT_MIN;

                en -= Y[i][length_pdf_] * Tau[j][i] * (FLOAT)log(Tau[j][i]);
            }
            else {
                Tau[j][i] = (FLOAT)1.0;
            }
        }
    }

    C = (int*)malloc(c * sizeof(int));

    Error = NULL == C; if (Error) goto E0;

    for (i = 0; i < c; i++) {
        C[i] = i; F[i] = T[i] = 0; ED[i] = (FLOAT)0.0; EN[i] = en;       
    }

    i = c;

    while (i > 1) {
        II = J = 0; ED[i - 2] = (FLOAT)0.0;
 
        for (ii = 0; ii < i - 1; ii++) {
            for (j = ii + 1; j < i; j++) {
                ed = (FLOAT)0.0;
 
                for (jj = 0; jj < k; jj++) {
                    ed -= Y[jj][length_pdf_] * (Tau[ii][jj] * (FLOAT)log(Tau[ii][jj]) + Tau[j][jj] * (FLOAT)log(Tau[j][jj]));
                
                    ed += Y[jj][length_pdf_] * (Tau[ii][jj] + Tau[j][jj]) * (FLOAT)log(Tau[ii][jj] + Tau[j][jj]);
                }

                if (ed > ED[i - 2]) {
                    ED[i - 2] = ed; II = ii; J = j;
                }
            }
        }

        F[i - 2] = C[J] + 1, T[i - 2] = C[II] + 1;

        for (j = 0; j < k; j++) {
            Tau[II][j] += Tau[J][j];
        }

        Tmpf = Tau[J]; Tmpi = C[J];

        for (j = J; j < c - 1; j++) {
            Tau[j] = Tau[j + 1]; C[j] = C[j + 1];
        }

        Tau[c - 1] = Tmpf; C[c - 1] = Tmpi;

        EN[i - 2] = (FLOAT)0.0;

        for (ii = 0; ii < i - 1; ii++) {
            for (jj = 0; jj < k; jj++) {
                EN[i - 2] -= Y[jj][length_pdf_] * Tau[ii][jj] * (FLOAT)log(Tau[ii][jj]);
            }
        }

        i--;
    }

E0: if (C) free(C);

    if (Tau) {
        for (i = 0; i < c; i++) {
            if (Tau[i]) free(Tau[i]);
        }

        free(Tau);
    }

    return Error;
} // CombineComponentsH 

// REBMIX algorithm for k-nearest neighbours.

int Rebmix::REBMIXKNN()
{
    FLOAT                **Y = NULL;
    FLOAT                *ymin = NULL, *ymax = NULL, *h = NULL;
    FLOAT                *R = NULL, *E = NULL, *Epsilon = NULL;
    FLOAT                *W = NULL;
    CompnentDistribution **RigidTheta = NULL, **LooseTheta = NULL; 
    FLOAT                **FirstM = NULL, **SecondM = NULL;
    int                  opt_length;    
    int                  *O = NULL, *opt_c = NULL;        
    FLOAT                *opt_IC = NULL;       
    FLOAT                *opt_logL = NULL;     
    FLOAT                *opt_D = NULL;
    int                  c = 0, i, I, j, J, l, m, M;
    FLOAT                Dmin, r, lognl, nl, elp, eln, epsilonlmax, logfl, fl, Dl, f, IC, logL, D;
    FLOAT                A = (FLOAT)0.0, ar;
    int                  Error = 0, Stop = 0, Found = 0, Outlier = 0;

    // Allocation and initialisation.

    summary_.IC = FLOAT_MAX;

    summary_.h = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == summary_.h; if (Error) goto E0;

    summary_.ymin = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == summary_.ymin; if (Error) goto E0;

    summary_.ymax = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == summary_.ymax; if (Error) goto E0;

    W_ = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));

    Error = NULL == W_; if (Error) goto E0;

    MixTheta_ = new CompnentDistribution* [(unsigned int)cmax_];

    Error = NULL == MixTheta_; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        MixTheta_[i] = new CompnentDistribution(this);

        Error = NULL == MixTheta_[i]; if (Error) goto E0;

        Error = MixTheta_[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;
    }

    opt_length_ = ItMax;

    opt_c_ = (int*)malloc(ItMax * sizeof(int));

    Error = NULL == opt_c_; if (Error) goto E0;

    opt_IC_ = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_IC_; if (Error) goto E0;

    opt_logL_ = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_logL_; if (Error) goto E0;

    opt_D_ = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_D_; if (Error) goto E0;

    all_length_ = K_[length_K_ - 1] - K_[0] + 1;

    all_I_ = (int*)calloc((size_t)all_length_, sizeof(int));

    Error = NULL == all_I_; if (Error) goto E0;

    all_K_ = (int*)calloc((size_t)all_length_, sizeof(int));

    Error = NULL == all_K_; if (Error) goto E0;

    for (i = 0; i < length_K_; i++) {
        all_K_[K_[i] - K_[0]] = K_[i];
    }

    additional_.Bracket = 1;

    all_IC_ = (FLOAT*)malloc(all_length_ * sizeof(FLOAT));

    Error = NULL == all_IC_; if (Error) goto E0;

    for (i = 0; i < all_length_; i++) {
        all_IC_[i] = FLOAT_MAX;
    }

    Y = (FLOAT**)malloc(n_ * sizeof(FLOAT*));

    Error = NULL == Y; if (Error) goto E0;

    for (i = 0; i < n_; i++) {
        Y[i] = (FLOAT*)malloc((length_pdf_ + 3) * sizeof(FLOAT));

        Error = NULL == Y[i]; if (Error) goto E0;

        for (j = 0; j < length_pdf_; j++) Y[i][j] = Y_[i][j];
    }

    ymin = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == ymin; if (Error) goto E0;

    if (ymin_ == NULL) {
        for (i = 0; i < length_pdf_; i++) {
            ymin[i] = Y_[0][i];

            for (j = 1; j < n_; j++) {
                if (Y_[j][i] < ymin[i]) ymin[i] = Y_[j][i];
            }
        }
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            ymin[i] = ymin_[i];
        }
    }

    ymax = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == ymax; if (Error) goto E0;

    if (ymax_ == NULL) {
        for (i = 0; i < length_pdf_; i++) {
            ymax[i] = Y_[0][i];

            for (j = 1; j < n_; j++) {
                if (Y_[j][i] > ymax[i]) ymax[i] = Y_[j][i];
            }
        }
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            ymax[i] = ymax_[i];
        }
    }

    h = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == h; if (Error) goto E0;

    for (i = 0; i < length_pdf_; i++) {
        h[i] = ymax[i] - ymin[i];
    }

    R = (FLOAT*)malloc(n_ * sizeof(FLOAT));

    Error = NULL == R; if (Error) goto E0;

    E = (FLOAT*)malloc(n_ * sizeof(FLOAT));

    Error = NULL == E; if (Error) goto E0;

    Epsilon = (FLOAT*)malloc(n_ * sizeof(FLOAT));

    Error = NULL == Epsilon; if (Error) goto E0;

    W = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));

    Error = NULL == W; if (Error) goto E0;

    RigidTheta = new CompnentDistribution* [(unsigned int)cmax_];

    Error = NULL == RigidTheta; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        RigidTheta[i] = new CompnentDistribution(this);

        Error = NULL == RigidTheta[i]; if (Error) goto E0;

        Error = RigidTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;

        Error = RigidTheta[i]->Memmove(IniTheta_);

        if (Error) goto E0;
    }

    LooseTheta = new CompnentDistribution* [(unsigned int)cmax_];

    Error = NULL == LooseTheta; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        LooseTheta[i] = new CompnentDistribution(this);

        Error = NULL == LooseTheta[i]; if (Error) goto E0;

        Error = LooseTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;

        Error = LooseTheta[i]->Memmove(IniTheta_);

        if (Error) goto E0;
    }

    FirstM = (FLOAT**)malloc(cmax_ * sizeof(FLOAT*));

    Error = NULL == FirstM; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        FirstM[i] = (FLOAT*)malloc(length_theta_[0] * sizeof(FLOAT));

        Error = NULL == FirstM[i]; if (Error) goto E0;
    }

    SecondM = (FLOAT**)malloc(cmax_ * sizeof(FLOAT*));

    Error = NULL == SecondM; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        SecondM[i] = (FLOAT*)malloc(length_theta_[1] * sizeof(FLOAT));

        Error = NULL == SecondM[i]; if (Error) goto E0;
    }


    opt_length = ItMax;

    opt_c = (int*)malloc(opt_length * sizeof(int));

    Error = NULL == opt_c; if (Error) goto E0;

    opt_IC = (FLOAT*)malloc(opt_length * sizeof(FLOAT));

    Error = NULL == opt_IC; if (Error) goto E0;

    opt_logL = (FLOAT*)malloc(opt_length * sizeof(FLOAT));

    Error = NULL == opt_logL; if (Error) goto E0;

    opt_D = (FLOAT*)malloc(opt_length * sizeof(FLOAT));

    Error = NULL == opt_D; if (Error) goto E0;

    O = (int*)malloc(n_ * sizeof(int));

    Error = NULL == O; if (Error) goto E0;

    Error = Initialize();

    if (Error) goto E0;

    do for (i = 0; i < all_length_; i++) if (all_K_[i] && (all_I_[i] == 0)) {
        // Preprocessing of observations.

        Error = PreprocessingKNN(all_K_[i], h, Y);

        if (Error) goto E0;

        all_I_[i] = 1;

        Found = 0; Dmin = (FLOAT)1.0; J = 1;

        // Outer loop.

        while (J <= ItMax) {
            l = 0; r = (FLOAT)n_; nl = (FLOAT)n_;

            memset(O, 0, n_ * sizeof(int));

            // Middle loop.

            while (nl / n_ > Dmin * l) {
                // Global mode detection.

                Error = GlobalModeKNN(&m, Y, h, O);

                if (Error) goto E0;

                I = 1; W[l] = nl / n_; memset(R, 0, n_ * sizeof(FLOAT));

                // Inner loop.

                while (I <= ItMax) {
                    // Rough component parameter estimation.

                    Error = RoughEstimationKNN(Y, all_K_[i], h, nl, m, RigidTheta[l], LooseTheta[l]);

                    if (Error) goto E0;

                    elp = eln = epsilonlmax = (FLOAT)0.0; lognl = (FLOAT)log(nl);

                    for (j = 0; j < n_; j++) {
                        E[j] = Epsilon[j] = (FLOAT)0.0;

                        if ((Y[j][length_pdf_] > FLOAT_MIN) || (R[j] > FLOAT_MIN)) {
                            Error = LogComponentDist(Y[j], LooseTheta[l], &logfl, NULL);

                            if (Error) goto E0;

                            E[j] = Y[j][length_pdf_] - (FLOAT)exp(lognl + logfl + Y[j][length_pdf_ + 1]) / all_K_[i];

                            if (E[j] > (FLOAT)0.0) {
                                Epsilon[j] = E[j] / Y[j][length_pdf_]; 
                                
                                if (Epsilon[j] > epsilonlmax) epsilonlmax = Epsilon[j]; 
                                
                                elp += E[j];
                            }
                            else {
                                if (E[j] < -R[j]) E[j] = -R[j]; 
                                
                                eln -= E[j];
                            }
                        }
                    }
                    
                    Dl = elp / nl;

                    if ((Dl <= Dmin / W[l]) || (I == ItMax) || (nl <= (FLOAT)1.0)) {
                        // Enhanced component parameter estimation.

                        EnhancedEstimationKNN(Y, nl, RigidTheta[l], LooseTheta[l]);

                        break;
                    }
                    else {
                        if (I == 1) {
                            A = ((FLOAT)1.0 - ar_) / ar_ / (Dl * W[l] - Dmin);
                        }

                        ar = (FLOAT)1.0 / (A * (Dl * W[l] - Dmin) + (FLOAT)1.0);

                        epsilonlmax *= (FLOAT)1.0 - ar;

                        for (j = 0; j < n_; j++) if (Epsilon[j] > epsilonlmax) {
                            Y[j][length_pdf_] -= E[j]; R[j] += E[j]; nl -= E[j];
                        }

                        if (eln > FLOAT_MIN) {
                            elp = elp / Dl - nl; if (eln > elp) f = elp / eln; else f = (FLOAT)1.0;

                            for (j = 0; j < n_; j++) if (E[j] < (FLOAT)0.0) {
                                E[j] *= f; Y[j][length_pdf_] -= E[j]; R[j] += E[j]; nl -= E[j];
                            }
                        }

                        W[l] = nl / n_;
                    }

                    I++;
                }

                // Outlier detection.

                for (j = 0; j < n_; j++) {
                    Error = ComponentDist(Y[j], LooseTheta[l], &fl, &Outlier);

                    if (Error) goto E0;
 
                    if (!Outlier) O[j] = 1;
                }

                // Moments calculation.

                Error = MomentsCalculation(LooseTheta[l], FirstM[l], SecondM[l]);

                if (Error) goto E0;

                c = ++l;

                r -= nl; nl = r; for (j = 0; j < n_; j++) Y[j][length_pdf_] = R[j];

                Stop = (c >= n_) || (c >= cmax_);

                if (Stop) break;
            }

            // Bayes classification of the remaining observations.
            
            Error = BayesClassificationKNN(Y, c, W, LooseTheta, FirstM, SecondM);

            if (Error) goto E0;

            for (j = 0; j < n_; j++) Y[j][length_pdf_] = (FLOAT)1.0;

            Error = InformationCriterionKNN(all_K_[i], Y, c, W, LooseTheta, &IC, &logL, &M, &D);
            
            if (Error) goto E0;

            if (IC < all_IC_[i]) all_IC_[i] = IC;

            if (IC < summary_.IC) {
                Found = 1;

                summary_.k = all_K_[i];

                memmove(summary_.h, h, length_pdf_ * sizeof(FLOAT));  

                memmove(summary_.ymin, ymin, length_pdf_ * sizeof(FLOAT)); 

                memmove(summary_.ymax, ymax, length_pdf_ * sizeof(FLOAT));  
                
                summary_.IC = IC; summary_.logL = logL; summary_.M = M; summary_.c = c; 

                memmove(W_, W, c * sizeof(FLOAT));  

                for (j = 0; j < c; j++) {
                    Error = MixTheta_[j]->Memmove(LooseTheta[j]);

                    if (Error) goto E0;
                }
            }

            j = J - 1; opt_c[j] = c; opt_IC[j] = IC; opt_logL[j] = logL; opt_D[j] = D;

            Dmin *= c / (c + (FLOAT)1.0); J++;

            if (Stop) break; 
        }

        opt_length = J - 1;

        if (Found) {
            opt_length_ = opt_length;

            memmove(opt_c_, opt_c, opt_length_ * sizeof(int));  
            memmove(opt_IC_, opt_IC, opt_length_ * sizeof(FLOAT));  
            memmove(opt_logL_, opt_logL, opt_length_ * sizeof(FLOAT));  
            memmove(opt_D_, opt_D, opt_length_ * sizeof(FLOAT));  
        }
    }
    while (!Golden());

E0: if (O) free(O);

    if (opt_D) free(opt_D);

    if (opt_logL) free(opt_logL);

    if (opt_IC) free(opt_IC);

    if (opt_c) free(opt_c);
    
    if (SecondM) {
        for (i = 0; i < cmax_; i++) {
            if (SecondM[i]) free(SecondM[i]);
        }
         
        free(SecondM);
    }

    if (FirstM) {
        for (i = 0; i < cmax_; i++) {
            if (FirstM[i]) free(FirstM[i]);
        }
         
        free(FirstM);
    }

    if (LooseTheta) {
        for (i = 0; i < cmax_; i++) {
            if (LooseTheta[i]) delete LooseTheta[i];
        }

        delete[] LooseTheta;
    }

    if (RigidTheta) {
        for (i = 0; i < cmax_; i++) {
            if (RigidTheta[i]) delete RigidTheta[i];
        }

        delete[] RigidTheta;
    }

    if (W) free(W);

    if (Epsilon) free(Epsilon);

    if (E) free(E);

    if (R) free(R);

    if (h) free(h);

    if (ymax) free(ymax);

    if (ymin) free(ymin);

    if (Y) {
        for (i = 0; i < n_; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    return Error;
} // REBMIXKNN

// REBMIX algorithm for Parzen window.

int Rebmix::REBMIXPW()
{
    FLOAT                **Y = NULL;
    FLOAT                *ymin = NULL, *ymax = NULL, *h = NULL;
    FLOAT                *R = NULL, *E = NULL, *Epsilon = NULL;
    FLOAT                *W = NULL;
    CompnentDistribution **RigidTheta = NULL, **LooseTheta = NULL; 
    FLOAT                **FirstM = NULL, **SecondM = NULL;
    int                  opt_length;
    int                  *O = NULL, *opt_c = NULL;
    FLOAT                *opt_IC = NULL;
    FLOAT                *opt_logL = NULL;
    FLOAT                *opt_D = NULL;
    int                  c = 0, i, I, j, J, l, m, M;
    FLOAT                logV, Dmin, r, lognl, nl, elp, eln, epsilonlmax, logfl, fl, Dl, f, IC, logL, D;
    FLOAT                A = (FLOAT)0.0, ar;
    int                  Error = 0, Stop = 0, Found = 0, Outlier = 0;

    // Allocation and initialisation.

    summary_.IC = FLOAT_MAX;

    summary_.h = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == summary_.h; if (Error) goto E0;

    summary_.ymin = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == summary_.ymin; if (Error) goto E0;

    summary_.ymax = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == summary_.ymax; if (Error) goto E0;

    W_ = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));

    Error = NULL == W_; if (Error) goto E0;

    MixTheta_ = new CompnentDistribution* [(unsigned int)cmax_];

    Error = NULL == MixTheta_; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        MixTheta_[i] = new CompnentDistribution(this);

        Error = NULL == MixTheta_[i]; if (Error) goto E0;

        Error = MixTheta_[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;
    }

    opt_length_ = ItMax;

    opt_c_ = (int*)malloc(ItMax * sizeof(int));

    Error = NULL == opt_c_; if (Error) goto E0;

    opt_IC_ = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_IC_; if (Error) goto E0;

    opt_logL_ = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_logL_; if (Error) goto E0;

    opt_D_ = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_D_; if (Error) goto E0;

    all_length_ = K_[length_K_ - 1] - K_[0] + 1;

    all_I_ = (int*)calloc((size_t)all_length_, sizeof(int));

    Error = NULL == all_I_; if (Error) goto E0;

    all_K_ = (int*)calloc((size_t)all_length_, sizeof(int));

    Error = NULL == all_K_; if (Error) goto E0;

    for (i = 0; i < length_K_; i++) {
        all_K_[K_[i] - K_[0]] = K_[i];
    }

    additional_.Bracket = 1;

    all_IC_ = (FLOAT*)malloc(all_length_ * sizeof(FLOAT));

    Error = NULL == all_IC_; if (Error) goto E0;

    for (i = 0; i < all_length_; i++) {
        all_IC_[i] = FLOAT_MAX;
    }

    Y = (FLOAT**)malloc(n_ * sizeof(FLOAT*));

    Error = NULL == Y; if (Error) goto E0;

    for (i = 0; i < n_; i++) {
        Y[i] = (FLOAT*)malloc((length_pdf_ + 2) * sizeof(FLOAT));

        Error = NULL == Y[i]; if (Error) goto E0;

        for (j = 0; j < length_pdf_; j++) Y[i][j] = Y_[i][j];
    }

    ymin = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == ymin; if (Error) goto E0;

    if (ymin_ == NULL) {
        for (i = 0; i < length_pdf_; i++) {
            ymin[i] = Y_[0][i];

            for (j = 1; j < n_; j++) {
                if (Y_[j][i] < ymin[i]) ymin[i] = Y_[j][i];
            }
        }
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            ymin[i] = ymin_[i];
        }
    }

    ymax = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == ymax; if (Error) goto E0;

    if (ymax_ == NULL) {
        for (i = 0; i < length_pdf_; i++) {
            ymax[i] = Y_[0][i];

            for (j = 1; j < n_; j++) {
                if (Y_[j][i] > ymax[i]) ymax[i] = Y_[j][i];
            }
        }
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            ymax[i] = ymax_[i];
        }
    }

    h = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == h; if (Error) goto E0;

    R = (FLOAT*)malloc(n_ * sizeof(FLOAT));

    Error = NULL == R; if (Error) goto E0;

    E = (FLOAT*)malloc(n_ * sizeof(FLOAT));

    Error = NULL == E; if (Error) goto E0;

    Epsilon = (FLOAT*)malloc(n_ * sizeof(FLOAT));

    Error = NULL == Epsilon; if (Error) goto E0;

    W = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));

    Error = NULL == W; if (Error) goto E0;

    RigidTheta = new CompnentDistribution* [(unsigned int)cmax_];

    Error = NULL == RigidTheta; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        RigidTheta[i] = new CompnentDistribution(this);

        Error = NULL == RigidTheta[i]; if (Error) goto E0;

        Error = RigidTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;

        Error = RigidTheta[i]->Memmove(IniTheta_);

        if (Error) goto E0;
    }

    LooseTheta = new CompnentDistribution* [(unsigned int)cmax_];

    Error = NULL == LooseTheta; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        LooseTheta[i] = new CompnentDistribution(this);

        Error = NULL == LooseTheta[i]; if (Error) goto E0;

        Error = LooseTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;

        Error = LooseTheta[i]->Memmove(IniTheta_);

        if (Error) goto E0;
    }

    FirstM = (FLOAT**)malloc(cmax_ * sizeof(FLOAT*));

    Error = NULL == FirstM; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        FirstM[i] = (FLOAT*)malloc(length_theta_[0] * sizeof(FLOAT));

        Error = NULL == FirstM[i]; if (Error) goto E0;
    }

    SecondM = (FLOAT**)malloc(cmax_ * sizeof(FLOAT*));

    Error = NULL == SecondM; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        SecondM[i] = (FLOAT*)malloc(length_theta_[1] * sizeof(FLOAT));

        Error = NULL == SecondM[i]; if (Error) goto E0;
    }

    opt_length = ItMax;

    opt_c = (int*)malloc(ItMax * sizeof(int));

    Error = NULL == opt_c; if (Error) goto E0;

    opt_IC = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_IC; if (Error) goto E0;

    opt_logL = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_logL; if (Error) goto E0;

    opt_D = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_D; if (Error) goto E0;

    O = (int*)malloc(n_ * sizeof(int));

    Error = NULL == O; if (Error) goto E0;

    Error = Initialize();

    if (Error) goto E0;

    do for (i = 0; i < all_length_; i++) if (all_K_[i] && (all_I_[i] == 0)) {
        // Preprocessing of observations.

        logV = (FLOAT)0.0;
        
        for (j = 0; j < length_pdf_; j++) {
            switch (Variables_[j]) {
            case vtContinuous:
                h[j] = (ymax[j] - ymin[j]) / all_K_[i]; logV += (FLOAT)log(h[j]); 

                break;
            case vtDiscrete:
                h[j] = (FLOAT)1.0;
            }
        }

        Error = PreprocessingPW(h, Y);

        if (Error) goto E0;

        all_I_[i] = 1;

        Found = 0; Dmin = (FLOAT)1.0; J = 1;

        // Outer loop.

        while (J <= ItMax) {
            l = 0; r = (FLOAT)n_; nl = (FLOAT)n_;

            memset(O, 0, n_ * sizeof(int));

            // Middle loop.

            while (nl / n_ > Dmin * l) {
                // Global mode detection.

                Error = GlobalModePW(&m, Y, h, O);

                if (Error) goto E0;

                I = 1; W[l] = nl / n_; memset(R, 0, n_ * sizeof(FLOAT));

                // Inner loop.

                while (I <= ItMax) {
                    // Rough component parameter estimation.

                    Error = RoughEstimationPW(Y, h, nl, m, RigidTheta[l], LooseTheta[l]);

                    if (Error) goto E0;

                    elp = eln = epsilonlmax = (FLOAT)0.0; lognl = (FLOAT)log(nl);

                    for (j = 0; j < n_; j++) {
                        E[j] = Epsilon[j] = (FLOAT)0.0;

                        if ((Y[j][length_pdf_] > FLOAT_MIN) || (R[j] > FLOAT_MIN)) {
                            Error = LogComponentDist(Y[j], LooseTheta[l], &logfl, NULL);

                            if (Error) goto E0;

                            E[j] = Y[j][length_pdf_] - (FLOAT)exp(lognl + logfl  + logV) / Y[j][length_pdf_ + 1];

                            if (E[j] > (FLOAT)0.0) {
                                Epsilon[j] = E[j] / Y[j][length_pdf_];

                                if (Epsilon[j] > epsilonlmax) epsilonlmax = Epsilon[j];

                                elp += E[j];
                            }
                            else {
                                if (E[j] < -R[j]) E[j] = -R[j]; 
                                
                                eln -= E[j];
                            }
                        }
                    }

                    Dl = elp / nl;

                    if ((Dl <= Dmin / W[l]) || (I == ItMax) || (nl <= (FLOAT)1.0)) {
                        // Enhanced component parameter estimation.

                        EnhancedEstimationPW(Y, nl, RigidTheta[l], LooseTheta[l]);

                        break;
                    }
                    else {
                        if (I == 1) {
                            A = ((FLOAT)1.0 - ar_) / ar_ / (Dl * W[l] - Dmin);
                        }

                        ar = (FLOAT)1.0 / (A * (Dl * W[l] - Dmin) + (FLOAT)1.0);

                        epsilonlmax *= (FLOAT)1.0 - ar;

                        for (j = 0; j < n_; j++) if (Epsilon[j] > epsilonlmax) {
                            Y[j][length_pdf_] -= E[j]; R[j] += E[j]; nl -= E[j];
                        }

                        if (eln > FLOAT_MIN) {
                            elp = elp / Dl - nl; if (eln > elp) f = elp / eln; else f = (FLOAT)1.0;

                            for (j = 0; j < n_; j++) if (E[j] < (FLOAT)0.0) {
                                E[j] *= f; Y[j][length_pdf_] -= E[j]; R[j] += E[j]; nl -= E[j];
                            }
                        }

                        W[l] = nl / n_;
                    }

                    I++;
                }

                // Outlier detection.

                for (j = 0; j < n_; j++) {
                    Error = ComponentDist(Y[j], LooseTheta[l], &fl, &Outlier);

                    if (Error) goto E0;
 
                    if (!Outlier) O[j] = 1;
                }

                // Moments calculation.

                Error = MomentsCalculation(LooseTheta[l], FirstM[l], SecondM[l]);

                if (Error) goto E0;

                c = ++l;

                r -= nl; nl = r; for (j = 0; j < n_; j++) Y[j][length_pdf_] = R[j];

                Stop = (c >= n_) || (c >= cmax_);

                if (Stop) break;
            }

            // Bayes classification of the remaining observations.

            Error = BayesClassificationPW(Y, c, W, LooseTheta, FirstM, SecondM);

            if (Error) goto E0;

            for (j = 0; j < n_; j++) Y[j][length_pdf_] = (FLOAT)1.0;

            Error = InformationCriterionPW(logV, Y, c, W, LooseTheta, &IC, &logL, &M, &D);
            
            if (Error) goto E0;

            if (IC < all_IC_[i]) all_IC_[i] = IC;

            if (IC < summary_.IC) {
                Found = 1;

                summary_.k = all_K_[i];

                memmove(summary_.h, h, length_pdf_ * sizeof(FLOAT));  

                memmove(summary_.ymin, ymin, length_pdf_ * sizeof(FLOAT)); 

                memmove(summary_.ymax, ymax, length_pdf_ * sizeof(FLOAT));
                
                summary_.IC = IC; summary_.logL = logL; summary_.M = M; summary_.c = c; 

                memmove(W_, W, c * sizeof(FLOAT));  

                for (j = 0; j < c; j++) {
                    Error = MixTheta_[j]->Memmove(LooseTheta[j]);

                    if (Error) goto E0;
                }
            }

            j = J - 1; opt_c[j] = c; opt_IC[j] = IC; opt_logL[j] = logL; opt_D[j] = D;

            Dmin *= c / (c + (FLOAT)1.0); J++; 

            if (Stop) break;
        }

        opt_length = J - 1;

        if (Found) {
            opt_length_ = opt_length;

            memmove(opt_c_, opt_c, opt_length_ * sizeof(int));  
            memmove(opt_IC_, opt_IC, opt_length_ * sizeof(FLOAT));  
            memmove(opt_logL_, opt_logL, opt_length_ * sizeof(FLOAT));  
            memmove(opt_D_, opt_D, opt_length_ * sizeof(FLOAT));  
        }
    }
    while (!Golden());

E0: if (O) free(O);

    if (opt_D) free(opt_D);

    if (opt_logL) free(opt_logL);

    if (opt_IC) free(opt_IC);

    if (opt_c) free(opt_c);
    
    if (SecondM) {
        for (i = 0; i < cmax_; i++) {
            if (SecondM[i]) free(SecondM[i]);
        }
         
        free(SecondM);
    }

    if (FirstM) {
        for (i = 0; i < cmax_; i++) {
            if (FirstM[i]) free(FirstM[i]);
        }
         
        free(FirstM);
    }

    if (LooseTheta) {
        for (i = 0; i < cmax_; i++) {
            if (LooseTheta[i]) delete LooseTheta[i];
        }

        delete[] LooseTheta;
    }

    if (RigidTheta) {
        for (i = 0; i < cmax_; i++) {
            if (RigidTheta[i]) delete RigidTheta[i];
        }

        delete[] RigidTheta;
    }

    if (W) free(W);

    if (Epsilon) free(Epsilon);

    if (E) free(E);

    if (R) free(R);

    if (ymax) free(ymax);

    if (ymin) free(ymin);

    if (h) free(h);

    if (Y) {
        for (i = 0; i < n_; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    return Error;
} // REBMIXPW

// REBMIX algorithm for histogram.

int Rebmix::REBMIXH()
{
    FLOAT                **Y = NULL;
    FLOAT                *y0 = NULL, *ymin = NULL, *ymax = NULL, *h = NULL;
    FLOAT                *R = NULL, *E = NULL, *Epsilon = NULL;
    FLOAT                *K = NULL;
    FLOAT                *W = NULL;
    CompnentDistribution **RigidTheta = NULL, **LooseTheta = NULL; 
    FLOAT                **FirstM = NULL, **SecondM = NULL;
    int                  opt_length;
    int                  *O = NULL, *opt_c = NULL;
    FLOAT                *opt_IC = NULL;       
    FLOAT                *opt_logL = NULL;     
    FLOAT                *opt_D = NULL;
    int                  c = 0, i, I, j, J, k, l, m, M;
    FLOAT                logV, Dmin, r, lognl, nl, elp, eln, epsilonlmax, logfl, fl, Dl, f, IC, logL, D;
    FLOAT                A = (FLOAT)0.0, ar;
    int                  Error = 0, Stop = 0, Found = 0, Outlier = 0, State = 0;

    // Allocation and initialisation.

    summary_.IC = FLOAT_MAX;

    summary_.h = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == summary_.h; if (Error) goto E0;

    summary_.y0 = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == summary_.y0; if (Error) goto E0;

    summary_.ymin = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == summary_.ymin; if (Error) goto E0;

    summary_.ymax = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == summary_.ymax; if (Error) goto E0;

    W_ = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));

    Error = NULL == W_; if (Error) goto E0;

    MixTheta_ = new CompnentDistribution* [(unsigned int)cmax_];

    Error = NULL == MixTheta_; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        MixTheta_[i] = new CompnentDistribution(this);

        Error = NULL == MixTheta_[i]; if (Error) goto E0;

        Error = MixTheta_[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;
    }

    opt_length_ = ItMax;

    opt_c_ = (int*)malloc(ItMax * sizeof(int));

    Error = NULL == opt_c_; if (Error) goto E0;

    opt_IC_ = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_IC_; if (Error) goto E0;

    opt_logL_ = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_logL_; if (Error) goto E0;

    opt_D_ = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_D_; if (Error) goto E0;

    all_length_ = K_[length_K_ - 1] - K_[0] + 1;

    all_I_ = (int*)calloc((size_t)all_length_, sizeof(int));

    Error = NULL == all_I_; if (Error) goto E0;

    all_K_ = (int*)calloc((size_t)all_length_, sizeof(int));

    Error = NULL == all_K_; if (Error) goto E0;

    for (i = 0; i < length_K_; i++) {
        all_K_[K_[i] - K_[0]] = K_[i];
    }

    additional_.Bracket = 1;

    all_IC_ = (FLOAT*)malloc(all_length_ * sizeof(FLOAT));

    Error = NULL == all_IC_; if (Error) goto E0;

    for (i = 0; i < all_length_; i++) {
        all_IC_[i] = FLOAT_MAX;
    }

    Y = (FLOAT**)malloc(n_ * sizeof(FLOAT*));

    Error = NULL == Y; if (Error) goto E0;

    for (i = 0; i < n_; i++) {
        Y[i] = (FLOAT*)malloc((length_pdf_ + 1) * sizeof(FLOAT));

        Error = NULL == Y[i]; if (Error) goto E0;
    }

    y0 = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == y0; if (Error) goto E0;

    ymin = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == ymin; if (Error) goto E0;

    if (ymin_ == NULL) {
        for (i = 0; i < length_pdf_; i++) {
            ymin[i] = Y_[0][i];

            for (j = 1; j < n_; j++) {
                if (Y_[j][i] < ymin[i]) ymin[i] = Y_[j][i];
            }

            if (Variables_[i] == vtContinuous) {
                ymin[i] -= Eps;
            }
        }
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            ymin[i] = ymin_[i];

            if (Variables_[i] == vtContinuous) {
                ymin[i] -= Eps;
            }
        }
    }

    ymax = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == ymax; if (Error) goto E0;

    if (ymax_ == NULL) {
        for (i = 0; i < length_pdf_; i++) {
            ymax[i] = Y_[0][i];

            for (j = 1; j < n_; j++) {
                if (Y_[j][i] > ymax[i]) ymax[i] = Y_[j][i];
            }

            if (Variables_[i] == vtContinuous) {
                ymax[i] += Eps;
            }
        }
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            ymax[i] = ymax_[i];

            if (Variables_[i] == vtContinuous) {
                ymax[i] += Eps;
            }
        }
    }

    h = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == h; if (Error) goto E0;

    R = (FLOAT*)malloc(n_ * sizeof(FLOAT));

    Error = NULL == R; if (Error) goto E0;

    E = (FLOAT*)malloc(n_ * sizeof(FLOAT));

    Error = NULL == E; if (Error) goto E0;

    Epsilon = (FLOAT*)malloc(n_ * sizeof(FLOAT));

    Error = NULL == Epsilon; if (Error) goto E0;

    K = (FLOAT*)malloc(n_ * sizeof(FLOAT));

    Error = NULL == K; if (Error) goto E0;

    W = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));

    Error = NULL == W; if (Error) goto E0;

    RigidTheta = new CompnentDistribution* [(unsigned int)cmax_];

    Error = NULL == RigidTheta; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        RigidTheta[i] = new CompnentDistribution(this);

        Error = NULL == RigidTheta[i]; if (Error) goto E0;

        Error = RigidTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;

        Error = RigidTheta[i]->Memmove(IniTheta_);

        if (Error) goto E0;
    }

    LooseTheta = new CompnentDistribution* [(unsigned int)cmax_];

    Error = NULL == LooseTheta; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        LooseTheta[i] = new CompnentDistribution(this);

        Error = NULL == LooseTheta[i]; if (Error) goto E0;

        Error = LooseTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;

        Error = LooseTheta[i]->Memmove(IniTheta_);

        if (Error) goto E0;
    }

    FirstM = (FLOAT**)malloc(cmax_ * sizeof(FLOAT*));

    Error = NULL == FirstM; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        FirstM[i] = (FLOAT*)malloc(length_theta_[0] * sizeof(FLOAT));

        Error = NULL == FirstM[i]; if (Error) goto E0;
    }

    SecondM = (FLOAT**)malloc(cmax_ * sizeof(FLOAT*));

    Error = NULL == SecondM; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        SecondM[i] = (FLOAT*)malloc(length_theta_[1] * sizeof(FLOAT));

        Error = NULL == SecondM[i]; if (Error) goto E0;
    }

    opt_length = ItMax;

    opt_c = (int*)malloc(ItMax * sizeof(int));

    Error = NULL == opt_c; if (Error) goto E0;

    opt_IC = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_IC; if (Error) goto E0;

    opt_logL = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_logL; if (Error) goto E0;

    opt_D = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == opt_D; if (Error) goto E0;

    O = (int*)malloc(n_ * sizeof(int));

    Error = NULL == O; if (Error) goto E0;

    Error = Initialize();

    if (Error) goto E0;

    do for (i = 0; i < all_length_; i++) if (all_K_[i] && (all_I_[i] == 0)) {
        // Preprocessing of observations.

        k = all_K_[i]; logV = (FLOAT)0.0; 
        
        for (j = 0; j < length_pdf_; j++) {
            switch (Variables_[j]) {
            case vtContinuous:
                h[j] = (ymax[j] - ymin[j]) / all_K_[i]; 
                
                if (y0_ == NULL) {
                    y0[j] = ymin[j] + (FLOAT)0.5 * h[j]; 
                }
                else {
                    y0[j] = y0_[j];
                }
                
                logV += (FLOAT)log(h[j]);

                break;
            case vtDiscrete:
                h[j] = (FLOAT)1.0; y0[j] = ymin[j];
            }
        }

        State = i > 0;

        Error = PreprocessingH(h, y0, &all_K_[i], Y, &State);

        if (Error) goto E0;

        all_I_[i] = 1;

        if (State == 2) {
            for (j = i; j < all_length_; j++) all_I_[j] = 2;
            
            goto E1;
        }

        for (j = 0; j < all_K_[i]; j++) K[j] = Y[j][length_pdf_];

        Found = 0; Dmin = (FLOAT)1.0; J = 1;

        // Outer loop.

        while (J <= ItMax) {
            l = 0; r = (FLOAT)n_; nl = (FLOAT)n_;

            memset(O, 0, all_K_[i] * sizeof(int));

            // Middle loop.
            
            while (nl / n_ > Dmin * l) {
                // Global mode detection.

                Error = GlobalModeH(&m, all_K_[i], Y, h, O);

                if (Error) goto E0;

                I = 1; W[l] = nl / n_; memset(R, 0, all_K_[i] * sizeof(FLOAT));

                // Inner loop.

                while (I <= ItMax) { 
                    // Rough component parameter estimation.

                    Error = RoughEstimationH(all_K_[i], Y, h, nl, m, RigidTheta[l], LooseTheta[l]);

                    if (Error) goto E0;

                    elp = eln = epsilonlmax = (FLOAT)0.0; lognl = (FLOAT)log(nl);

                    for (j = 0; j < all_K_[i]; j++) {
                        E[j] = Epsilon[j] = (FLOAT)0.0;

                        if ((Y[j][length_pdf_] > FLOAT_MIN) || (R[j] > FLOAT_MIN)) {
                            Error = LogComponentDist(Y[j], LooseTheta[l], &logfl, NULL);

                            if (Error) goto E0;

                            E[j] = Y[j][length_pdf_] - (FLOAT)exp(lognl + logfl + logV);

                            if (E[j] > (FLOAT)0.0) {
                                Epsilon[j] = E[j] / Y[j][length_pdf_]; 
                                
                                if (Epsilon[j] > epsilonlmax) epsilonlmax = Epsilon[j];

                                elp += E[j];
                            }
                            else {
                                if (E[j] < -R[j]) E[j] = -R[j]; 
                                
                                eln -= E[j];
                            }
                        }
                    }

                    Dl = elp / nl; 

                    if ((Dl <= Dmin / W[l]) || (I == ItMax) || (nl <= (FLOAT)1.0)) {
                        // Enhanced component parameter estimation.

                        EnhancedEstimationH(all_K_[i], Y, nl, RigidTheta[l], LooseTheta[l]);

                        break;
                    }
                    else {
                        if (I == 1) {
                            A = ((FLOAT)1.0 - ar_) / ar_ / (Dl * W[l] - Dmin);
                        }

                        ar = (FLOAT)1.0 / (A * (Dl * W[l] - Dmin) + (FLOAT)1.0);

                        epsilonlmax *= (FLOAT)1.0 - ar;

                        for (j = 0; j < all_K_[i]; j++) if (Epsilon[j] > epsilonlmax) {
                            Y[j][length_pdf_] -= E[j]; R[j] += E[j]; nl -= E[j];
                        }

                        if (eln > FLOAT_MIN) {
                            elp = elp / Dl - nl; if (eln > elp) f = elp / eln; else f = (FLOAT)1.0;

                            for (j = 0; j < all_K_[i]; j++) if (E[j] < (FLOAT)0.0) {
                                E[j] *= f; Y[j][length_pdf_] -= E[j]; R[j] += E[j]; nl -= E[j];
                            }
                        }

                        W[l] = nl / n_;
                    }

                    I++;
                }

                // Outlier detection.

                for (j = 0; j < all_K_[i]; j++) {
                    Error = ComponentDist(Y[j], LooseTheta[l], &fl, &Outlier);

                    if (Error) goto E0;
 
                    if (!Outlier) O[j] = 1;
                }

                // Moments calculation.

                Error = MomentsCalculation(LooseTheta[l], FirstM[l], SecondM[l]);

                if (Error) goto E0;

                c = ++l;

                r -= nl; nl = r; for (j = 0; j < all_K_[i]; j++) Y[j][length_pdf_] = R[j];

                Stop = (c >= all_K_[i]) || (c >= cmax_);

                if (Stop) break;
            }

            // Bayes classification of the remaining observations.

            Error = BayesClassificationH(all_K_[i], Y, c, W, LooseTheta, FirstM, SecondM);

            if (Error) goto E0;

            for (j = 0; j < all_K_[i]; j++) Y[j][length_pdf_] = K[j];

            Error = InformationCriterionH(logV, all_K_[i], Y, c, W, LooseTheta, &IC, &logL, &M, &D);
            
            if (Error) goto E0;

            if (IC < all_IC_[i]) all_IC_[i] = IC;

            if (IC < summary_.IC) {
                Found = 1; 

                summary_.k = k;

                memmove(summary_.h, h, length_pdf_ * sizeof(FLOAT));

                memmove(summary_.y0, y0, length_pdf_ * sizeof(FLOAT));

                memmove(summary_.ymin, ymin, length_pdf_ * sizeof(FLOAT)); 

                memmove(summary_.ymax, ymax, length_pdf_ * sizeof(FLOAT));  

                summary_.IC = IC; summary_.logL = logL; summary_.M = M; summary_.c = c; 

                memmove(W_, W, c * sizeof(FLOAT));  

                for (j = 0; j < c; j++) {
                    Error = MixTheta_[j]->Memmove(LooseTheta[j]);

                    if (Error) goto E0;
                }
            }

            j = J - 1; opt_c[j] = c; opt_IC[j] = IC; opt_logL[j] = logL; opt_D[j] = D;

            Dmin *= c / (c + (FLOAT)1.0); J++;

            if (Stop) break;
        }

        opt_length = J - 1;

        if (Found) {
            opt_length_ = opt_length;

            memmove(opt_c_, opt_c, opt_length_ * sizeof(int));  
            memmove(opt_IC_, opt_IC, opt_length_ * sizeof(FLOAT));  
            memmove(opt_logL_, opt_logL, opt_length_ * sizeof(FLOAT));  
            memmove(opt_D_, opt_D, opt_length_ * sizeof(FLOAT));
        }

E1:     all_K_[i] = k;
    }
    while (!Golden());

E0: if (O) free(O);

    if (opt_D) free(opt_D);

    if (opt_logL) free(opt_logL);

    if (opt_IC) free(opt_IC);

    if (opt_c) free(opt_c);
    
    if (SecondM) {
        for (i = 0; i < cmax_; i++) {
            if (SecondM[i]) free(SecondM[i]);
        }
         
        free(SecondM);
    }

    if (FirstM) {
        for (i = 0; i < cmax_; i++) {
            if (FirstM[i]) free(FirstM[i]);
        }
         
        free(FirstM);
    }

    if (LooseTheta) {
        for (i = 0; i < cmax_; i++) {
            if (LooseTheta[i]) delete LooseTheta[i];
        }

        delete[] LooseTheta;
    }

    if (RigidTheta) {
        for (i = 0; i < cmax_; i++) {
            if (RigidTheta[i]) delete RigidTheta[i];
        }

        delete[] RigidTheta;
    }

    if (W) free(W);

    if (K) free(K);

    if (Epsilon) free(Epsilon);

    if (E) free(E);

    if (R) free(R);

    if (h) free(h);

    if (ymax) free(ymax);

    if (ymin) free(ymin);

    if (y0) free(y0);

    if (Y) {
        for (i = 0; i < n_; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    return Error;
} // REBMIXH

// Reads data file.

int Rebmix::ReadDataFile()
{
    char line[65536];
    char *pchar = NULL;
    int  i, j, BufSize = 0;
    FILE *fp = NULL;
    int  Error = 0;

    if ((fp = fopen(curr_, "r")) == NULL) {
        Error = 1; goto E0;
    }

    Y_ = (FLOAT**)malloc((BufSize + BufInc) * sizeof(FLOAT*));

    Error = NULL == Y_; if (Error) goto E0;

    for (i = BufSize; i < BufSize + BufInc; i++) {
        Y_[i] = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

        Error = NULL == Y_[i]; if (Error) goto E0;
    }

    BufSize += BufInc;

    n_ = 0;

S0: while (fgets(line, 2048, fp) != NULL) {
        pchar = strtok(line, "\n");

        if (!pchar) goto S0;

        j = 0;
        for (i = 0; i < (int)strlen(pchar); i++) {
            if (pchar[i] == ',') {
                line[j] = '.'; j++;
            }
            else
            if (pchar[i] != ' ') {
                line[j] = pchar[i]; j++;
            }
        }

        line[j] = '\0';

        if (!j) goto S0;

        pchar = strtok(pchar, "\t");

        if (n_ == BufSize) {
            Y_ = (FLOAT**)realloc(Y_, (BufSize + BufInc) * sizeof(FLOAT*));

            Error = NULL == Y_; if (Error) goto E0;

            for (i = BufSize; i < BufSize + BufInc; i++) {
                Y_[i] = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

                Error = NULL == Y_[i]; if (Error) goto E0;
            }

            BufSize += BufInc;
        }

        i = 0;
        while (pchar) {
            Y_[n_][i] = (FLOAT)atof(pchar); pchar = strtok(NULL, "\t"); i++;
        }

        n_++;
    }

    for (i = n_; i < BufSize; i++) {
        if (Y_[i]) free(Y_[i]); 
    }

    Y_ = (FLOAT**)realloc(Y_, n_ * sizeof(FLOAT*));

    Error = NULL == Y_; if (Error) goto E0;

    X_ = (FLOAT**)malloc(n_ * sizeof(FLOAT*));

    Error = NULL == X_; if (Error) goto E0;

    for (i = 0; i < n_; i++) {
        X_[i] = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

        Error = NULL == X_[i]; if (Error) goto E0;
    }

E0: if (fp) fclose(fp);

    return Error;
} // ReadDataFile

// Writes data file.

int Rebmix::WriteDataFile()
{
    int  i, j, k;
    char line[65536];
    char mode[2];
    char path[FILENAME_MAX];
    char ext[FILENAME_MAX];
    char *pchar = NULL;
    FILE *fp0 = NULL, *fp1 = NULL;
    int  Error = 0;

    if (curr_ == open_[0])
        strcpy(mode, "w");
    else
        strcpy(mode, "a");

    strcpy(path, save_); 
        
    pchar = strrchr(path, '.'); 
        
    if (pchar) {
        strcpy(ext, pchar); pchar[0] = '\0';
    }
    else {
        strcpy(ext, "");
    }
        
    sprintf(line, "%s%s%s", path, "_1", ext);

    if ((fp0 = fopen(line, mode)) == NULL) {
        Error = 1; goto E0;
    }

    strcpy(path, save_); 
        
    pchar = strrchr(path, '.'); 
        
    if (pchar) {
        strcpy(ext, pchar); pchar[0] = '\0';
    }
    else {
        strcpy(ext, "");
    }
        
    sprintf(line, "%s%s%s", path, "_2", ext);

    if ((fp1 = fopen(line, mode)) == NULL) {
        Error = 1; goto E0;
    }

    if (!strcmp(mode, "w")) {
        fprintf(fp0, "%s\t%s\t%s\t%s\t%s\t%s\t%s", "Dataset",
                                                   "Preprocessing",
                                                   "cmax",
                                                   "Criterion",
                                                   "ar",
                                                   "Restraints", 
                                                   "c");

        switch (Preprocessing_) {
        case poHistogram:
            fprintf(fp0, "\t%s", "k");

            for (i = 0; i < length_pdf_; i++) {
                if (length_pdf_ == 1)
                    sprintf(line, "%s", "y0");
                else
                    sprintf(line, "%s%d", "y0", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }

            for (i = 0; i < length_pdf_; i++) {
                if (length_pdf_ == 1)
                    sprintf(line, "%s", "h");
                else
                    sprintf(line, "%s%d", "h", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }

            break;
        case poParzenWindow:
            fprintf(fp0, "\t%s", "k");

            for (i = 0; i < length_pdf_; i++) {
                if (length_pdf_ == 1)
                    sprintf(line, "%s", "h");
                else
                    sprintf(line, "%s%d", "h", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }

            break;
        case poKNearestNeighbour:
            fprintf(fp0, "\t%s", "k");

            for (i = 0; i < length_pdf_; i++) {
                if (length_pdf_ == 1)
                    sprintf(line, "%s", "h");
                else
                    sprintf(line, "%s%d", "h", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }
        }

        for (i = 0; i < length_pdf_; i++) {
            if (length_pdf_ == 1)
                sprintf(line, "%s", "ymin");
            else
                sprintf(line, "%s%d", "ymin", i + 1);
                  
            fprintf(fp0, "\t%s", line);
        }

        for (i = 0; i < length_pdf_; i++) {
            if (length_pdf_ == 1)
                sprintf(line, "%s", "ymax");
            else
                sprintf(line, "%s%d", "ymax", i + 1);
                  
            fprintf(fp0, "\t%s", line);
        }

        fprintf(fp0, "\t%s\t%s\n", "IC",
                                   "logL");

        fprintf(fp1, "%s\t%s", "Dataset",
                               "w");

        for (i = 0; i < length_pdf_; i++) {
            if (length_pdf_ == 1)
                fprintf(fp1, "\t%s", "pdf");
            else
                fprintf(fp1, "\t%s%d", "pdf", i + 1);
        }

        for (j = 0; j < Min(2, length_Theta_); j++) {
            for (k = 0; k < length_theta_[j]; k++) {
                fprintf(fp1, "\t%s%d%s%d", "theta", j + 1, ".", k + 1);
            }
        }

        fprintf(fp1, "\n");
    }

    strcpy(path, curr_); 

    pchar = strrchr(path, '\\');

    if (!pchar) {
        pchar = strrchr(path, '/');
    }

    if (pchar) {
        strcpy(path, pchar + 1);
    }

    pchar = strrchr(path, '.'); 
        
    if (pchar) pchar[0] = '\0';

    fprintf(fp0, "%s", path);

    switch (Preprocessing_) {
    case poHistogram: 
        strcpy(line, "histogram");

        break;
    case poParzenWindow: 
        strcpy(line, "Parzen window");

        break;
    case poKNearestNeighbour:
        strcpy(line, "k-nearest neighbour");
    }

    fprintf(fp0, "\t%s", line);

    fprintf(fp0, "\t%d", cmax_);

    switch (Criterion_) {
    case icAIC:
        strcpy(line, "AIC");

        break;
    case icAIC3:
        strcpy(line, "AIC3");

        break;
    case icAIC4:
        strcpy(line, "AIC4");

        break;
    case icAICc:
        strcpy(line, "AICc");

        break;
    case icBIC:
        strcpy(line, "BIC");

        break;
    case icCAIC:
        strcpy(line, "CAIC");

        break;
    case icHQC:
        strcpy(line, "HQC");

        break;
    case icMDL2:
        strcpy(line, "MDL2");

        break;
    case icMDL5:
        strcpy(line, "MDL5");

        break;
    case icAWE:
        strcpy(line, "AWE");

        break;
    case icCLC:
        strcpy(line, "CLC");

        break;
    case icICL:
        strcpy(line, "ICL");

        break;
    case icPC:
        strcpy(line, "PC");

        break;
    case icICLBIC:
        strcpy(line, "ICL-BIC");

        break;
    case icD:
        strcpy(line, "D");
        
        break;
    case icSSE:
        strcpy(line, "SSE");
    }

    fprintf(fp0, "\t%s", line);

    fprintf(fp0, "\t%E", ar_);

    switch (Restraints_) {
    case rtRigid:
        strcpy(line, "rigid");

        break;
    case rtLoose:
        strcpy(line, "loose");
    }

    fprintf(fp0, "\t%s", line);

    fprintf(fp0, "\t%d", summary_.c);

    switch (Preprocessing_) {
    case poHistogram:
        fprintf(fp0, "\t%d", summary_.k);

        for (i = 0; i < length_pdf_; i++) {
            fprintf(fp0, "\t%E", summary_.y0[i]);
        }

        for (i = 0; i < length_pdf_; i++) {
            fprintf(fp0, "\t%E", summary_.h[i]);
        }

        break;
    case poParzenWindow:
        fprintf(fp0, "\t%d", summary_.k);

        for (i = 0; i < length_pdf_; i++) {
            fprintf(fp0, "\t%E", summary_.h[i]);
        }

        break;
    case poKNearestNeighbour:
        fprintf(fp0, "\t%d", summary_.k);

        for (i = 0; i < length_pdf_; i++) {
            fprintf(fp0, "\t%E", summary_.h[i]);
        }
    }

    for (i = 0; i < length_pdf_; i++) {
        fprintf(fp0, "\t%E", summary_.ymin[i]);
    }

    for (i = 0; i < length_pdf_; i++) {
        fprintf(fp0, "\t%E", summary_.ymax[i]);
    }

    fprintf(fp0, "\t%E\t%E\n", summary_.IC,
                               summary_.logL);

    for (i = 0; i < summary_.c; i++) {
        fprintf(fp1, "%s\t%E", path, W_[i]);

        for (j = 0; j < length_pdf_; j++) switch (MixTheta_[i]->pdf_[j]) {
            case pfNormal:
                fprintf(fp1, "\t%s", "normal");

                break;
            case pfLognormal:
                fprintf(fp1, "\t%s", "lognormal");

                break;
            case pfWeibull:
                fprintf(fp1, "\t%s", "Weibull");

                break;
            case pfGamma:
                fprintf(fp1, "\t%s", "gamma");

                break;
            case pfvonMises:
                fprintf(fp1, "\t%s", "vonMises");

                break;
            case pfBinomial:
                fprintf(fp1, "\t%s", "binomial");

                break;
            case pfPoisson:
                fprintf(fp1, "\t%s", "Poisson");

                break;
            case pfDirac:
                fprintf(fp1, "\t%s", "Dirac");
            
                break;
            case pfUniform:
                fprintf(fp1, "\t%s", "uniform");
        }

        for (j = 0; j < Min(2, length_Theta_); j++) {
            for (k = 0; k < length_theta_[j]; k++) {
                fprintf(fp1, "\t%E", MixTheta_[i]->Theta_[j][k]);
            }
        }

        fprintf(fp1, "\n");
    }

E0: if (fp0) fclose(fp0);
    if (fp1) fclose(fp1);

    if (all_IC_) {
        free(all_IC_); all_IC_ = NULL;
    }

    if (all_K_) {
        free(all_K_); all_K_ = NULL;
    }

    if (opt_D_) {
        free(opt_D_); opt_D_ = NULL;
    }

    if (opt_logL_) {
        free(opt_logL_); opt_logL_ = NULL;
    }

    if (opt_IC_) {
        free(opt_IC_); opt_IC_ = NULL;
    }

    if (opt_c_) {
        free(opt_c_); opt_c_ = NULL;
    }

    if (summary_.h) {
        free(summary_.h); summary_.h = NULL;
    }

    if (summary_.y0) {
        free(summary_.y0); summary_.y0 = NULL;
    }

    if (MixTheta_) {
        for (i = 0; i < cmax_; i++) {
            if (MixTheta_[i]) delete MixTheta_[i];
        }

        delete[] MixTheta_; MixTheta_ = NULL;
    }

    if (W_) {
        free(W_); W_ = NULL;
    }

    if (X_) {
        for (i = 0; i < n_; i++) {
            if (X_[i]) free(X_[i]);
        }

        free(X_); X_ = NULL;
    }

    if (Y_) {
        for (i = 0; i < n_; i++) {
            if (Y_[i]) free(Y_[i]);
        }

        free(Y_); Y_ = NULL;
    }

    return Error;
} // WriteDataFile

// Runs template file.

int Rebmix::RunTemplateFile(char *file)
{
    char  line[65536], ident[65536], list[65536];
    char  *pchar = NULL, *rchar = NULL;
    int   i, imin, imax, iinc, j, k, isI;
    FLOAT isF;
    FILE  *fp = NULL;
    int   Error = 0;

    #if (_REBMIXEXE)
    printf("REBMIX Version 2.10.2\n");
    #endif

    if ((fp = fopen(file, "r")) == NULL) {
        Error = 1; goto E0;
    }

S0: while (fgets(line, 2048, fp) != NULL) {
        pchar = strtok(line, "\n"); 
        
        pchar = strtok(pchar, "=");

        if (!pchar) goto S0;

        j = 0;

        for (i = 0; i < (int)strlen(pchar); i++) {
            if (pchar[i] != ' ') {
                ident[j] = (char)toupper(pchar[i]); j++;
            }
        }

        ident[j] = '\0';

        j = 0; list[j] = '\0'; imin = 0; imax = 0;

        while((pchar = strtok(NULL, ",")) != NULL) {
            if (!strcmp(ident, "DATASET") || !strcmp(ident, "SAVE")) {
                for (i = 0; pchar[i] != '\0'; i++) {
                    list[j] = pchar[i]; j++;

                    if (pchar[i] != ' ') {
                        imax = i; if (!imin) imin = i;
                    }
                }

                j = imax + 1 - imin;

                for (i = 0; i < j; i++) {
                    list[i] = list[imin + i];
                }
            }
            else {
                for (i = 0; pchar[i] != '\0'; i++) if (pchar[i] != '[' && pchar[i] != ']') {
                    if (pchar[i] == ' ') {
                    }
                    else {
                        list[j] = (char)toupper(pchar[i]); j++;
                    }
                }
            }

            list[j] = '\t'; j++;
        }

        if (!j) goto S0; else list[j - 1] = '\0';

        pchar = strtok(list, "\t");

        if (!strcmp(ident, "RUN")) {
            for (k = 0; k < o_; k++) {
                curr_ = open_[k];

                Error = ReadDataFile();

                if (Error) goto E0;

                #if (_REBMIXEXE)
                printf("Dataset = %s\n", curr_);
                #endif

                Error = REBMIX();

                if (Error) goto E0;

                Error = WriteDataFile();

                if (Error) goto E0;
            }
        }
        else
        if (!strcmp(ident, "DATASET")) {
            open_ = (char**)realloc(open_, (o_ + 1) * sizeof(char*));

            Error = NULL == open_; if (Error) goto E0;

            open_[o_] = (char*)malloc((strlen(pchar) + 1) * sizeof(char));

            Error = NULL == open_[o_]; if (Error) goto E0;

            strcpy(open_[o_], pchar); o_++;
        }
        else
        if (!strcmp(ident, "PREPROCESSING")) {
            if (!strcmp(pchar, "HISTOGRAM"))
                Preprocessing_ = poHistogram;
            else
            if (!strcmp(pchar, "PARZENWINDOW"))
                Preprocessing_ = poParzenWindow;
            else
            if (!strcmp(pchar, "K-NEARESTNEIGHBOUR"))
                Preprocessing_ = poKNearestNeighbour;
            else {
                Error = 1; goto E0;
            }
        } 
        else
        if (!strcmp(ident, "CMAX")) {
            cmax_ = isI = (int)atol(pchar);

            Error = isI <= 0; if (Error) goto E0;
        } 
        else
        if (!strcmp(ident, "CRITERION")) {
            if (!strcmp(pchar, "AIC"))
                Criterion_ = icAIC;
            else
            if (!strcmp(pchar, "AIC3"))
                Criterion_ = icAIC3;
            else
            if (!strcmp(pchar, "AIC4"))
                Criterion_ = icAIC4;
            else
            if (!strcmp(pchar, "AICC"))
                Criterion_ = icAICc;
            else
            if (!strcmp(pchar, "BIC"))
                Criterion_ = icBIC;
            else
            if (!strcmp(pchar, "CAIC"))
                Criterion_ = icCAIC;
            else
            if (!strcmp(pchar, "HQC"))
                Criterion_ = icHQC;
            else
            if (!strcmp(pchar, "MDL2"))
                Criterion_ = icMDL2;
            else
            if (!strcmp(pchar, "MDL5"))
                Criterion_ = icMDL5;
            else
            if (!strcmp(pchar, "AWE"))
                Criterion_ = icAWE;
            else
            if (!strcmp(pchar, "CLC"))
                Criterion_ = icCLC;
            else
            if (!strcmp(pchar, "ICL"))
                Criterion_ = icICL;
            else
            if (!strcmp(pchar, "PC"))
                Criterion_ = icPC;
            else
            if (!strcmp(pchar, "ICL-BIC"))
                Criterion_ = icICLBIC;
            else
            if (!strcmp(pchar, "D"))
                Criterion_ = icD;
            else
            if (!strcmp(pchar, "SSE"))
                Criterion_ = icSSE;
            else {
                Error = 1; goto E0;
            }
        }
        else
        if (!strcmp(ident, "VARIABLES")) {
            i = 0;

            while (pchar) {
                Variables_ = (VariablesType_e*)realloc(Variables_, (i + 1) * sizeof(VariablesType_e));

                Error = NULL == Variables_; if (Error) goto E0;

                if (!strcmp(pchar, "CONTINUOUS"))
                    Variables_[i] = vtContinuous;
                else
                if (!strcmp(pchar, "DISCRETE"))
                    Variables_[i] = vtDiscrete;
                else {
                    Error = 1; goto E0;
                }
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((length_pdf_ > 0) && (length_pdf_ != i)) {
                Error = 1; goto E0;
            }
            else {
                length_pdf_ = i;
            }
        } 
        else
        if (!strcmp(ident, "LENGTHPDF")) {
            length_pdf_ = isI = (int)atol(pchar);

            Error = isI < 1; if (Error) goto E0;
        }
        else
        if (!strcmp(ident, "LENGTHTHETA")) {
            i = 0;

            while (pchar) {
                length_theta_ = (int*)realloc(length_theta_, (i + 1) * sizeof(int));

                Error = NULL == length_theta_; if (Error) goto E0;

                length_theta_[i] = isI = (int)atol(pchar);

                Error = isI == 0; if (Error) goto E0;

                pchar = strtok(NULL, "\t"); ++i;
            }

            length_Theta_ = i;

            IniTheta_ = new CompnentDistribution(this);

            Error = NULL == IniTheta_; if (Error) goto E0;

            Error = IniTheta_->Realloc(length_pdf_, length_Theta_, length_theta_);

            if (Error) goto E0;
        }
        else
        if (!strcmp(ident, "PDF")) {
            i = 0;

            while (pchar) {
                if (!strcmp(pchar, "NORMAL"))
                    IniTheta_->pdf_[i] = pfNormal;
                else
                if (!strcmp(pchar, "LOGNORMAL"))
                    IniTheta_->pdf_[i] = pfLognormal;
                else
                if (!strcmp(pchar, "WEIBULL"))
                    IniTheta_->pdf_[i] = pfWeibull;
                else
                if (!strcmp(pchar, "GAMMA"))
                    IniTheta_->pdf_[i] = pfGamma;
                else
                if (!strcmp(pchar, "VONMISES"))
                    IniTheta_->pdf_[i] = pfvonMises;
                else
                if (!strcmp(pchar, "BINOMIAL"))
                    IniTheta_->pdf_[i] = pfBinomial;
                else
                if (!strcmp(pchar, "POISSON"))
                    IniTheta_->pdf_[i] = pfPoisson;
                else
                if (!strcmp(pchar, "DIRAC"))
                    IniTheta_->pdf_[i] = pfDirac;
                else
                if (!strcmp(pchar, "UNIFORM"))
                    IniTheta_->pdf_[i] = pfUniform;
                else {
                    Error = 1; goto E0;
                }
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((length_pdf_ > 0) && (length_pdf_ != i)) {
                Error = 1; goto E0;
            }
        } 
        else
        if (!strcmp(ident, "THETA1")) {
            i = 0;

            while (pchar) {
                IniTheta_->Theta_[0] = (FLOAT*)realloc(IniTheta_->Theta_[0], (i + 1) * sizeof(FLOAT));

                Error = NULL == IniTheta_->Theta_[0]; if (Error) goto E0;

                IniTheta_->Theta_[0][i] = (FLOAT)atof(pchar);
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((length_theta_[0] > 0) && (length_theta_[0] != i)) {
                Error = 1; goto E0;
            }
        } 
        else
        if (!strcmp(ident, "THETA2")) {
            i = 0;

            while (pchar) {
                IniTheta_->Theta_[1] = (FLOAT*)realloc(IniTheta_->Theta_[1], (i + 1) * sizeof(FLOAT));

                Error = NULL == IniTheta_->Theta_[1]; if (Error) goto E0;

                IniTheta_->Theta_[1][i] = (FLOAT)atof(pchar);

                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((length_theta_[1] > 0) && (length_theta_[1] != i)) {
                Error = 1; goto E0;
            }
        } 
        else
        if (!strcmp(ident, "K")) {
            i = 0;

            while (pchar != NULL) {
                if ((rchar = strrchr(pchar, '-')) != NULL) {
                    imin = (int)atol(pchar); imax = (int)atol(rchar + 1);

                    if (imin > imax) {
                        j = imin; imin = imax; imax = j;
                    }

                    if ((rchar = strrchr(pchar, '+')) != NULL)
                        iinc = (int)atol(rchar + 1);
                    else
                        iinc = 1;

                    K_ = (int*)realloc(K_, (i + (imax - imin) / iinc + 1) * sizeof(int));

                    Error = NULL == K_; if (Error) goto E0;

                    for (j = imin; j <= imax; j += iinc) {
                        K_[i] = isI = j; 

                        Error = isI <= 0; if (Error) goto E0;
                        
                        i++;
                    }
                
                    length_K_ = i;
                }
                else {
                    K_ = (int*)realloc(K_, (i + 1) * sizeof(int));

                    Error = NULL == K_; if (Error) goto E0;

                    K_[i] = isI = (int)atol(pchar);

                    Error = isI <= 0; if (Error) goto E0;
                
                    length_K_ = ++i;
                }

                pchar = strtok(NULL, "\t"); 
            }
        } 
        else
        if (!strcmp(ident, "Y0")) {
            i = 0;

            while (pchar) {
                y0_ = (FLOAT*)realloc(y0_, (i + 1) * sizeof(FLOAT));

                Error = NULL == y0_; if (Error) goto E0;

                y0_[i] = (FLOAT)atof(pchar);
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((length_pdf_ > 0) && (length_pdf_ != i)) {
                Error = 1; goto E0;
            }
            else {
                length_pdf_ = i;
            }
        } 
        else
        if (!strcmp(ident, "YMIN")) {
            i = 0;

            while (pchar) {
                ymin_ = (FLOAT*)realloc(ymin_, (i + 1) * sizeof(FLOAT));

                Error = NULL == ymin_; if (Error) goto E0;

                ymin_[i] = (FLOAT)atof(pchar);
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((length_pdf_ > 0) && (length_pdf_ != i)) {
                Error = 1; goto E0;
            }
            else {
                length_pdf_ = i;
            }
        } 
        else
        if (!strcmp(ident, "YMAX")) {
            i = 0;

            while (pchar) {
                ymax_ = (FLOAT*)realloc(ymax_, (i + 1) * sizeof(FLOAT));

                Error = NULL == ymax_; if (Error) goto E0;

                ymax_[i] = (FLOAT)atof(pchar);
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((length_pdf_ > 0) && (length_pdf_ != i)) {
                Error = 1; goto E0;
            }
            else {
                length_pdf_ = i;
            }
        } 
        else
        if (!strcmp(ident, "AR")) {
            ar_ = isF = (FLOAT)atof(pchar);

            Error = (isF <= (FLOAT)0.0) || (isF > (FLOAT)1.0); if (Error) goto E0;

        } 
        else
        if (!strcmp(ident, "RESTRAINTS")) {
            if (!strcmp(pchar, "RIGID"))
                Restraints_ = rtRigid;
            else
            if (!strcmp(pchar, "LOOSE"))
                Restraints_ = rtLoose;
            else {
                Error = 1; goto E0;
            }
        } 
        else
        if (!strcmp(ident, "SAVE")) {
            save_ = (char*)realloc(save_, (strlen(pchar) + 1) * sizeof(char));

            Error = NULL == save_; if (Error) goto E0;

            strcpy(save_, pchar);
        }
    }

E0: if (fp) fclose(fp);

    return Error;
} // RunTemplateFile

// REBMIX algorithm.

int Rebmix::REBMIX()
{
    int  Error = 0;

    switch (Preprocessing_) {
    case poHistogram:
        Error = REBMIXH();

        if (Error) goto E0;

        break;
    case poParzenWindow:
        Error = REBMIXPW();

        if (Error) goto E0;

        break;
    case poKNearestNeighbour:
        Error = REBMIXKNN();

        if (Error) goto E0;
    }

E0: return Error;
} // REBMIX