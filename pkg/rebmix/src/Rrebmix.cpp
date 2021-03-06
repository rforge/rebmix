#include "rngmixf.h"

#include <math.h>

extern "C" {

// Runs RNGMIX in R.

void RRNGMIX(int    *IDum,         // Random seed.
             int    *d,            // Number of independent random variables.
             int    *c,            // Number of components.
             int    *N,            // Numbers of observations.
             int    *length_pdf,   // Length of pdf.
             int    *length_Theta, // Length of Theta.
             int    *length_theta, // Length of Theta[i].
             char   **pdf,         // Parametric family types.
             double *Theta,        // Component parameters.
             int    *n,            // Number of observations.
             double *Y,            // Dataset.
             int    *Z,            // Component membership.
             int    *Error)        // Error code.
{
    Rngmix *rngmix = NULL;
    int    i, j, k, l;

    rngmix = new Rngmix;

    *Error = NULL == rngmix; if (*Error) goto E0;

    rngmix->IDum_ = *IDum;
    rngmix->length_pdf_ = *d;
    rngmix->c_ = *c;

    rngmix->N_ = (int*)malloc(rngmix->c_ * sizeof(int));

    *Error = NULL == rngmix->N_; if (*Error) goto E0;

    for (i = 0; i < rngmix->c_; i++) rngmix->N_[i] = N[i];

    rngmix->IniTheta_ = new CompnentDistribution(rngmix);

    *Error = NULL == rngmix->IniTheta_; if (*Error) goto E0;

    rngmix->length_pdf_ = *length_pdf;

    rngmix->length_Theta_ = *length_Theta;

    rngmix->length_theta_ = (int*)malloc(rngmix->length_Theta_ * sizeof(int));

    *Error = NULL == rngmix->length_theta_; if (*Error) goto E0;

    *Error = rngmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    if (*Error) goto E0;

    for (i = 0; i < rngmix->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rngmix->IniTheta_->pdf_[i] = pfNormal;
        }
        else
        if (!strcmp(pdf[i], "lognormal")) {
            rngmix->IniTheta_->pdf_[i] = pfLognormal;
        }
        else
        if (!strcmp(pdf[i], "Weibull")) {
            rngmix->IniTheta_->pdf_[i] = pfWeibull;
        }
        else
        if (!strcmp(pdf[i], "gamma")) {
            rngmix->IniTheta_->pdf_[i] = pfGamma;
        }
        else
        if (!strcmp(pdf[i], "Gumbel")) {
            rngmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rngmix->IniTheta_->pdf_[i] = pfvonMises;
        }
        else
        if (!strcmp(pdf[i], "binomial")) {
            rngmix->IniTheta_->pdf_[i] = pfBinomial;
        }
        else
        if (!strcmp(pdf[i], "Poisson")) {
            rngmix->IniTheta_->pdf_[i] = pfPoisson;
        }
        else
        if (!strcmp(pdf[i], "Dirac")) {
            rngmix->IniTheta_->pdf_[i] = pfDirac;
        }
        else
        if (!strcmp(pdf[i], "uniform")) {
            rngmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    rngmix->MixTheta_ = new CompnentDistribution* [(unsigned int)rngmix->c_];

    *Error = NULL == rngmix->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < rngmix->c_; i++) {
        rngmix->MixTheta_[i] = new CompnentDistribution(rngmix);

        *Error = NULL == rngmix->MixTheta_[i]; if (*Error) goto E0;

        *Error = rngmix->MixTheta_[i]->Realloc(rngmix->length_pdf_, rngmix->length_Theta_, rngmix->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < rngmix->c_; i++) {
        for (j = 0; j < rngmix->length_pdf_; j++) {
            rngmix->MixTheta_[i]->pdf_[j] = rngmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rngmix->length_Theta_; j++) if (rngmix->IniTheta_->Theta_[j]) {
        for (k = 0; k < rngmix->c_; k++) {
            for (l = 0; l < rngmix->length_theta_[j]; l++) {
                rngmix->MixTheta_[k]->Theta_[j][l] = Theta[i];

                i++;
            }
        }
    }

    *Error = rngmix->RNGMIX();

    if (*Error) goto E0;

    *n = rngmix->n_; i = 0;

    for (j = 0; j < rngmix->length_pdf_; j++) {
        for (k = 0; k < rngmix->n_; k++) {
            Y[i] = rngmix->Y_[j][k]; i++;
        }
    }

    for (i = 0; i < rngmix->n_; i++) {
        Z[i] = rngmix->Z_[i];
    }

E0: if (rngmix) delete rngmix;
} // RRNGMIX

// Runs REBMIX in R.

void RREBMIX(char   **Preprocessing, // Preprocessing type.
             int    *cmax,           // Maximum number of components.
             int    *cmin,           // Minimum number of components.
             char   **Criterion,     // Information criterion type.
             int    *d,              // Number of independent random variables.
             char   **Variables,     // Types of variables.
             int    *length_pdf,     // Length of pdf.
             char   **pdf,           // Parametric family types.
             int    *length_Theta,   // Length of Theta.
             int    *length_theta,   // Length of Theta[i].
             double *Theta,          // Component parameters.
             int    *length_K,       // Length of K.
             int    *K,              // Numbers of bins v or numbers of nearest neighbours k.
             int    *length_y0,      // Length of y0.
             double *y0,             // Origins.
             int    *length_ymin,    // Length of ymin.
             double *ymin,           // Minimum observations.
             int    *length_ymax,    // Length of ymax.
             double *ymax,           // Maximum observations.
             double *ar,             // Acceleration rate.
             char   **Restraints,    // Restraints type.
             int    *n,              // Number of observations.
             double *Y,              // Dataset.
/// Panic Branislav.
             char   **EMStrategy,       // Strategy for EM algorithm.
             char   **EMVariant,        // EM algorithm variant.
             char   **EMAcceleration,   // Acceleration for the standard EM algorithm.
             double *EMTolerance,       // Tolerance for EM algortihm.
             double *EMAccelerationMul, // Acceleration rate for Em algorithm.
             int    *EMMaxIter,         // Maximum number of iterations in EM algorithm.
	         int    *EMK,               // Number of bins for histogram EM algorithm.
             int    *n_iter,            // Number of iterations for optimal case.
             int    *n_iter_sum,        // Number of iterations in whole run.
/// End
             int    *summary_k,      // Optimal v or optimal k.
             double *summary_h,      // Optimal class widths of length d.
             double *summary_y0,     // Optimal origins of length d.
             double *summary_ymin,   // Optimal minimum observations of length d.
             double *summary_ymax,   // Optimal maximum observations of length d.
             double *summary_IC,     // Optimal information criterion.
             double *summary_logL,   // Log-likelihood.
             int    *summary_M,      // Degrees of freedom.
             int    *summary_c,      // Optimal number of components.
             double *W,              // Component weights.
             double *theta1,         // Component parameters.
             double *theta2,         // Component parameters.
             int    *opt_length,     // Length of opt_c, opt_IC, opt_logL and opt_D.
             int    *opt_c,          // Numbers of components for optimal v or for optimal k.
             double *opt_IC,         // Information criteria for optimal v or for optimal k.
             double *opt_logL,       // Log-likelihoods for optimal v or for optimal k.
             double *opt_D,          // Totals of positive relative deviations for optimal v or for optimal k.
             int    *all_length,     // Length of all_K and all_IC.
             int    *all_K,          // All processed numbers of bins v or all processed numbers of nearest neighbours k.
             double *all_IC,         // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.
             int    *Error)          // Error code.
{
    Rebmix *rebmix = NULL;
    int    i, j, k, l;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    if (!strcmp(Preprocessing[0], "histogram")) {
        rebmix->Preprocessing_ = poHistogram;
    }
    else
    if (!strcmp(Preprocessing[0], "kernel density estimation")) {
        rebmix->Preprocessing_ = poKDE;
    }
    else
    if (!strcmp(Preprocessing[0], "k-nearest neighbour")) {
        rebmix->Preprocessing_ = poKNearestNeighbour;
    }
    else {
        *Error = 1; goto E0;
    }

    rebmix->cmax_ = *cmax;

    rebmix->cmin_ = *cmin;

    if (!strcmp(Criterion[0], "AIC"))
        rebmix->Criterion_ = icAIC;
    else
    if (!strcmp(Criterion[0], "AIC3"))
        rebmix->Criterion_ = icAIC3;
    else
    if (!strcmp(Criterion[0], "AIC4"))
        rebmix->Criterion_ = icAIC4;
    else
    if (!strcmp(Criterion[0], "AICc"))
        rebmix->Criterion_ = icAICc;
    else
    if (!strcmp(Criterion[0], "BIC"))
        rebmix->Criterion_ = icBIC;
    else
    if (!strcmp(Criterion[0], "CAIC"))
        rebmix->Criterion_ = icCAIC;
    else
    if (!strcmp(Criterion[0], "HQC"))
        rebmix->Criterion_ = icHQC;
    else
    if (!strcmp(Criterion[0], "MDL2"))
        rebmix->Criterion_ = icMDL2;
    else
    if (!strcmp(Criterion[0], "MDL5"))
        rebmix->Criterion_ = icMDL5;
    else
    if (!strcmp(Criterion[0], "AWE"))
        rebmix->Criterion_ = icAWE;
    else
    if (!strcmp(Criterion[0], "CLC"))
        rebmix->Criterion_ = icCLC;
    else
    if (!strcmp(Criterion[0], "ICL"))
        rebmix->Criterion_ = icICL;
    else
    if (!strcmp(Criterion[0], "PC"))
        rebmix->Criterion_ = icPC;
    else
    if (!strcmp(Criterion[0], "ICL-BIC"))
        rebmix->Criterion_ = icICLBIC;
    else
    if (!strcmp(Criterion[0], "D"))
        rebmix->Criterion_ = icD;
    else
    if (!strcmp(Criterion[0], "SSE"))
        rebmix->Criterion_ = icSSE;
    else {
        *Error = 1; goto E0;
    }

    rebmix->length_pdf_ = *d;

    rebmix->Variables_ = (VariablesType_e*)malloc(rebmix->length_pdf_ * sizeof(VariablesType_e));

    *Error = NULL == rebmix->Variables_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        if (!strcmp(Variables[i], "continuous")) {
            rebmix->Variables_[i] = vtContinuous;
        }
        else
        if (!strcmp(Variables[i], "discrete")) {
            rebmix->Variables_[i] = vtDiscrete;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    *Error = NULL == rebmix->IniTheta_; if (*Error) goto E0;

    rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (int*)malloc(rebmix->length_Theta_ * sizeof(int));

    *Error = NULL == rebmix->length_theta_; if (*Error) goto E0;

    *Error = rebmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmix->IniTheta_->pdf_[i] = pfNormal;
        }
        else
        if (!strcmp(pdf[i], "lognormal")) {
            rebmix->IniTheta_->pdf_[i] = pfLognormal;
        }
        else
        if (!strcmp(pdf[i], "Weibull")) {
            rebmix->IniTheta_->pdf_[i] = pfWeibull;
        }
        else
        if (!strcmp(pdf[i], "gamma")) {
            rebmix->IniTheta_->pdf_[i] = pfGamma;
        }
        else
        if (!strcmp(pdf[i], "Gumbel")) {
            rebmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rebmix->IniTheta_->pdf_[i] = pfvonMises;
        }
        else
        if (!strcmp(pdf[i], "binomial")) {
            rebmix->IniTheta_->pdf_[i] = pfBinomial;
        }
        else
        if (!strcmp(pdf[i], "Poisson")) {
            rebmix->IniTheta_->pdf_[i] = pfPoisson;
        }
        else
        if (!strcmp(pdf[i], "Dirac")) {
            rebmix->IniTheta_->pdf_[i] = pfDirac;
        }
        else
        if (!strcmp(pdf[i], "uniform")) {
            rebmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) if (rebmix->IniTheta_->Theta_[j]) {
        for (k = 0; k < rebmix->length_theta_[j]; k++) {
            rebmix->IniTheta_->Theta_[j][k] = Theta[i];

            i++;
        }
    }

    rebmix->length_K_ = *length_K;

    rebmix->K_ = (int*)malloc(rebmix->length_K_ * rebmix->length_pdf_ * sizeof(int));

    *Error = NULL == rebmix->K_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_K_ * rebmix->length_pdf_; i++) {
        rebmix->K_[i] = K[i];
    }

    if (*length_y0 > 0) {
        rebmix->y0_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

        *Error = NULL == rebmix->y0_; if (*Error) goto E0;

        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->y0_[i] = y0[i];
        }
    }
    else {
        rebmix->y0_ = NULL;
    }

    if (*length_ymin > 0) {
        rebmix->ymin_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

        *Error = NULL == rebmix->ymin_; if (*Error) goto E0;

        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymin_[i] = ymin[i];
        }
    }
    else {
        rebmix->ymin_ = NULL;
    }

    if (*length_ymax > 0) {
        rebmix->ymax_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

        *Error = NULL == rebmix->ymax_; if (*Error) goto E0;

        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymax_[i] = ymax[i];
        }
    }
    else {
        rebmix->ymax_ = NULL;
    }

    rebmix->ar_ = *ar;

    if (!strcmp(Restraints[0], "rigid")) {
        rebmix->Restraints_ = rtRigid;
    }
    else
    if (!strcmp(Restraints[0], "loose")) {
        rebmix->Restraints_ = rtLoose;
    }
    else {
        *Error = 1; goto E0;
    }

/// Panic Branislav.

    if (!strcmp(EMStrategy[0], "exhaustive")) {
        rebmix->EM_strategy_ = strategy_exhaustive;
    }
    else
    if (!strcmp(EMStrategy[0], "best")) {
        rebmix->EM_strategy_ = strategy_best;
    }
    else
    if (!strcmp(EMStrategy[0], "single")) {
        rebmix->EM_strategy_ = strategy_single;
    }
    else{
        rebmix->EM_strategy_ = strategy_none;
    }

    if (!strcmp(EMVariant[0], "EM")) {
        rebmix->EM_variant_ = varEM;
    }
    else
    if (!strcmp(EMVariant[0], "ECM")) {
        rebmix->EM_variant_ = varECM;
    }
    else {
        if (rebmix->EM_strategy_ != strategy_none) {
            *Error = 1; goto E0;
        }

        rebmix->EM_variant_ = varEM;
    }

    if (!strcmp(EMAcceleration[0], "fixed")) {
        rebmix->EM_accel_ = acc_fixed;
    }
    else
    if (!strcmp(EMAcceleration[0], "line")) {
        rebmix->EM_accel_ = acc_line;
    }
    else
    if (!strcmp(EMAcceleration[0], "golden")) {
        rebmix->EM_accel_ = acc_golden;
    }
    else {
        if (rebmix->EM_strategy_ != strategy_none) {
            *Error = 1; goto E0;
        }

        rebmix->EM_accel_ = acc_fixed;
    }

    rebmix->EM_TOL_ = *EMTolerance;

    rebmix->EM_am_ = *EMAccelerationMul;

    rebmix->EM_max_iter_ = *EMMaxIter;

	rebmix->EM_K_ = *EMK;

/// End

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    rebmix->X_ = (FLOAT**)malloc(rebmix->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->X_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->X_[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == rebmix->X_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[j][l] = Y[i]; i++;
        }
    }

    *Error = rebmix->REBMIX();

    if (*Error) goto E0;

/// Panic Branislav.

    *n_iter = rebmix->n_iter_;

    *n_iter_sum = rebmix->n_iter_sum_;

/// End

    *summary_k = rebmix->summary_.k;

    if (rebmix->summary_.h) for (i = 0; i < rebmix->length_pdf_; i++) {
        summary_h[i] = rebmix->summary_.h[i];
    }

    if (rebmix->summary_.y0) for (i = 0; i < rebmix->length_pdf_; i++) {
        summary_y0[i] = rebmix->summary_.y0[i];
    }

    if (rebmix->summary_.ymin) for (i = 0; i < rebmix->length_pdf_; i++) {
        summary_ymin[i] = rebmix->summary_.ymin[i];
    }

    if (rebmix->summary_.ymax) for (i = 0; i < rebmix->length_pdf_; i++) {
        summary_ymax[i] = rebmix->summary_.ymax[i];
    }

    *summary_IC = rebmix->summary_.IC;
    *summary_logL = rebmix->summary_.logL;
    *summary_M = rebmix->summary_.M;
    *summary_c = rebmix->summary_.c;

    for (j = 0; j < rebmix->summary_.c; j++) {
        W[j] = rebmix->W_[j];
    }

    i = 0;

    for (j = 0; j < rebmix->summary_.c; j++) {
        for (l = 0; l < rebmix->length_theta_[0]; l++) {
            theta1[i] = rebmix->MixTheta_[j]->Theta_[0][l];

            i++;
        }
    }

    i = 0;

    for (j = 0; j < rebmix->summary_.c; j++) {
        for (l = 0; l < rebmix->length_theta_[1]; l++) {
            theta2[i] = rebmix->MixTheta_[j]->Theta_[1][l];

            i++;
        }
    }

    *opt_length = rebmix->opt_length_;

    for (i = 0; i < rebmix->opt_length_; i++) {
        opt_c[i] = rebmix->opt_c_[i];
        opt_IC[i] = rebmix->opt_IC_[i];
        opt_logL[i] = rebmix->opt_logL_[i];
        opt_D[i] = rebmix->opt_D_[i];
    }

    i = 0;

    for (j = 0; j < rebmix->all_length_; j++) if (rebmix->all_K_[j]) {
        all_K[i] = rebmix->all_K_[j]; all_IC[i] = rebmix->all_IC_[j];

        i++;
    }

    *all_length = i;

E0: if (rebmix) delete rebmix;
} // RREBMIX

// Returns k-nearest neighbour empirical densities in R.

void RdensKNearestNeighbourXY(int    *n,     // Total number of independent observations.
                              double *x,     // Pointer to the input array x.
                              double *y,     // Pointer to the input array y.
                              double *p,     // Pointer to the output array p.
                              int    *k,     // k-nearest neighbours.
                              double *hx,    // Normalizing vector.
                              double *hy,    // Normalizing vector.
                              int    *Error) // Error code.
{
    FLOAT *Dk = NULL;
    FLOAT Dc, R, C;
    int   i, j, K, l, m, q;

    *Error = *n < 1; if (*Error) return;

    K = *k; if (K > 1) K -= 1; else K = 1;

    Dk = (FLOAT*)malloc(K * sizeof(FLOAT));

    *Error = NULL == Dk; if (*Error) goto E0;

    C = (*k) / ((*n) * Pi * (*hx) * (*hy));

    for (i = 0; i < *n; i++) {
        Dk[0] = FLOAT_MAX; q = 0;

        for (j = 0; j < *n; j++) if (i != j) {
            R = (x[i] - x[j]) / (*hx); Dc = R * R;
            R = (y[i] - y[j]) / (*hy); Dc += R * R;

            q += Dc <= FLOAT_MIN;

            for (l = 0; l < K; l++) {
                if (Dc < Dk[l]) {
                    for (m = K - 1; m > l; m--) Dk[m] = Dk[m - 1];

                    if ((Dc > FLOAT_MIN) || (l != K - 1)) Dk[l] = Dc;

                    break;
                }
            }
        }

        R = (FLOAT)sqrt(Dk[K - 1]);

        if (q >= K) R *= (FLOAT)sqrt((K + (FLOAT)1.0) / (q + (FLOAT)2.0));

        p[i] = C / (R * R);
    }

E0: if (Dk) free(Dk);
} // RdensKNearestNeighbourXY

// Returns kernel density estimation empirical densities in R.

void RdensKDEXY(int    *n,     // Total number of independent observations.
                double *x,     // Pointer to the input array x.
                double *y,     // Pointer to the input array y.
                double *p,     // Pointer to the output array p.
                double *hx,    // Side of the hypersquare.
                double *hy,    // Side of the hypersquare.
                int    *Error) // Error code.
{
    int   i, j;
    FLOAT C, rx, ry;

    *Error = *n < 1; if (*Error) return;

    C = (FLOAT)1.0 / (*hx) / (*hy) / (*n); rx = (FLOAT)0.5 * (*hx); ry = (FLOAT)0.5 * (*hy);

    for (i = 0; i < *n; i++) {
        p[i] = (FLOAT)0.0;
    }

    for (i = 0; i < *n; i++) {
        for (j = i; j < *n; j++) {
            if (((FLOAT)fabs(x[j] - x[i]) <= rx) && ((FLOAT)fabs(y[j] - y[i]) <= ry)) {
                p[i] += C; if (i != j) p[j] += C;
            }
        }
    }
} // RdensKDEXY

// Returns histogram empirical densities in R.

void RdensHistogramXY(int    *k,     // Total number of bins.
                      int    *n,     // Total number of independent observations.
                      double *x,     // Pointer to the input array x.
                      double *y,     // Pointer to the input array y.
                      double *p,     // Pointer to the output array p.
                      double *x0,    // Origin.
                      double *xmin,  // Minimum observation.
                      double *xmax,  // Maximim observation.
                      double *y0,    // Origin.
                      double *ymin,  // Minimum observation.
                      double *ymax,  // Maximim observation.
                      double *hx,    // Side of the hypersquare.
                      double *hy,    // Side of the hypersquare.
                      char   **px,   // Parametric family type.
                      char   **py,   // Parametric family type.
                      int    *Error) // Error code.
{
    int                    i, j;
    ParametricFamilyType_e pdfx, pdfy;
    FLOAT                  C, rx, ry;

    *Error = *n < 1; if (*Error) return;

    if (!strcmp(px[0], "normal")) {
        pdfx = pfNormal;
    }
    else
    if (!strcmp(px[0], "lognormal")) {
        pdfx = pfLognormal;
    }
    else
    if (!strcmp(px[0], "Weibull")) {
        pdfx = pfWeibull;
    }
    else
    if (!strcmp(px[0], "gamma")) {
        pdfx = pfGamma;
    }
    else
    if (!strcmp(px[0], "Gumbel")) {
        pdfx = pfGumbel;
    }
    else
    if (!strcmp(px[0], "vonMises")) {
        pdfx = pfvonMises;
    }
    else
    if (!strcmp(px[0], "binomial")) {
        pdfx = pfBinomial;
    }
    else
    if (!strcmp(px[0], "Poisson")) {
        pdfx = pfPoisson;
    }
    else
    if (!strcmp(px[0], "Dirac")) {
        pdfx = pfDirac;
    }
    else
    if (!strcmp(px[0], "uniform")) {
        pdfx = pfUniform;
    }
    else {
      *Error = 1; return;
    }

    if (!strcmp(py[0], "normal")) {
        pdfy = pfNormal;
    }
    else
    if (!strcmp(py[0], "lognormal")) {
        pdfy = pfLognormal;
    }
    else
    if (!strcmp(py[0], "Weibull")) {
        pdfy = pfWeibull;
    }
    else
    if (!strcmp(py[0], "gamma")) {
        pdfy = pfGamma;
    }
    else
    if (!strcmp(py[0], "Gumbel")) {
        pdfy = pfGumbel;
    }
    else
    if (!strcmp(py[0], "vonMises")) {
        pdfy = pfvonMises;
    }
    else
    if (!strcmp(py[0], "binomial")) {
        pdfy = pfBinomial;
    }
    else
    if (!strcmp(py[0], "Poisson")) {
        pdfy = pfPoisson;
    }
    else
    if (!strcmp(py[0], "Dirac")) {
        pdfy = pfDirac;
    }
    else
    if (!strcmp(py[0], "uniform")) {
        pdfy = pfUniform;
    }
    else {
      *Error = 1; return;
    }

    C = (FLOAT)1.0 / (*hx) / (*hy) / (*n); rx = (FLOAT)0.5 * (*hx); ry = (FLOAT)0.5 * (*hy);

    *k = 0;

    for (i = 0; i < *n; i++) {
        j = (int)floor((x[i] - (*x0)) / (*hx) + (FLOAT)0.5);

        x[*k] = (*x0) + j * (*hx);

        if (x[*k] < *xmin) {
            x[*k] += (*hx);
        }
        else
        if (x[*k] > *xmax) {
            x[*k] -= (*hx);
        }

        switch (pdfx) {
        case pfNormal: case pfTNormal: case pfvonMises: case pfBinomial: case pfPoisson: case pfDirac: case pfUniform: case pfGumbel: default:
            break;
        case pfLognormal: case pfWeibull: case pfGamma:
            if (x[*k] <= FLOAT_MIN) x[*k] += (*hx);
        }

        j = (int)floor((y[i] - (*y0)) / (*hy) + (FLOAT)0.5);

        y[*k] = (*y0) + j * (*hy);

        if (y[*k] < *ymin) {
            y[*k] += (*hy);
        }
        else
        if (y[*k] > *ymax) {
            y[*k] -= (*hy);
        }

        switch (pdfy) {
        case pfNormal: case pfTNormal: case pfvonMises: case pfBinomial: case pfPoisson: case pfDirac: case pfUniform: case pfGumbel: default:
            break;
        case pfLognormal: case pfWeibull: case pfGamma:
            if (y[*k] <= FLOAT_MIN) y[*k] += (*hy);
        }

        for (j = 0; j < *k; j++) {
            if (((FLOAT)fabs(x[j] - x[*k]) > rx) || ((FLOAT)fabs(y[j] - y[*k]) > ry)) goto S0;

            p[j] += C; goto S1;
S0:;    }

        p[*k] = C; (*k)++;
S1:;}
} // RdensHistogramXY

// Returns k-nearest neighbour empirical densities in R.

void RdensKNearestNeighbourX(int    *n,     // Total number of independent observations.
                             double *x,     // Pointer to the input array x.
                             double *p,     // Pointer to the output array p.
                             int    *k,     // k-nearest neighbours.
                             double *hx,    // Normalizing vector.
                             int    *Error) // Error code.
{
    FLOAT *Dk = NULL;
    FLOAT Dc, R, C;
    int   i, j, K, l, m, q;

    *Error = *n < 1; if (*Error) return;

    K = *k; if (K > 1) K -= 1; else K = 1;

    Dk = (FLOAT*)malloc(K * sizeof(FLOAT));

    *Error = NULL == Dk; if (*Error) goto E0;

    C = (*k) / ((*n) * (FLOAT)2.0 * (*hx));

    for (i = 0; i < *n; i++) {
        Dk[0] = FLOAT_MAX; q = 0;

        for (j = 0; j < *n; j++) if (i != j) {
            R = (FLOAT)fabs((x[i] - x[j]) / (*hx)); Dc = R;

            q += Dc <= FLOAT_MIN;

            for (l = 0; l < K; l++) {
                if (Dc < Dk[l]) {
                    for (m = K - 1; m > l; m--) Dk[m] = Dk[m - 1];

                    if ((Dc > FLOAT_MIN) || (l != K - 1)) Dk[l] = Dc;

                    break;
                }
            }
        }

        R = Dk[K - 1];

        if (q >= K) R *= (FLOAT)((K + (FLOAT)1.0) / (q + (FLOAT)2.0));

        p[i] = C / R;
    }

E0: if (Dk) free(Dk);
} // RdensKNearestNeighbourX

// Returns kernel density estimation empirical densities in R.

void RdensKDEX(int    *n,     // Total number of independent observations.
               double *x,     // Pointer to the input array x.
               double *p,     // Pointer to the output array p.
               double *hx,    // Side of the hypersquare.
               int    *Error) // Error code.
{
    int   i, j;
    FLOAT C, rx;

    *Error = *n < 1; if (*Error) return;

    C = (FLOAT)1.0 / (*hx) / (*n); rx = (FLOAT)0.5 * (*hx);

    for (i = 0; i < *n; i++) {
        p[i] = (FLOAT)0.0;
    }

    for (i = 0; i < *n; i++) {
        for (j = i; j < *n; j++) {
            if ((FLOAT)fabs(x[j] - x[i]) <= rx) {
                p[i] += C; if (i != j) p[j] += C;
            }
        }
    }
} // RdensKDEX

// Returns histogram empirical densities in R.

void RdensHistogramX(int    *k,     // Total number of bins.
                     int    *n,     // Total number of independent observations.
                     double *x,     // Pointer to the input array x.
                     double *p,     // Pointer to the output array p.
                     double *x0,    // Origin.
                     double *xmin,  // Minimum observation.
                     double *xmax,  // Maximum observation.
                     double *hx,    // Side of the hypersquare.
                     char   **px,   // Parametric family type.
                     int    *Error) // Error code.
{
    int                    i, j;
    ParametricFamilyType_e pdfx;
    FLOAT                  C, rx;

    *Error = *n < 1; if (*Error) return;

    if (!strcmp(px[0], "normal")) {
        pdfx = pfNormal;
    }
    else
    if (!strcmp(px[0], "lognormal")) {
        pdfx = pfLognormal;
    }
    else
    if (!strcmp(px[0], "Weibull")) {
        pdfx = pfWeibull;
    }
    else
    if (!strcmp(px[0], "gamma")) {
        pdfx = pfGamma;
    }
    else
    if (!strcmp(px[0], "Gumbel")) {
        pdfx = pfGumbel;
    }
    else
    if (!strcmp(px[0], "vonMises")) {
        pdfx = pfvonMises;
    }
    else
    if (!strcmp(px[0], "binomial")) {
        pdfx = pfBinomial;
    }
    else
    if (!strcmp(px[0], "Poisson")) {
        pdfx = pfPoisson;
    }
    else
    if (!strcmp(px[0], "Dirac")) {
        pdfx = pfDirac;
    }
    else
    if (!strcmp(px[0], "uniform")) {
        pdfx = pfUniform;
    }
    else {
      *Error = 1; return;
    }

    C = (FLOAT)1.0 / (*hx) / (*n); rx = (FLOAT)0.5 * (*hx);

    *k = 0;

    for (i = 0; i < *n; i++) {
        j = (int)floor((x[i] - (*x0)) / (*hx) + (FLOAT)0.5);

        x[*k] = (*x0) + j * (*hx);

        if (x[*k] < *xmin) {
            x[*k] += (*hx);
        }
        else
        if (x[*k] > *xmax) {
            x[*k] -= (*hx);
        }

        switch (pdfx) {
        case pfNormal: case pfTNormal: case pfvonMises: case pfBinomial: case pfPoisson: case pfDirac: case pfUniform: case pfGumbel: default:
            break;
        case pfLognormal: case pfWeibull: case pfGamma:
            if (x[*k] <= FLOAT_MIN) x[*k] += (*hx);
        }

        for (j = 0; j < *k; j++) {
            if ((FLOAT)fabs(x[j] - x[*k]) > rx) goto S0;

            p[j] += C; goto S1;
S0:;    }

        p[*k] = C; (*k)++;
S1:;}
} // RdensHistogramX

// Returns classified observations in R.

void RCLSMIX(int    *n,      // Total number of independent observations.
             double *X,      // Pointer to the input array X.
             int    *s,      // Number of classes.
             int    *o,      // Number of input REBMIX objects.
             int    *d,      // Number of independent random variables in REBMIX objects.
             int    *c,      // Number of components.
             double *W,      // Component weights.
             char   **pdf,   // Component parameters.
             double *theta1, // Component parameters.
             double *theta2, // Component parameters.
             double *P,      // Prior probabilities.
             int    *Z,      // Pointer to the output array Z.
             int    *Error)  // Error code.
{
    Rebmix               *rebmix = NULL;
    int                  **C = NULL;
    int                  A[2];
    FLOAT                ***Q = NULL;
    FLOAT                **Y = NULL;
    CompnentDistribution ****Theta = NULL;
    FLOAT                CmpDist, MixDist, MaxMixDist;
    int                  i, j, k, l, m, dmax = 0;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    C = (int**)malloc(*s * sizeof(int*));

    *Error = NULL == C; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < *s; j++) {
        C[j] = (int*)malloc(*o * sizeof(int));

        *Error = NULL == C[j]; if (*Error) goto E0;

        for (k = 0; k < *o; k++) {
            C[j][k] = c[i]; i++;
        }
    }

    Q = (FLOAT***)malloc(*s * sizeof(FLOAT**));

    *Error = NULL == Q; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < *s; j++) {
        Q[j] = (FLOAT**)malloc(*o * sizeof(FLOAT*));

        *Error = NULL == Q[j]; if (*Error) goto E0;

        for (k = 0; k < *o; k++) {
            Q[j][k] = (FLOAT*)malloc(C[j][k] * sizeof(FLOAT));

            *Error = NULL == Q[j][k]; if (*Error) goto E0;

            for (l = 0; l < C[j][k]; l++) {
                Q[j][k][l] = W[i]; i++;
            }
        }
    }

    Theta = new CompnentDistribution*** [(unsigned int)(*s)];

    *Error = NULL == Theta; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < *s; j++) {
        Theta[j] = new CompnentDistribution** [(unsigned int)(*o)];

        *Error = NULL == Theta[j]; if (*Error) goto E0;

        for (k = 0; k < *o; k++) {
            Theta[j][k] = new CompnentDistribution* [(unsigned int)C[j][k]];

            *Error = NULL == Theta[j][k]; if (*Error) goto E0;

            for (l = 0; l < C[j][k]; l++) {
                Theta[j][k][l] = new CompnentDistribution(rebmix);

                *Error = NULL == Theta[j][k][l]; if (*Error) goto E0;

                A[0] = A[1] = d[k];

                *Error = Theta[j][k][l]->Realloc(d[k], 2, A);

                if (*Error) goto E0;

                for (m = 0; m < d[k]; m++) {
                    if (!strcmp(pdf[i], "normal")) {
                        Theta[j][k][l]->pdf_[m] = pfNormal;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "lognormal")) {
                        Theta[j][k][l]->pdf_[m] = pfLognormal;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "Weibull")) {
                        Theta[j][k][l]->pdf_[m] = pfWeibull;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "gamma")) {
                        Theta[j][k][l]->pdf_[m] = pfGamma;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "Gumbel")) {
                        Theta[j][k][l]->pdf_[m] = pfGumbel;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "vonMises")) {
                        Theta[j][k][l]->pdf_[m] = pfvonMises;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "binomial")) {
                        Theta[j][k][l]->pdf_[m] = pfBinomial;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "Poisson")) {
                        Theta[j][k][l]->pdf_[m] = pfPoisson;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                    }
                    else
                    if (!strcmp(pdf[i], "Dirac")) {
                        Theta[j][k][l]->pdf_[m] = pfDirac;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "uniform")) {
                        Theta[j][k][l]->pdf_[m] = pfUniform;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else {
                        *Error = 1; goto E0;
                    }

                    i++;
                }
            }
        }
    }

    dmax = d[0]; for (i = 1; i < *o; i++) if (d[i] > dmax) dmax = d[i];

    Y = (FLOAT**)malloc(dmax * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < dmax; i++) {
        Y[i] = (FLOAT*)malloc(sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    for (i = 0; i < *n; i++) {
        Z[i] = 1; MaxMixDist = (FLOAT)0.0;

        for (j = 0; j < *s; j++) {
            k = 0; MixDist = (FLOAT)1.0;

            for (l = 0; l < *o; l++) {
                for (m = 0; m < d[l]; m++) {
                    Y[m][0] = X[i + (*n) * (m + k)];
                }

                *Error = rebmix->MixtureDist(0, Y, C[j][l], Q[j][l], Theta[j][l], &CmpDist);

                if (*Error) goto E0;

                k += d[l]; MixDist *= CmpDist;
            }

            MixDist *= P[j];

            if (MixDist > MaxMixDist) {
                Z[i] = j + 1; MaxMixDist = MixDist;
            }
        }
    }

E0: if (Y) {
        for (i = 0; i < dmax; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (Theta) {
        for (i = 0; i < *s; i++) {
            if (Theta[i]) {
                for (j = 0; j < *o; j++) {
                    if (Theta[i][j]) {
                        for (k = 0; k < C[i][j]; k++) {
                            if (Theta[i][j][k]) delete Theta[i][j][k];
                        }

                        delete[] Theta[i][j];
                    }
                }

                delete[] Theta[i];
            }
        }

        delete[] Theta;
    }

    if (Q) {
        for (i = 0; i < *s; i++) {
            if (Q[i]) {
                for (j = 0; j < *o; j++) if (Q[i][j]) free(Q[i][j]);

                free(Q[i]);
            }
        }

        free(Q);
    }

    if (C) {
        for (i = 0; i < *s; i++) {
            if (C[i]) free(C[i]);
        }

        free(C);
    }

    if (rebmix) delete rebmix;
} // RCLSMIX

// Returns clustered observations in R.

void RCLRMIX(int    *n,      // Total number of independent observations.
             double *X,      // Pointer to the input array X.
             int    *d,      // Number of independent random variables.
             int    *c,      // Number of components.
             double *W,      // Component weights.
             char   **pdf,   // Component parameters.
             double *theta1, // Component parameters.
             double *theta2, // Component parameters.
             int    *Z,      // Pointer to the output array Z.
             int    *Error)  // Error code.
{
    Rebmix               *rebmix = NULL;
    FLOAT                **Y = NULL;
    int                  A[2];
    CompnentDistribution **Theta = NULL;
    FLOAT                CmpDist, MaxCmpDist;
    int                  i, j, k;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->length_pdf_ = *d;

    Theta = new CompnentDistribution* [(unsigned int)(*c)];

    *Error = NULL == Theta; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < *c; j++) {
        Theta[j] = new CompnentDistribution(rebmix);

        *Error = NULL == Theta[j]; if (*Error) goto E0;

        A[0] = A[1] = *d;

        *Error = Theta[j]->Realloc(*d, 2, A);

        if (*Error) goto E0;

        for (k = 0; k < *d; k++) {
            if (!strcmp(pdf[i], "normal")) {
                Theta[j]->pdf_[k] = pfNormal;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "lognormal")) {
                Theta[j]->pdf_[k] = pfLognormal;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "Weibull")) {
                Theta[j]->pdf_[k] = pfWeibull;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "gamma")) {
                Theta[j]->pdf_[k] = pfGamma;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "Gumbel")) {
                Theta[j]->pdf_[k] = pfGumbel;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "vonMises")) {
                Theta[j]->pdf_[k] = pfvonMises;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "binomial")) {
                Theta[j]->pdf_[k] = pfBinomial;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "Poisson")) {
                Theta[j]->pdf_[k] = pfPoisson;

                Theta[j]->Theta_[0][k] = theta1[i];
            }
            else
            if (!strcmp(pdf[i], "Dirac")) {
                Theta[j]->pdf_[k] = pfDirac;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "uniform")) {
                Theta[j]->pdf_[k] = pfUniform;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else {
                *Error = 1; goto E0;
            }

            i++;
        }
    }

    Y = (FLOAT**)malloc(*d * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < *d; i++) {
        Y[i] = (FLOAT*)malloc(sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    for (i = 0; i < *n; i++) {
        for (j = 0; j < *d; j++) {
            Y[j][0] = X[i + (*n) * j];
        }

        Z[i] = 1; MaxCmpDist = (FLOAT)0.0;

        for (j = 0; j < *c; j++) {
            *Error = rebmix->ComponentDist(0, Y, Theta[j], &CmpDist, NULL);

            if (*Error) goto E0;

            CmpDist *= W[j];

            if (CmpDist > MaxCmpDist) {
                Z[i] = j + 1; MaxCmpDist = CmpDist;
            }
        }
    }

E0: if (Y) {
        for (i = 0; i < *d; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (Theta) {
        for (i = 0; i < *c; i++) {
            if (Theta[i]) delete Theta[i];
        }

        delete[] Theta;
    }

    if (rebmix) delete rebmix;
} // RCLRMIX

void RPreprocessingKNNMIX(int    *k,     // k-nearest neighbours.
                          double *h,     // Normalizing vector.
                          int    *n,     // Total number of independent observations.
                          int    *d,     // Number of independent random variables.
                          double *x,     // Pointer to the input array x.
                          double *y,     // Pointer to the output array y.
                          int    *Error) // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    int    i, j, l;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->n_ = *n;
    rebmix->length_pdf_ = *d;

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 3) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_ + 3; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            Y[j][l] = x[i]; i++;
        }
    }

    *Error = rebmix->PreprocessingKNN(*k, h, Y);

    if (*Error) goto E0;

    i = 0;

    for (j = 0; j < rebmix->length_pdf_ + 3; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

E0: if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 3; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;
} // RPreprocessingKNNMIX

void RPreprocessingKDEMIX(double *h,     // Sides of the hypersquare.
                          int    *n,     // Total number of independent observations.
                          int    *d,     // Number of independent random variables.
                          double *x,     // Pointer to the input array x.
                          double *y,     // Pointer to the output array y.
                          int    *Error) // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    int    i, j, l;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->n_ = *n;
    rebmix->length_pdf_ = *d;

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 2) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_ + 2; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            Y[j][l] = x[i]; i++;
        }
    }

    *Error = rebmix->PreprocessingKDE(h, Y);

    if (*Error) goto E0;

    i = 0;

    for (j = 0; j < rebmix->length_pdf_ + 2; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

E0: if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 2; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;
} // RPreprocessingKDEMIX

void RPreprocessingHMIX(double *h,          // Sides of the hypersquare.
                        double *y0,         // Origins.
                        double *ymin,       // Minimum observations.
                        double *ymax,       // Maximum observations.
                        int    *k,          // Total number of bins.
                        int    *n,          // Total number of independent observations.
                        int    *d,          // Number of independent random variables.
                        double *x,          // Pointer to the input array x.
                        double *y,          // Pointer to the output array y.
                        int    *Error)      // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    int    i, j, l;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->n_ = *n;
    rebmix->length_pdf_ = *d;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 1) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    *Error = rebmix->PreprocessingH(h, y0, ymin, ymax, k, Y);

    if (*Error) goto E0;

    i = 0;

    for (j = 0; j < rebmix->length_pdf_ + 1; j++) {
        for (l = 0; l < *k; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

E0: if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;
} // RPreprocessingHMIX

void RInformationCriterionKNNMIX(double *h,            // Sides of the hypersquare.
                                 int    *k,            // k-nearest neighbours.
                                 char   **Criterion,   // Information criterion type.
                                 int    *c,            // Number of components.
                                 double *W,            // Component weights.
                                 int    *length_pdf,   // Length of pdf.
                                 int    *length_Theta, // Length of Theta.
                                 int    *length_theta, // Length of Theta[i].
                                 char   **pdf,         // Parametric family types.
                                 double *Theta,        // Component parameters.
                                 int    *n,            // Number of observations.
                                 double *x,            // Dataset.
                                 double *IC,           // Information criterion.
                                 double *logL,         // log-likelihood.
                                 int    *M,            // Degrees of freedom.
                                 double *D,            // Total of positive relative deviations.
                                 int    *Error)        // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    int    i, j, l, m;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    if (!strcmp(Criterion[0], "AIC"))
        rebmix->Criterion_ = icAIC;
    else
    if (!strcmp(Criterion[0], "AIC3"))
        rebmix->Criterion_ = icAIC3;
    else
    if (!strcmp(Criterion[0], "AIC4"))
        rebmix->Criterion_ = icAIC4;
    else
    if (!strcmp(Criterion[0], "AICc"))
        rebmix->Criterion_ = icAICc;
    else
    if (!strcmp(Criterion[0], "BIC"))
        rebmix->Criterion_ = icBIC;
    else
    if (!strcmp(Criterion[0], "CAIC"))
        rebmix->Criterion_ = icCAIC;
    else
    if (!strcmp(Criterion[0], "HQC"))
        rebmix->Criterion_ = icHQC;
    else
    if (!strcmp(Criterion[0], "MDL2"))
        rebmix->Criterion_ = icMDL2;
    else
    if (!strcmp(Criterion[0], "MDL5"))
        rebmix->Criterion_ = icMDL5;
    else
    if (!strcmp(Criterion[0], "AWE"))
        rebmix->Criterion_ = icAWE;
    else
    if (!strcmp(Criterion[0], "CLC"))
        rebmix->Criterion_ = icCLC;
    else
    if (!strcmp(Criterion[0], "ICL"))
        rebmix->Criterion_ = icICL;
    else
    if (!strcmp(Criterion[0], "PC"))
        rebmix->Criterion_ = icPC;
    else
    if (!strcmp(Criterion[0], "ICL-BIC"))
        rebmix->Criterion_ = icICLBIC;
    else
    if (!strcmp(Criterion[0], "D"))
        rebmix->Criterion_ = icD;
    else
    if (!strcmp(Criterion[0], "SSE"))
        rebmix->Criterion_ = icSSE;
    else {
        *Error = 1; goto E0;
    }

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == rebmix->W_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    *Error = NULL == rebmix->IniTheta_; if (*Error) goto E0;

    rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (int*)malloc(rebmix->length_Theta_ * sizeof(int));

    *Error = NULL == rebmix->length_theta_; if (*Error) goto E0;

    *Error = rebmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmix->IniTheta_->pdf_[i] = pfNormal;
        }
        else
        if (!strcmp(pdf[i], "lognormal")) {
            rebmix->IniTheta_->pdf_[i] = pfLognormal;
        }
        else
        if (!strcmp(pdf[i], "Weibull")) {
            rebmix->IniTheta_->pdf_[i] = pfWeibull;
        }
        else
        if (!strcmp(pdf[i], "gamma")) {
            rebmix->IniTheta_->pdf_[i] = pfGamma;
        }
        else
        if (!strcmp(pdf[i], "Gumbel")) {
            rebmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rebmix->IniTheta_->pdf_[i] = pfvonMises;
        }
        else
        if (!strcmp(pdf[i], "binomial")) {
            rebmix->IniTheta_->pdf_[i] = pfBinomial;
        }
        else
        if (!strcmp(pdf[i], "Poisson")) {
            rebmix->IniTheta_->pdf_[i] = pfPoisson;
        }
        else
        if (!strcmp(pdf[i], "Dirac")) {
            rebmix->IniTheta_->pdf_[i] = pfDirac;
        }
        else
        if (!strcmp(pdf[i], "uniform")) {
            rebmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    rebmix->MixTheta_ = new CompnentDistribution*[(unsigned int)(*c)];

    *Error = NULL == rebmix->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        *Error = NULL == rebmix->MixTheta_[i]; if (*Error) goto E0;

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) if (rebmix->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 3) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_ + 3; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;

        if (i < rebmix->length_pdf_) {
            for (j = 0; j < rebmix->n_; j++) Y[i][j] = rebmix->Y_[i][j];
        }
    }

    *Error = rebmix->PreprocessingKNN(*k, h, Y);

    if (*Error) goto E0;

    rebmix->cmax_ = *c;

    *Error = rebmix->InformationCriterionKNN(*k,
                                             Y,
                                             *c,
                                             rebmix->W_,
                                             rebmix->MixTheta_,
                                             IC,
                                             logL,
                                             M,
                                             D);

    if (*Error) goto E0;

E0: if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 3; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;
} // RInformationCriterionKNNMIX

void RInformationCriterionKDEMIX(double *h,            // Sides of the hypersquare.
                                 char   **Criterion,   // Information criterion type.
                                 int    *c,            // Number of components.
                                 double *W,            // Component weights.
                                 int    *length_pdf,   // Length of pdf.
                                 int    *length_Theta, // Length of Theta.
                                 int    *length_theta, // Length of Theta[i].
                                 char   **pdf,         // Parametric family types.
                                 double *Theta,        // Component parameters.
                                 int    *n,            // Number of observations.
                                 double *x,            // Dataset.
                                 double *IC,           // Information criterion.
                                 double *logL,         // log-likelihood.
                                 int    *M,            // Degrees of freedom.
                                 double *D,            // Total of positive relative deviations.
                                 int    *Error)        // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    FLOAT  logV;
    int    i, j, l, m;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    if (!strcmp(Criterion[0], "AIC"))
        rebmix->Criterion_ = icAIC;
    else
    if (!strcmp(Criterion[0], "AIC3"))
        rebmix->Criterion_ = icAIC3;
    else
    if (!strcmp(Criterion[0], "AIC4"))
        rebmix->Criterion_ = icAIC4;
    else
    if (!strcmp(Criterion[0], "AICc"))
        rebmix->Criterion_ = icAICc;
    else
    if (!strcmp(Criterion[0], "BIC"))
        rebmix->Criterion_ = icBIC;
    else
    if (!strcmp(Criterion[0], "CAIC"))
        rebmix->Criterion_ = icCAIC;
    else
    if (!strcmp(Criterion[0], "HQC"))
        rebmix->Criterion_ = icHQC;
    else
    if (!strcmp(Criterion[0], "MDL2"))
        rebmix->Criterion_ = icMDL2;
    else
    if (!strcmp(Criterion[0], "MDL5"))
        rebmix->Criterion_ = icMDL5;
    else
    if (!strcmp(Criterion[0], "AWE"))
        rebmix->Criterion_ = icAWE;
    else
    if (!strcmp(Criterion[0], "CLC"))
        rebmix->Criterion_ = icCLC;
    else
    if (!strcmp(Criterion[0], "ICL"))
        rebmix->Criterion_ = icICL;
    else
    if (!strcmp(Criterion[0], "PC"))
        rebmix->Criterion_ = icPC;
    else
    if (!strcmp(Criterion[0], "ICL-BIC"))
        rebmix->Criterion_ = icICLBIC;
    else
    if (!strcmp(Criterion[0], "D"))
        rebmix->Criterion_ = icD;
    else
    if (!strcmp(Criterion[0], "SSE"))
        rebmix->Criterion_ = icSSE;
    else {
        *Error = 1; goto E0;
    }

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == rebmix->W_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    *Error = NULL == rebmix->IniTheta_; if (*Error) goto E0;

    rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (int*)malloc(rebmix->length_Theta_ * sizeof(int));

    *Error = NULL == rebmix->length_theta_; if (*Error) goto E0;

    *Error = rebmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmix->IniTheta_->pdf_[i] = pfNormal;
        }
        else
        if (!strcmp(pdf[i], "lognormal")) {
            rebmix->IniTheta_->pdf_[i] = pfLognormal;
        }
        else
        if (!strcmp(pdf[i], "Weibull")) {
            rebmix->IniTheta_->pdf_[i] = pfWeibull;
        }
        else
        if (!strcmp(pdf[i], "gamma")) {
            rebmix->IniTheta_->pdf_[i] = pfGamma;
        }
        else
        if (!strcmp(pdf[i], "Gumbel")) {
            rebmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rebmix->IniTheta_->pdf_[i] = pfvonMises;
        }
        else
        if (!strcmp(pdf[i], "binomial")) {
            rebmix->IniTheta_->pdf_[i] = pfBinomial;
        }
        else
        if (!strcmp(pdf[i], "Poisson")) {
            rebmix->IniTheta_->pdf_[i] = pfPoisson;
        }
        else
        if (!strcmp(pdf[i], "Dirac")) {
            rebmix->IniTheta_->pdf_[i] = pfDirac;
        }
        else
        if (!strcmp(pdf[i], "uniform")) {
            rebmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    rebmix->MixTheta_ = new CompnentDistribution*[(unsigned int)(*c)];

    *Error = NULL == rebmix->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        *Error = NULL == rebmix->MixTheta_[i]; if (*Error) goto E0;

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) if (rebmix->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 2) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_ + 2; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;

        if (i < rebmix->length_pdf_) {
            for (j = 0; j < rebmix->n_; j++) Y[i][j] = rebmix->Y_[i][j];
        }
    }

    *Error = rebmix->PreprocessingKDE(h, Y);

    if (*Error) goto E0;

    rebmix->cmax_ = *c;

    logV = (FLOAT)0.0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        logV += (FLOAT)log(h[i]);
    }

    *Error = rebmix->InformationCriterionKDE(logV,
                                             Y,
                                             *c,
                                             rebmix->W_,
                                             rebmix->MixTheta_,
                                             IC,
                                             logL,
                                             M,
                                             D);

    if (*Error) goto E0;

E0: if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 2; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;
} // RInformationCriterionKDEMIX

void RInformationCriterionHMIX(double *h,            // Sides of the hypersquare.
                               double *y0,           // Origins.
                               double *ymin,         // Minimum observations.
                               double *ymax,         // Maximum observations.
                               int    *k,            // Total number of bins.
                               char   **Criterion,   // Information criterion type.
                               int    *c,            // Number of components.
                               double *W,            // Component weights.
                               int    *length_pdf,   // Length of pdf.
                               int    *length_Theta, // Length of Theta.
                               int    *length_theta, // Length of Theta[i].
                               char   **pdf,         // Parametric family types.
                               double *Theta,        // Component parameters.
                               int    *n,            // Number of observations.
                               double *x,            // Dataset.
                               double *IC,           // Information criterion.
                               double *logL,         // log-likelihood.
                               int    *M,            // Degrees of freedom.
                               double *D,            // Total of positive relative deviations.
                               int    *Error)        // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    FLOAT  logV;
    int    i, j, l, m;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    if (!strcmp(Criterion[0], "AIC"))
        rebmix->Criterion_ = icAIC;
    else
    if (!strcmp(Criterion[0], "AIC3"))
        rebmix->Criterion_ = icAIC3;
    else
    if (!strcmp(Criterion[0], "AIC4"))
        rebmix->Criterion_ = icAIC4;
    else
    if (!strcmp(Criterion[0], "AICc"))
        rebmix->Criterion_ = icAICc;
    else
    if (!strcmp(Criterion[0], "BIC"))
        rebmix->Criterion_ = icBIC;
    else
    if (!strcmp(Criterion[0], "CAIC"))
        rebmix->Criterion_ = icCAIC;
    else
    if (!strcmp(Criterion[0], "HQC"))
        rebmix->Criterion_ = icHQC;
    else
    if (!strcmp(Criterion[0], "MDL2"))
        rebmix->Criterion_ = icMDL2;
    else
    if (!strcmp(Criterion[0], "MDL5"))
        rebmix->Criterion_ = icMDL5;
    else
    if (!strcmp(Criterion[0], "AWE"))
        rebmix->Criterion_ = icAWE;
    else
    if (!strcmp(Criterion[0], "CLC"))
        rebmix->Criterion_ = icCLC;
    else
    if (!strcmp(Criterion[0], "ICL"))
        rebmix->Criterion_ = icICL;
    else
    if (!strcmp(Criterion[0], "PC"))
        rebmix->Criterion_ = icPC;
    else
    if (!strcmp(Criterion[0], "ICL-BIC"))
        rebmix->Criterion_ = icICLBIC;
    else
    if (!strcmp(Criterion[0], "D"))
        rebmix->Criterion_ = icD;
    else
    if (!strcmp(Criterion[0], "SSE"))
        rebmix->Criterion_ = icSSE;
    else {
        *Error = 1; goto E0;
    }

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == rebmix->W_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    *Error = NULL == rebmix->IniTheta_; if (*Error) goto E0;

    rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (int*)malloc(rebmix->length_Theta_ * sizeof(int));

    *Error = NULL == rebmix->length_theta_; if (*Error) goto E0;

    *Error = rebmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmix->IniTheta_->pdf_[i] = pfNormal;
        }
        else
        if (!strcmp(pdf[i], "lognormal")) {
            rebmix->IniTheta_->pdf_[i] = pfLognormal;
        }
        else
        if (!strcmp(pdf[i], "Weibull")) {
            rebmix->IniTheta_->pdf_[i] = pfWeibull;
        }
        else
        if (!strcmp(pdf[i], "gamma")) {
            rebmix->IniTheta_->pdf_[i] = pfGamma;
        }
        else
        if (!strcmp(pdf[i], "Gumbel")) {
            rebmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rebmix->IniTheta_->pdf_[i] = pfvonMises;
        }
        else
        if (!strcmp(pdf[i], "binomial")) {
            rebmix->IniTheta_->pdf_[i] = pfBinomial;
        }
        else
        if (!strcmp(pdf[i], "Poisson")) {
            rebmix->IniTheta_->pdf_[i] = pfPoisson;
        }
        else
        if (!strcmp(pdf[i], "Dirac")) {
            rebmix->IniTheta_->pdf_[i] = pfDirac;
        }
        else
        if (!strcmp(pdf[i], "uniform")) {
            rebmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    rebmix->MixTheta_ = new CompnentDistribution* [(unsigned int)(*c)];

    *Error = NULL == rebmix->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        *Error = NULL == rebmix->MixTheta_[i]; if (*Error) goto E0;

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) if (rebmix->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 1) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    *Error = rebmix->PreprocessingH(h, y0, ymin, ymax, k, Y);

    if (*Error) goto E0;

    rebmix->cmax_ = *c;

    logV = (FLOAT)0.0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        logV += (FLOAT)log(h[i]);
    }

    *Error = rebmix->InformationCriterionH(logV,
                                           *k,
                                           Y,
                                           *c,
                                           rebmix->W_,
                                           rebmix->MixTheta_,
                                           IC,
                                           logL,
                                           M,
                                           D);

    if (*Error) goto E0;

E0: if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;
} // RInformationCriterionHMIX

void RCombineComponentsMIX(int    *c,            // Number of components.
                           double *W,            // Component weights.
                           int    *length_pdf,   // Length of pdf.
                           int    *length_Theta, // Length of Theta.
                           int    *length_theta, // Length of Theta[i].
                           char   **pdf,         // Parametric family types.
                           double *Theta,        // Component parameters.
                           int    *n,            // Number of observations.
                           double *x,            // Dataset.
                           double *tau,          // Conditional probabilities.
                           int    *F,            // From components.
                           int    *T,            // To components.
                           double *EN,           // Entropy.
                           double *ED,           // Entropy decrease.
                           int    *Error)        // Error code.
{
    Rebmix *rebmix = NULL;
    int    i, j, l, m;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == rebmix->W_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    *Error = NULL == rebmix->IniTheta_; if (*Error) goto E0;

    rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (int*)malloc(rebmix->length_Theta_ * sizeof(int));

    *Error = NULL == rebmix->length_theta_; if (*Error) goto E0;

    *Error = rebmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmix->IniTheta_->pdf_[i] = pfNormal;
        }
        else
        if (!strcmp(pdf[i], "lognormal")) {
            rebmix->IniTheta_->pdf_[i] = pfLognormal;
        }
        else
        if (!strcmp(pdf[i], "Weibull")) {
            rebmix->IniTheta_->pdf_[i] = pfWeibull;
        }
        else
        if (!strcmp(pdf[i], "gamma")) {
            rebmix->IniTheta_->pdf_[i] = pfGamma;
        }
        else
        if (!strcmp(pdf[i], "Gumbel")) {
            rebmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rebmix->IniTheta_->pdf_[i] = pfvonMises;
        }
        else
        if (!strcmp(pdf[i], "binomial")) {
            rebmix->IniTheta_->pdf_[i] = pfBinomial;
        }
        else
        if (!strcmp(pdf[i], "Poisson")) {
            rebmix->IniTheta_->pdf_[i] = pfPoisson;
        }
        else
        if (!strcmp(pdf[i], "Dirac")) {
            rebmix->IniTheta_->pdf_[i] = pfDirac;
        }
        else
        if (!strcmp(pdf[i], "uniform")) {
            rebmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    rebmix->MixTheta_ = new CompnentDistribution* [(unsigned int)(*c)];

    *Error = NULL == rebmix->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        *Error = NULL == rebmix->MixTheta_[i]; if (*Error) goto E0;

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) if (rebmix->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    rebmix->cmax_ = *c;

    *Error = rebmix->CombineComponents(*c,
                                       rebmix->W_,
                                       rebmix->MixTheta_,
                                       tau,
                                       F,
                                       T,
                                       EN,
                                       ED);

    if (*Error) goto E0;

E0: if (rebmix) delete rebmix;
} // RCombineComponentsMIX

void RvonMisesPdf(int *n, double *y, double *Mean, double *Kappa, double *f)
{
    FLOAT A;
    int   i;

    A = Pi2 * BesselI0(*Kappa);

    for (i = 0; i < *n; i++) {
        if (y[i] > Pi2) {
            f[i] = (FLOAT)0.0;
        }
        else
        if (y[i] < (FLOAT)0.0) {
            f[i] = (FLOAT)0.0;
        }
        else {
            f[i] = (FLOAT)exp(*Kappa * (FLOAT)cos(y[i] - *Mean)) / A;
        }
    }
} // RvonMisesPdf

void RvonMisesCdf(int *n, double *y, double *Mean, double *Kappa, double *F)
{
    FLOAT A[3], Io, In, I0, I1;
    int   i, j, Error;

    I0 = BesselI0(*Kappa); I1 = BesselI1(*Kappa);

    A[0] = (FLOAT)1.0 / Pi2; A[1] = (FLOAT)2.0 * A[0] / I0;

    for (i = 0; i < *n; i++) {
        if (y[i] > Pi2) {
            F[i] = (FLOAT)1.0;
        }
        else
        if (y[i] < (FLOAT)0.0) {
            F[i] = (FLOAT)0.0;
        }
        else {
            Io = I0; In = I1; j = 1; F[i] = A[0] * y[i]; Error = 1;
            while ((j <= ItMax) && Error) {
                F[i] += A[1] * In * ((FLOAT)sin(j * (y[i] - *Mean)) + (FLOAT)sin(j * *Mean)) / j;

                if (In < Eps) Error = 0;

                A[2] = Io - (FLOAT)2.0 * j * In / *Kappa; Io = In; In = A[2];

                j++;
            }

            if (F[i] > (FLOAT)1.0) {
                F[i] = (FLOAT)1.0;
            }
            else
            if (F[i] < (FLOAT)0.0) {
                F[i] = (FLOAT)0.0;
            }
        }
    }
} // RvonMisesCdf

void RGumbelPdf(int *n, double *y, double *Mean, double *Beta, double *f)
{
    FLOAT A;
    int   i;

    for (i = 0; i < *n; i++) {
		A = -(y[i] - *Mean) / (*Beta);

        f[i] = (FLOAT)exp(A - (FLOAT)exp(A)) / (*Beta);
    }
} // RGumbelPdf

void RGumbelCdf(int *n, double *y, double *Mean, double *Beta, double *F)
{
    FLOAT A;
	int   i;

	for (i = 0; i < *n; i++) {
	    A = -(y[i] - *Mean) / (*Beta);
      
		F[i] = (FLOAT)exp(-(FLOAT)exp(A));
    }
} // RGumbelCdf

/// Panic Branislav & Marko Nagode.

/// Optimal number of bins calculation.

void Roptbins(int    *d,           // Number of independent random variables.
              int    *n,           // Number of observations.
              double *x,           // Dataset.
              char   **Rule,       // Rule.
              int    *length_y0,   // Length of y0.
              double *y0,          // Origins.
              int    *length_ymin, // Length of ymin.
              double *ymin,        // Minimum observations.
              int    *length_ymax, // Length of ymax.
              double *ymax,        // Maximum observations.
              int    *kmin,        // Minimum number of bins.
              int    *kmax,        // Maximum number of bins.
              int    *opt_k,       // Optimal number of bins.
              int    *Error)       // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    FLOAT  *h = NULL;
    int    *opt_k_tmp = NULL;
    int    i, j, l;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->length_pdf_ = *d;

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    rebmix->y0_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    *Error = NULL == rebmix->y0_; if (*Error) goto E0;

    if (*length_y0 > 0) {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->y0_[i] = y0[i];
        }
    }

    rebmix->ymin_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    *Error = NULL == rebmix->ymin_; if (*Error) goto E0;

    if (*length_ymin > 0) {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymin_[i] = ymin[i];
        }
    }
    else {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymin_[i] = rebmix->Y_[i][0];

            for (j = 1; j < rebmix->n_; j++) {
                if (rebmix->Y_[i][j] < rebmix->ymin_[i]) rebmix->ymin_[i] = rebmix->Y_[i][j];
            }
        }
    }

    rebmix->ymax_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    *Error = NULL == rebmix->ymax_; if (*Error) goto E0;

    if (*length_ymax > 0) {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymax_[i] = ymax[i];
        }
    }
    else {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymax_[i] = rebmix->Y_[i][0];

            for (j = 1; j < rebmix->n_; j++) {
                if (rebmix->Y_[i][j] > rebmix->ymax_[i]) rebmix->ymax_[i] = rebmix->Y_[i][j];
            }
        }
    }

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 1) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    h = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    *Error = NULL == h; if (*Error) goto E0;

    rebmix->Initialize();

    if (!strcmp(Rule[0], "Sturges")) {
        opt_k[0] = (int)ceil((FLOAT)1.0 + (FLOAT)log2((FLOAT)rebmix->n_));

        for (j = 1; j < rebmix->length_pdf_; j++) {
            opt_k[j] = opt_k[0];
        }
    }
    else
    if (!strcmp(Rule[0], "Log10")) {
        opt_k[0] = (int)ceil((FLOAT)10.0 * (FLOAT)log10((FLOAT)rebmix->n_));

        for (j = 1; j < rebmix->length_pdf_; j++) {
            opt_k[j] = opt_k[0];
        }
    }
    else
    if (!strcmp(Rule[0], "RootN")) {
        opt_k[0] = (int)ceil((FLOAT)2.0 * (FLOAT)sqrt((FLOAT)rebmix->n_));

        for (j = 1; j < rebmix->length_pdf_; j++) {
            opt_k[j] = opt_k[0];
        }
    }
    else
    if (!strcmp(Rule[0], "Knuth equal")) {
        FLOAT logp, logpopt, M;
        int   k;

        logpopt = -FLOAT_MAX;

        for (l = *kmin; l < *kmax + 1; l++) {
            for (j = 0; j < rebmix->length_pdf_; j++) {
                h[j] = (rebmix->ymax_[j] - rebmix->ymin_[j]) / l;

                if (*length_y0 == 0) {
                    rebmix->y0_[j] = rebmix->ymin_[j] + (FLOAT)0.5 * h[j];
                }
            }

            *Error = rebmix->PreprocessingH(h, rebmix->y0_, rebmix->ymin_, rebmix->ymax_, &k, Y);

            if (*Error) goto E0;

            if (k > rebmix->kmax_) {
                break;
            }

            M = (FLOAT)1.0;

            for (j = 0; j < rebmix->length_pdf_; j++) {
                M *= (FLOAT)l;
            }

            logp = (FLOAT)rebmix->n_ * (FLOAT)log(M) + (FLOAT)Gammaln((FLOAT)0.5 * M) - M * (FLOAT)Gammaln((FLOAT)0.5) - (FLOAT)Gammaln((FLOAT)rebmix->n_ + (FLOAT)0.5 * M);

            for (j = 0; j < k; j++) {
                logp += Gammaln(Y[rebmix->length_pdf_][j] + (FLOAT)0.5);
            }

            logp += Max(M - (FLOAT)k, (FLOAT)0.0) * Gammaln((FLOAT)0.5);

            if (logp > logpopt) {
                logpopt = logp; opt_k[0] = l;
            }
        }

        for (j = 1; j < rebmix->length_pdf_; j++) {
            opt_k[j] = opt_k[0];
        }
    }
    else
    if (!strcmp(Rule[0], "Knuth unequal")) {
        FLOAT logp, logpopt, M, M_broken;
        int   k, num_e, lim_broken, div, inc, Converged = 1;
        int   num_r, MAX_ITER_NUM = 131072;

        logpopt = -FLOAT_MAX;

        opt_k_tmp = (int*)malloc(rebmix->length_pdf_ * sizeof(int));

        *Error = NULL == opt_k_tmp; if (*Error) goto E0;

        for (j = 0; j < rebmix->length_pdf_; j++) {
            opt_k_tmp[j] = 1;
        }

        while (Converged) {
            for (j = 0; j < rebmix->length_pdf_; j++) {
                for (i = 0; i < rebmix->length_pdf_; i++) {
                    if (i != j) {
                        h[i] = (rebmix->ymax_[i] - rebmix->ymin_[i]) / opt_k_tmp[i];

                        if (*length_y0 == 0) {
                            rebmix->y0_[i] = rebmix->ymin_[i] + (FLOAT)0.5 * h[i];
                        }
                    }
                }

                opt_k[j] = opt_k_tmp[j];

                for (l = *kmin; l < *kmax + 1; l++) {
                    h[j] = (rebmix->ymax_[j] - rebmix->ymin_[j]) / l;

                    if (*length_y0 == 0) {
                        rebmix->y0_[j] = rebmix->ymin_[j] + (FLOAT)0.5 * h[j];
                    }

                    *Error = rebmix->PreprocessingH(h, rebmix->y0_, rebmix->ymin_, rebmix->ymax_, &k, Y);

                    if (*Error) goto E0;

                    if (k > rebmix->kmax_ && opt_k_tmp[j] != 1) {
                        break;
                    }

                    M = (FLOAT)1.0;

                    for (i = 0; i < rebmix->length_pdf_; i++) {
                        if (i != j) {
                            M *= (FLOAT)opt_k_tmp[i];
                        }
                        else{
                            M *= (FLOAT)l;
                        }
                    }

                    logp = (FLOAT)rebmix->n_ * (FLOAT)log(M) + (FLOAT)Gammaln((FLOAT)0.5 * M) - M * (FLOAT)Gammaln((FLOAT)0.5) - (FLOAT)Gammaln((FLOAT)rebmix->n_ + (FLOAT)0.5 * M);

                    for (i = 0; i < k; i++) {
                        logp += Gammaln(Y[rebmix->length_pdf_][i] + (FLOAT)0.5);
                    }

                    logp += Max(M - (FLOAT)k, (FLOAT)0.0) * Gammaln((FLOAT)0.5);

                    if (logp > logpopt || opt_k_tmp[j] == 1) {
                        logpopt = logp; opt_k_tmp[j] = l;
                    }
                }

                if (j == 0 && opt_k_tmp[j] == opt_k[j]) {
                    for (i = 0; i < rebmix->length_pdf_; i++) opt_k[i] = opt_k_tmp[i];

                    Converged = 0; break;
                }
            }
        }

        *kmin = *kmax = opt_k[0];

        for (j = 1; j < rebmix->length_pdf_; j++) {
            if (opt_k[j] < *kmin) {
                *kmin = opt_k[j];
            }
            else
            if (opt_k[j] > *kmax) {
                *kmax = opt_k[j];
            }
        }

        num_e = *kmax - *kmin + 1;

        num_r = (int)floor(pow((FLOAT)num_e, rebmix->length_pdf_));

        if (num_r > MAX_ITER_NUM || num_r <= 0) {
            for (j = *kmax - 1; j >= *kmin; j--) {
                *kmax = j;

                num_e = *kmax - *kmin + 1;

                num_r = (int)floor(pow((FLOAT)num_e, rebmix->length_pdf_));

                if (num_r <= MAX_ITER_NUM && num_r > 0) {
                    break;
                }
            }
        }

        lim_broken = 0; M_broken = (FLOAT)0.0;

        if (*kmax == *kmin) goto E0;

        for (i = 0; i < num_r; i++) {
            for (j = rebmix->length_pdf_ - 1; j >= 0; j--) {
                div = (int)floor(pow((FLOAT)num_e, j));

                inc = (i / div) % num_e;

                opt_k_tmp[j] = *kmin + inc;
            }

            for (j = 0; j < rebmix->length_pdf_; j++) {
                h[j] = (rebmix->ymax_[j] - rebmix->ymin_[j]) / opt_k_tmp[j];

                if (*length_y0 == 0) {
                    rebmix->y0_[j] = rebmix->ymin_[j] + (FLOAT)0.5 * h[j];
                }
            }

            M = (FLOAT)1.0;

            for (j = 0; j < rebmix->length_pdf_; j++) {
                M *= (FLOAT)opt_k_tmp[j];
            }

            if (lim_broken) {
                if (M > M_broken) continue;
            }

            *Error = rebmix->PreprocessingH(h, rebmix->y0_, rebmix->ymin_, rebmix->ymax_, &k, Y);

            if (*Error) goto E0;

            if (k > rebmix->kmax_) {
                lim_broken = 1; M_broken = M; continue;
            }

            logp = (FLOAT)rebmix->n_ * (FLOAT)log(M) + (FLOAT)Gammaln((FLOAT)0.5 * M) - M * (FLOAT)Gammaln((FLOAT)0.5) - (FLOAT)Gammaln((FLOAT)rebmix->n_ + (FLOAT)0.5 * M);

            for (j = 0; j < k; j++) {
                logp += Gammaln(Y[rebmix->length_pdf_][j] + (FLOAT)0.5);
            }

            logp += Max(M - (FLOAT)k, (FLOAT)0.0) * Gammaln((FLOAT)0.5);

            if (logp > logpopt) {
                logpopt = logp;

                for (j = 0; j < rebmix->length_pdf_; j++) opt_k[j] = opt_k_tmp[j];
            }
        }
    }

E0: if (opt_k_tmp) free(opt_k_tmp);
    
    if (h) free(h);
    
    if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;
} // Roptbins

/// Binning of data.

void Rbins(int    *d,           // Number of independent random variables.
           int    *n,           // Number of observations.
           double *x,           // Dataset.
           int    *length_y0,   // Length of y0.
           double *y0,          // Origins.
           int    *length_ymin, // Length of ymin.
           double *ymin,        // Minimum observations.
           int    *length_ymax, // Length of ymax.
           double *ymax,        // Maximum observations.
           int    *k,           // Number of bins.
           int    *length_y,    // Length of y.
           double *y,           // Binned dataset.
           int    *Error)       // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    FLOAT  *h = NULL;
    int    i, j, l;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->length_pdf_ = *d;

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    rebmix->y0_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    *Error = NULL == rebmix->y0_; if (*Error) goto E0;

    if (*length_y0 > 0) {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->y0_[i] = y0[i];
        }
    }

    rebmix->ymin_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    *Error = NULL == rebmix->ymin_; if (*Error) goto E0;

    if (*length_ymin > 0) {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymin_[i] = ymin[i];
        }
    }
    else {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymin_[i] = rebmix->Y_[i][0];

            for (j = 1; j < rebmix->n_; j++) {
                if (rebmix->Y_[i][j] < rebmix->ymin_[i]) rebmix->ymin_[i] = rebmix->Y_[i][j];
            }
        }
    }

    rebmix->ymax_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    *Error = NULL == rebmix->ymax_; if (*Error) goto E0;

    if (*length_ymax > 0) {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymax_[i] = ymax[i];
        }
    }
    else {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymax_[i] = rebmix->Y_[i][0];

            for (j = 1; j < rebmix->n_; j++) {
                if (rebmix->Y_[i][j] > rebmix->ymax_[i]) rebmix->ymax_[i] = rebmix->Y_[i][j];
            }
        }
    }

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 1) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->n_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    h = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    *Error = NULL == h; if (*Error) goto E0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        h[j] = (rebmix->ymax_[j] - rebmix->ymin_[j]) / k[j];

        if (*length_y0 == 0) {
            rebmix->y0_[j] = rebmix->ymin_[j] + (FLOAT)0.5 * h[j];
        }
    }

    *Error = rebmix->PreprocessingH(h, rebmix->y0_, rebmix->ymin_, rebmix->ymax_, length_y, Y);

    if (*Error) goto E0;

    i = 0;

    for (j = 0; j < rebmix->length_pdf_ + 1; j++) {
        for (l = 0; l < *length_y; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

E0: if (h) free(h);

    if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;
} // Rbins

/// End

}
