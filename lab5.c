#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double m, v0, alpha_deg, alpha_rad, c, g;
double v0x, v0y, k, gv;

void fara_rezistenta(double *xM, double *yM, double *range)
{
    double tM, tA;
    tM = v0 * sin(alpha_rad) / g;
    *xM = v0 * cos(alpha_rad) * tM;
    *yM = v0 * sin(alpha_rad) * tM - g * tM * tM / 2.0;
    tA = 2.0 * v0 * sin(alpha_rad) / g;
    *range = v0 * cos(alpha_rad) * tA;
}

double x_rez(double t)
{
    return (v0x / k) * (1.0 - exp(-k * t));
}

double y_rez(double t)
{
    return ((v0y + gv) / k) * (1.0 - exp(-k * t)) - gv * t;
}

void cu_rezistenta_M(double *xM, double *yM, double *tM)
{
    *tM = (1.0 / k) * log((v0y + gv) / gv);
    *xM = x_rez(*tM);
    *yM = y_rez(*tM);
}

double timp_bataie_rezistenta(double eps)
{
    double t1 = 0.0, t2, tm;
    double y1, y2, ym;

    t2 = 2.0 * v0y / g;
    if (t2 <= 0.0) t2 = 1.0;
    while (y_rez(t2) > 0.0)
        t2 *= 2.0;

    y1 = y_rez(t1);
    y2 = y_rez(t2);

    while (t2 - t1 > eps)
    {
        tm = (t1 + t2) / 2.0;
        ym = y_rez(tm);

        if (y1 * ym <= 0.0)
        {
            t2 = tm;
            y2 = ym;
        }
        else
        {
            t1 = tm;
            y1 = ym;
        }
    }

    return (t1 + t2) / 2.0;
}

int main(void)
{
    double xM0, yM0, range0;
    double xM1, yM1, tM1, range1;
    double eps;

    printf("m (kg) = ");
    scanf("%lf", &m);
    printf("v0 (m/s) = ");
    scanf("%lf", &v0);
    printf("alpha (grade) = ");
    scanf("%lf", &alpha_deg);
    printf("c (coef. rezistenta, N*s/m) = ");
    scanf("%lf", &c);
    printf("g (m/s^2), de obicei 9.81 = ");
    scanf("%lf", &g);
    printf("eps (precizie pentru metoda injumatatirii) = ");
    scanf("%lf", &eps);

    alpha_rad = alpha_deg * M_PI / 180.0;
    v0x = v0 * cos(alpha_rad);
    v0y = v0 * sin(alpha_rad);

    if (c != 0.0)
    {
        k  = c / m;
        gv = g / k;
    }

    fara_rezistenta(&xM0, &yM0, &range0);

    printf("\nFARA rezistenta aerului:\n");
    printf("Punctul M (inaltime maxima):\n");
    printf("  xM0 = %lf m\n", xM0);
    printf("  yM0 = %lf m\n", yM0);
    printf("Bataia OA (raza de actiune):\n");
    printf("  OA0 = %lf m\n", range0);

    if (c > 0.0)
    {
        cu_rezistenta_M(&xM1, &yM1, &tM1);
        printf("\nCU rezistenta aerului (forta -c*v):\n");
        printf("Punctul M (inaltime maxima):\n");
        printf("  tM1 = %lf s\n", tM1);
        printf("  xM1 = %lf m\n", xM1);
        printf("  yM1 = %lf m\n", yM1);

        {
            double tA1 = timp_bataie_rezistenta(eps);
            range1 = x_rez(tA1);
            printf("\nBataia OA cu rezistenta:\n");
            printf("  tA1 = %lf s\n", tA1);
            printf("  OA1 = %lf m\n", range1);
        }
    }
    else
    {
        printf("\nCoeficientul de rezistenta c este 0 -> miscarea cu rezistenta nu se mai calculeaza.\n");
    }

    return 0;
}