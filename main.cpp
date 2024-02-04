/**
 * IMS 
 * T9: Spojitý model z oblasti fyziky a biologie
 * Téma: Hemolýza erytrocytů v důsledku ozvučení ultrazvukem
 * 
 * @file main.cpp
 * @brief Soubor se simulaci
 * @date 10.12.2022
 * 
 * @author Tomáš Janečka (xjanec35)
 * @author Jan Novák (xnovak3i)
 */





#include "simlib.h"
#include "cmath"


#define KONSTANTA_UMRTI 90
#define KONSTANTA_STARNUTI 0.8
#define OBSAH_HEMOGLOBINU 30 // pocet miligramu v 10^12 erytrocytech
#define KONSTANTA_ROZPADU 0.0083 // pri frekvenci 1, intenzite 0.1 a citlivosti erytrocytu 50 by melo dojit k uplne hemolyze za 285 minut 0.0083

// hodnoty pro dospeleho muze
double koncentraceErytrocytu = Normal(5, 0.8); // 5+-0.8 * 10^12 erytrocytu na litr krve => obsah hemoglobinu je 5+-0.7 * 30 miligramu
double koncentraceHemoglobinu = Uniform(40,80); // (40; 80) miligramu volneho hemoglobinu na litr krve

double rozpadErytrocytu;

double citlivostErytrocytu = 1;

double frekvence = 1; // 1
/*
    Nizka:   0.2-1    MHz
    Stredni: 1-10     MHz
    Vysoka:  10-100   MHz
*/ // 1 - 3 

double intenzita = 0.1; // 0.1
/*
    Nizka:   0-1  mW/cm^2
    Stredni: 1-3  mW/cm^2
    Vysoka:  3-10 mW/cm^2
    Nastup kavitace: > 1000 mW/cm^2 
*/ // 0 - 2

double cas_zareni = 285;
double cas_cekani = 76;


// STATS:
double poskozeniStarim = 0;
double poskozeniZarenim = 0;
double poskozeniKavitaci = 0;




class KoncentraceHemoglobinu: public Process
{
    double dHdt;
    double ubytekErytrocytu = rozpadErytrocytu;
    
    double k1 = -1*OBSAH_HEMOGLOBINU; // 10^12 erytrocytu obshaju 30 mg hemoglobinu
    void Behavior() 
    {  
        dHdt = ubytekErytrocytu * k1;
        if(dHdt < 0)
        {
            dHdt = dHdt * (-1);
        }
        koncentraceHemoglobinu = koncentraceHemoglobinu + dHdt;
    }
};


class Kavitace: public Process
{
    double nabeh = Normal(Uniform(1, 1000)/(frekvence*intenzita*60000), 500/(frekvence*intenzita*60000));
    const double rho = 1060; // hustota krve
    const double P0 = 1.06;    // tlak v kapalině
    const double Pv = 0.5;    // tlak nasycené páry v kapalině
    const double sigma = 0.072; // povrchové napětí kapaliny
    const double R0 = 1e-6;   // rovnovážný poloměr bubliny
    const double Pinf = 1.0;  // tlak kapaliny ve větší vzdálenosti od bubliny

    double R = R0;
    double dRdt = 0.0;
    double dt = 0.0001*frekvence*intenzita; // krok casu
    double cas = 0;
    double dopad;


    double RayleighPlesset(double R, double dRdt) 
    {
        return 1.5 * std::pow(dRdt, 2) + R * std::pow(dRdt, 2) - (P0 - Pv - 2 * sigma / R - Pinf) / rho
           - ((2 * sigma / R0 - P0) * std::pow(R0, 3) / (rho * std::pow(R, 3)));
    }


    void RungeKutta(double& R, double& dRdt, double dt) 
    {
        double k1, k2, k3, k4;

        k1 = dt * dRdt;
        k2 = dt * (dRdt + 0.5 * k1);
        k3 = dt * (dRdt + 0.5 * k2);
        k4 = dt * (dRdt + k3);

        R = R + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        dRdt = dRdt + RayleighPlesset(R, dRdt) * dt;
    }


    void Behavior()
    {
        while(cas <= nabeh)
        {
            RungeKutta(R, dRdt, dt);
            Wait(0.0001); // 6ms
            cas = cas + dt;
        }

        dopad = R*koncentraceErytrocytu*10*citlivostErytrocytu;
        poskozeniKavitaci = poskozeniKavitaci + dopad;
        if(koncentraceErytrocytu - dopad > 0)
        {
            koncentraceErytrocytu = koncentraceErytrocytu - dopad;
            rozpadErytrocytu = dopad;
            (new KoncentraceHemoglobinu)->Activate();
        }
        else
        {
            poskozeniKavitaci = poskozeniKavitaci + koncentraceErytrocytu;
            koncentraceErytrocytu = 0;
        }
    }
};




class Zareni: public Process
{
    //double treciSila;
    //double tlakovySpad;

    double dopadZareni = 1;
    double dEdt;
    double vznikKavitaceHranice;
    double vznik;
    void Behavior() 
    {
           while(Time <= cas_zareni + cas_cekani)
           {
                if(koncentraceErytrocytu > 0)
                {
                    dopadZareni = (frekvence*intenzita*citlivostErytrocytu)/koncentraceErytrocytu;
                    dEdt = (-dopadZareni)*KONSTANTA_ROZPADU;
                    koncentraceErytrocytu = koncentraceErytrocytu + dEdt;
                    poskozeniZarenim = poskozeniZarenim - dEdt;
                    if(koncentraceErytrocytu <= 0)
                    {
                        poskozeniZarenim = poskozeniZarenim + koncentraceErytrocytu;
                        koncentraceErytrocytu = 0;
                    }
                    else
                    {
                        rozpadErytrocytu = dEdt;  
                        (new KoncentraceHemoglobinu)->Activate();
                    }
                }
                else
                {
                    koncentraceErytrocytu = 0;
                }
                    
                vznikKavitaceHranice = (-1/(intenzita*intenzita+1))+1;
                vznik = Random();
                if(vznik <= vznikKavitaceHranice)
                {
                    (new Kavitace)->Activate();
                }
                Wait(1);
           }
    }
};


class Cekani: public Process
{
    /*
        Vlastnosti erytrocytu ovlivnenene starim, ktere maji vliv na miru hemolyzy ultrazvukem:
            Flexibilita
            Membranove lipidy
            Enzymalni aktivita
        Data jsou reprezentovany v %
    */
   double flexibilita;
   double membrana;
   double enzymy;
   double konstanta;
   double vek = 0;
   double dEdt;

    double var = 0;
    void Behavior() //override
    {
        while(Time < cas_cekani)
        {
            flexibilita = Uniform(vek, 100.000001);
            membrana = Uniform(vek, 100.000001);
            enzymy = Uniform(vek, 100.000001);
            if(vek < 100)
            {
                vek = vek + KONSTANTA_STARNUTI;
                if(vek > 100)
                {
                    vek = 100;
                }
            }
            konstanta = (flexibilita + membrana + enzymy)/3;

            if(citlivostErytrocytu > vek*konstanta/100)
            {
                citlivostErytrocytu = vek*konstanta/100 + (citlivostErytrocytu - vek*konstanta/100);
            }
            else
            {
                citlivostErytrocytu = vek*konstanta/100;
            }
            if(citlivostErytrocytu <= 0)
            {
                citlivostErytrocytu = 0.000001;
            }
            if(citlivostErytrocytu >= 100)
            {
                citlivostErytrocytu = 99.999999;
            }

            if(konstanta > KONSTANTA_UMRTI || flexibilita == 100 || membrana == 100 || enzymy == 100)
            {
                if(konstanta <= 95)
                {
                    konstanta = konstanta/10000;
                }
                if(konstanta > 95 && konstanta <= 97.5)
                {
                    konstanta = konstanta/1000;
                }
                if(konstanta > 97.5 && konstanta <= 99)
                {
                    konstanta = konstanta/500;
                }
                if(konstanta > 99)
                {
                    konstanta = konstanta/100;
                }
                dEdt = -konstanta*koncentraceErytrocytu;
                koncentraceErytrocytu = koncentraceErytrocytu + dEdt;
                rozpadErytrocytu = dEdt;
                poskozeniStarim = poskozeniStarim - dEdt; 
                if(koncentraceErytrocytu < 0)
                {
                    poskozeniStarim = poskozeniStarim + koncentraceErytrocytu;
                    koncentraceErytrocytu = 0;
                }
                (new KoncentraceHemoglobinu)->Activate();
            }

            Wait(1);
        }

        (new Zareni)->Activate();

    }
};


int main(int argc, char *argv[]) 
{

    cas_cekani = std::stod(argv[1]);
    cas_zareni = std::stod(argv[2]);
    frekvence = std::stod(argv[3]);
    intenzita = std::stod(argv[4]);
    if(frekvence <= 0 || intenzita <= 0)
    {
        printf("Frekvence a intenzita musi byt vesti nez 0\n");
        return 1;
    }

    printf("-------------------------ZACATEK-------------------------\n");

    printf("Pocatecni koncetrace erytrocytu:   %f\n", koncentraceErytrocytu);
    printf("Pocatecni koncentrace hemoglobinu: %f\n", koncentraceHemoglobinu);
    printf("Doba cekani:                       %f\n", cas_cekani);
    printf("Doba ozvuceni:                     %f\n", cas_zareni);
    printf("Celkova doba:                      %f\n", cas_cekani+cas_zareni);
    printf("Pocatecni citlivost erytrocytu:    %f\n", citlivostErytrocytu);
    printf("Frekvence ultrazvuku:              %f\n", frekvence);
    printf("Intenzita ultrazvuku:              %f\n", intenzita);

    printf("-------------------------ZACATEK-------------------------\n\n");

    Init(0, cas_cekani + cas_zareni + 1);

    (new Cekani)->Activate();

    Run();

    printf("--------------------------KONEC--------------------------\n");

    printf("Konecna koncetrace erytrocytu:     %f\n", koncentraceErytrocytu);
    printf("Konecna koncentrace hemoglobinu:   %f\n", koncentraceHemoglobinu);
    printf("Konecna citlivost erytrocytu:      %f\n", citlivostErytrocytu);
    printf("Poskozeni starim:                  %f\n", poskozeniStarim);
    printf("Poskozeni ozvucenim:               %f\n", poskozeniZarenim);
    printf("Poskozeni kavitaci:                %f\n", poskozeniKavitaci);

    printf("--------------------------KONEC--------------------------\n");

    return 0;
}





