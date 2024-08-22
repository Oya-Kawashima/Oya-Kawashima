#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

double dz=0.2e-3;//mm
double Icur=2.5;//[A]
double sigma=5.67e-8;//Stefan Boltzmann Constant [W/m2*K4]
double T0=300;//Ambient temp [K]


double RichardsonConstant=1.20e6;// [A/m2*K2]
double BoltzmannnConstant=1.38e-23;// [J/K]

double functionF(double T, double d, double e, double alpha, double rho, double k);
double functionforiginal(double T, double d, double e, double alpha, double rho, double k);
double RTcontact(double TA, double TB, double alphaA, double alphaB, double kA, double kB, double rhoA, double rhoB, double REcontact);

int main()
{
    FILE *GPLT;
    FILE *fp;
    
    GPLT=popen("gnuplot -persist","w");
    fp=fopen("T.dat","w");
    
    double totalpower=0.;
    double totalcurrent=0.;
    double fil_power=0.;
    double cable_power=0.;
    double SUS_power=0.;
    
    
    double REcontact_SUS_fil_circuit = 0.23;
    double REcontact_cable_SUS_fil_circuit = 1.10;
    
    double k_fil=147;//filament thermal conductivity[W/m*K]
    double e_fil=0.3;//filament emissisvity
    double A_fil;//filament Emissive surface area[m^2]
    double S_fil;//filament cross section area[m^2]
    double d_fil=0.125e-3;//filament diameter[m]
    double l_fil=4.e-3;//filament half length[m]
    double alpha_fil = 7e-3;//temp coefficient of resistivity [Ohm*m/K]
    double rho_fil = 5e-8;//volume resistivity [Ohm*m]
    double workfunction_fil = 5.0e-19;//[J] 2.0eV for Y2O3
    
    double k_SUS=16.7;//SUS304 thermal conductivity[W/m*K]
    double e_SUS=0.3;//SUS304 emissisvity
    double A_SUS;//SUS304 Emissive surface area[m^2]
    double S_SUS;//SUS304 cross section area[m^2]
    double d_SUS = 1.6e-3;//SUS304 diameter [m]
    double l_SUS = 10.5e-3;//SUS304 length[m]
    double alpha_SUS = 0.07;//temp coefficient of resistivity [Ohm*m/K]
    double rho_SUS = 71.e-8;//volume resistivity [Ohm*m]
    
    double k_cable=65.;//Tinned soft copper wire[W/m*K]
    double e_cable=0.1;//cable emissisvity
    double A_cable;//cable Emissive surface area[m^2]
    double S_cable;//cable Emissive surface area[m^2]
    double d_cable = 0.5e-3;//cable diameter
    double l_cable = 20.e-3;//cable length[m]
    double alpha_cable = 0.007;//temp coefficient of resistivity [Ohm*m/K]
    double rho_cable = 1.55e-8;//volume resistivity [Ohm*m]
    
    
    A_cable = M_PI * d_cable * dz;
    A_SUS = M_PI * d_SUS * dz;
    A_fil = M_PI * d_fil * dz;
    
    S_cable = M_PI * (d_cable/2.) * (d_cable/2.);
    S_SUS = M_PI * (d_SUS/2.) * (d_SUS/2.);
    S_fil = M_PI * (d_fil/2.) * (d_fil/2.);
    
    double REcontact_SUS_fil;
    double REcontact_cable_SUS;
    
    REcontact_SUS_fil = (REcontact_SUS_fil_circuit - 2.*rho_fil*l_fil/S_fil -2.*rho_SUS*l_SUS/S_SUS)/2.; //0.0363
    
    REcontact_cable_SUS = (REcontact_cable_SUS_fil_circuit - REcontact_SUS_fil_circuit -2.*rho_cable*l_cable/S_cable)/2.; //0.03705
    printf("REcontact_cable_SUS %f\n", REcontact_cable_SUS);
    printf("REcontact_SUS_fil %f\n", REcontact_SUS_fil);
    
    int cable_ini;
    int cable_end;
    int SUS_ini;
    int SUS_end;
    int fil_ini;
    int fil_end;
    int DIM;
    
    cable_ini=0;
    cable_end=(int)((l_cable)/dz);
    SUS_ini=(int)((l_cable)/dz)+1;
    SUS_end=(int)((l_cable+l_SUS)/dz)+1;
    fil_ini=(int)((l_cable+l_SUS)/dz)+2;
    fil_end=(int)((l_cable+l_SUS+l_fil)/dz)+2;
    
    DIM=fil_end+1;
    
    printf("Cable:  %d  to  %d\n",cable_ini, cable_end);
    printf("SUS304:  %d  to  %d\n",SUS_ini, SUS_end);
    printf("Filament:  %d  to  %d\n",fil_ini, fil_end);
    printf("DIM/all    %d \n",DIM);
    
    double error=0.;
    double preerror=0.;
    int count=0;
       
    int i, j, m, n, dim, maxline;
    double temp, reserve, sum1, sum2;
    
    double *b;
    double *Temperature;
    double *DeltaTemperature;
    double *y;
    b=(double*)calloc(DIM,sizeof(double));
    Temperature=(double*)calloc(DIM,sizeof(double));
    DeltaTemperature=(double*)calloc(DIM,sizeof(double));
    y=(double*)calloc(DIM,sizeof(double));
    
    double **A;
    A=(double**)calloc(DIM,sizeof(double*));
    for(i=0; i<DIM; i++){
        A[i]=(double*)calloc(3,sizeof(double));
    }
    
    double **U;
    double **L;
    A=(double**)calloc(DIM,sizeof(double*));
    for(int i=0; i<DIM; i++){
        A[i]=(double*)calloc(5,sizeof(double));
    }
    U=(double**)calloc(DIM,sizeof(double*));
    L=(double**)calloc(DIM,sizeof(double*));
    for(int i=0; i<DIM; i++){
        U[i]=(double*)calloc(DIM,sizeof(double));
        L[i]=(double*)calloc(DIM,sizeof(double));
    }
    
    
    
    double DeltaTemperaturenorm=0.;
    double Breakvalue=1.e-1;
    
    double RTcontact_cable_SUS;
    double RTcontact_SUS_fil;
    
    for(i=cable_ini; i<=cable_end; i++){
        Temperature[i]=300.;
    }
    for(i=SUS_ini; i<=SUS_end; i++){
        Temperature[i]=300.;//+3.*(i-SUS_ini);
    }
    for(i=fil_ini; i<=fil_end; i++){
        Temperature[i]=300.;
    }
    
    int LOOPcount = 0;
    int LOOPABORT = 20000;
    double TEMPABORT = 3000;
    double CURRENTABORT = 100e-3;
    int displaycount = 1000;
    
    //  Endless loop til convergence
    while(1){
        LOOPcount += 1;
        
        RTcontact_cable_SUS = RTcontact(Temperature[cable_end], Temperature[SUS_ini], alpha_cable, alpha_SUS, k_cable, k_SUS, rho_cable, rho_SUS, REcontact_cable_SUS);
        RTcontact_SUS_fil = RTcontact(Temperature[SUS_end], Temperature[fil_ini], alpha_SUS, alpha_fil, k_SUS, k_fil, rho_SUS, rho_fil, REcontact_SUS_fil);
        
        if(LOOPcount % displaycount == 0){
            printf("RTcontact_cable_SUS %f\n", RTcontact_cable_SUS);
            printf("RTcontact_SUS_fil %f\n", RTcontact_SUS_fil);
        }
        
        
    //  input matrix-------------------------------------------
    //  A
        for(i = 0; i < DIM; i++){
            for(j = 0; j < 2; j++){
                
                if(i == cable_ini){
                    if(j == 0){
                        A[i][j]= 0.;
                    }
                    if(j == 1){
                        A[i][j]= -2.+dz*dz*
                        functionF(Temperature[i], d_cable, e_cable, alpha_cable, rho_cable, k_cable);
                    }
                    if(j == 2){
                        A[i][j]= 1.;
                    }
                } else if(cable_ini < i && cable_end > i){
                    if(j == 0){
                        A[i][j]=1.;
                    }
                    if(j == 1){
                        A[i][j]= -2.+dz*dz*
                        functionF(Temperature[i], d_cable, e_cable, alpha_cable, rho_cable, k_cable);
                    }
                    if(j == 2){
                        A[i][j]=1.;
                    }
                } else if(i == cable_end){
                    if(j == 0){
                        A[i][j]= k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz;
                    }
                    if(j == 1){
                        A[i][j]= -1./RTcontact_cable_SUS - k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz;
                    }
                    if(j == 2){
                        A[i][j]= 1./RTcontact_cable_SUS;
                    }
                }
                
                else if(i == SUS_ini){
                    if(j == 0){
                        A[i][j]= 1./RTcontact_cable_SUS;
                    }
                    if(j == 1){
                        A[i][j]= -1./RTcontact_cable_SUS - k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz;
                    }
                    if(j == 2){
                        A[i][j]= k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz;
                    }
                } else if(SUS_ini < i && SUS_end > i){
                    if(j == 0){
                        A[i][j]=1.;
                    }
                    if(j == 1){
                        A[i][j]= -2.+dz*dz*
                        functionF(Temperature[i], d_SUS, e_SUS, alpha_SUS, rho_SUS, k_SUS);
                    }
                    if(j == 2){
                        A[i][j]=1.;
                    }
                } else if(i == SUS_end){
                    if(j == 0){
                        A[i][j]= k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz;
                    }
                    if(j == 1){
                        A[i][j]= -1./RTcontact_SUS_fil - k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz;
                    }
                    if(j == 2){
                        A[i][j]= 1./RTcontact_SUS_fil;
                    }
                }
                
                else if(i == fil_ini){
                    if(j == 0){
                        A[i][j]= 1./RTcontact_SUS_fil;
                    }
                    if(j == 1){
                        A[i][j]= -1./RTcontact_SUS_fil - k_fil*M_PI*(d_fil/2.)*(d_fil/2.)/dz;
                    }
                    if(j == 2){
                        A[i][j]= k_fil*M_PI*(d_fil/2.)*(d_fil/2.)/dz;
                    }
                }
                else if(fil_ini < i && fil_end > i){
                    if(j == 0){
                        A[i][j]=1.;
                    }
                    if(j == 1){
                        A[i][j]= -2.+dz*dz*
                        functionF(Temperature[i], d_fil, e_fil, alpha_fil, rho_fil, k_fil);
                    }
                    if(j == 2){
                        A[i][j]=1.;
                    }
                }
                else if(i == fil_end){
                    if(j == 0){
                        A[i][j]=1.;
                    }
                    if(j == 1){
                        A[i][j]= -1.+dz*dz*
                        functionF(Temperature[i], d_fil, e_fil, alpha_fil, rho_fil, k_fil);
                    }
                    if(j == 2){
                        A[i][j]=0.;
                    }
                }
            }
        }
        
        for(i = 0; i < DIM; i++){
            
            if(i == cable_ini){
                b[i]= -(Temperature[i-1] -2.*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_cable, e_cable, alpha_cable, rho_cable, k_cable)+300.);
            } else if(cable_ini < i && cable_end > i){
                b[i]= -(Temperature[i-1] -2.*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_cable, e_cable, alpha_cable, rho_cable, k_cable));
            } else if(i == cable_end){
                b[i]= -(Temperature[i-1]*k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz -Temperature[i]/RTcontact_cable_SUS -Temperature[i]*k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz +Temperature[i+1]/RTcontact_cable_SUS);
            }
            
            else if(i == SUS_ini){
                b[i]= -(Temperature[i-1]/RTcontact_cable_SUS -Temperature[i]/RTcontact_cable_SUS -Temperature[i]*k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz +Temperature[i+1]*k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz);
            } else if(SUS_ini < i && SUS_end > i){
                b[i]= -(Temperature[i-1] -2.*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_SUS, e_SUS, alpha_SUS, rho_SUS, k_SUS));
            } else if(i == SUS_end){
                b[i]= -(Temperature[i-1]*k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz -Temperature[i]/RTcontact_SUS_fil -Temperature[i]*k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz +Temperature[i+1]/RTcontact_SUS_fil);
            }
            
            else if(i == fil_ini){
                b[i]= -(Temperature[i-1]/RTcontact_SUS_fil -Temperature[i]/RTcontact_SUS_fil -Temperature[i]*k_fil*M_PI*(d_fil/2.)*(d_fil/2.)/dz +Temperature[i+1]*k_fil*M_PI*(d_fil/2.)*(d_fil/2.)/dz);
            } else if(fil_ini < i && fil_end > i){
                b[i]= -(Temperature[i-1] -2.*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_fil, e_fil, alpha_fil, rho_fil, k_fil));
            } else if(i == fil_end){
                b[i]= -(Temperature[i-1] -Temperature[i] + dz*dz*functionforiginal(Temperature[i], d_fil, e_fil, alpha_fil, rho_fil, k_fil));
            }
        }
        //  -------------------------------------------------------
        //  (A=LU)x = b ----------------------------------------------------
            
        for(i = 0; i < DIM; i++){
            L[i][i] = 1.0;
        }
        
        for(i = 0; i < DIM; i++){
            for(j = 0; j < DIM; j++){
                if( i <= j){
                    sum1 = 0.0;
                    for(m = 0; m <= i-1; m++){
                        sum1 += L[i][m] * U[m][j];
                    }
                    
                    if(j == i){
                        U[i][j] = A[i][1] - sum1;
                    }
                    else if(j == i+1){
                        U[i][j] = A[i][2] - sum1;
                    }
                    else{
                        U[i][j] = - sum1;
                    }
                }
                else if(i > j){
                    sum2 = 0.0;
                    for(n = 0; n <= j-1; n++){
                        sum2 += L[i][n] * U[n][j];
                    }
                    
                    if(j == i-1){
                        L[i][j] = (A[i][0] - sum2) / U[j][j];
                        
                        if ( isnan(L[i][j] ) ) {
                            printf("A Not a Number. : WARNING ERROR %d %d\n", i, j);
                            break;
                        }
                    }
                    else{
                        L[i][j] = - sum2 / U[j][j];
                        if ( isnan(L[i][j]) ) {
                            printf("B Not a Number. : WARNING ERROR %d %d\n", i, j);
                            break;
                        }
                    }
                }
            }
        }
        
        //  -------------------------------------------------------
        //  solve the problem--------------------------------------
        //  Ly = b-------------------------------------------------
        for(i = 0; i < DIM; i++){
            y[i] = b[i];
        }

        for(j = 0; j < DIM-1; j++){
            for(i = j+1; i < DIM; i++){
                y[i] -= y[j] * L[i][j];
            }
        }

        //  Ux = y
        for(i = 0; i < DIM; i++){
            DeltaTemperature[i] = y[i];
        }

        for(j = DIM-1; j >= 0; j--){
            DeltaTemperature[j] /= U[j][j];
            for (i = 0; i <= j-1; i++){
                DeltaTemperature[i] -= U[i][j] * DeltaTemperature[j];
            }
        }
        
        //LOOP END judge
        for(i=0; i<DIM; i++){
            DeltaTemperaturenorm += DeltaTemperature[i] * DeltaTemperature[i];
        }
        
        if ( LOOPcount > LOOPABORT ) {
            printf("LOOPcount over\n\n");
            break;
        }
        
        if ( DeltaTemperaturenorm < Breakvalue ) {
            printf("LOOP %d: norm %1.2f ,filtop temp %1.2f\n\n", LOOPcount, DeltaTemperaturenorm, Temperature[DIM-1]);
            printf("Convergence\n\n");
            break;
        }
        else{
            if(LOOPcount % displaycount == 0){
                printf("LOOP %d: norm %1.2f ,filtop temp %1.2f\n\n", LOOPcount, DeltaTemperaturenorm, Temperature[DIM-1]);
            }
            DeltaTemperaturenorm = 0;
            for(i=0; i<DIM; i++){
                Temperature[i] += DeltaTemperature[i];
            }
        }
        
        if ( Temperature[DIM-1] > TEMPABORT) {
            printf("TEMP over\n\n");
            break;
        }
        
        totalcurrent = 0.;
        for(i = fil_ini; i <= fil_end; i++){
            totalcurrent += A_fil * RichardsonConstant * Temperature[i] * Temperature[i] * exp(-workfunction_fil/Temperature[i]/BoltzmannnConstant);
        }
        if ( totalcurrent > CURRENTABORT) {
            printf("CURRENT over\n\n");
            break;
        }
        
    }//  Endless loop til convergence
    
    
    for(i = 0; i < DIM; i++){
        if(SUS_ini <= i && i <= SUS_end){
            SUS_power += pow(Temperature[i],4.)*e_SUS*sigma*A_SUS;
        }
        else if(fil_ini <= i && i <= fil_end){
            fil_power += pow(Temperature[i],4.)*e_fil*sigma*A_fil;
        }
        else if(cable_ini <= i && i <= cable_end){
            cable_power += pow(Temperature[i],4.)*e_cable*sigma*A_cable;
        }
    }
    totalpower = 2.*(SUS_power + fil_power + cable_power + (Temperature[0]-300.)*k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz);
    
    fprintf(GPLT,"set title font 'Arial,24'\n");
    fprintf(GPLT,"set tics font 'Arial,14'\n");
    fprintf(GPLT,"set xlabel font 'Arial,14'\n");
    fprintf(GPLT,"set ylabel font 'Arial,14'\n");
    
    fprintf(GPLT,"unset key\n");
    fprintf(GPLT,"set xr [0.:%f]\n", DIM*dz*1e3/*2*DIM*dz*1e3*/);
    fprintf(GPLT,"set yr [200.:%f]\n", Temperature[DIM-1]+100.);
    fprintf(GPLT,"set xlabel 'length [mm]'\n");
    fprintf(GPLT,"set ylabel 'temperature [K]'\n");
    fprintf(GPLT,"set title 'Re-W current %1.2f[A]      TOTAL %1.2f[W]'\n",Icur,totalpower);
    
    
    fprintf(GPLT,"p '-' w p pt 7 ps 1 lc 'black'\n");
    
    
    
    
    for(i = cable_ini; i <= cable_end; i++){
        fprintf(GPLT,"%f %f\n",i*dz*1e3,Temperature[i]);
        //fprintf(GPLT,"%f %f\n",(2*DIM - i)*dz*1e3,Temperature[i]);
    }
    for(i = SUS_ini; i <= SUS_end; i++){
        fprintf(GPLT,"%f %f\n",(i-1)*dz*1e3,Temperature[i]);
        //fprintf(GPLT,"%f %f\n",(2*DIM - i + 1)*dz*1e3,Temperature[i]);
    }
    for(i = fil_ini; i <= fil_end; i++){
        fprintf(GPLT,"%f %f\n",(i-2)*dz*1e3,Temperature[i]);
        //fprintf(GPLT,"%f %f\n",(2*DIM - i + 2)*dz*1e3,Temperature[i]);
    }
    
    
    
    
     printf("\n\nTEMPERATURE PROFILE\n\n");
     
     for(i = 0; i < DIM; i++){
        //fprintf(GPLT,"%f %f\n",i*dz,Temperature[i]);
        if(fil_ini <= i && i <= fil_end){
            printf("%f %f FIL\n",i*dz,Temperature[i]);
        }else{
            printf("%f %f\n",i*dz,Temperature[i]);
        }
    }
    fprintf(GPLT,"e\n");
    
    
    printf("consumption TOTAL %1.2f[W]  Electronemission %1.2f[mA]\n",totalpower , totalcurrent*1e3);
    printf("              fil %1.2f[W]  SUS %1.2f[W] cable %1.2f[W] cableTHR %1.2f[W]\n", 2.*fil_power, 2.*SUS_power, 2.*cable_power, (Temperature[0]-300.)*k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz);
    
}
            

double functionF(double T, double d, double e, double alpha, double rho, double k){
    double r;
    double Ans;
    
    r = d/2.;
    Ans = 1./k *(alpha * rho * Icur * Icur/(M_PI*r*r)/(M_PI*r*r) - 8. * sigma * e * T*T*T / r );
    
    return Ans;
}

double functionforiginal(double T, double d, double e, double alpha, double rho, double k){
    double r;
    double Ans;
    
    r = d/2.;
    Ans = 1./k *(rho * (1. + alpha * (T-T0)) * Icur * Icur/(M_PI*r*r)/(M_PI*r*r) - 2. * sigma * e * (pow(T,4.) - pow(T0,4.))/ r );
    
    return Ans;
}

double RTcontact(double TA, double TB, double alphaA, double alphaB, double kA, double kB, double rhoA, double rhoB, double REcontact){
    double Ans;
    
    Ans = REcontact/( (kA+kB)*(rhoA*(1+alphaA*TA) + rhoB*(1+alphaB*TB) ) );
    
    return Ans;
}
