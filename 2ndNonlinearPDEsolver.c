#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

double dz=0.2e-3;//m
double dt=1.;//s
double Icur=2.898;//[A]
double sigma=5.67e-8;//Stefan Boltzmann Constant [W/m2*K4]
double T0=300;//Ambient temp [K]

int Mode = 2;
//Mode1 = Static solver
//Mode2 = Time solver

double RichardsonConstant=1.20e6;// [A/m2*K2]
double BoltzmannnConstant=1.38e-23;// [J/K]

double functionF(double T, double d, double e, double alpha, double rho, double k);
double functionforiginal(double T, double d, double e, double alpha, double rho, double k);
double RTcontact(double TA, double TB, double alphaA, double alphaB, double kA, double kB, double rhoA, double rhoB, double REcontact);

int main()
{
    FILE *GPLT;
    FILE *fp;
    FILE *fp2;
    FILE *fp3;
    FILE *fp_time;
    
    GPLT=popen("gnuplot -persist","w");
    fp=fopen("circuitT.dat","w");
    fp2=fopen("coatingT.dat","w");
    fp3=fopen("Filtop_with_Time.dat","w");
    
    double totalpower=0.;
    double totalcurrent=0.;
    double fil_power=0.;
    double cable_power=0.;
    double SUS_power=0.;
    double sumRI2=0.;
    
    double EmitcalcT=0.;
    double REcontact_SUS_fil_circuit = 0.30;
    double REcontact_cable_SUS_fil_circuit = 1.00;
    
    double k_fil=147;//filament thermal conductivity[W/m*K]
    double e_fil=0.33;//filament emissisvity
    double A_fil;//filament Emissive surface area[m^2]
    double S_fil;//filament cross section area[m^2]
    double d_fil=0.125e-3;//0.075e-3;//filament diameter[m]
    double l_fil=4.e-3;//filament half length[m]
    double alpha_fil = 5.1e-3;//temp coefficient of resistivity [Ohm*m/K]
    double rho_fil = 5.3e-8;//volume resistivity [Ohm*m]
    double rhoMASS_fil = 22.5e3;//density [kg/m^3]
    double SpecificHeat_fil = 130;//Specific Heat [J/kg*K]
    
    double workfunction_fil = 4.48e-19;//[J] 2.8eV
    double coating_thickness = 0.05e-3;
    double coating_half_length = 2.e-3;
    
    double k_SUS=10.7;//SUS304 thermal conductivity[W/m*K]
    double e_SUS=0.4;//SUS304 emissisvity
    double A_SUS;//SUS304 Emissive surface area[m^2]
    double S_SUS;//SUS304 cross section area[m^2]
    double d_SUS = 1.6e-3;//SUS304 diameter [m]
    double l_SUS = 10e-3;//SUS304 length[m]
    double alpha_SUS = 5.0e-3;//temp coefficient of resistivity [Ohm*m/K]
    double rho_SUS = 74.e-8;//volume resistivity [Ohm*m]
    double rhoMASS_SUS = 7.93e3;//density [kg/m^3]
    double SpecificHeat_SUS = 500;//Specific Heat [J/kg*K]
    
    double k_cable=65.;//Tinned soft copper wire[W/m*K]
    double e_cable=0.;//cable emissisvity
    double A_cable;//cable Emissive surface area[m^2]
    double S_cable;//cable Emissive surface area[m^2]
    double d_cable = 0.5e-3;//cable diameter
    double l_cable = 20.e-3;//cable length[m]
    double alpha_cable = 4.7e-3;//temp coefficient of resistivity [Ohm*m/K]
    double rho_cable = 1.5e-8;//volume resistivity [Ohm*m]
    double rhoMASS_cable = 8.96e3;//density [kg/m^3]
    double SpecificHeat_cable = 385;//Specific Heat [J/kg*K]
    
    double omegafil;
    double omegaSUS;
    double omegacable;
    double Breakvalue;
    
    double Time_delta;
    double Time_Breakvalue;
    int timel_MAX;
    
    if(Mode == 1){
        printf("Static Solver\n");
        omegafil = 0.;
        omegaSUS = 0.;
        omegacable = 0.;
        timel_MAX=1;
        Breakvalue=1.e-6;
        Time_Breakvalue=0.;
    }else if(Mode == 2){
        printf("Time Solver\n");
        omegafil = dz*dz*rhoMASS_fil*SpecificHeat_fil/dt/k_fil;
        omegaSUS = dz*dz*rhoMASS_SUS*SpecificHeat_SUS/dt/k_SUS;
        omegacable = dz*dz*rhoMASS_cable*SpecificHeat_cable/dt/k_cable;
        timel_MAX=1000;
        Breakvalue=1.e-6;
        Time_Breakvalue=1.e-6;
    }
    
    A_cable = M_PI * d_cable * dz;
    A_SUS = M_PI * d_SUS * dz;
    A_fil = M_PI * d_fil * dz;
    
    S_cable = M_PI * (d_cable/2.) * (d_cable/2.);
    S_SUS = M_PI * (d_SUS/2.) * (d_SUS/2.);
    S_fil = M_PI * (d_fil/2.) * (d_fil/2.);
    
    double REcontact_SUS_fil;
    double REcontact_cable_SUS;
    
    REcontact_SUS_fil = 0.12;//(REcontact_SUS_fil_circuit - 2.*rho_fil*l_fil/S_fil -2.*rho_SUS*l_SUS/S_SUS)/2.;
    
    REcontact_cable_SUS = 0.21;//(REcontact_cable_SUS_fil_circuit - REcontact_SUS_fil_circuit -2.*rho_cable*l_cable/S_cable)/2.;
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
    printf("Filament coating:  %d  to  %d\n",fil_end - (int)(coating_half_length/dz), fil_end);
    printf("DIM/all    %d \n",DIM);
    
    double error=0.;
    double preerror=0.;
    int count=0;
       
    int i, j, m, n, timel, dim, maxline;
    double temp, reserve, sum1, sum2;
    
    double *b;
    double *Temperature;
    double *Temperature_before;
    double *DeltaTemperature;
    double *y;
    b=(double*)calloc(DIM,sizeof(double));
    Temperature=(double*)calloc(DIM,sizeof(double));
    Temperature_before=(double*)calloc(DIM,sizeof(double));
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
    
    double RTcontact_cable_SUS;
    double RTcontact_SUS_fil;
    
    for(i=cable_ini; i<=cable_end; i++){
        Temperature[i]=300.;
        Temperature_before[i]=300.;
    }
    for(i=SUS_ini; i<=SUS_end; i++){
        Temperature[i]=300.;
        Temperature_before[i]=300.;
    }
    for(i=fil_ini; i<=fil_end; i++){
        Temperature[i]=300.;
        Temperature_before[i]=300.;
    }
    
    int LOOPcount = 0;
    int LOOPABORT = 20000;
    double TEMPABORT = 10000;
    double CURRENTABORT = 1e3;
    int displaycount = 1000;
    char fname[30];
    
    int Time_output=1;
    
    
    fprintf(fp3,"Time Filtop_Temperature\n");
    fprintf(fp3,"%f %f\n",(timel)*dt,Temperature[DIM-1]);
    sprintf(fname, "Timestep_%d.dat", 0);
    fp_time = fopen(fname, "w");
    for(i = cable_ini; i <= cable_end; i++){
        fprintf(fp_time,"%f %f\n",i*dz*1e3,Temperature[i]);
    }
    for(i = SUS_ini; i <= SUS_end; i++){
        fprintf(fp_time,"%f %f\n",(i-1)*dz*1e3,Temperature[i]);
    }
    for(i = fil_ini; i <= fil_end; i++){
        fprintf(fp_time,"%f %f\n",(i-2)*dz*1e3,Temperature[i]);
    }
    for(i = fil_end; i >= fil_ini; i--){
        fprintf(fp_time,"%f %f\n",(2.*fil_end-i-1)*dz*1e3,Temperature[i]);
    }
    for(i = SUS_end; i >= SUS_ini; i--){
        fprintf(fp_time,"%f %f\n",(2.*fil_end-i-2)*dz*1e3,Temperature[i]);
    }
    for(i = cable_end; i >= cable_ini; i--){
        fprintf(fp_time,"%f %f\n",(2.*fil_end-i-3)*dz*1e3,Temperature[i]);
    }
    fclose(fp_time);
    
    for(int timel = 0; timel < timel_MAX ; timel++){// Time Loop
        
        Time_delta = 0;
        printf("\n\n *********************************************\n");
        printf("Time Loop  %1.1f s",(timel+1)*dt);
        printf("\n *********************************************\n\n");
        
        sprintf(fname, "Timestep_%d.dat", timel+1);
        fp_time = fopen(fname, "w");
        LOOPcount = 0;
        
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
                            A[i][j]= -2.-omegacable+dz*dz*
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
                            A[i][j]= -2.-omegacable+dz*dz*
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
                            A[i][j]= -2.-omegaSUS+dz*dz*
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
                            A[i][j]= -2.-omegafil+dz*dz*
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
                            A[i][j]= -1.-omegafil+dz*dz*
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
                    b[i]= -( (-2.-omegacable)*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_cable, e_cable, alpha_cable, rho_cable, k_cable)+300.+omegacable*Temperature_before[i]);
                } else if(cable_ini < i && cable_end > i){
                    b[i]= -(Temperature[i-1] +(-2.-omegacable)*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_cable, e_cable, alpha_cable, rho_cable, k_cable) +omegacable*Temperature_before[i]);
                } else if(i == cable_end){
                    b[i]= -(Temperature[i-1]*k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz -Temperature[i]/RTcontact_cable_SUS -Temperature[i]*k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz +Temperature[i+1]/RTcontact_cable_SUS);
                }
                
                else if(i == SUS_ini){
                    b[i]= -(Temperature[i-1]/RTcontact_cable_SUS -Temperature[i]/RTcontact_cable_SUS -Temperature[i]*k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz +Temperature[i+1]*k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz);
                } else if(SUS_ini < i && SUS_end > i){
                    b[i]= -(Temperature[i-1] +(-2.-omegaSUS)*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_SUS, e_SUS, alpha_SUS, rho_SUS, k_SUS)+omegaSUS*Temperature_before[i]);
                } else if(i == SUS_end){
                    b[i]= -(Temperature[i-1]*k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz -Temperature[i]/RTcontact_SUS_fil -Temperature[i]*k_SUS*M_PI*(d_SUS/2.)*(d_SUS/2.)/dz +Temperature[i+1]/RTcontact_SUS_fil);
                }
                
                else if(i == fil_ini){
                    b[i]= -(Temperature[i-1]/RTcontact_SUS_fil -Temperature[i]/RTcontact_SUS_fil -Temperature[i]*k_fil*M_PI*(d_fil/2.)*(d_fil/2.)/dz +Temperature[i+1]*k_fil*M_PI*(d_fil/2.)*(d_fil/2.)/dz);
                } else if(fil_ini < i && fil_end > i){
                    b[i]= -(Temperature[i-1] +(-2.-omegafil)*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_fil, e_fil, alpha_fil, rho_fil, k_fil)+omegafil*Temperature_before[i]);
                } else if(i == fil_end){
                    b[i]= -(Temperature[i-1] +(-1.-omegafil)*Temperature[i] + dz*dz*functionforiginal(Temperature[i], d_fil, e_fil, alpha_fil, rho_fil, k_fil) +omegafil*Temperature_before[i]);
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
                printf("LOOP %d: norm %1.7f ,filtop temp %1.2f\n\n", LOOPcount, DeltaTemperaturenorm, Temperature[DIM-1]);
                printf("Convergence\n\n");
                break;
            }
            else{
                if(LOOPcount % displaycount == 0){
                    printf("LOOP %d: norm %1.7f ,filtop temp %1.2f\n\n", LOOPcount, DeltaTemperaturenorm, Temperature[DIM-1]);
                }
                DeltaTemperaturenorm = 0;
                for(i=0; i<DIM; i++){
                    Temperature[i] += DeltaTemperature[i];
                }
            }
            
            if ( Temperature[DIM-1] > TEMPABORT) {
                printf("TEMP over\n\n");
                printf("filtop temp %1.2f\n\n", Temperature[DIM-1]);
                break;
            }
            
            totalcurrent = 0.;
            for(i = fil_end; i >=  fil_end - (int)(coating_half_length/dz) ; i--){
                EmitcalcT = pow(e_fil * d_fil/(d_fil+coating_thickness) ,1./4.) * Temperature[i];
                
                totalcurrent += 2.* A_fil * RichardsonConstant * EmitcalcT * EmitcalcT * exp(-workfunction_fil/EmitcalcT/BoltzmannnConstant);
            }
            if ( totalcurrent > CURRENTABORT) {
                printf("CURRENT over\n\n");
                break;
            }
            
        }//  Endless loop til convergence
        
        fprintf(fp3,"%f %f\n",(timel+1)*dt,Temperature[DIM-1]);
        
        for(i=0; i<DIM; i++){
            Time_delta += (Temperature_before[i] - Temperature[i])*(Temperature_before[i] - Temperature[i]);
            Temperature_before[i] = Temperature[i];
        }
        
        if(timel % Time_output == 0){
            for(i = cable_ini; i <= cable_end; i++){
                fprintf(fp_time,"%f %f\n",i*dz*1e3,Temperature[i]);
            }
            for(i = SUS_ini; i <= SUS_end; i++){
                fprintf(fp_time,"%f %f\n",(i-1)*dz*1e3,Temperature[i]);
            }
            for(i = fil_ini; i <= fil_end; i++){
                fprintf(fp_time,"%f %f\n",(i-2)*dz*1e3,Temperature[i]);
            }
            for(i = fil_end; i >= fil_ini; i--){
                fprintf(fp_time,"%f %f\n",(2.*fil_end-i-1)*dz*1e3,Temperature[i]);
            }
            for(i = SUS_end; i >= SUS_ini; i--){
                fprintf(fp_time,"%f %f\n",(2.*fil_end-i-2)*dz*1e3,Temperature[i]);
            }
            for(i = cable_end; i >= cable_ini; i--){
                fprintf(fp_time,"%f %f\n",(2.*fil_end-i-3)*dz*1e3,Temperature[i]);
            }
        }
        
        fclose(fp_time);
        
        if(Time_delta < Time_Breakvalue) break;
        
    }// Time Loop
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
    fprintf(GPLT,"set xr [0.:%f]\n", 2*DIM*dz*1e3/*2*DIM*dz*1e3*/);
    fprintf(GPLT,"set yr [300.:%f]\n", 3000./*Temperature[DIM-1]+100.*/);
    fprintf(GPLT,"set xlabel 'length [mm]'\n");
    fprintf(GPLT,"set ylabel 'temperature [K]'\n");
    fprintf(GPLT,"set title 'Circuit current %1.2f[A]    Power %1.2f[W]'\n",Icur,totalpower);
    
    
    fprintf(GPLT,"p '-' w p pt 5 ps 1 lc 'black', '-' w p pt 5 ps 1 lc 'red'\n");
    
    
    for(i = cable_ini; i <= cable_end; i++){
        fprintf(GPLT,"%f %f\n",i*dz*1e3,Temperature[i]);
        fprintf(fp,"%f %f\n",i*dz*1e3,Temperature[i]);
    }
    for(i = SUS_ini; i <= SUS_end; i++){
        fprintf(GPLT,"%f %f\n",(i-1)*dz*1e3,Temperature[i]);
        fprintf(fp,"%f %f\n",(i-1)*dz*1e3,Temperature[i]);
    }
    for(i = fil_ini; i <= fil_end; i++){
        fprintf(GPLT,"%f %f\n",(i-2)*dz*1e3,Temperature[i]);
        fprintf(fp,"%f %f\n",(i-2)*dz*1e3,Temperature[i]);
    }
    
    for(i = fil_end; i >= fil_ini; i--){
        fprintf(GPLT,"%f %f\n",(2.*fil_end-i-1)*dz*1e3,Temperature[i]);
        fprintf(fp,"%f %f\n",(2.*fil_end-i-1)*dz*1e3,Temperature[i]);
    }
    for(i = SUS_end; i >= SUS_ini; i--){
        fprintf(GPLT,"%f %f\n",(2.*fil_end-i-2)*dz*1e3,Temperature[i]);
        fprintf(fp,"%f %f\n",(2.*fil_end-i-2)*dz*1e3,Temperature[i]);
    }
    for(i = cable_end; i >= cable_ini; i--){
        fprintf(GPLT,"%f %f\n",(2.*fil_end-i-3)*dz*1e3,Temperature[i]);
        fprintf(fp,"%f %f\n",(2.*fil_end-i-3)*dz*1e3,Temperature[i]);
    }
    fprintf(GPLT,"e\n");
    
    
    
    //printf("Coating %d\n", fil_end - (int)(coating_half_length/dz) );
    totalcurrent = 0.;
    for(i = fil_end - (int)(coating_half_length/dz); i <= fil_end ; i++){
        EmitcalcT = pow(e_fil * d_fil/(d_fil+coating_thickness) ,1./4.) * Temperature[i];
        
        printf("%d   EmitcalcT:%1.2f\n", i, EmitcalcT);
        totalcurrent += 2.* A_fil * RichardsonConstant * EmitcalcT * EmitcalcT * exp(-workfunction_fil/EmitcalcT/BoltzmannnConstant);
        
        fprintf(GPLT,"%f %f\n",(i-2)*dz*1e3,EmitcalcT);
        fprintf(fp2,"%f %f\n",(i-2)*dz*1e3,EmitcalcT);
    }
    for(i = fil_end; i >= fil_end - (int)(coating_half_length/dz); i--){
        EmitcalcT = pow(e_fil * d_fil/(d_fil+coating_thickness) ,1./4.) * Temperature[i];
        
        fprintf(GPLT,"%f %f\n",(2*fil_end-i-1)*dz*1e3,EmitcalcT);
        fprintf(fp2,"%f %f\n",(2*fil_end-i-1)*dz*1e3,EmitcalcT);
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
    
    for(i = 0; i < DIM; i++){
       if(fil_ini <= i && i <= fil_end){
           sumRI2 += dz * rho_fil * (1. + alpha_fil * (Temperature[i]-T0)) * Icur * Icur/(M_PI*d_fil/2.*d_fil/2.);
       }else if (SUS_ini <= i && i <= SUS_end){
           sumRI2 += dz * rho_SUS * (1. + alpha_SUS * (Temperature[i]-T0)) * Icur * Icur/(M_PI*d_SUS/2.*d_SUS/2.);
       }
       else if (cable_ini <= i && i <= cable_end){
           sumRI2 += dz * rho_cable * (1. + alpha_cable * (Temperature[i]-T0)) * Icur * Icur/(M_PI*d_cable/2.*d_cable/2.);
       }
   }
    
    
    printf("consumption TOTAL %1.2f[W]  RI^2 %1.2f[W]  Electronemission %1.2f[mA]\n",totalpower , sumRI2, totalcurrent*1e3);
    printf("              fil %1.2f[W]  SUS %1.2f[W] cable %1.2f[W] cableTHR %1.2f[W]\n", 2.*fil_power, 2.*SUS_power, 2.*cable_power, 2.*(Temperature[0]-300.)*k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz/*, filTHR %1.2f[W] 2.*(Temperature[DIM-1]-Temperature[DIM-2])*k_fil*M_PI*(d_fil/2.)*(d_fil/2.)/dz*/);
    
    fclose(GPLT);
    fclose(fp);
    fclose(fp2);
    fclose(fp3);
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
