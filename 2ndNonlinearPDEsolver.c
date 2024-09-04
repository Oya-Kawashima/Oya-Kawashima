#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

double dz=0.2e-3;//m
double dt=1.;//s
double sigma=5.67e-8;//Stefan Boltzmann Constant [W/m2*K4]
double RichardsonConstant=1.20e6;// [A/m2*K2]
double BoltzmannnConstant=1.38e-23;// [J/K]

double T0=300;//Ambient temp [K]

double functionF(double T, double d, double e, double alpha, double rho, double k, double Current);
double functionforiginal(double T, double d, double e, double alpha, double rho, double k, double Current);
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
    
    FILE *readfile;
    char *readfilename = "./CalcParameter.csv";
    int ret;
    char buf[31][40];
    int Mode;
    double data[30];

    readfile = fopen( readfilename, "r" );
    if( readfile == NULL ){
        printf( "Failed to open Parameter file -%s-\n", readfilename );
        return -1;
    }
    printf("\n");

    fscanf(readfile, "%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%s", buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7], buf[8], buf[9], buf[10], buf[11], buf[12], buf[13], buf[14], buf[15], buf[16], buf[17], buf[18], buf[19], buf[20], buf[21], buf[22], buf[23], buf[24], buf[25], buf[26], buf[27], buf[28], buf[29], buf[30]);
    fscanf(readfile, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &Mode, &data[0], &data[1], &data[2], &data[3], &data[4], &data[5], &data[6], &data[7], &data[8], &data[9], &data[10], &data[11], &data[12], &data[13], &data[14], &data[15], &data[16], &data[17], &data[18], &data[19], &data[20], &data[21], &data[22], &data[23], &data[24], &data[25], &data[26], &data[27], &data[28], &data[29]);
    
    printf("%s %d\n",buf[0], Mode);
    for(int n=1; n<31; n++){
        printf("%s %1.2e\n",buf[n], data[n-1]);
    }

    printf("\n");
    fclose(readfile);
    
    
    double Icur;                Icur = data[0];//[A]
    double k_fil;               k_fil = data[1];//[W/m*K]
    double d_fil;               d_fil = data[2];//[m]
    double e_fil;               e_fil = data[3];//
    double l_half_fil;               l_half_fil = data[4]/2.;//[m]
    double rho_fil;             rho_fil = data[5];//[Ohm*m]
    double alpha_fil;           alpha_fil = data[6];//[Ohm*m/K]
    double rhoMASS_fil;         rhoMASS_fil = data[7];//[kg/m^3]
    double SpecificHeat_fil;    SpecificHeat_fil = data[8];//[J/kg*K]
    double workfunction_fil;    workfunction_fil = data[9];//[J]
    double coating_thickness;   coating_thickness = data[10];//[m]
    double coating_half_length; coating_half_length = data[11]/2.;//[m]
    
    double k_Pin;               k_Pin = data[12];//[W/m*K]
    double d_Pin;               d_Pin = data[13];//[m]
    double e_Pin;               e_Pin = data[14];//
    double l_Pin;               l_Pin = data[15];//[m]
    double rho_Pin;             rho_Pin = data[16];//[Ohm*m]
    double alpha_Pin;           alpha_Pin = data[17];//[Ohm*m/K]
    double rhoMASS_Pin;         rhoMASS_Pin = data[18];//[kg/m^3]
    double SpecificHeat_Pin;    SpecificHeat_Pin = data[19];//[J/kg*K]
    
    double k_cable;             k_cable = data[20];//[W/m*K]
    double d_cable;             d_cable = data[21];//[m]
    double e_cable;             e_cable = data[22];//
    double l_cable;             l_cable = data[23];//[m]
    double rho_cable;           rho_cable = data[24];//[Ohm*m]
    double alpha_cable;         alpha_cable = data[25];//[Ohm*m/K]
    double rhoMASS_cable;       rhoMASS_cable = data[26];//[kg/m^3]
    double SpecificHeat_cable;  SpecificHeat_cable = data[27];//[J/kg*K]
    
    double A_fil;//filament Emissive surface area[m^2]
    double S_fil;//filament cross section area[m^2]
    double A_Pin;//Pin Emissive surface area[m^2]
    double S_Pin;//Pin cross section area[m^2]
    double A_cable;//cable Emissive surface area[m^2]
    double S_cable;//cable Emissive surface area[m^2]
    A_cable = M_PI * d_cable * dz;
    A_Pin = M_PI * d_Pin * dz;
    A_fil = M_PI * d_fil * dz;
    S_cable = M_PI * (d_cable/2.) * (d_cable/2.);
    S_Pin = M_PI * (d_Pin/2.) * (d_Pin/2.);
    S_fil = M_PI * (d_fil/2.) * (d_fil/2.);
    
    
    double REcontact_Pin_fil_circuit;
    REcontact_Pin_fil_circuit = data[28];
    double REcontact_cable_Pin_fil_circuit = 1.00;
    REcontact_cable_Pin_fil_circuit = data[29];
    
    double REcontact_Pin_fil;
    double REcontact_cable_Pin;
    REcontact_Pin_fil = (REcontact_Pin_fil_circuit - 2.*rho_fil*l_half_fil/S_fil -2.*rho_Pin*l_Pin/S_Pin)/2.;
    REcontact_cable_Pin = (REcontact_cable_Pin_fil_circuit - REcontact_Pin_fil_circuit -2.*rho_cable*l_cable/S_cable)/2.;
    printf("REcontact_cable_Pin %f\n", REcontact_cable_Pin);
    printf("REcontact_Pin_fil %f\n", REcontact_Pin_fil);
    
    double totalpower=0.;
    double totalcurrent=0.;
    double fil_power=0.;
    double cable_power=0.;
    double Pin_power=0.;
    double sumRI2=0.;
    double EmitcalcT=0.;
    
    
    double omegafil;
    double omegaPin;
    double omegacable;
    double Breakvalue;
    double Time_delta;
    double Time_Breakvalue;
    int timel_MAX;
    
    if(Mode == 1){
        printf("Steady State\n");
        omegafil = 0.;
        omegaPin = 0.;
        omegacable = 0.;
        timel_MAX=1;
        Breakvalue=1.e-6;
        Time_Breakvalue=0.;
    }else if(Mode == 2){
        printf("Time evolution\n");
        omegafil = dz*dz*rhoMASS_fil*SpecificHeat_fil/dt/k_fil;
        omegaPin = dz*dz*rhoMASS_Pin*SpecificHeat_Pin/dt/k_Pin;
        omegacable = dz*dz*rhoMASS_cable*SpecificHeat_cable/dt/k_cable;
        timel_MAX=1000;
        Breakvalue=1.e-6;
        Time_Breakvalue=1.e-6;
    }
    
    int cable_ini;
    int cable_end;
    int Pin_ini;
    int Pin_end;
    int fil_ini;
    int fil_end;
    int DIM;
    
    cable_ini=0;
    cable_end=(int)((l_cable)/dz);
    Pin_ini=(int)((l_cable)/dz)+1;
    Pin_end=(int)((l_cable+l_Pin)/dz)+1;
    fil_ini=(int)((l_cable+l_Pin)/dz)+2;
    fil_end=(int)((l_cable+l_Pin+l_half_fil)/dz)+2;
    DIM=fil_end+1;
    
    printf("Cable:  %d  to  %d\n",cable_ini, cable_end);
    printf("Pin304:  %d  to  %d\n",Pin_ini, Pin_end);
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
    
    double RTcontact_cable_Pin;
    double RTcontact_Pin_fil;
    
    //Initial Temperature Condition
    for(i=cable_ini; i<=cable_end; i++){
        Temperature[i]=300.;
        Temperature_before[i]=300.;
    }
    for(i=Pin_ini; i<=Pin_end; i++){
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
    for(i = Pin_ini; i <= Pin_end; i++){
        fprintf(fp_time,"%f %f\n",(i-1)*dz*1e3,Temperature[i]);
    }
    for(i = fil_ini; i <= fil_end; i++){
        fprintf(fp_time,"%f %f\n",(i-2)*dz*1e3,Temperature[i]);
    }
    for(i = fil_end; i >= fil_ini; i--){
        fprintf(fp_time,"%f %f\n",(2.*fil_end-i-1)*dz*1e3,Temperature[i]);
    }
    for(i = Pin_end; i >= Pin_ini; i--){
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
            
            RTcontact_cable_Pin = RTcontact(Temperature[cable_end], Temperature[Pin_ini], alpha_cable, alpha_Pin, k_cable, k_Pin, rho_cable, rho_Pin, REcontact_cable_Pin);
            RTcontact_Pin_fil = RTcontact(Temperature[Pin_end], Temperature[fil_ini], alpha_Pin, alpha_fil, k_Pin, k_fil, rho_Pin, rho_fil, REcontact_Pin_fil);
            
            if(LOOPcount % displaycount == 0){
                printf("RTcontact_cable_Pin %f\n", RTcontact_cable_Pin);
                printf("RTcontact_Pin_fil %f\n", RTcontact_Pin_fil);
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
                            functionF(Temperature[i], d_cable, e_cable, alpha_cable, rho_cable, k_cable, Icur);
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
                            functionF(Temperature[i], d_cable, e_cable, alpha_cable, rho_cable, k_cable, Icur);
                        }
                        if(j == 2){
                            A[i][j]=1.;
                        }
                    } else if(i == cable_end){
                        if(j == 0){
                            A[i][j]= k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz;
                        }
                        if(j == 1){
                            A[i][j]= -1./RTcontact_cable_Pin - k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz;
                        }
                        if(j == 2){
                            A[i][j]= 1./RTcontact_cable_Pin;
                        }
                    }
                    
                    else if(i == Pin_ini){
                        if(j == 0){
                            A[i][j]= 1./RTcontact_cable_Pin;
                        }
                        if(j == 1){
                            A[i][j]= -1./RTcontact_cable_Pin - k_Pin*M_PI*(d_Pin/2.)*(d_Pin/2.)/dz;
                        }
                        if(j == 2){
                            A[i][j]= k_Pin*M_PI*(d_Pin/2.)*(d_Pin/2.)/dz;
                        }
                    } else if(Pin_ini < i && Pin_end > i){
                        if(j == 0){
                            A[i][j]=1.;
                        }
                        if(j == 1){
                            A[i][j]= -2.-omegaPin+dz*dz*
                            functionF(Temperature[i], d_Pin, e_Pin, alpha_Pin, rho_Pin, k_Pin, Icur);
                        }
                        if(j == 2){
                            A[i][j]=1.;
                        }
                    } else if(i == Pin_end){
                        if(j == 0){
                            A[i][j]= k_Pin*M_PI*(d_Pin/2.)*(d_Pin/2.)/dz;
                        }
                        if(j == 1){
                            A[i][j]= -1./RTcontact_Pin_fil - k_Pin*M_PI*(d_Pin/2.)*(d_Pin/2.)/dz;
                        }
                        if(j == 2){
                            A[i][j]= 1./RTcontact_Pin_fil;
                        }
                    }
                    
                    else if(i == fil_ini){
                        if(j == 0){
                            A[i][j]= 1./RTcontact_Pin_fil;
                        }
                        if(j == 1){
                            A[i][j]= -1./RTcontact_Pin_fil - k_fil*M_PI*(d_fil/2.)*(d_fil/2.)/dz;
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
                            functionF(Temperature[i], d_fil, e_fil, alpha_fil, rho_fil, k_fil, Icur);
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
                            functionF(Temperature[i], d_fil, e_fil, alpha_fil, rho_fil, k_fil, Icur);
                        }
                        if(j == 2){
                            A[i][j]=0.;
                        }
                    }
                }
            }
            
            for(i = 0; i < DIM; i++){
                
                if(i == cable_ini){
                    b[i]= -( (-2.-omegacable)*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_cable, e_cable, alpha_cable, rho_cable, k_cable, Icur)+300.+omegacable*Temperature_before[i]);
                } else if(cable_ini < i && cable_end > i){
                    b[i]= -(Temperature[i-1] +(-2.-omegacable)*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_cable, e_cable, alpha_cable, rho_cable, k_cable, Icur) +omegacable*Temperature_before[i]);
                } else if(i == cable_end){
                    b[i]= -(Temperature[i-1]*k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz -Temperature[i]/RTcontact_cable_Pin -Temperature[i]*k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz +Temperature[i+1]/RTcontact_cable_Pin);
                }
                
                else if(i == Pin_ini){
                    b[i]= -(Temperature[i-1]/RTcontact_cable_Pin -Temperature[i]/RTcontact_cable_Pin -Temperature[i]*k_Pin*M_PI*(d_Pin/2.)*(d_Pin/2.)/dz +Temperature[i+1]*k_Pin*M_PI*(d_Pin/2.)*(d_Pin/2.)/dz);
                } else if(Pin_ini < i && Pin_end > i){
                    b[i]= -(Temperature[i-1] +(-2.-omegaPin)*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_Pin, e_Pin, alpha_Pin, rho_Pin, k_Pin, Icur)+omegaPin*Temperature_before[i]);
                } else if(i == Pin_end){
                    b[i]= -(Temperature[i-1]*k_Pin*M_PI*(d_Pin/2.)*(d_Pin/2.)/dz -Temperature[i]/RTcontact_Pin_fil -Temperature[i]*k_Pin*M_PI*(d_Pin/2.)*(d_Pin/2.)/dz +Temperature[i+1]/RTcontact_Pin_fil);
                }
                
                else if(i == fil_ini){
                    b[i]= -(Temperature[i-1]/RTcontact_Pin_fil -Temperature[i]/RTcontact_Pin_fil -Temperature[i]*k_fil*M_PI*(d_fil/2.)*(d_fil/2.)/dz +Temperature[i+1]*k_fil*M_PI*(d_fil/2.)*(d_fil/2.)/dz);
                } else if(fil_ini < i && fil_end > i){
                    b[i]= -(Temperature[i-1] +(-2.-omegafil)*Temperature[i] + Temperature[i+1] + dz*dz*functionforiginal(Temperature[i], d_fil, e_fil, alpha_fil, rho_fil, k_fil, Icur)+omegafil*Temperature_before[i]);
                } else if(i == fil_end){
                    b[i]= -(Temperature[i-1] +(-1.-omegafil)*Temperature[i] + dz*dz*functionforiginal(Temperature[i], d_fil, e_fil, alpha_fil, rho_fil, k_fil, Icur) +omegafil*Temperature_before[i]);
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
            for(i = Pin_ini; i <= Pin_end; i++){
                fprintf(fp_time,"%f %f\n",(i-1)*dz*1e3,Temperature[i]);
            }
            for(i = fil_ini; i <= fil_end; i++){
                fprintf(fp_time,"%f %f\n",(i-2)*dz*1e3,Temperature[i]);
            }
            for(i = fil_end; i >= fil_ini; i--){
                fprintf(fp_time,"%f %f\n",(2.*fil_end-i-1)*dz*1e3,Temperature[i]);
            }
            for(i = Pin_end; i >= Pin_ini; i--){
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
        if(Pin_ini <= i && i <= Pin_end){
            Pin_power += pow(Temperature[i],4.)*e_Pin*sigma*A_Pin;
        }
        else if(fil_ini <= i && i <= fil_end){
            fil_power += pow(Temperature[i],4.)*e_fil*sigma*A_fil;
        }
        else if(cable_ini <= i && i <= cable_end){
            cable_power += pow(Temperature[i],4.)*e_cable*sigma*A_cable;
        }
    }
    totalpower = 2.*(Pin_power + fil_power + cable_power + (Temperature[0]-300.)*k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz);
    
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
    for(i = Pin_ini; i <= Pin_end; i++){
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
    for(i = Pin_end; i >= Pin_ini; i--){
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
       }else if (Pin_ini <= i && i <= Pin_end){
           sumRI2 += dz * rho_Pin * (1. + alpha_Pin * (Temperature[i]-T0)) * Icur * Icur/(M_PI*d_Pin/2.*d_Pin/2.);
       }
       else if (cable_ini <= i && i <= cable_end){
           sumRI2 += dz * rho_cable * (1. + alpha_cable * (Temperature[i]-T0)) * Icur * Icur/(M_PI*d_cable/2.*d_cable/2.);
       }
   }
    
    
    printf("consumption TOTAL %1.2f[W]  RI^2 %1.2f[W]  Electronemission %1.2f[mA]\n",totalpower , sumRI2, totalcurrent*1e3);
    printf("              fil %1.2f[W]  Pin %1.2f[W] cable %1.2f[W] cableTHR %1.2f[W]\n", 2.*fil_power, 2.*Pin_power, 2.*cable_power, 2.*(Temperature[0]-300.)*k_cable*M_PI*(d_cable/2.)*(d_cable/2.)/dz/*, filTHR %1.2f[W] 2.*(Temperature[DIM-1]-Temperature[DIM-2])*k_fil*M_PI*(d_fil/2.)*(d_fil/2.)/dz*/);
    
    fclose(GPLT);
    fclose(fp);
    fclose(fp2);
    fclose(fp3);
}
            



double functionF(double T, double d, double e, double alpha, double rho, double k, double Current){
    double r;
    double Ans;
    
    r = d/2.;
    Ans = 1./k *(alpha * rho * Current * Current/(M_PI*r*r)/(M_PI*r*r) - 8. * sigma * e * T*T*T / r );
    
    return Ans;
}

double functionforiginal(double T, double d, double e, double alpha, double rho, double k, double Current){
    double r;
    double Ans;
    
    r = d/2.;
    Ans = 1./k *(rho * (1. + alpha * (T-T0)) * Current * Current/(M_PI*r*r)/(M_PI*r*r) - 2. * sigma * e * (pow(T,4.) - pow(T0,4.))/ r );
    
    return Ans;
}

double RTcontact(double TA, double TB, double alphaA, double alphaB, double kA, double kB, double rhoA, double rhoB, double REcontact){
    double Ans;
    
    Ans = REcontact/( (kA+kB)*(rhoA*(1+alphaA*TA) + rhoB*(1+alphaB*TB) ) );
    
    return Ans;
}
