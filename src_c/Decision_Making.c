/* TODO: continue with around line 100
 * use pointer to switch between slow/fast circuit,
 * and don't forget to update the reference line number.
 */

#include <stdio.h>
#include <stdlib.h>

#define BoundNe 8
#define BoundNeA1 9   ///
#define ShiftNe 14
#define InhibitionNe 1
#define FMNe 1
#define BaseFrequencyNe 1
#define CoupledNe 2
#define TotalNe 27
#define TotalNeA1 28  /// 
#define DoubTotalNe 54 /// 
#define Vr 10
#define Vth 130
#define Tau1 10
#define Tau2 5
#define LenFiring 100000 // SimulationTime / DeltaT
const int SimulationTime = 2000; // ms
const double DeltaT = 0.0100000000000; // ms
const int NoiseStrengthBase = 0;
int Velocity[2] = {1, 8}; // Length subject to change
int VelocityLen = 2;
//const double SaSpeed[4] = {1, 1.5, 2.5, 4};
const double SaSpeed[4] = {4, 2.5, 1.5, 1};
int SaSpeedLen = 4;
int SaPosition[4] = {1, 1, 1, 1};
int SaPositionLen = 4;
const double SaBorder=1.7308337;
const double StimulusStrength = 15.5;
int StimulusNeuron=1;
const double StimulationOnset= 10000; 
const double StimulationOffset= 13000;
const int DirectionSlow[3] = {0,4,-4};
const int DirectionFast[3] = {0,7,-7};
// const double TimeWindow = 15.0 / DeltaT;

const int FrameInterval=3;
// const int SpeedChecking= FrameInterval * 16.0 /DeltaT;
// BoundNe+ShiftNe+InhibitionNe+FMNe+BaseFrequencyNe+CoupledNe;

// % Bump,Shift,Inh,Couple,Base,FM
const double a[6]={0.04,0.04,0.04,0.02,0.02,0.02}, b[6]={0.1,0.1,0.1,0.2,0.2,0.2}, c[6]={100-39.5,100-52,100-57,100-45,100-39.5,100-57}, dFast[6]={0.1,0.1,0.1,0.1,0.1,0.1}, dSlow[6]={0.1,0.1,0.1,0.1,0.1,6}, IBFast[6]={6,-35,-2,16,10,-11}, IBSlow[6]={9.4,2,-2,16,10,-11};
void funca(double fa[TotalNe], const double v[TotalNe],const double u[TotalNe],const double b[TotalNe],const double I[TotalNe],const double dT, const double dtt, const double arg1[TotalNe], const double arg2[TotalNe]);
void funcb(double fb[TotalNe], const double a[TotalNe],const double b[TotalNe],const double v[TotalNe],const double u[TotalNe],const double dT, const double dtt, const double arg1[TotalNe], const double arg2[TotalNe]);

double MYABS(double x) {
	if(x < 0) return -x;
	else return x;
}

int main()
{
    int ModulationCurrent=DirectionSlow[1];
    //int Speed=Velocity[1];
    int debug_ = 2;
    int i, j, k;
    double Speed = 0;
    double mychart_x[29]={0}, mychart_y[29]={0};
    double tmpVal, tmpVal_1;
    double UpdateStimulus = 3000;
    // line 50 init A,B,C,D
    double A[TotalNe], B[TotalNe], C[TotalNe], DSlow[TotalNe], DFast[TotalNe], IbiasSlow[TotalNe], IbiasFast[TotalNe];
    //double A[TotalNe], B[TotalNe], C[TotalNe], DSlow[TotalNe], DFast[TotalNe], IbiasSlow[TotalNe], IbiasFast[TotalNe], ExternalI[TotalNe] = {0}, Inoise[TotalNe] = {0}, v[TotalNe], u[TotalNe], S1[TotalNe] = {0}, S2[TotalNe] = {0}, Potential[TotalNeA1] = {0}, SynapticCurrent1[TotalNeA1] = {0}, SynapticCurrent2[TotalNeA1] = {0};
    double ExternalI[TotalNe] = {0}, Inoise[TotalNe] = {0}, v[TotalNe], u[TotalNe], S1[TotalNe] = {0}, S2[TotalNe] = {0}, Potential[TotalNeA1] = {0}, SynapticCurrent1[TotalNeA1] = {0}, SynapticCurrent2[TotalNeA1] = {0};
    double *D, *Ibias;
    for(i = 0; i < BoundNe; i++) {  A[i] = a[0];    B[i] = b[0];    C[i] = c[0];    DSlow[i] = dSlow[0];    DFast[i] = dFast[0];    IbiasSlow[i] = IBSlow[0];    IbiasFast[i] = IBFast[0]; }
    for(i = BoundNe; i < ShiftNe + BoundNe; i++) {  A[i] = a[1];    B[i] = b[1];    C[i] = c[1];    DSlow[i] = dSlow[1];    DFast[i] = dFast[1];    IbiasSlow[i] = IBSlow[1];    IbiasFast[i] = IBFast[1];  }
    for(i = BoundNe + ShiftNe; i < ShiftNe + BoundNe + InhibitionNe; i++) {  A[i] = a[2];    B[i] = b[2];    C[i] = c[2];    DSlow[i] = dSlow[2];    DFast[i] = dFast[2];    IbiasSlow[i] = IBSlow[2];    IbiasFast[i] = IBFast[2];  }
    for(i = BoundNe + ShiftNe + InhibitionNe; i < ShiftNe + BoundNe + InhibitionNe + CoupledNe; i++) {  A[i] = a[3];    B[i] = b[3];    C[i] = c[3];    DSlow[i] = dSlow[3];    DFast[i] = dFast[3];    IbiasSlow[i] = IBSlow[3];    IbiasFast[i] = IBFast[3];  }
    for(i = BoundNe + ShiftNe + InhibitionNe + CoupledNe; i < ShiftNe + BoundNe + InhibitionNe + CoupledNe + BaseFrequencyNe; i++) {  A[i] = a[4];    B[i] = b[4];    C[i] = c[4];    DSlow[i] = dSlow[4];    DFast[i] = dFast[4];    IbiasSlow[i] = IBSlow[4];    IbiasFast[i] = IBFast[4];  }
    for(i = BoundNe + ShiftNe + InhibitionNe + CoupledNe + BaseFrequencyNe; i < ShiftNe + BoundNe + InhibitionNe + CoupledNe + BaseFrequencyNe + FMNe; i++) {  A[i] = a[5];    B[i] = b[5];    C[i] = c[5];    DSlow[i] = dSlow[5];    DFast[i] = dFast[5];    IbiasSlow[i] = IBSlow[5];    IbiasFast[i] = IBFast[5];  }
    D = DSlow;
    Ibias = IbiasSlow;
    for(i = 0; i < TotalNe; i++) { v[i] = Vr; u[i] = 10.0; }
    v[TotalNe - 2] = 70;

    for(i = 0; i < TotalNe; i++) { Potential[i+1] = v[i]; SynapticCurrent1[i+1] = S1[i]; SynapticCurrent2[i+1] = S2[i]; }
    
    // line 76 g1 g2
    double g1_Slow[TotalNe][TotalNe] = {0}, g1_Fast[TotalNe][TotalNe] = {0}, g2_Slow[TotalNe][TotalNe] = {0}, g2_Fast[TotalNe][TotalNe] = {0};  // Assume dimension same as TotalNe
    double (*g1)[TotalNe] = g1_Slow, (*g2)[TotalNe] = g2_Slow;
    FILE *fptr; 
    fptr = fopen("../Connection_Table_Slow.txt", "r");
    while(fscanf(fptr, "%d %d %lf", &i, &j, &tmpVal) == 3) { g1_Slow[j-1][i-1] = tmpVal; }
    fclose(fptr);
    fptr = fopen("../Connection_Table_Slow_short.txt", "r");
    while(fscanf(fptr, "%d %d %lf", &i, &j, &tmpVal) == 3) { g2_Slow[j-1][i-1] = tmpVal; }
    fclose(fptr);
    fptr = fopen("../Connection_Table_Fast.txt", "r");
    while(fscanf(fptr, "%d %d %lf", &i, &j, &tmpVal) == 3) { g1_Fast[j-1][i-1] = tmpVal; }
    fclose(fptr);
    fptr = fopen("../Connection_Table_Fast_short.txt", "r");
    while(fscanf(fptr, "%d %d %lf", &i, &j, &tmpVal) == 3) { g2_Fast[j-1][i-1] = tmpVal; }
    fclose(fptr);
    i = 0;
    fptr = fopen("../speedchart.csv", "r");
    while(fscanf(fptr, " %lf,%lf", &tmpVal, &tmpVal_1) == 2) {
        mychart_x[i] = tmpVal;
        mychart_y[i] = tmpVal_1;
        i++;
    }
    fclose(fptr);
    // VelocityLen = 1;
    int l = 1;
    // Line 84 - 92
    Speed = SaSpeed[0];
    for(i = 0; i < TotalNe; i++) { v[i] = Vr; u[i] = 10.0; S1[i] = 0; S2[i] = 0;}
    v[TotalNe - 2] = 70;
    Potential[0] = 0, SynapticCurrent1[0] = 0, SynapticCurrent2[0] = 0;
    for(i = 0; i < TotalNe; i++) { Potential[i+1] = v[i]; SynapticCurrent1[i+1] = S1[i]; SynapticCurrent2[i+1] = S2[i]; }
    // line 96
    double firings[DoubTotalNe] = {0}, FirRate[BoundNeA1] = {0}, I[TotalNe] = {0};
    int n = 1, t;
    double Position[1000][3] = {0};
    int PositionIdx = 0, x = 0;
    FILE* tempfptr = NULL;
    if (debug_ == 2) {
        tempfptr = fopen("plot.txt", "w");
        // fprintf(tempfptr,"Simulation Time: %d, DeltaT: %f\n\n", SimulationTime, DeltaT);
    }

    double SaSpeed_temp = 0;
    int m = 0;
    int t_updateTime = UpdateStimulus;
    double tmpi, tmpvi, tmpui, minchar, tmpsu;
    for(t = 1; t <= SimulationTime/DeltaT; t++) 
    // for(t = 1; t < 3638; t++) 
    {
        if((t % 50000) == 5000) {
            SaSpeed_temp = SaSpeed[m];
            t_updateTime = 0;
            StimulusNeuron = SaPosition[m];
            if(MYABS(SaSpeed_temp) <= SaBorder) {
                if(SaSpeed_temp < 0) {
                    ModulationCurrent = DirectionSlow[2];
                }
                else {
                    ModulationCurrent = DirectionSlow[1];
                }
                UpdateStimulus = 5000;
                //mode = 0;
                D = DSlow;
                Ibias = IbiasSlow;
                g1 = g1_Slow;
                g2 = g2_Slow;
            }
            else {
                if(SaSpeed_temp < 0) {
                    ModulationCurrent = DirectionFast[2];
                }
                else {
                    ModulationCurrent = DirectionFast[1];
                }
                UpdateStimulus = 3000;
                //mode = 1;
                D = DFast;
                Ibias = IbiasFast;
                g1 = g1_Fast;
                g2 = g2_Fast;
            }
            tmpi = MYABS(SaSpeed_temp);
            minchar = MYABS(mychart_y[0]-tmpi);
            tmpsu = mychart_x[0];
            for(i = 1; i < 29; i++) {
                tmpvi = MYABS(mychart_y[i]-tmpi);
                if(tmpvi < minchar) {
                    minchar = tmpvi;
                    tmpsu = mychart_x[i];
                }
            }
            Speed = tmpsu;
            if(m < SaSpeedLen) m++;
        }

        //if (t > StimulationOnset && t <= StimulationOffset) { ExternalI[StimulusNeuron-1] = StimulusStrength[0]; }
        //else { ExternalI[StimulusNeuron - 1] = 0; }
        if(t_updateTime < UpdateStimulus) {
            ExternalI[StimulusNeuron-1] = StimulusStrength;
            t_updateTime++;
        }
        else {
            ExternalI[StimulusNeuron-1] = 0;
        }

        int fired1[TotalNe] = {0}, fired2[TotalNe] = {0}, fired3[TotalNe] = {0}, fired[TotalNe] = {0};
        int fired1Num = 0, fired2Num = 0, fired3Num = 0;
        for(i = 0, j = 0; i < BoundNe + ShiftNe + InhibitionNe; i++) { if(v[i] >= Vth) {fired1[j] = i; fired3[j] = i; fired[j] = i; v[fired[j]] = C[fired[j]];  u[fired[j]] += D[fired[j]]; j++;}} 
        fired1Num = j;
        for(i = BoundNe + ShiftNe + InhibitionNe, j = 0; i < TotalNe; i++) { 
            if(v[i] >= Vth) {
                fired2[j] = i; 
                fired3[j + fired1Num] = i;
                fired[j + fired1Num] = i; 
                v[fired[j+fired1Num]] = C[fired[j+fired1Num]]; 
                u[fired[j+fired1Num]] += D[fired[j+fired1Num]]; j++;
            }
        } 
        fired2Num = j; fired3Num = fired2Num + fired1Num;
        for(i = 0; i < fired3Num; i++) { firings[i] = t * DeltaT + fired3[i]; firings[i+fired3Num] = fired3[i]; }
        
        int tem = {0};
        double temp[BoundNeA1] = {0};
        temp[0] = t * DeltaT;
        // line 138 ~ 150 , 152 153 ???
        // if (t <= TimeWindow) {
        //     // line 139
        //     for(i = 0, j = 0; i <= t; i++) { if(firings[i][0] <= t) { tem[j] = i; j++; }
        //     for(i = 0; i < BoundNe; i++) temp[i+1] = 1000 * 

        // } else {

        // }
        if (debug_ == 1) {
            char flnm[30];
            snprintf(flnm, 30, "../couts/t%d.txt", t);
            printf("%s", flnm);
            tempfptr = fopen(flnm, "w");
            fprintf(tempfptr,"before fired2 = ");
        }
        
        // for(j = 0; j < fired2Num; j++) fprintf(tempfptr,"%3.4f ", (fired2[j] + 1) * 1.0);  // + 1 to match matlab 
        // fprintf(tempfptr,"\n");
        // line 145 144
        for(i = 0; i < TotalNe; i++) {
            double tmpSum = 0;
            for(j = 0; j < fired1Num; j++) tmpSum += g1[i][fired1[j]];
            S1[i] = S1[i] + tmpSum - (S1[i]/Tau1)*DeltaT;
            
            tmpSum = 0;
            for(j = 0; j < fired2Num; j++) tmpSum += g2[i][fired2[j]]; 
            if (debug_ == 1) {fprintf(tempfptr,"%3.4f ", tmpSum); }
            S2[i] = S2[i] + tmpSum - (S2[i]/Tau2)*DeltaT;
            
        }
        if (debug_ == 1) { fprintf(tempfptr,"\n"); }

        

        //line 160
        tmpvi = Speed;
        if(MYABS(SaSpeed_temp) <= SaBorder) {
            ExternalI[23] = Speed;
            ExternalI[24] = Speed;
            tmpvi = 0;
        }
        tmpui = MYABS(ModulationCurrent) + tmpvi;
        if(ModulationCurrent > 0) {
            // 4 
            for(i = BoundNe; i < ShiftNe/2+BoundNe; i++) ExternalI[i] = tmpui;
        } else if(ModulationCurrent < 0) {
            for(i = BoundNe+ShiftNe/2; i < ShiftNe+BoundNe; i++) ExternalI[i] = tmpui;
        } else {
            for(i = BoundNe; i < ShiftNe+BoundNe; i++) ExternalI[i] = 0;
        }
                                    
        for(i = 0; i < TotalNe; i++) I[i] = ExternalI[i] + S1[i] + S2[i] + Inoise[i] + Ibias[i];

        // Line 172
        double tmpf[TotalNe] = {0};
        double fa1[TotalNe], fa2[TotalNe], fa3[TotalNe], fa4[TotalNe], fb1[TotalNe], fb2[TotalNe], fb3[TotalNe], fb4[TotalNe];
        if (debug_ == 1) {
            fprintf(tempfptr,"ExternalI = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", ExternalI[i]); }
            fprintf(tempfptr,"\n");

            fprintf(tempfptr,"before S1 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", S1[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"before S2 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", S2[i]); }
            fprintf(tempfptr,"\n");

            fprintf(tempfptr,"before u = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", u[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"before v = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", v[i]); }
            fprintf(tempfptr,"\n");
        }
        

        funca(fa1, v,u,B,I,DeltaT, 0, tmpf, tmpf);
        funcb(fb1, A,B,v,u,DeltaT, 0, tmpf, tmpf);
        
        funca(fa2, v, u, B, I, DeltaT, DeltaT*0.5, fa1, fb1);
        funcb(fb2, A, B, v, u, DeltaT, DeltaT*0.50000, fa1, fb1);

        funca(fa3, v, u, B, I, DeltaT, DeltaT*0.5, fa2, fb2);
        funcb(fb3, A, B, v, u, DeltaT, DeltaT*0.5, fa2, fb2);

        funca(fa4, v, u, B, I, DeltaT, DeltaT, fa3, fb3);
        funcb(fb4, A, B, v, u, DeltaT, DeltaT, fa3, fb3);

        for(i = 0; i < TotalNe; i++) {
            v[i] += (DeltaT / 6.0) * (fa1[i] + 2 * fa2[i] + 2 * fa3[i] + fa4[i]);
            u[i] += (DeltaT / 6.0) * (fb1[i] + 2 * fb2[i] + 2 * fb3[i] + fb4[i]);
        }
        if(debug_ == 1) {
            fprintf(tempfptr,"Ibias = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", Ibias[i]);}
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"fa1 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", fa1[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"fa2 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", fa2[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"fa3 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", fa3[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"fa4 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", fa4[i]); }
            fprintf(tempfptr,"\n");

            fprintf(tempfptr,"fb1 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", fb1[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"fb2 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", fb2[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"fb3 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", fb3[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"fb4 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", fb4[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"I = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", I[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"S1 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", S1[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"S2 = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", S2[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"u = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", u[i]); }
            fprintf(tempfptr,"\n");
            fprintf(tempfptr,"v = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", v[i]); }
            fprintf(tempfptr,"\n");
            fclose(tempfptr);
            printf("\n");
        } else if (debug_ == 2) {
            fprintf(tempfptr,"v = ");
            for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", v[i]); }
            fprintf(tempfptr,"\n");
        }
        

    }
    


// }


}

void funca(double fa[TotalNe], const double v[TotalNe],const double u[TotalNe],const double b[TotalNe],const double I[TotalNe],const double dT, const double dtt, const double arg1[TotalNe], const double arg2[TotalNe])
{ 
    double tmpv, tmpu;
    for(int ii = 0; ii < TotalNe; ii++) {
        tmpv = v[ii] + dtt * arg1[ii];
        tmpu = u[ii] + dtt * arg2[ii];
        fa[ii] = 0.04 * tmpv * tmpv - tmpv * 3 + 40 + 100 * b[ii] - tmpu + I[ii];
    }
}

void funcb(double fb[TotalNe], const double a[TotalNe],const double b[TotalNe],const double v[TotalNe],const double u[TotalNe],const double dT, const double dtt, const double arg1[TotalNe], const double arg2[TotalNe])
{
    double tmpv, tmpu;
    for(int ii = 0; ii < TotalNe; ii++){
        tmpv = v[ii] + dtt * arg1[ii];
        tmpu = u[ii] + dtt * arg2[ii];
        fb[ii] = a[ii] * ( b[ii] * tmpv - tmpu);
    } 
}
