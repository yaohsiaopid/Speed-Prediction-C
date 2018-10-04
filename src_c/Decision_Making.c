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
const int SimulationTime = 1000; // ms
const double DeltaT = 0.0100000000000; // ms
const int NoiseStrengthBase = 0;
int Velocity[2] = {1, 8}; // Length subject to change
int VelocityLen = 2;
const double StimulusStrength[1] = {5.5};
const int StimulusNeuron=1;
const double StimulationOnset= 0; 
const double StimulationOffset= 5000;
const int Direction[3] = {0,4,-4};
// const double TimeWindow = 15.0 / DeltaT;

const int FrameInterval=3;
// const int SpeedChecking= FrameInterval * 16.0 /DeltaT;
// BoundNe+ShiftNe+InhibitionNe+FMNe+BaseFrequencyNe+CoupledNe;

// % Bump,Shift,Inh,Couple,Base,FM
double a[6]={0.04,0.04,0.04,0.02,0.02,0.02}, b[6]={0.1,0.1,0.1,0.2,0.2,0.2}, c[6]={100-39.5,100-52,100-57,100-45,100-39.5,100-57}, d[6]={0.1,0.1,0.1,0.1,0.1,0.1},IB[6]={9.4,2,-2,16,10,-11};
void funca(double fa[TotalNe], const double v[TotalNe],const double u[TotalNe],const double b[TotalNe],const double I[TotalNe],const double dT, const double dtt, const double arg1[TotalNe], const double arg2[TotalNe]);
void funcb(double fb[TotalNe], const double a[TotalNe],const double b[TotalNe],const double v[TotalNe],const double u[TotalNe],const double dT, const double dtt, const double arg1[TotalNe], const double arg2[TotalNe]);
int main()
{
    const int ModulationCurrent=Direction[1];
    int Speed=Velocity[1];
    
    int i, j, k;
    double tmpVal;
    // line 50 init A,B,C,D
    double A[TotalNe], B[TotalNe], C[TotalNe], D[TotalNe], Ibias[TotalNe], ExternalI[TotalNe] = {0}, Inoise[TotalNe] = {0}, v[TotalNe], u[TotalNe], S1[TotalNe] = {0}, S2[TotalNe] = {0}, Potential[TotalNeA1] = {0}, SynapticCurrent1[TotalNeA1] = {0}, SynapticCurrent2[TotalNeA1] = {0};
    for(i = 0; i < BoundNe; i++) {  A[i] = a[0];    B[i] = b[0];    C[i] = c[0];    D[i] = d[0];    Ibias[i] = IB[0];  }
    for(i = BoundNe; i < ShiftNe + BoundNe; i++) {  A[i] = a[1];    B[i] = b[1];    C[i] = c[1];    D[i] = d[1];    Ibias[i] = IB[1];  }
    for(i = BoundNe + ShiftNe; i < ShiftNe + BoundNe + InhibitionNe; i++) {  A[i] = a[2];    B[i] = b[2];    C[i] = c[2];    D[i] = d[2];    Ibias[i] = IB[2];  }
    for(i = BoundNe + ShiftNe + InhibitionNe; i < ShiftNe + BoundNe + InhibitionNe + CoupledNe; i++) {  A[i] = a[3];    B[i] = b[3];    C[i] = c[3];    D[i] = d[3];    Ibias[i] = IB[3];  }
    for(i = BoundNe + ShiftNe + InhibitionNe + CoupledNe; i < ShiftNe + BoundNe + InhibitionNe + CoupledNe + BaseFrequencyNe; i++) {  A[i] = a[4];    B[i] = b[4];    C[i] = c[4];    D[i] = d[4];    Ibias[i] = IB[4];  }
    for(i = BoundNe + ShiftNe + InhibitionNe + CoupledNe + BaseFrequencyNe; i < ShiftNe + BoundNe + InhibitionNe + CoupledNe + BaseFrequencyNe + FMNe; i++) {  A[i] = a[5];    B[i] = b[5];    C[i] = c[5];    D[i] = d[5];    Ibias[i] = IB[5];  }
    
    for(i = 0; i < TotalNe; i++) { v[i] = Vr; u[i] = 10.0; }
    v[TotalNe - 2] = 70;

    for(i = 0; i < TotalNe; i++) { Potential[i+1] = v[i]; SynapticCurrent1[i+1] = S1[i]; SynapticCurrent2[i+1] = S2[i]; }
    
    // line 76 g1 g2
    double g1[TotalNe][TotalNe] = {0} , g2[TotalNe][TotalNe] = {0}; // Assume dimension same as TotalNe
    FILE *fptr; 
    fptr = fopen("../Connection_Table_temp.txt", "r");
    while(fscanf(fptr, "%d %d %lf", &i, &j, &tmpVal) == 3) {
        // printf("tmp %f\n", tmpVal);
        g1[j-1][i-1] = tmpVal; 
    }
    // for(i = 0; i < TotalNe; i++) {
    //     for(j = 0; j < TotalNe; j++) {
    //         printf("%5.5f   ", g1[i][j] );
    //     }
    //     printf("\n");
    // }
    // printf("------g1-----\n");
    fclose(fptr);

    fptr = fopen("../Connection_Table_temp_short.txt", "r");
    while(fscanf(fptr, "%d %d %lf", &i, &j, &tmpVal) == 3) {
        g2[j-1][i-1] = tmpVal; 
    }
    // for(i = 0; i < TotalNe; i++) {
    //     for(j = 0; j < TotalNe; j++) {
    //         printf("%5.5f   ", g2[i][j] );
    //     }
    //     printf("\n");
    // }
    // printf("------g2-----\n");
    fclose(fptr);
    // VelocityLen = 1;
    int l = 1;

    // Line 84 - 92
    Speed = Velocity[l];
    for(i = 0; i < TotalNe; i++) { v[i] = Vr; u[i] = 10.0; S1[i] = 0; S2[i] = 0;}
    v[TotalNe - 2] = 70;
    Potential[0] = 0, SynapticCurrent1[0] = 0, SynapticCurrent2[0] = 0;
    for(i = 0; i < TotalNe; i++) { Potential[i+1] = v[i]; SynapticCurrent1[i+1] = S1[i]; SynapticCurrent2[i+1] = S2[i]; }
    // line 96
    double firings[DoubTotalNe] = {0}, FirRate[BoundNeA1] = {0}, I[TotalNe] = {0};
    int n = 1, t;
    double Position[1000][3] = {0};
    int PositionIdx = 0, x = 0;
    // for(t = 1; t < SimulationTime / DeltaT; t++) 
    for(t = 1; t < 200; t++) 
    {
        // printf("%d,", t);
        if (t > StimulationOnset && t <= StimulationOffset) { ExternalI[StimulusNeuron-1] = StimulusStrength[0]; }
        else { ExternalI[StimulusNeuron - 1] = 0; }

        // lien 114
        // if (t % SpeedChecking == 0) {
        //     x = 0;
        //     for(int m = 1; m < 9; m++) {
        //         if(FirRate[t-1][m] > 950) { x = m; Position[PositionIdx][0] = t; Position[PositionIdx][1] = (t/16)*DeltaT; Position[PositionIdx][2] = x;  PositionIdx++; }
        //     }
        //     if(x == 8)
        //         break;
        // }

        int fired1[TotalNe] = {0}, fired2[TotalNe] = {0}, fired3[TotalNe] = {0}, fired[TotalNe] = {0};
        int fired1Num = 0, fired2Num = 0, fired3Num = 0;
        for(i = 0, j = 0; i < BoundNe + ShiftNe + InhibitionNe; i++) { if(v[i] >= Vth) {fired1[j] = i; fired3[j] = i; fired[j] = i; v[fired[j]] = C[fired[j]];  u[fired[j]] += D[fired[j]]; j++;}} 
        fired1Num = j;
        for(i = BoundNe + ShiftNe + InhibitionNe, j = 0; i < TotalNe; i++) { if(v[i] >= Vth) {fired2[j] = i; fired3[j + fired1Num] = i; fired[j + fired1Num] = i; v[fired[j+fired1Num]] = C[fired[j+fired1Num]]; u[fired[j+fired1Num]] += D[fired[j+fired1Num]]; j++;}} 
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

        // line 145 144
        for(i = 0; i < TotalNe; i++) {
            int tmpSum = 0;
            for(j = 0; j < fired1Num; j++) tmpSum += g1[i][fired1[j]];
            S1[i] = S1[i] + tmpSum - (S1[i]/Tau1)*DeltaT;
            
            tmpSum = 0;
            for(j = 0; j < fired2Num; j++) tmpSum += g2[i][fired2[j]];
            S2[i] = S2[i] + tmpSum - (S2[i]/Tau2)*DeltaT;
            
        }
        
        //line 160
        if(ModulationCurrent > 0) {
            // 4 
            for(i = BoundNe; i < ShiftNe/2+BoundNe; i++) ExternalI[i] = ModulationCurrent;
            for(i = ShiftNe+BoundNe+InhibitionNe; i < ShiftNe+BoundNe+1+CoupledNe; i++) ExternalI[i] = Speed;
            
        } else if(ModulationCurrent < 0) {
            for(i = ShiftNe; i < ShiftNe+BoundNe; i++) ExternalI[i] = -ModulationCurrent;
            for(i = ShiftNe+BoundNe+InhibitionNe; i < ShiftNe+BoundNe+1+CoupledNe; i++) ExternalI[i] = Speed;
        } else {
            for(i = BoundNe; i < ShiftNe+BoundNe; i++) ExternalI[i] = 0;
        }
                                    
        for(i = 0; i < TotalNe; i++) I[i] = ExternalI[i] + S1[i] + S2[i] + Inoise[i] + Ibias[i];

        // Line 172
        double tmpf[TotalNe] = {0};
        double fa1[TotalNe], fa2[TotalNe], fa3[TotalNe], fa4[TotalNe], fb1[TotalNe], fb2[TotalNe], fb3[TotalNe], fb4[TotalNe];
        char flnm[30];
        snprintf(flnm, 30, "../couts/t%d.txt", t);
        printf("%s", flnm);
        FILE* tempfptr = fopen(flnm, "w");
        fprintf(tempfptr,"ExternalI = ");
        for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", ExternalI[i]); }
        fprintf(tempfptr,"\n");
        fprintf(tempfptr,"before u = ");
        for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", u[i]); }
        fprintf(tempfptr,"\n");
        fprintf(tempfptr,"before v = ");
        for(i = 0; i < TotalNe; i++) { fprintf(tempfptr,"%3.4f ", v[i]); }
        fprintf(tempfptr,"\n");

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