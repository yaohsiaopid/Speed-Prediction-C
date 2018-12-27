% Created by Eugene M. Izhikevich, February 25, 2003
% Modified by Chen-Fu Yeh, 2018
% Excitatory neurons    Inhibitory neurons
% TODO:
% Use firing rate detection to output the bump position.

%{
ver distcomp        % show the version of Parallel Computing Toolbox
parpool('local',8)  % parallel with 4 cores. (i7-7700 4C8T -> 4)
                    % please adjust to the actual processor core number.
%}

%tic
clc
clear
SimulationTime=2000;        % ms
DeltaT=0.01;                % ms
CameraFps=60;               % Hz of motion_estimation camera
DeltaC=1000/CameraFps;      % ms
Vr=10;                      % membrane potential initial value
Vth=130;                    % membrane potential threshold value
NoiseStrengthBase=0;        % set nonzero value to add noises
Velocity=0:0.5:8;           % target Speed current(nA).
SaPosition=[1,1,1,1];       % saliency position, update every 10 frames(167ms)
%SaPosition=5;
SaSpeed=[1,1.5,2.5,4];      % saliency speed, update every 200ms
%SaSpeed=-4;
SaBorder=1.7308337;         % border to switching between slow/fast

% set parameters to produce first bump
%StimuluStrength=5.5;
StimuluStrength=15.5;
StimulusNeuron=1;
StimulationOnset=100;               % start at 100 ms
StimulationOffset=130;
StimulationOnset=StimulationOnset/DeltaT;
StimulationOffset=StimulationOffset/DeltaT;
DirectionSlow=[0,4,-4];             % direction=[no,right,left]
DirectionFast=[0,7,-7];             % same as above. Fast circuit needs higher current
ModulationCurrent=DirectionSlow(2); % set direction
%ModulationCurrent=DirectionFast(2); % set direction

% Neuron index
BoundNe=1:8;        % main bump neurons
RightShiftNe=9:15;
LeftShiftNe=16:22;
ShiftNe=9:22;
InhibitionNe=23;
CoupledNe=24:25;
BaseFrequencyNe=26;
FMNe=27;            % frequency modulation neuron
TotalNe=27;

%{
a=[0.04,0.04,0.04,0.02,0.02,0.02];
b=[0.1,0.1,0.1,0.2,0.2,0.2];
c=[-39.5,-52,-57,-45,-39.5,-60];
d=[0.1,0.1,0.1,0.1,0.1,0.1];
IB=[9.5,2,-2,16,10,-11];
c=c+100;ExternalI
%}

% set parameters for Izhikevich model equation
% Bump,Shift,Inh,Couple,Base,FM
a=[0.04,0.04,0.04,0.02,0.02,0.02];
b=[0.1,0.1,0.1,0.2,0.2,0.2];
c=[-39.5,-52,-57,-45,-39.5,-57];
dFast=[0.1,0.1,0.1,0.1,0.1,0.1];
%dSlow=[0.1,0.1,0.1,0.1,0.1,10];
dSlow=[0.1,0.1,0.1,0.1,0.1,6];
IBslow=[9.4,2,-2,16,10,-11];
IBfast=[6,-35,-2,16,10,-11];
c=c+100;                            % change offset to make all voltage positive

A=ones(TotalNe,1);
A(BoundNe)=a(1);
A(ShiftNe)=a(2);
A(InhibitionNe)=a(3);
A(CoupledNe)=a(4);
A(BaseFrequencyNe)=a(5);
A(FMNe)=a(6);

B=ones(TotalNe,1);
B(BoundNe)=b(1);
B(ShiftNe)=b(2);
B(InhibitionNe)=b(3);
B(CoupledNe)=b(4);
B(BaseFrequencyNe)=b(5);
B(FMNe)=b(6);

C=ones(TotalNe,1);
C(BoundNe)=c(1);
C(ShiftNe)=c(2);
C(InhibitionNe)=c(3);
C(CoupledNe)=c(4);
C(BaseFrequencyNe)=c(5);
C(FMNe)=c(6);

D_Slow=ones(TotalNe,1);
D_Slow(BoundNe)=dSlow(1);
D_Slow(ShiftNe)=dSlow(2);
D_Slow(InhibitionNe)=dSlow(3);
D_Slow(CoupledNe)=dSlow(4);
D_Slow(BaseFrequencyNe)=dSlow(5);
D_Slow(FMNe)=dSlow(6);
D_Fast=ones(TotalNe,1);
D_Fast(BoundNe)=dFast(1);
D_Fast(ShiftNe)=dFast(2);
D_Fast(InhibitionNe)=dFast(3);
D_Fast(CoupledNe)=dFast(4);
D_Fast(BaseFrequencyNe)=dFast(5);
D_Fast(FMNe)=dFast(6);
D=D_Slow;
%D=D_Fast;

Ibias_Slow=ones(TotalNe,1);
Ibias_Slow(BoundNe)=IBslow(1);
Ibias_Slow(ShiftNe)=IBslow(2);
Ibias_Slow(InhibitionNe)=IBslow(3);
Ibias_Slow(CoupledNe)=IBslow(4);
Ibias_Slow(BaseFrequencyNe)=IBslow(5);
Ibias_Slow(FMNe)=IBslow(6);
Ibias_Fast=ones(TotalNe,1);
Ibias_Fast(BoundNe)=IBfast(1);
Ibias_Fast(ShiftNe)=IBfast(2);
Ibias_Fast(InhibitionNe)=IBfast(3);
Ibias_Fast(CoupledNe)=IBfast(4);
Ibias_Fast(BaseFrequencyNe)=IBfast(5);
Ibias_Fast(FMNe)=IBfast(6);
Ibias=Ibias_Slow;
%Ibias=Ibias_Fast;

Inoise=NoiseStrengthBase*normrnd(0,1,[TotalNe,1]);  % add some noises
ExternalI=0*ones(TotalNe,1);

v=Vr.*ones(TotalNe,1);
v(TotalNe-1)=70;
u=10*ones(TotalNe,1);
S1=0*ones(TotalNe,1);
S2=0*ones(TotalNe,1);
Tau1=10;
Tau2=5;
SynapticCurrent1=transpose([0;S1]);
SynapticCurrent2=transpose([0;S2]);

% read neuron weights
g1_Slow=transpose(full(spconvert(dlmread('Connection_Table_Slow.txt'))));
g2_Slow=transpose(full(spconvert(dlmread('Connection_Table_Slow_short.txt'))));
g1_Fast=transpose(full(spconvert(dlmread('Connection_Table_Fast.txt'))));
g2_Fast=transpose(full(spconvert(dlmread('Connection_Table_Fast_short.txt'))));
g1=g1_Slow;
g2=g2_Slow;
%g1=g1_Fast;
%g2=g2_Fast;

% read speed chart
chart = csvread('speedchart.csv');

% simulating for each speeds
for l=1
%for l=1:length(Velocity)
%parfor l=1:length(Velocity)

    %Speed=Velocity(l);
    Speed=0;
    
    % do some initializations
    ExternalI=0*ones(TotalNe,1);
    v=Vr.*ones(TotalNe,1);
    v(BaseFrequencyNe)=70;
    u=10*ones(TotalNe,1);
    S1=0*ones(TotalNe,1);
    S2=0*ones(TotalNe,1);
    Tau1=10;
    Tau2=5;
    Potential=zeros(SimulationTime/DeltaT+1,TotalNe+1);
    Potential(1,:)=transpose([0;v]);
    SynapticCurrent1=transpose([0;S1]);
    SynapticCurrent2=transpose([0;S2]);
    I=0*ones(TotalNe,1);
    TotalCurrent=zeros(SimulationTime/DeltaT+1,TotalNe+1);
    TotalCurrent(1,:)=transpose([0;I]);
    SaSpeed_temp=0;
    UpdateStimulus=StimulationOffset-StimulationOnset;
    t_updateTime = UpdateStimulus;
    %{
    BumpIsAt=zeros(SimulationTime/DeltaT+1,2);
    BumpIsAt_temp=0;
    Prediction=zeros(round((SimulationTime-StimulationOnset)/DeltaC),2);
    CameraI=1;
    %}

    % main simulation loop
    tic
    m=1;
    for t=1:SimulationTime/DeltaT
        %disp(t);
        %Speed=SaSpeed(m);
        
        
        % update saliency info
        %if mod(mod(t,500/DeltaT),167/DeltaT) == 50/DeltaT
        if mod(t,500/DeltaT) == 50/DeltaT
            SaSpeed_temp = SaSpeed(m);
            t_updateTime = 0;
            if abs(SaSpeed_temp) <= SaBorder    % slow circuit
                if SaSpeed_temp < 0
                    ModulationCurrent=DirectionSlow(3); % left
                    StimulusNeuron=SaPosition(m);
                else
                    ModulationCurrent=DirectionSlow(2); % right
                    StimulusNeuron=SaPosition(m);
                end
                UpdateStimulus=50/DeltaT;
                D=D_Slow;
                Ibias=Ibias_Slow;
                g1=g1_Slow;
                g2=g2_Slow;
            else                                % fast circuit
                if SaSpeed_temp < 0
                    ModulationCurrent=DirectionFast(3); % left
                    StimulusNeuron=SaPosition(m);
                else
                    ModulationCurrent=DirectionFast(2); % right
                    StimulusNeuron=SaPosition(m);
                end
                UpdateStimulus=30/DeltaT;
                D=D_Fast;
                Ibias=Ibias_Fast;
                g1=g1_Fast;
                g2=g2_Fast;
            end
            
            % set speed current
            [err, index]=min(abs(chart(:,2)-abs(SaSpeed_temp)));
            Speed=chart(index,1);
            if m < length(SaSpeed)
                m = m + 1;
            end
        end
        

        Inoise=NoiseStrengthBase*normrnd(0,1,[TotalNe,1]);
        
        %{
        % start the bump
        if t>StimulationOnset && t<=StimulationOffset
            ExternalI(StimulusNeuron)=StimuluStrength;
        else
            ExternalI(StimulusNeuron)=0;
        end
        %}
        
        
        % update(overwrite) the bump
        if t_updateTime < UpdateStimulus
            ExternalI(StimulusNeuron) = StimuluStrength;
            t_updateTime = t_updateTime + 1;
        else
            ExternalI(StimulusNeuron) = 0;
        end
        
        
        fired1=find(v(1:InhibitionNe)>=Vth);
        fired2=find(v(InhibitionNe+1:end)>=Vth);
        fired2=fired2+InhibitionNe;
        fired3=find(v(BoundNe)>=Vth);
        fired=[fired1;fired2];
        v(fired)=C(fired);
        u(fired)=u(fired)+D(fired);

        %{
        % output where the bump is, using camera time stamp (currently 60fps)
        BumpIsAt(t,1)=t-1;
        if fired3
            BumpIsAt_temp=fired3(1);
        end
        %disp(BumpIsAt_temp);
        BumpIsAt(t,2)=BumpIsAt_temp;
        if (mod(t*DeltaT,DeltaC) < 0.9) && (mod(t*DeltaT,1) == 0) && (t > StimulationOnset)
            Prediction(CameraI,1)=CameraI;
            Prediction(CameraI,2)=BumpIsAt_temp;
            CameraI = CameraI+1;
            disp(t);
        end
        %}

        S1=S1+sum(1*g1(:,fired1),2)-(S1/Tau1)*DeltaT;
        S2=S2+sum(1*g2(:,fired2),2)-(S2/Tau2)*DeltaT;
        
        
        % set shift neuron bias current & coupled neuron firing rate
        ExternalI(ShiftNe)=0;
        if ModulationCurrent>0                  % right
            if abs(SaSpeed_temp) > SaBorder     % fast
                ExternalI(RightShiftNe)=ModulationCurrent+Speed;
            else                                % slow
                ExternalI(RightShiftNe)=ModulationCurrent;
                ExternalI(CoupledNe)=Speed;
            end
        elseif ModulationCurrent<0              % left
            if abs(SaSpeed_temp) > SaBorder     % fast
                ExternalI(LeftShiftNe)=abs(ModulationCurrent)+Speed;
            else                                % slow
                ExternalI(CoupledNe)=Speed;
                ExternalI(LeftShiftNe)=abs(ModulationCurrent);
            end
        else
	        ExternalI(ShiftNe)=0;
        end
        
        
        %{
        if ModulationCurrent>0
	        ExternalI(RightShiftNe)=ModulationCurrent;
            %ExternalI(RightShiftNe)=ModulationCurrent+Speed;
	        ExternalI(CoupledNe)=Speed;
        elseif ModulationCurrent<0
	        ExternalI(LeftShiftNe)=abs(ModulationCurrent);
            %ExternalI(LeftShiftNe)=ModulationCurrent+Speed;
	        ExternalI(CoupledNe)=Speed;
        else
	        ExternalI(ShiftNe)=0;
        end
        %}
        
        I=ExternalI+S1+S2+Inoise+Ibias;
  
        % Izhikevich model implemented by Runge-Kutta 4 method
        fa1=funca(v,u,B,I,DeltaT);
        fb1=funcb(A,B,v,u,DeltaT);
        fa2=funca(v+DeltaT*0.5*fa1,u+DeltaT*0.5*fb1,B,I,DeltaT);
        fb2=funcb(A,B,v+DeltaT*0.5*fa1,u+DeltaT*0.5*fb1,DeltaT);
        fa3=funca(v+DeltaT*0.5*fa2,u+DeltaT*0.5*fb2,B,I,DeltaT);
        fb3=funcb(A,B,v+DeltaT*0.5*fa2,u+DeltaT*0.5*fb2,DeltaT);
        fa4=funca(v+DeltaT*fa3,u+DeltaT*fb3,B,I,DeltaT);
        fb4=funcb(A,B,v+DeltaT*fa3,u+DeltaT*fb3,DeltaT);
        v= v + (DeltaT/6)*(fa1 + 2*fa2 + 2*fa3 + fa4);
        u= u + (DeltaT/6)*(fb1 + 2*fb2 + 2*fb3 + fb4);
  
        temp=transpose([t*DeltaT;v]);
        Potential(t+1,:)=temp;
        temp=transpose([t*DeltaT;I]);
        TotalCurrent(t+1,:)=temp;
                
    end
    toc
    
    %foldername=num2str(Velocity(l));
    %foldername=num2str(SaSpeed(l));
    foldername=num2str(l);
    mkdir (foldername);

    
    % save some datas
    for N=1:TotalNe

        % save neuron potential
        fig=plot(Potential(:,1),Potential(:,N+1),'b-');
        saveas(fig,strcat(num2str(l),'NeuronV_',num2str(N),'.png'));
        copyfile(strcat(num2str(l),'NeuronV_',num2str(N),'.png'),foldername);
        %{
        m=matfile(sprintf('%dNeuronV_%d.mat',l,N),'writable',true);
        m.x=Potential(:,1);
        m.y=Potential(:,N+1);
        copyfile(strcat(num2str(l),'NeuronV_',num2str(N),'.mat'),foldername);
        %}
        
        % save neuron input current
        fig=plot(TotalCurrent(:,1),TotalCurrent(:,N+1),'b-');
        saveas(fig,strcat(num2str(l),'CurrentV_',num2str(N),'.png'));
        copyfile(strcat(num2str(l),'CurrentV_',num2str(N),'.png'),foldername);
        n=matfile(sprintf('%dCurrentV_%d.mat',l,N),'writable',true);
        n.x=TotalCurrent(:,1);
        n.y=TotalCurrent(:,N+1);
        copyfile(strcat(num2str(l),'CurrentV_',num2str(N),'.mat'),foldername);
        
        % save neuron triggering histories
        formatSpec = '%dNeuronV_%d.txt';
        filename = sprintf(formatSpec,l,N);
        fileID = fopen(filename,'w');
        threshold_record(Potential(:,N+1),DeltaT,42.5,fileID);
        fclose(fileID);
        copyfile(filename,foldername);
                        
    end
    
    
    %{
    % save where the bump is at
    m=matfile(sprintf('%dBump.mat',l),'writable',true);
    m.x=BumpIsAt(:,1);
    m.y=BumpIsAt(:,2);
    copyfile(strcat(num2str(l),'Bump.mat'),foldername);
    
    % save prediction result
    formatSpec = '%dPrediction.txt';
    filename = sprintf(formatSpec,l);
    fileID = fopen(filename,'w');
    fprintf(fileID,'%d %d\n',transpose(Prediction));
    copyfile(filename,foldername);
    %}

end

% clear unnecessary files
delete *NeuronV_*.*
delete *CurrentV_*.*
delete *Prediction.txt
delete *Bump.mat

%toc
%delete(gcp('nocreate'))