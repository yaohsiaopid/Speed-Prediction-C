% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons    Inhibitory neurons

ver distcomp        % show the version of Parallel Computing Toolbox
parpool('local',4)  % parallel with 4 cores.
                    % adjust to the actual processor core number.

clc
clear
tic
SimulationTime=2000;    % ms
DeltaT=0.01;            % ms
Vr=10;
Vth=130;
NoiseStrengthBase=0;
%Velocity=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7];
Velocity=[0:0.5:9.5];   % Speed current(nA) of interest
StimuluStrength=5.5;
StimulusNeuron=1;
StimulationOnset=300;   % start at 300ms
StimulationOffset=305;
StimulationOnset=StimulationOnset/DeltaT;
StimulationOffset=StimulationOffset/DeltaT;
Direction=[0,4,-4.5];
ModulationCurrent=Direction(2);
%FMCurrent=14.6:0.1:15.7;
%InhibitionCurrent=0;
ShiftCurrent=4.5:0.1:4.8;
%Speed=1.5;

% Neuron index
BoundNe=1:8;        % main bump neurons
RightShiftNe=9:15;
LeftShiftNe=16:22;
ShiftNe=9:22;
InhibitionNe=23;
CoupledNe=24:25;
BaseFrequencyNe=26;
FMNe=27;            % frequency modulation neuron
TotalNe=FMNe;

%{
a=[0.04,0.04,0.04,0.02,0.02,0.02];
b=[0.1,0.1,0.1,0.2,0.2,0.2];
c=[-39.5,-52,-57,-45,-39.5,-60];
d=[0.1,0.1,0.1,0.1,0.1,0.1];
IB=[9.5,2,-2,16,10,-11];
c=c+100;
%}

% Bump,Shift,Inh,Couple,Base,FM
a=[0.04,0.04,0.04,0.02,0.02,0.02];
b=[0.1,0.1,0.1,0.2,0.2,0.2];
c=[-39.5,-52,-57,-45,-39.5,-57];
d=[0.1,0.1,0.1,0.1,0.1,0.1];
IB=[9.4,2,-2,16,10,-11];
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

D=ones(TotalNe,1);
D(BoundNe)=d(1);
D(ShiftNe)=d(2);
D(InhibitionNe)=d(3);
D(CoupledNe)=d(4);
D(BaseFrequencyNe)=d(5);
D(FMNe)=d(6);

Ibias=ones(TotalNe,1);
Ibias(BoundNe)=IB(1);
Ibias(ShiftNe)=IB(2);
Ibias(InhibitionNe)=IB(3);
Ibias(CoupledNe)=IB(4);
Ibias(BaseFrequencyNe)=IB(5);
Ibias(FMNe)=IB(6);

%{
A=[a(1)*ones(BoundNe,1);a(2)*ones(ShiftNe,1);a(3)*ones(InhibitionNe,1);a(4)*ones(CoupledNe,1);a(5)*ones(BaseFrequencyNe,1);a(6)*ones(FMNe,1)];
B=[b(1)*ones(BoundNe,1);b(2)*ones(ShiftNe,1);b(3)*ones(InhibitionNe,1);b(4)*ones(CoupledNe,1);b(5)*ones(BaseFrequencyNe,1);b(6)*ones(FMNe,1)];
C=[c(1)*ones(BoundNe,1);c(2)*ones(ShiftNe,1);c(3)*ones(InhibitionNe,1);c(4)*ones(CoupledNe,1);c(5)*ones(BaseFrequencyNe,1);c(6)*ones(FMNe,1)];
D=[d(1)*ones(BoundNe,1);d(2)*ones(ShiftNe,1);d(3)*ones(InhibitionNe,1);d(4)*ones(CoupledNe,1);d(5)*ones(BaseFrequencyNe,1);d(6)*ones(FMNe,1)];
Ibias=[IB(1)*ones(BoundNe,1);IB(2)*ones(ShiftNe,1);IB(3)*ones(InhibitionNe,1);IB(4)*ones(CoupledNe,1);IB(5)*ones(BaseFrequencyNe,1);IB(6)*ones(FMNe,1)];
%}

Inoise=NoiseStrengthBase*normrnd(0,1,[TotalNe,1]);
ExternalI=0*ones(TotalNe,1);

v=Vr.*ones(TotalNe,1);
v(TotalNe-1)=70;
u=10*ones(TotalNe,1);
S1=0*ones(TotalNe,1);
S2=0*ones(TotalNe,1);
Tau1=10;
Tau2=5;

%Potential=zeros(SimulationTime/DeltaT+1,TotalNe+1);
%Potential(1,:)=transpose([0;v]);
%Potential=transpose([0;v]);
SynapticCurrent1=transpose([0;S1]);
SynapticCurrent2=transpose([0;S2]);
g1=transpose(full(spconvert(dlmread('Connection_Table_temp.txt'))));
g2=transpose(full(spconvert(dlmread('Connection_Table_temp_short.txt'))));

%for l=1:length(ShiftCurrent)
parfor l=1:length(ShiftCurrent) % parallel for
                                % change back to "for" loop if encounter problems
    %Speed=Velocity(l);
    Speed=Velocity(5);
    ExternalI=0*ones(TotalNe,1);
    v=Vr.*ones(TotalNe,1);
    v(TotalNe-1)=70;
    u=10*ones(TotalNe,1);
    S1=0*ones(TotalNe,1);
    S2=0*ones(TotalNe,1);
    Tau1=10;
    Tau2=5;
    Potential=zeros(SimulationTime/DeltaT+1,TotalNe+1);
    Potential(1,:)=transpose([0;v]);
    %Potential=transpose([0;v]);
    SynapticCurrent1=transpose([0;S1]);
    SynapticCurrent2=transpose([0;S2]);
    I=0*ones(TotalNe,1);
    TotalCurrent=zeros(SimulationTime/DeltaT+1,TotalNe+1);
    TotalCurrent(1,:)=transpose([0;I]);
    %TotalCurrent= transpose([0;I]);

    for t=1:SimulationTime/DeltaT   % simulation time ms
        disp(ShiftCurrent(l));
        disp(t);
        Inoise=NoiseStrengthBase*normrnd(0,1,[TotalNe,1]);
 
        if t>StimulationOnset && t<=StimulationOffset
            ExternalI(StimulusNeuron)=StimuluStrength;
        else
            ExternalI(StimulusNeuron)=0;
        end
        
        fired1=find(v(1:InhibitionNe)>=Vth);
        fired2=find(v(InhibitionNe+1:end)>=Vth);
        fired2=fired2+InhibitionNe;
        fired=[fired1;fired2];
        v(fired)=C(fired);
        u(fired)=u(fired)+D(fired);

        S1=S1+sum(1*g1(:,fired1),2)-(S1/Tau1)*DeltaT;
        S2=S2+sum(1*g2(:,fired2),2)-(S2/Tau2)*DeltaT;
        if ModulationCurrent>0
	        ExternalI(RightShiftNe)=ModulationCurrent;
	        ExternalI(CoupledNe)=Speed;
        elseif ModulationCurrent<0
	        ExternalI(LeftShiftNe)=abs(ModulationCurrent);
	        ExternalI(CoupledNe)=Speed;
        else
	        ExternalI(ShiftNe)=0;
        end
        
        ExternalI(RightShiftNe)=ShiftCurrent(l);
        %ExternalI(InhibitionNe)=InhibitionCurrent(l);
 
        I=ExternalI+S1+S2+Inoise+Ibias;
  
        % Izhikevich model implemented by Runge-Kutta methods
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
  
        %v=v+funca(v,u,B,I,DeltaT)*DeltaT;
        %u=u+funcb(A,B,v,u,DeltaT)*DeltaT;
        %xn+1 = xn + h 
        %yn+1 = yn + (h/2) (f(xn, yn) + f(xn + h, yn +  h f(xn, yn)))

        %{
        VK1=funca(v,u,B,I,DeltaT);
        UK1=funcb(A,B,v,u,DeltaT);
        VK2=funca(v+VK1*DeltaT,u+UK1*DeltaT,B,I,DeltaT);
        UK2=funcb(A,B,v+VK1*DeltaT,u+UK1*DeltaT,DeltaT);
  
        v=v+(DeltaT/2)*(VK1+VK2);
        u=u+(DeltaT/2)*(UK1+UK2);
        %}
  
        temp=transpose([t*DeltaT;v]);
        Potential(t+1,:)=temp;
        %Potential=cat(1,Potential,temp);
        %temp=transpose([t*DeltaT;S1]);
        %SynapticCurrent1=cat(1,SynapticCurrent1,temp);
        %temp=transpose([t*DeltaT;S2]);
        %SynapticCurrent2=cat(1,SynapticCurrent2,temp);
        temp=transpose([t*DeltaT;I]);
        TotalCurrent(t+1,:)=temp;
        %TotalCurrent=cat(1,TotalCurrent,temp);
                
    end
    foldername=num2str(ShiftCurrent(l));
    mkdir (foldername);

    % save Potential pictures, datas, and triggering records
    % save SynapticCurrent2 pictures and datas
    for N=1:TotalNe     % save Potential pictures, datas, and triggering records
        fig=plot(Potential(:,1),Potential(:,N+1),'b-');
        saveas(fig,strcat(num2str(l),'NeuronV_',num2str(N),'.png'));
        copyfile(strcat(num2str(l),'NeuronV_',num2str(N),'.png'),foldername);
        m=matfile(sprintf('%dNeuronV_%d.mat',l,N),'writable',true);
        m.x=Potential(:,1);
        m.y=Potential(:,N+1);
        copyfile(strcat(num2str(l),'NeuronV_',num2str(N),'.mat'),foldername);
        
        
        fig=plot(TotalCurrent(:,1),TotalCurrent(:,N+1),'b-');
        saveas(fig,strcat(num2str(l),'CurrentV_',num2str(N),'.png'));
        copyfile(strcat(num2str(l),'CurrentV_',num2str(N),'.png'),foldername);
        n=matfile(sprintf('%dCurrentV_%d.mat',l,N),'writable',true);
        n.x=TotalCurrent(:,1);
        n.y=TotalCurrent(:,N+1);
        copyfile(strcat(num2str(l),'CurrentV_',num2str(N),'.mat'),foldername);
        
        
        %fig=plot(SynapticCurrent1(:,1),SynapticCurrent1(:,N+1),'b-');
        %saveas(fig,strcat('NeuronC1','_',num2str(N),'.png'));
        %fig=plot(SynapticCurrent2(:,1),SynapticCurrent2(:,N+1),'b-');
        %saveas(fig,strcat('NeuronC2','_',num2str(N),'.png'));
        formatSpec = '%dNeuronV_%d.txt';
        filename = sprintf(formatSpec,l,N);
        fileID = fopen(filename,'w');
        threshold_record(Potential(:,N+1),DeltaT,40,fileID);    % use 40mV as threshold voltage
        fclose(fileID);
        copyfile(filename,foldername);
        
        
        formatSpec = '%dCurrentV_%d.txt';
        filename = sprintf(formatSpec,l,N);
        fileID = fopen(filename,'w');
        median_record(TotalCurrent(:,N+1),fileID);  % record the current median for comparing
        fclose(fileID);
        copyfile(filename,foldername);
        
    end
end

delete *NeuronV_*.*         % clear unnecessary files
delete *CurrentV_*.*
toc

delete(gcp('nocreate'));


