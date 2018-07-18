% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons    Inhibitory neurons

ver distcomp
%parpool('local',4) % parallel with 4 cores.
                   % adjust to the actual processor core number.

clc
clear
tic
SimulationTime=2000; %ms;
DeltaT=0.01; %ms;
Vr=10;
Vth=130;
NoiseStrengthBase=0;
Velocity=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7];
StimuluStrength=[5.5];
StimulusNeuron=1;
StimulationOnset=300;
StimulationOffset=305;
StimulationOnset=StimulationOnset/DeltaT;
StimulationOffset=StimulationOffset/DeltaT;
Direction=[0,4,-4.5];
ModulationCurrent=Direction(2);
Speed=Velocity(4);


BoundNe=8;
ShiftNe=14;
InhibitionNe=1;
FMNe=1;
BaseFrequencyNe=1;
CoupledNe=2;
TotalNe=BoundNe+ShiftNe+InhibitionNe+FMNe+BaseFrequencyNe+CoupledNe;
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
c=c+100;

A=[a(1)*ones(BoundNe,1);a(2)*ones(ShiftNe,1);a(3)*ones(InhibitionNe,1);a(4)*ones(CoupledNe,1);a(5)*ones(BaseFrequencyNe,1);a(6)*ones(FMNe,1)];
B=[b(1)*ones(BoundNe,1);b(2)*ones(ShiftNe,1);b(3)*ones(InhibitionNe,1);b(4)*ones(CoupledNe,1);b(5)*ones(BaseFrequencyNe,1);b(6)*ones(FMNe,1)];
C=[c(1)*ones(BoundNe,1);c(2)*ones(ShiftNe,1);c(3)*ones(InhibitionNe,1);c(4)*ones(CoupledNe,1);c(5)*ones(BaseFrequencyNe,1);c(6)*ones(FMNe,1)];
D=[d(1)*ones(BoundNe,1);d(2)*ones(ShiftNe,1);d(3)*ones(InhibitionNe,1);d(4)*ones(CoupledNe,1);d(5)*ones(BaseFrequencyNe,1);d(6)*ones(FMNe,1)];

Ibias=[IB(1)*ones(BoundNe,1);IB(2)*ones(ShiftNe,1);IB(3)*ones(InhibitionNe,1);IB(4)*ones(CoupledNe,1);IB(5)*ones(BaseFrequencyNe,1);IB(6)*ones(FMNe,1)];
Inoise=NoiseStrengthBase*normrnd(0,1,[TotalNe,1]);
ExternalI=0*ones(TotalNe,1);

v=Vr.*ones(TotalNe,1);
v(TotalNe-1)=70;
u=10*ones(TotalNe,1);
S1=0*ones(TotalNe,1);
S2=0*ones(TotalNe,1);
Tau1=10;
Tau2=5;

Potential=transpose([0;v]);
SynapticCurrent1=transpose([0;S1]);
SynapticCurrent2=transpose([0;S2]);
g1=transpose(full(spconvert(dlmread('Connection_Table_temp.txt'))));
g2=transpose(full(spconvert(dlmread('Connection_Table_temp_short.txt'))));
for l=1:7
%parfor l=1:7

Speed=Velocity(l);
ExternalI=0*ones(TotalNe,1);
v=Vr.*ones(TotalNe,1);
v(TotalNe-1)=70;
u=10*ones(TotalNe,1);
S1=0*ones(TotalNe,1);
S2=0*ones(TotalNe,1);
Tau1=10;
Tau2=5;
Potential=transpose([0;v]);
SynapticCurrent1=transpose([0;S1]);
SynapticCurrent2=transpose([0;S2]);
I=0*ones(TotalNe,1);

for t=1:SimulationTime/DeltaT        % simulation time ms
  disp(l);
  disp(t);
  Inoise=NoiseStrengthBase*normrnd(0,1,[TotalNe,1]);
 
  if t>StimulationOnset && t<=StimulationOffset
  ExternalI(StimulusNeuron)=StimuluStrength(1);
  else
  ExternalI(StimulusNeuron)=0;
  end

  

  fired1=find(v(1:BoundNe+ShiftNe+InhibitionNe)>=Vth);
  fired2=find(v(BoundNe+ShiftNe+InhibitionNe+1:end)>=Vth);
  fired2=fired2+(BoundNe+ShiftNe+InhibitionNe);
  fired=[fired1;fired2];
  v(fired)=C(fired);
  u(fired)=u(fired)+D(fired);

  S1=S1+sum(1*g1(:,fired1),2)-(S1/Tau1)*DeltaT;
  S2=S2+sum(1*g2(:,fired2),2)-(S2/Tau2)*DeltaT;
  if ModulationCurrent>0
      %ExternalI = [zeros(BoundNe,1);ModulationCurrent*ones(ShiftNe/2,1);zeros(ShiftNe/2+InhibitionNe,1);Speed*ones(CoupledNe-InhibitionNe+1,1);zeros(TotalNe-BoundNe-ShiftNe-CoupledNe-1,1)];
	  ExternalI(BoundNe+1:ShiftNe/2+BoundNe)=ModulationCurrent;
	  ExternalI(ShiftNe+BoundNe+InhibitionNe+1:ShiftNe+BoundNe+1+CoupledNe)=Speed;
  elseif ModulationCurrent<0
      %ExternalI = [zeros(ShiftNe,1);abs(ModulationCurrent)*ones(BoundNe,1);zeros(InhibitionNe,1);Speed*ones(CoupledNe-InhibitionNe+1,1);zeros(TotalNe-ShiftNe-BoundNe-CoupledNe-1,1)];
	  ExternalI(ShiftNe+1:ShiftNe+BoundNe)=abs(ModulationCurrent);
	  ExternalI(ShiftNe+BoundNe+InhibitionNe+1:ShiftNe+BoundNe+1+CoupledNe)=Speed;
  else
	  ExternalI(BoundNe+1:ShiftNe+BoundNe)=0;
  end
 
  I=ExternalI+S1+S2+Inoise+Ibias;
  
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
  Potential=cat(1,Potential,temp);
  %temp=transpose([t*DeltaT;S1]);
  %SynapticCurrent1=cat(1,SynapticCurrent1,temp);
  %temp=transpose([t*DeltaT;S2]);
  %SynapticCurrent2=cat(1,SynapticCurrent2,temp);

end
foldername=int2str(l);
mkdir (foldername);


for N=1:TotalNe
    fig=plot(Potential(:,1),Potential(:,N+1),'b-');
    saveas(fig,strcat('NeuronV','_',num2str(N),'.png'));
    copyfile(strcat('NeuronV','_',num2str(N),'.png'),foldername);
    %fig=plot(SynapticCurrent1(:,1),SynapticCurrent1(:,N+1),'b-');
    %saveas(fig,strcat('NeuronC1','_',num2str(N),'.png'));
    %fig=plot(SynapticCurrent2(:,1),SynapticCurrent2(:,N+1),'b-');
    %saveas(fig,strcat('NeuronC2','_',num2str(N),'.png'));
end




end

toc

%parpool close
%delete(gcp('nocreate'))



















