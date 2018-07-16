#!/usr/bin/python2
import numpy as np 
from pylab import*
import scipy as sp




if __name__ == '__main__':
	NeuronNumber=8
	ShiftingNeuron=7
	InhibitionNeuron=NeuronNumber+ShiftingNeuron*2+1
	CoupledNe=InhibitionNeuron+2
	BaseFrequencyNe=CoupledNe+1
	FMNe=BaseFrequencyNe+1
	
	
	
	output=open('Connection_Table_temp.txt','w')
	output2=open('Connection_Table_temp_short.txt','w')
	#output3=open('Connection_LeftShfiting_Table_temp.txt','w')
	#output4=open('Connection_Inhibition_Table_temp.txt','w')
	
	
	
	ExcitationToInhibition=1
	GlobalInhibition=-5.8
	BombToShift=0.35
	ShiftToBomb=0.3
	#FMEToBump=4
	FMEToShift=4	
	#BaseToFM=2
	CoupledInter=0.5
	CoupleToFM=BaseToFM=1.5
	'''
	ExcitationToInhibition=0
	GlobalInhibition=0
	BombToShift=0
	ShiftToBomb=0
	FMEToBump=0
	FMEToShift=0	
	BaseToFM=0
	CoupledInter=0
	CoupleToFM=0
	'''
	n=1
	for i in range(1,NeuronNumber+1):
		
		#print>>output2,i,i+1,InterExcitation
		print>>output,i,InhibitionNeuron,ExcitationToInhibition
		print>>output,InhibitionNeuron,i,GlobalInhibition
		#print>>output2,FMNe,i,FMEToBump

		if i==1:
			
			print>>output,i,1+ShiftingNeuron+i,BombToShift
			print>>output,1+ShiftingNeuron+i,i+1,ShiftToBomb
		elif i==NeuronNumber:
			
			print>>output,i,i+(2*ShiftingNeuron),BombToShift
			print>>output,i+(2*ShiftingNeuron),i-1,ShiftToBomb

		else:
			
			print>>output,i,1+ShiftingNeuron+i,BombToShift
			print>>output,1+ShiftingNeuron+i,i+1,ShiftToBomb
	
			print>>output,i,i+(2*ShiftingNeuron),BombToShift
			print>>output,i+(2*ShiftingNeuron),i-1,ShiftToBomb

	print>>output,FMNe,FMNe,0

	for i in range(NeuronNumber+1,NeuronNumber+ShiftingNeuron*2+1):
		print>>output2,FMNe,i,FMEToShift
	print>>output2,BaseFrequencyNe,FMNe,BaseToFM
	print>>output2,CoupledNe,FMNe,CoupleToFM
	print>>output2,CoupledNe,CoupledNe-1,CoupledInter
	print>>output2,CoupledNe-1,CoupledNe,CoupledInter
	
	









	
