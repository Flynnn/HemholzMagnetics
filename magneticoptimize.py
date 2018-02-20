from scipy.constants import *
import numpy as np
import random
import math 

def DetermineRegionOfUniformity(Distance, Width, Height, I, n, segments, toleranceChecks, tolerance):
	centerField = BFromCoilQuad(np.array([0,0,0]), np.array([1,0,0]), np.array([0,1,0]), Width, Height, Distance, I, n, segments, np.array([0,0,0]))
	print(centerField)
	ROU = np.array([Distance, Width, Height])
	for i in range(0, toleranceChecks):
		RP = np.array([random.uniform(-ROU[0],ROU[0]), random.uniform(-ROU[1],ROU[1]), random.uniform(-ROU[2],ROU[2])])
		sampleField = BFromCoilQuad(np.array([0,0,0]), np.array([1,0,0]), np.array([0,1,0]), Width, Height, Distance, I, n, segments, RP)
		diff = centerField - sampleField
		error = np.linalg.norm(diff)/np.linalg.norm(centerField)
		if (error > tolerance):#We are outside our tolerances. Reduce the ROU guess the least amount possible.
			reduced = False
			if (abs(RP[0]) < ROU[0] and abs(RP[0]) > max(abs(RP[1]), abs(RP[2]))):
				reduced = True
				ROU[0] = abs(RP[0])
			if (abs(RP[1]) < ROU[1] and abs(RP[1]) > max(abs(RP[0]), abs(RP[2]))):
				reduced = True
				ROU[1] = abs(RP[1])
			if (abs(RP[2]) < ROU[2] and abs(RP[2]) > max(abs(RP[1]), abs(RP[0]))):
				reduced = True
				ROU[2] = abs(RP[2])
			assert(reduced)
		print(ROU)
		print('.', end='', flush=True)
	return ROU
		
def BFromCoilQuad(Center, Axis, Horizontal, Width, Height, Distance, I, n, segments, POI):
	widthVector = Width * Horizontal
	heightVector = Height * np.cross(Horizontal, Axis)
	return BFromSquareCoil(Center + Axis * Distance/2, widthVector, heightVector, I, n, segments/4, POI) + BFromSquareCoil(Center - Axis * Distance/2, widthVector, heightVector, I, n, segments/4, POI) + BFromSquareCoil(Center + Axis * Distance/4, widthVector, heightVector, I, n, segments/4, POI) + BFromSquareCoil(Center - Axis * Distance/4, widthVector, heightVector, I, n, segments/4, POI)
		
def BFromCoilPair(Center, Axis, Horizontal, Width, Height, Distance, I, n, segments, POI):
	widthVector = Width * Horizontal
	heightVector = Height * np.cross(Horizontal, Axis)
	return BFromSquareCoil(Center + Axis * Distance/2, widthVector, heightVector, I, n, segments/2, POI) + BFromSquareCoil(Center - Axis * Distance/2, widthVector, heightVector, I, n, segments/2, POI)

def BFromSquareCoil(CoilCenter, CoilWidth, CoilHeight, I, n, segments, POI):
	A = CoilCenter + CoilWidth/2 + CoilHeight/2
	B = CoilCenter - CoilWidth/2 + CoilHeight/2
	C = CoilCenter - CoilWidth/2 - CoilHeight/2
	D = CoilCenter + CoilWidth/2 - CoilHeight/2
	return BFromLine(A,B,I*n,segments/4,POI) + BFromLine(B,C,I*n,segments/4,POI) + BFromLine(C,D,I*n,segments/4,POI) + BFromLine(D,A,I*n,segments/4,POI)

#Field returns are kg/(s^2 A)
#Inputs are meters, amps.
def BFromLine(A, B, I, segments, POI):
	segments = int(segments)
	dL = (A - B)/segments
	total = np.array([0.0,0.0,0.0])
	for i in range(0, segments):
		p = float(i+0.5)/float(segments)
		point = (B - A) * p + A
		rvector = point - POI
		rmag = np.linalg.norm(rvector)
		total += mu_0 * I * np.cross(dL, rvector / rmag) / (4 * pi * rmag * rmag)
	return total
		
#2ft by 2ft region.
#4.97A
#54 turns each coil
#0.508m between coils
#print(BFromCoilPair(np.array([0,0,0]), np.array([1,0,0]), np.array([0,1,0]), 0.6096, 0.6096, 0.508, 4.97, 54, 200, np.array([0,0,0])))
print(DetermineRegionOfUniformity(0.50, 0.6096, 0.6096, 4.97, 27, 100, 2000, 0.01))