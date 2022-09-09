#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import math
from equivalentJointLoad import equivalentJointLoad 


# In[ ]:


numberOfNodes = int(input ("Enter number of nodes: "))
numberOfBeams = numberOfNodes - 1


# In[ ]:


# getNodeDetails
coordinatesOfNodes=[]
for i in range(0,numberOfNodes):
    print('Enter the coordinates for node number ', i,' in (x,y) format:')
    x = input()
    x = tuple(int(a) for a in x.split(","))
    coordinatesOfNodes.append(x)
# print(coordinatesOfNodes)


# In[ ]:


# Get beam start and end nodes
detailsOfBeam=[]
equivalentJLoadsOfBeam = []
for i in range(0,numberOfBeams):
    print('Enter the Young\'s Modulus, Moment of Intertia, load being applied, no. of elements in the beam) for beam number ', i,':')
    #x = tuple(int(a) for a in x.split(",") if type(a)!=float else float(a))
    E=float(input("Enter the Young\'s Modulus( In Scientific Notation E ): "))
    I=float(input("Enter Moment of Intertia( In Scientific Notation E ): "))
    L=float(input("Enter load being applied: "))
    div=int(input("No. of elements in the beam: "))
    x=E,I,L,div
    
    eq=equivalentJointLoad(x[2])
    equivalentJLoadsOfBeam.append(eq)    
    
    detailsOfBeam.append(x)
# print(nodesOfBeam)


# In[ ]:


# Get beam start and end nodes
nodesOfBeam=[]
for i in range(0,numberOfBeams):
    print('Enter start and end nodes for beam number ', i,' in (start node number,end node number) format:')
    x = input()
    x = tuple(int(a) for a in x.split(","))
    nodesOfBeam.append(x)
    
    d=detailsOfBeam[i]
    d=list(d)
    d.append(x)
    detailsOfBeam[i]=tuple(d)
# print(nodesOfBeam)


# In[ ]:


#Making the connectivity matrix.
connectivityMatrix=[]

for i in range(0,numberOfNodes):
    config={}
    
    for j in range(0,numberOfBeams):
        bothNodesOfABeam = detailsOfBeam[j][4]
        
        if i == (bothNodesOfABeam[0]-1):
            config["i"]=j
        
        elif i == (bothNodesOfABeam[1]-1):
            config["j"]=j
        
    if "i" not in config.keys():
        config["i"]=None

    if "j" not in config.keys():
        config["j"]=None
    
    connectivityMatrix.append(config)


# In[ ]:


# identify known/Unknown displacement/force Of Nodes
knownUnknownOfNodes=[]
print('\nEnter details in (k,u) format (k:known, u:unknown):', i,':\n')
for i in range(0,numberOfNodes):
    print('Enter Known/Unknown displacement for node number ', i,':')
    print('Enter Known/Unknown rotation for node number ', i,':')
    # print('Enter known/Unknown force for node number ', i,':')
    # print('Enter known/Unknown moment for node number ', i,':')
    x = input()
    x = tuple(a for a in x.split(","))
    force = 'u' if x[0]=='k' else 'k'
    moment = 'u' if x[1]=='k' else 'k'
    x = x + tuple(force) + tuple(moment)
    knownUnknownOfNodes.append(x)
# print(DoFOfNodes)


# In[ ]:


# identify boundary conditions (Degree of Freedom)
DoFOfNodes=[]
print('\nUnknown details for a node are not asked.\n')
for i in range(0,numberOfNodes):
    knownUnknownForNodes = knownUnknownOfNodes[i]
    ithBeamForPresentNode = connectivityMatrix[i]["i"]
    jthBeamForPresentNode = connectivityMatrix[i]["j"]
    
    if ithBeamForPresentNode == None and jthBeamForPresentNode != None:
        equivalentForceToApply = equivalentJLoadsOfBeam[jthBeamForPresentNode][3]+0
    elif ithBeamForPresentNode != None and jthBeamForPresentNode == None:
        equivalentForceToApply = 0+equivalentJLoadsOfBeam[ithBeamForPresentNode][2]
    elif ithBeamForPresentNode != None and jthBeamForPresentNode != None:
        equivalentForceToApply = equivalentJLoadsOfBeam[jthBeamForPresentNode][3] + equivalentJLoadsOfBeam[ithBeamForPresentNode][2]
    else:
        print("The connectivity of beams is wrong.\n")
    
    if ithBeamForPresentNode == None and jthBeamForPresentNode != None:
        equivalentMomentToApply = equivalentJLoadsOfBeam[jthBeamForPresentNode][1]+0
    elif ithBeamForPresentNode != None and jthBeamForPresentNode == None:
        equivalentMomentToApply = 0-equivalentJLoadsOfBeam[ithBeamForPresentNode][0]
    elif ithBeamForPresentNode != None and jthBeamForPresentNode != None:
        equivalentMomentToApply = equivalentJLoadsOfBeam[jthBeamForPresentNode][1] - equivalentJLoadsOfBeam[ithBeamForPresentNode][0]
    else:
        print("The connectivity of beams is wrong.\n")
    
    
    if knownUnknownForNodes[0]=='k':
        print('Enter the known displacement for node number ', i,':')
        displacement = float(input())
    else:
        displacement = 1
    
    if knownUnknownForNodes[1]=='k':
        print('Enter the known rotation for node number ', i,':')
        rotation = float(input())
    else:
        rotation = 1
    
    if knownUnknownForNodes[2]=='k':
        print('Enter the known force for node number ', i,':')
        #force = float(input())
        extJointLoadApplied = float(input("Enter external point load applied at node[If none enter 0]: "))
        if extJointLoadApplied!=0:
            force = extJointLoadApplied-equivalentForceToApply
        else:
            force = -equivalentForceToApply
    else:
        force = 1
    
    if knownUnknownForNodes[3]=='k':
        print('Enter the known moment for node number ', i,':')
        #moment = float(input())
        extJointMomentApplied = float(input("Enter external concentrated moment applied at node[If none enter 0]: "))
        if extJointMomentApplied!=0:
            moment = extJointMomentApplied+equivalentMomentToApply
        else:
            moment = equivalentMomentToApply
    else:
        moment = 1
    
    x = (displacement,rotation,force,moment)
    DoFOfNodes.append(x)
# print(DoFOfNodes)


# In[ ]:


# Calculate the length of Beam
lengthOfBeam=[]
for i in range(1,numberOfNodes):
    x1 = coordinatesOfNodes[nodesOfBeam[i-1][0]-1][0]
    y1 = coordinatesOfNodes[nodesOfBeam[i-1][0]-1][1]
    x2 = coordinatesOfNodes[nodesOfBeam[i-1][1]-1][0]
    y2 = coordinatesOfNodes[nodesOfBeam[i-1][1]-1][1]
    #print(x1,y1,x2,y2)
    lengthOfBeam.append(math.sqrt(((x2-x1)**2) + ((y2-y1)**2)))
    print('length of Beam ',i,' is ',lengthOfBeam[i-1])


# In[ ]:


# Divide the beam in number of divisions
# Storing details in 2D array
detailsOfElements=[]
elementNumber=0
globalCountOfElements=0

# loop to itterate through all the beams
for i in range(0,numberOfBeams):
    firstCoordinateofBeam=coordinatesOfNodes[i]
    lastCoordinateofBeam=coordinatesOfNodes[i+1]
    # beamNumber=i
    
    numberOfElements=detailsOfBeam[i][3]
    lengthOfElement=lengthOfBeam[i]/numberOfElements
    print("Length Of Element:",lengthOfElement)
    
    # loop to itterate through all the elements
    for j in range(0,numberOfElements):
        
        YoungsModulus=detailsOfBeam[i][0]
        MomentofIntertia=detailsOfBeam[i][1]
        loadBeingApplied=detailsOfBeam[i][2]
        
        # assuming slope = 0
        unknownX=lengthOfElement+firstCoordinateofBeam[0] if j==0 else lastCoordinateofBeam[0]-lengthOfElement
        unknownY=firstCoordinateofBeam[1] if i==(numberOfBeams-1) else lastCoordinateofBeam[1]
        
        ithNode = firstCoordinateofBeam if j==0 else detailsOfElements[j-1][4]
        jthNode = lastCoordinateofBeam if j==(numberOfElements-1) else (unknownX,unknownY)
        
        # details=(elementNumber,YoungsModulus,MomentofIntertia,loadBeingApplied,ithNode,jthNode)
        details=(YoungsModulus,MomentofIntertia,lengthOfElement,loadBeingApplied,ithNode,jthNode)
        #elementDetail=[beamNumber,details]
        #detailsOfElements.append(elementDetail)        
        detailsOfElements.append(details)        
        
        # globalCountOfElements=elementNumber    
        elementNumber=elementNumber+1
        globalCountOfElements=elementNumber
        


# In[ ]:


def stiffnessMat(youngsModulus,momentOfInertia,elementLength):
        L=elementLength
        k=np.array([[12,6*L,-12,6*L],[6*L,4*L*L,-6*L,2*L*L],[-12,-6*L,12,-6*L],[6*L,2*L*L,-6*L,4*L*L]])
        K=((youngsModulus*momentOfInertia)/(elementLength**3))*k
        return K


# In[ ]:


# calculate element freedom table
d=0
r=1
globalNoOfNodes=0
countOfElements=0
elementFreedomTable=[]

# loop to itterate through all the beams
for i in range(0,numberOfBeams):
    # beamNumber=i
    numberOfElements=detailsOfBeam[i][3]
    # loop to itterate through all the elements
    for j in range(0,numberOfElements):
        # tempEFT=[d,r,d+2,r+2]
        # beamEFT=[beamNumber,tempEFT]
        elementEFT=[d,r,d+2,r+2]
        elementFreedomTable.append(elementEFT)
        d=d+2
        r=r+2
        countOfElements=countOfElements+1
globalNoOfNodes=countOfElements+1


# In[ ]:


# calculate stiffness matrix for all elements
# loop to itterate through all the elements
elementStiffnessMatrix=[]
for i in range(0,globalCountOfElements):
    stiffnessMatrix=stiffnessMat(detailsOfElements[i][0],detailsOfElements[i][1],detailsOfElements[i][2])
    elementStiffnessMatrix.append(stiffnessMatrix)


# In[ ]:


originalglobalStiffnessMatrix=np.zeros((globalNoOfNodes*2,globalNoOfNodes*2))
for i in range(0,globalCountOfElements):
        eFT = elementFreedomTable[i]
        k = elementStiffnessMatrix[i]
        
        for a in range(0,4):
            for b in range(0,4):
                originalglobalStiffnessMatrix[eFT[a],eFT[b]]=originalglobalStiffnessMatrix[eFT[a],eFT[b]]+k[a,b]


# In[ ]:


originalDisplacementMatrix=np.zeros((globalNoOfNodes*2,1))
#originalDisplacementMatrix=np.zeros((globalNoOfNodes*2))

firstDof=0
secondDof=1

#displacementDofNodeIndexes=[]

beamNumber=0
# numberOfBeams
for i in range(0,globalNoOfNodes):
    if i == 0:
        actualNodeNumber=0
        originalDisplacementMatrix[firstDof][0]=DoFOfNodes[actualNodeNumber][0]
        originalDisplacementMatrix[secondDof][0]=DoFOfNodes[actualNodeNumber][1]
        nextMainNode=i+detailsOfBeam[beamNumber][3]
        actualNodeNumber=actualNodeNumber+1
    elif i == nextMainNode:
        beamNumber=beamNumber+1
        originalDisplacementMatrix[firstDof][0]=DoFOfNodes[actualNodeNumber][0]
        originalDisplacementMatrix[secondDof][0]=DoFOfNodes[actualNodeNumber][1]
        if beamNumber!=numberOfBeams:
            nextMainNode=i+detailsOfBeam[beamNumber][3]
        else:
            pass
        actualNodeNumber=actualNodeNumber+1
    else:
        originalDisplacementMatrix[firstDof][0]=1 #  intermediate nodes 
        originalDisplacementMatrix[secondDof][0]=1 #  intermediate nodes
    
    
    #displacementDofNodeIndexes.append((firstIndex,secondINdex))
    
    firstDof=firstDof+2
    secondDof=secondDof+2
    #n+=1


# In[ ]:


originalForceMatrix=np.zeros((globalNoOfNodes*2,1))
#originalForceMatrix=np.zeros((globalNoOfNodes*2))

firstDof=0
secondDof=1

#forceDofNodeIndexes=[]

beamNumber=0
# numberOfBeams
for i in range(0,globalNoOfNodes):
    if i == 0:
        actualNodeNumber=0
        originalForceMatrix[firstDof][0]=DoFOfNodes[actualNodeNumber][2]
        originalForceMatrix[secondDof][0]=DoFOfNodes[actualNodeNumber][3]
        nextMainNode=i+detailsOfBeam[beamNumber][3]
        actualNodeNumber=actualNodeNumber+1
    elif i == nextMainNode:
        beamNumber=beamNumber+1
        originalForceMatrix[firstDof][0]=DoFOfNodes[actualNodeNumber][2]
        originalForceMatrix[secondDof][0]=DoFOfNodes[actualNodeNumber][3]
        if beamNumber!=numberOfBeams:
            nextMainNode=i+detailsOfBeam[beamNumber][3]
        else:
            pass
        actualNodeNumber=actualNodeNumber+1
    else:
        originalForceMatrix[firstDof][0]=0 # intermediate nodes
        originalForceMatrix[secondDof][0]=0 #  intermediate nodes

    #forceDofNodeIndexes.append((firstIndex,secondIndex))
    
    firstDof=firstDof+2
    secondDof=secondDof+2


# In[ ]:


#Creating the modified displacement matrix
u=originalDisplacementMatrix

knownDisplacementIndexes=[]
unknownDisplacementIndexes=[]
modifiedDisplacementIndexes=[]

knownDisplacementValues=[]
unknownDisplacementValues=[]

for i in range(0,len(originalDisplacementMatrix)):
    if u[i][0]==0:
        knownDisplacementIndexes.append(i)
        knownDisplacementValues.append(u[i][0])#what are the data types that are being appended to the knownDisplacementValues. np.array
    elif u[i][0]==1:
        unknownDisplacementIndexes.append(i)
        unknownDisplacementValues.append(u[i][0])
    else:
        knownDisplacementIndexes.append(i)
        knownDisplacementValues.append(u[i][0])
        
modifiedDisplacementIndexes=unknownDisplacementIndexes+knownDisplacementIndexes


knownDisplacementValuesMatrix=np.zeros((len(knownDisplacementValues),1))

for i in range(0,len(knownDisplacementValues)):
    
    knownDisplacementValuesMatrix[i][0] = knownDisplacementValues[i]

#knownDisplacementValuesMatrix=np.array(knownDisplacementValues)
#unknownDisplacementValues are mainly None value and needs to be replaced so IDK what to do to them


# In[ ]:


#Creating the modified force matrix
f=originalForceMatrix

knownForceIndexes=[]
unknownForceIndexes=[]
modifiedForceIndexes=[]

knownForceValues=[]
unknownForceValues=[]

for i in range(0,len(originalDisplacementMatrix)):
    if f[i][0]==0:
        knownForceIndexes.append(i)
        #print("Known value being added: ",u[i][0]," FOr index ",i)
        knownForceValues.append(f[i][0])
    elif f[i][0]==1:
        unknownForceIndexes.append(i)
        unknownForceValues.append(f[i][0])
    else:
        knownForceIndexes.append(i)
        knownForceValues.append(f[i][0])
        
modifiedForceIndexes=knownForceIndexes+unknownForceIndexes


knownForceValuesMatrix=np.zeros((len(knownForceValues),1))

for i in range(0,len(knownForceValues)):
    
    knownForceValuesMatrix[i][0] = knownForceValues[i]
#unknownForceValues are just None values and need to be replaced.


# In[ ]:


#Creating modified global Stiffness Matrix

modifiedglobalStiffnessMatrix=np.zeros((globalNoOfNodes*2,globalNoOfNodes*2))
for i in range(0,globalNoOfNodes*2):
    for j in range(0,globalNoOfNodes*2):
        #print("FOR I ",i," We have J ",j)
        #modifiedglobalStiffnessMatrix[modifiedForceIndexes[i],modifiedDisplacementIndexes[i]] = originalglobalStiffnessMatrix[i][j]
        modifiedglobalStiffnessMatrix[i][j] = originalglobalStiffnessMatrix[modifiedForceIndexes[i]][modifiedDisplacementIndexes[j]]


# In[ ]:


#Solving the matrix equations.

K1=modifiedglobalStiffnessMatrix[0:len(knownForceIndexes),0:len(unknownDisplacementIndexes)]

K2=modifiedglobalStiffnessMatrix[0:len(knownForceIndexes),len(unknownDisplacementIndexes):len(modifiedDisplacementIndexes)]

K3=modifiedglobalStiffnessMatrix[len(knownForceIndexes):len(modifiedForceIndexes),0:len(unknownDisplacementIndexes)]

K4=modifiedglobalStiffnessMatrix[len(knownForceIndexes):len(modifiedForceIndexes),len(unknownDisplacementIndexes):len(modifiedDisplacementIndexes)]


#Equation 1
#unknownDisplacementValuesMatrix = np.matmul((knownForceValuesMatrix.transpose() - np.matmul(K2,knownDisplacementValuesMatrix.transpose())),np.linalg.inv(K1))
#unknownDisplacementValuesMatrix = np.matmul(np.linalg.inv(K1),(knownForceValuesMatrix - np.matmul(K2,knownDisplacementValuesMatrix)))

#Solving Equation 1
firstMult=np.matmul(K2,knownDisplacementValuesMatrix)
firstSub=np.subtract(knownForceValuesMatrix, firstMult)
finalMult1=np.matmul(np.linalg.inv(K1),firstSub)

unknownDisplacementValuesMatrix=finalMult1



#Equation 2
#unknownForceValuesMatrix = np.add(np.matmul(K3,unknownDisplacementValuesMatrix),np.matmul(K4,knownDisplacementValuesMatrix.transpose()))
#unknownForceValuesMatrix = np.add(np.matmul(K3,unknownDisplacementValuesMatrix),np.matmul(K4,knownDisplacementValuesMatrix))

#Solving Equation 2
firstMult2=np.matmul(K4,knownDisplacementValuesMatrix)
secondMult2=np.matmul(K3,unknownDisplacementValuesMatrix)
finalAdd2=np.add(firstMult2,secondMult2)

unknownForceValuesMatrix=finalAdd2


# In[ ]:


for i in range(0,len(unknownDisplacementValuesMatrix)):
    originalDisplacementMatrix[unknownDisplacementIndexes[i]]=unknownDisplacementValuesMatrix[i][0]
    
for i in range(0,len(unknownForceValuesMatrix)):
    originalForceMatrix[unknownForceIndexes[i]]=unknownForceValuesMatrix[i][0]


# In[ ]:


print("THE FINAL ANSWER IS:\n\n")
print(originalDisplacementMatrix)
print("\n")
print(originalForceMatrix)


# In[ ]:




