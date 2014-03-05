import math, random, struct, numpy
import os

#Calculates distance between two points. Takes x,y,z coordinates of each point as arguments
def d(x1,y1,z1,x2,y2,z2):
    intermediate = math.pow((x1-x2),2) + math.pow((y1-y2),2) + math.pow((z1-z2),2)
    return math.sqrt(intermediate)
#Calculates the mean minimum distance between 2 fibers. f1 and f2 are arrays
def meanD(f1,f2):
    #sum1 and counter are used to calculate the average
    sum1=0
    counter=0
    #for each point in f1
    for x in range(len(f1)):
        #set the minimum distance high
        minD=999
        #Loop through each point in f2
        for y in range(len(f2)):
            #if the distance between x in f1 and y in f2 is less than minD
            if(d(f1[x][0],f1[x][1],f1[x][2],f2[y][0],f2[y][1],f2[y][2]) < minD):
                #change the value of minD to that distance
                minD=d(f1[x][0],f1[x][1],f1[x][2],f2[y][0],f2[y][1],f2[y][2])
        #Add the minD for each point x to sum1
        sum1+=minD
        counter+=1
        #print minD, "x, y = ", x, y
        #counter+=1
    #after finding the minD for each point in f1, the average is calculated and returned
    if(counter!=0):
        #print counter
        #print sum1
        return sum1/counter
    else:    
        return 1000

#Mean D between 2 fibers is not necessarily symmetric. This method averages the meanD of 2 fibers.
#This is the metric used in the Corouge paper.
def realMeanD(f1,f2):
    return (meanD(f1,f2)+meanD(f2,f1))/2
    
#Crappy Distance metric that returns the minimum distance between the two fibers
def minD(f1,f2):
    currentMin=1000
    #for each point in f1
    for x in range(len(f1)):
        #Loop through each point in f2
        for y in range(len(f2)):
            #if the distance between x in f1 and y in f2 is less than minD
            if(d(f1[x][0],f1[x][1],f1[x][2],f2[y][0],f2[y][1],f2[y][2]) < currentMin):
                #change the value of currentMin to that distance
                currentMin=d(f1[x][0],f1[x][1],f1[x][2],f2[y][0],f2[y][1],f2[y][2])
    return currentMin

#the Hausdorf distance is the maximum distance between a point in f1 and the closest point in f2
def hausD(f1,f2):
    #set the maximum distance low
    maxD=0
    #for each point in f1
    for x in range(len(f1)):
        currentD=1000
        #Loop through each point in f2
        for y in range(len(f2)):
            #if the distance between x in f1 and y in f2 is greater than maxD
            if(d(f1[x][0],f1[x][1],f1[x][2],f2[y][0],f2[y][1],f2[y][2]) < currentD):
                #change the value of maxD to that distance
                currentD=d(f1[x][0],f1[x][1],f1[x][2],f2[y][0],f2[y][1],f2[y][2])
        if(currentD>maxD):
            maxD=currentD
    return maxD

#returns the legit Hausdorf Distance
def maxHausD(f1,f2):
    if(hausD(f1,f2)>hausD(f2,f1)):
        return hausD(f1,f2)
    else:
        hausD(f2,f1)
    

#This method takes an array of fibers as input. It returns a matrix with entries
#corresponding to the mean distance between pairs of fibers
def makeClusterMatrix(fibers):
    #empty square matrix with numRows equal to the number of fibers
    clusterMatrix=numpy.zeros(shape=(len(fibers),len(fibers)))
    #for each pair of fibers, fill the corresponding entry in the matrix with something
    for i in range(len(fibers)):
        for j in range(len(fibers)):
            #if the fibers are the same, the distance between them is 0
            if(i==j):
                clusterMatrix[i][j]=0
            #if the fibers are different, fill the entry with the mean minimum distance
            else:
                clusterMatrix[i][j]=realMeanD(fibers[i],fibers[j])
    #print the similarity matrix
    #print clusterMatrix
    return clusterMatrix
    
class Tree(object):
    def __init__(self):
        self.left = None
        self.right = None
        self.data = None
        self.parent = None
 
 
#This method performs single linkage agglomerative clustering on the similarity matrix
def agglomCluster(fibers):
    #nodes contains all the leaves in our tree. The data for each node is the index of the fibers it contains
    nodes=[Tree() for i in range(len(fibers))]
    for i in range(len(fibers)):
        nodes[i].data=[i]
    #Making the similarity matrix
    clusterMatrix=makeClusterMatrix(fibers)
    #n is used to determine the stopping point
    n=1
    #while there is more than one cluster
    while(n<len(fibers)):
        #find the minimum non-zero entry in the matrix
        minEntry=1000
        x = 0
        y = 0
        for i in range(len(clusterMatrix)):
            for j in range(len(clusterMatrix)):
                if(clusterMatrix[i][j]<minEntry and clusterMatrix[i][j]!=0):
                    minEntry=clusterMatrix[i][j]
                    x = i
                    y = j
        #set the x column equal to the min of the x and y columns
        #set the y column equal to 0
        #set the xth entry of each column equal to the min of the xth and yth entry of that column
        #set the yth entry of each column to 0
        for i in range(len(clusterMatrix)):
            clusterMatrix[x][i]=min(clusterMatrix[x][i],clusterMatrix[y][i])
            clusterMatrix[y][i]=0
            clusterMatrix[i][x]=min(clusterMatrix[i][x],clusterMatrix[i][y])
            clusterMatrix[i][y]=0
        #newNode contains our new cluster
        newNode=Tree()
        
        #find the node that represents the cluster containing x
        currNode1=nodes[x]
        while(currNode1.parent!=None):
            currNode1=currNode1.parent
        #this node becomes the left child of our new cluster.
        newNode.left=currNode1
        currNode1.parent=newNode
        
        #repeat for y
        currNode2=nodes[y]
        while(currNode2.parent!=None):
            currNode2=currNode2.parent
        newNode.right=currNode2
        currNode2.parent=newNode
        
        #the cluster's data contains the data for both of its subclusters
        newData=[0 for i in range(len(currNode1.data)+len(currNode2.data))]
        for i in range(len(currNode1.data)):
            newData[i]=currNode1.data[i]
        for i in range(len(currNode2.data)):
            newData[i+len(currNode1.data)]=currNode2.data[i]
        newNode.data=newData
        
        exceptions=[]
        for i in range(len(nodes)):
            first=True
            for j in range(len(exceptions)):
                if(i==exceptions[j]):
                    first=False
            if(i!=x and i!=y and first):
                currentNode=nodes[i]
                newParent=Tree()
                while(currentNode.parent!=None):
                    currentNode=currentNode.parent
                if(currentNode.data!=newNode.data):
                    currentNode.parent=newParent
                    newParent.left=currentNode
                    newParent.data=currentNode.data
                    for j in range(len(newParent.data)):
                        exceptions.append(newParent.data[j])
        
        #print newNode.data
        n=n+1
    #finds the root when the clustering is complete
    root=nodes[0]
    while(root.parent!=None):
            root=root.parent
    #print root.data
    #the root is returned. We may want to return the array nodes, which contains the leaves instead
    #I'm not really sure how making "cuts" works
    return nodes
    
def findKthLevel(nodes,k):
    exceptions=[]
    kLevel=[0 for i in range(len(nodes)-k)]
    next=0
    for i in range(len(nodes)):
        first=True
        for j in range(len(exceptions)):
            if(i==exceptions[j]):
                first=False
        if(first):
            currentNode=nodes[i]
            n=0
            while(n<k):
                currentNode=currentNode.parent
                #print currentNode.data
                n=n+1
            kLevel[next]=currentNode.data
            next=next+1
            for j in range(len(currentNode.data)):
                exceptions.append(currentNode.data[j])
    print kLevel
    return kLevel
        
def readinFibersUtil():
    f = open('/Accounts/lynnz/Desktop/cs318_project_relevant/dti.trk', mode='rb')
    binContent = f.read()
    #Skip the bytes in the header
    offset=1000
    #Number of fibers in the data set (on our current set, it is about 250000)
    numF=struct.unpack("i", binContent[988:992])[0]#Number of fibers in the set
    memoryAddresses=[[0,0] for i in range(numF)]
    for i in range(numF):
        memoryAddresses[i][0]=offset
        tempLength=struct.unpack("i", binContent[offset:offset+4])[0]
        memoryAddresses[i][1]=tempLength
        offset=12*tempLength+offset+4
    
    f.close()
    return memoryAddresses

def readInFibers():
    memoryAddresses=readinFibersUtil()
    f = open('/Accounts/lynnz/Desktop/cs318_project_relevant/dti.trk', mode='rb')
    binContent = f.read()
        
    #print memoryAddresses[i]
    numF=100
    #make the blank fiberArray
    fiberArray=[0 for i in range(numF)]
    for i in range(numF):
        #find the number of points in each fiber
        #numPoints=struct.unpack("i", binContent[offset:offset+4])[0]
        #offset=offset+4
        numPoints=memoryAddresses[i][1]
        offset=memoryAddresses[i][0]+4
        #print numPoints
        #set the ith entry of fiberArray equal to a blank array with the correct number of points
        fiberArray[i]=[[0 for x in range(3)] for y in range(numPoints)]
        #fill the array corresponding to a fiber with points
        for x in range(numPoints):
            for y in range(3):
                fiberArray[i][x][y]=struct.unpack("f", binContent[offset:offset+4])[0]
                print fiberArray[i][x][y]
                offset=offset+4
    #print fiberArray
    return fiberArray

def selectRandomFibers(numFibers=100, mode='default', lengthThreshold=25):
    memoryAddresses=readinFibersUtil()
    f = open('/Accounts/lynnz/Desktop/cs318_project_relevant/dti.trk', mode='rb')
    binContent = f.read()
    fibersUsed=[]
    while(len(fibersUsed)!=numFibers):
        x=random.randint(0,len(memoryAddresses))
        if (fibersUsed.count(x)!=0):
            pass
        elif(memoryAddresses[x][1]<lengthThreshold):
            pass
        else:
            fibersUsed.append(x)
    for i in range(numFibers):
        newFiber=fibersUsed[i]
        fibersUsed[i]=memoryAddresses[newFiber]
    memoryAddresses=fibersUsed
    return memoryAddresses

def readWriteTest():
    f = open('/Accounts/lynnz/Desktop/cs318_project_relevant/dti.trk', mode='rb')
    binContent = f.read()
    directoryPath=os.getcwd()
    newFile = open(directoryPath + "/newTrack1.trk", "w")
    newFile.write(binContent[:])
    
#this method writes a track file for the randomized set of fibers (below the
#length threshold). 
def writeTrkFile(memoryAddresses):
    numF=len(memoryAddresses)
    f = open('/Accounts/lynnz/Desktop/cs318_project_relevant/dti.trk', mode='rb')
    binContent = f.read()
    directoryPath=os.getcwd()
    newFile = open(directoryPath + "/newTrack.trk", "w")
    header=binContent[0:988]
    newFile.write(header)
    newFile.write(struct.pack('<i',numF))
    endHeader=binContent[992:1000]
    newFile.write(endHeader)

    for i in range(numF):
        offset=memoryAddresses[i][0]
        length=memoryAddresses[i][1]
        newFile.write(binContent[offset:offset+12*length+4])
            

def readFibers(memoryAddresses):
    f = open('/Accounts/lynnz/Desktop/cs318_project_relevant/dti.trk', mode='rb')
    binContent = f.read()
    fiberArray=[0 for i in range(len(memoryAddresses))]
    for i in range(len(memoryAddresses)):
        #find the number of points in each fiber
        #numPoints=struct.unpack("i", binContent[offset:offset+4])[0]
        #offset=offset+4
        numPoints=memoryAddresses[i][1]
        offset=memoryAddresses[i][0]+4
        #print numPoints
        #set the ith entry of fiberArray equal to a blank array with the correct number of points
        fiberArray[i]=[[0 for x in range(3)] for y in range(numPoints)]
        #fill the array corresponding to a fiber with points
        for x in range(numPoints):
            for y in range(3):
                fiberArray[i][x][y]=struct.unpack("f", binContent[offset:offset+4])[0]
                offset=offset+4
    #print fiberArray
    return fiberArray
    

def main():
    #readWriteTest()
    memoryAddresses=selectRandomFibers(500)
    print memoryAddresses
    fibers = readFibers(memoryAddresses)
    writeTrkFile(memoryAddresses)
    #makeClusterMatrix(fibers)
    #nodes=agglomCluster(fibers)
    #for i in range(len(fibers)):
        #findKthLevel(nodes,i)

if __name__ == "__main__":

    main()