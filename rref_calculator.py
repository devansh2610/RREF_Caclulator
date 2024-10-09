#INPUT FROM FILE
with open('matrix.txt') as myfile:
    for line in myfile:
        if 'Enter number of rows:' in line:
            rows=int(line.strip('Enter number of rows: '))
        if 'Enter number of columns' in line:
            columns=int(line.strip('Enter number of columns: '))
    A=[[0 for c in range(columns)] for r in range(rows)]
    myfile.seek(0,0)
    for line in myfile:
        if 'Enter entries in A' in line:
            colonindex=line.find(":")
            line=line.strip('\n')
            line=line.strip('Enter entries in A:')
            if ", " in line:
                listofentries=line.split(", ")
            else:
                listofentries=line.split(" ")

    if len(listofentries)!=(rows*columns):
        print("Please enter correct number of rows and columns in the input file")
    else:
        for box in A:
            for i in range(len(box)):
                box[i]=float(listofentries[i])
            for j in range(columns):    
                listofentries.pop(0)
                
        #FUNCTIONS:
        def interchange(i,j,M=A):
            M[i],M[j]=M[j],M[i]
            return M

        def scalarmult(i,k,M=A):
            for j in range(len(M[i])):
                M[i][j] = k*(M[i][j])
            return M

        def subtrow(j,i,c,M=A):
            for k in range(len(M[i])):
                M[j][k]= M[j][k] - (c*(M[i][k]))
            return M
        
        def transpose(M=A):
            return [[M[p][q] for p in range(len(M))] for q in range(len(M[0]))]

        def negposzero(M=A):
            for r in range(len(M)):
                for p in range(len(M[r])):
                    if M[r][p]== -0.0:
                        M[r][p] = 0.0

        #BRINGING PIVOT TO A11
        if A[0][0]==0:
            for box in A:
                if box[0]!=0:
                    interchange(A.index(box),0)
                    break

        #CREATING ECHELON FORM
        i=1
        while True:
            if i>(columns-1):
                break
            for row in range(i,rows):
                if A[row][i-1]==0:
                    pass
                else:
                    if (A[i-1][i-1])==0:
                        interchange(row,i-1)
                    else:
                        subtrow(row,i-1,(A[(row)][i-1]/A[i-1][i-1]))
            i=i+1

        #CONVERTING PIVOTS INTO 1 (DIVIDING ROW BY ITS LEADING ENTRY)
        if A[0][0]!=0:
            scalarmult(0,(1/A[0][0]))
        for box in A:
            if A.index(box)==0:
                pass
            else:
                for elt in box:
                    if elt!=0:
                        scalarmult(A.index(box),(1/elt))
                        break

        #CREATING 0'S ABOVE 1'S
        j=1
        while True:
            if j>(columns-1):
                break
            for row in range(j,rows):
                if A[row][j]==1:
                    for k in range(row):
                        if A[k][j]!=0:
                            subtrow(k,row,(A[k][j]/A[row][j]))
            j=j+1

        #SWAPPING FOR CREATING 1'S AS PIVOT
        m=1
        while True:
            if m>(columns-1):
                break
            for row in range(m,rows):
                if A[row][m]==1 and A[row-1][m-1]!=1:
                    interchange(row,(row-1))
            m=m+1

        #CREATING 0 ABOVE A PIVOT IF ANY LEFT
        n=1
        while True:
            if n>(columns-1):
                break
            for row in range(n,rows):
                if A[row][n]!=1:
                    if 1.0 in A[row]:
                        index=A[row].index(1.0)
                        for f in range(row):
                            if A[f][index]!=0:
                                subtrow(f,row,(A[f][index]/A[row][index]))
            n=n+1

        #MAKING -0.0 INTO +0.0
        negposzero(A)

        #IF A ROW HAS ALL ZEROES AND ONE PIVOT WHICH IS NOT 1
        for box in A:
            if box.count(0)==(columns-1):
                for t in range(len(box)):
                    if box[t]!=0:
                        scalarmult(A.index(box),(1/box[t]))
                        for row in range(A.index(box)):
                            subtrow(row,A.index(box),A[row][t])

        #PRINTING RREF
        for box in A:
            for u in range(len(box)):
                box[u]=round(box[u],5)
        print("The matrix A in RREF is: ")
        for box in A:
            print(box)
        print()

        #IDENTIFYING BASIC AND FREE VARIABLES
        basicvar=set()
        freevar = set()
        X=transpose(A)
        for p in range(len(X)):
            if X[p]==[0 for x in range(rows)]:
                freevar.add(p)
        for p in range(len(X)):
            for q in range(len(X[p])):
                if X[p][q]>1 or X[p][q]<0 or (X[p][q]<1 and X[p][q]>0):
                    freevar.add(p)
        for p in range(len(X)):
            if p not in freevar:
                basicvar.add(p)

        freevar=[x for x in freevar]
        basicvar=[x for x in basicvar]

        #CREATING A DICTIONARY THAT CONTAINS RREF OF EQUATIONS (FOR EG: x0 = -6*x1 -2*x3)
        somedict={}
        for b in basicvar:
            somedict[f'x{b}'] = []
        for t in freevar:
            somedict[f'x{t}'] = [f'x{t}']
        for b in basicvar:
            index=X[b].index(1)
            for t in freevar:
                somedict[f'x{b}'].append((-X[t][index], f'x{t}'))
        vardict={}
        for key in sorted(somedict):
            vardict[key] = somedict[key]

        #CREATING A DUMMY MATRIX FOR CONTAINING COEFFICIENT VECTORS OF FREE VARIABLES
        Y=[[0.0 for x in range(len(freevar))] for y in range(columns)]
        freevarindex={}
        for p in range(len(freevar)):
            freevarindex[f'x{freevar[p]}']= p
        
        #CREATING COEFFICIENT MATRICES OF FREE VARIABLES
        for key in vardict.keys():
            if key in freevarindex.keys():
                for k,v in vardict.items():
                    if type(v[0])==tuple:
                        for p in range(len(v)):
                            Y[int(k.split('x')[1])][freevarindex[v[p][1]]]= v[p][0]
                    if v[0]==k:
                        Y[int(k.split('x')[1])][freevarindex[k]]=1.0
        
        #PRINTING GENERAL SOLUTION IN PARAMETRIC FORM
        negposzero(Y)
        vectorX = transpose(Y)
        print("The general solution of AX=0 in parametric form is: ")
        if len(vectorX)==0:
            print(f'{[0 for x in range(columns)]}')
        else:
            for key in freevarindex.keys():
                if freevarindex[key]!=freevar.index(freevar[-1]):
                    print(f'{key}*{vectorX[freevarindex[key]]}',end=" + ")
                else:
                    print(f'{key}*{vectorX[freevarindex[key]]}',end=" ")
                    



    
    


    

