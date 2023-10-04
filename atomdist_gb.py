#################################
######### Hope it helps #########
######### Gourab Bhattacharje ###
######### For public use ########
import numpy as np

def distance_calculation(Coordinates1, Coordinates2): ######################Euclidean Distance Calculation
       
    d1=Coordinates1[0]-Coordinates2[0]
    d2=Coordinates1[1]-Coordinates2[1]
    d3=Coordinates1[2]-Coordinates2[2]

    dist=d1*d1+d2*d2+d3*d3
    dist=np.sqrt(dist)
    return (dist)


PDB_input=input("Enter PDB ID\n")
print ("Entered PDB ID is: "+ PDB_input + "\n")

Chain_input_1=input("Enter First Chain ID\n")
print ("Entered First Chain ID is: "+ Chain_input_1 + "\n")

Chain_input_2=input("Enter Second Chain ID\n")
print ("Entered Second Chain ID is: "+ Chain_input_2 + "\n")

PDB_input_file=(PDB_input+".pdb")
print(PDB_input_file)

#PDB_input=6ifc.pdb
#Chain_input_1=A
#Chain_input_2=B

source=open(PDB_input_file, 'r')     ########### Writing protein coordinates and RNA coordinates separately in two files


out1=open('temporary1.pdb','w')
out2=open('temporary2.pdb','w')

for line in source:
    if (line.startswith('ATOM')):
        if(line[21]==Chain_input_1):
            out1.write(line)
    if (line.startswith('ATOM')):
        if (line[21]==Chain_input_2):
            out2.write(line)
source.close()
out1.close()
out2.close()

Distance_cutoff_input=input("Enter cut-off distance in Angstrom please :P\n")
print ("Entered cut-off distance is: "+ Distance_cutoff_input + " Angstrom\n")

        
source1 = open('temporary1.pdb', 'r')      ###############Reading the coordinates from the protein file and the RNA file
source2 = open('temporary2.pdb', 'r')


        
protein_1=[]
protein_2=[]
for line in source1:
    splitline=line.split()
    protein_1.append([splitline[1], splitline[2], splitline[3], splitline[4], splitline[5],float(splitline[6]), float(splitline[7]), float (splitline[8])])
source1.close()

for line in source2:
    splitline=line.split()
    protein_2.append([splitline[1], splitline[2], splitline[3], splitline[4], splitline[5],float(splitline[6]), float(splitline[7]), float (splitline[8])])
source2.close()

Tot_protein_1= len(protein_1)
Tot_protein_2= len(protein_2)

s=str(Distance_cutoff_input)
print(s)

filname1="distance_less_than_"+s+"angstrom.out"

out3 = open(filname1, 'w')

                                         ###########Calculating the distance between the atoms, checking if they are less than 4.5 Ang and writing the output
count=0
for i in range(Tot_protein_1):
    int_atoms=''
    primary_atom=''
    for j in range(Tot_protein_2):
        a=protein_1[i][5:]
        b=protein_2[j][5:]
        d=distance_calculation(a,b)
        if (d<float(Distance_cutoff_input)):
            s='' + protein_2[j][4]+ '_' + protein_2[j][3]+ '_' + protein_2[j][2]+ '_'+ protein_2[j][1]+ ' ,'
            int_atoms=int_atoms+s
            count+=1
            print(count, d)
    if (int_atoms!=''):
        primary_atom = protein_1[i][4]+ '_' + protein_1[i][3]+ '_' + protein_1[i][2]+ '_' + protein_1[i][1]+ ' ' + 'interacts with:'+ '\t' + int_atoms + '\n'
        out3.write(primary_atom)

out3.close()



