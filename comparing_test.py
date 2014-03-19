#!/usr/bin/env python

from reshape_outputs import * #Import the functions to get the list of the files, and the annotated genes
import os

path="" #Path of where the files are
os.chdir(path)

# Begin is the first position of the gene and end the last one
typ_options=["begin","end"]

#This block creates the list for the begining of the gene predictions of the following files:
typ=typ_options[0]

augustus_b=augustus("augustus.gtf", typ) #Results from Augustus
genemark_b=genemark("genemark.hmm", typ) #Results from genemark
easygene_b=easygene("easygene.gff", typ) #Results from easygene
prodigal_b=prodigal("prodigal.gff", typ) #Results from prodigal
glimmer_b=glimmer("glimmer.gff", typ) #Results from glimmer
annotated_b=features("Annotated_genes.txt", typ) # The reference of the test data extracted from NCBI

#This block creates the list for the end of the gene predictions of the following files:
typ=typ_options[1]

augustus_e=augustus("augustus.gtf", typ) #Results from Augustus
genemark_e=genemark("genemark.hmm", typ) #Results from genemark
easygene_e=easygene("easygene.gff", typ) #Results from easygene
prodigal_e=prodigal("prodigal.gff", typ) #Results from prodigal
glimmer_e=glimmer("glimmer.gff", typ) #Results from glimmer
annotated_e=features("Annotated_genes.txt", typ) # The reference of the test data extracted from NCBI

#Print the length of each predictions
print len(augustus_b)
print len(genemark_b)
print len(easygene_b)
print len(prodigal_b)
print len(glimmer_b)

# A function to check between each how many data is shared
def SharedData(list1, list2, list3, list4, list5, remove=False):
	"""Given 4 lists check if the first has the same data than  the others. If add is true return the ones not shared if false, return the list of shared data."""
	sets=list1[:]
	least_1=[]

	for position in list1:
		#Check if position is inside any other object.
		if position in list2: 
			least_1.append(position)
			sets.remove(position)
		if position in list3:
			try:	
				sets.remove(position)
				least_1.append(position)
			except:
				pass #In python 2.7 is complusory the except clause... pass is when some other is already delted. 
		if position in list4:
			try:
				sets.remove(position)
				least_1.append(position)
			except:
				pass
		if position in list5:
			try:
				sets.remove(position)
				least_1.append(position)
			except:
				pass

	print len(sets) / float(len(list1)) * 100, "% of genes start are different between the predictors and genemark"
	if remove==True:
		return sets
	elif remove==False:
		return least_1


SharedData(genemark_b, easygene_b, prodigal_b, glimmer_b, augustus_b)

#Some results
## Results of the above code just changing the order (The order of the 2nd, 3rd, 4th, 5th objects doesn't affect the result)
#3.23 % of genes start are different between the predictors and augustus
#71.14 % of genes start are different between the predictors and glimmer
#36.66 % of genes start are different between the predictors and easygene
#22.23 % of genes start are different between the predictors and prodigal
#14.04 % of genes start are different between the predictors and genemark


## Check if the predictors overlap each other.
def Comp(ORP1, stop1, ORP2, stop2):
	"""The ORP1 (Open reading position) and stop1 should come from the same file and program, and have the same length. Given that it checks if any gene of file1 is overlap with other gene"""

	#Defining some variables
	sets=[]
	Shared=0
	Soon_end=0
	Late_end=0
	Soon_begin=0
	Late_begin=0
	Includes=0
	Overlap_begin=0
	Overlap_end=0
	Between=0
	Different=0

	for a in range(0,len(ORP1)):
		# Checking for existing genes even if they are not in the same order!
		# The genes start in the same position as other in the second list
		if ORP1[a] in ORP2: 
			try:
				c=stop2[ORP2.index(ORP1[a])] #Get the stop position in the second list according to the index where it was found the same opening. It ensure that is the same gene what is compared.
			except ValueError:
				pass # This should be addressed in the following block of code
			else:
				if stop1[a]==c:
					Shared+=1
					#print "The gene predicted between positions", ORP1[a], "and", stop1[a], "is shared."
				elif stop1[a]>c:
					Soon_end+=1
					#print "The gene predicted between positions", ORP1[a], "and", stop1[a], "begins at the same position but ends sooner in the second program", stop2[a]
				elif stop1[a]<c:
					Late_end+=1
					#print "The gene predicted between positions", ORP1[a], "and", stop1[a], "begins at the same position but ends later in the second program", stop2[a]

		# If they end in the same position
		elif stop1[a] in stop2:
			try:
				c=ORP2[stop2.index(stop1[a])] #Get the ORP position in the second list according to the index where it was found the same stop. It ensure that is the same gene what is compared.
			except:
				pass # This should 
			else:
				if ORP1[a]<c: # It begins later in the list2
					Soon_begin+=1
					#print "The gene predicted between positions", ORP1[a], "and", stop1[a], "ends at the same position but begins later in", c
				elif ORP1[a]>c: # It begins before in the list2
					Late_begin+=1
					#print "The gene predicted between positions", ORP1[a], "and", stop1[a], "ends at the same position but begins before in", c
				# There is not any equal because it is already solved in the previous block
			
		elif ORP1[a] not in ORP2 and stop1[a] not in stop2: #Not matching any other parameters
			Different+=1

#This was for further classification, but I didn't succed yet
#			for i in range(ORP1[a], stop1[a]): #Check for overlaping intragenic positions
#				if i in ORP2:
#					c=True
#				if i in stop2:
#					e=True
#				try:
#					if c!=None and e!=None:
#						Includes+=1
#						#print "This gene includes one gene" #It could be that one genes starts there and other ends befor and it is not the whole gene included.
#						break
#					elif c==True:
#						Overlap_begin+=1
#						#print "The gene overlaps partially with other gene begining"
#						break
#					elif e==True:
#						Overlap_end+=1
#						#print "The gene overlaps partially with other gene end"
#						break
#				except:
#					pass
#				
#
#			for b in range(1, len(ORP2)):
#				for i in range(stop2[b-1], ORP2[b]): # Check if it is between the end and the new ORF
#					if i==ORP1[a]:
#						e=True
#					if i==stop1[a]:
#						g=True
#					try:
#						if e==True and g==True:
#							Between+=1
#							#print "The gene is between two genes in the second program."							
#							break
#					except:
#						pass

	# Summary of the results
	print "First program predictions", len(ORP1)
	print "Second program predictions", len(ORP2)
	print Shared, "predictions shared"
	print Soon_begin, "predictions begins sooner in the second prediction"
	print Soon_end, "predictions ends sooner in the second prediction"
	print Late_begin, "predictions begin later than the second prediction"
	print Late_end, "predictions end later in the second prediction"
	print Different, "predictions are different...."

	#print Includes, "predictions include predictions of the second program"
	#print Overlap_begin, "predictions overlap with the begining of a prediction in the second program"
	#print Overlap_end, "predictions overlap with the end of a prediction in the second program"
	#print Between, "predictions are between the predictions in the second program."

# Apply function to the glimmer and the reference data.
Comp(glimmer_b, glimmer_e, annotated_b, annotated_e)
# The otput really interesting is the shared between programs or/and data, the ones that are sligthly different, and the ones that are completely different.
# Note that it is not symmetric, and changing the order of the elements to compare affects the result.
# (mainly the Different data, but also in some cases the similar ones -I don't understand why-) 
