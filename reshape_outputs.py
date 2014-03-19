#!/usr/bin/env python

# This code is somewhat similar to a parser, but it just when the user got already the page of the results in the computer.
# It get two list for each wraper, the begin and the end 

import re



#Function working with glimmer! But it is not completely refined, it needs further improvements of the regular expression to use the raw glimmer file. 
def glimmer(files, typ):
	"""Given an output of glimmer it finds the gene positions begin or end.
Requires just the table. Further improvements must be done to read raw glimmer files"""
	
	# Read, open and get into more readable format
	data=open(files)
	info=data.readlines()
	begin=[]
	end=[]
	regex=re.compile('\s+')

	for line in info:
		genes=re.match("orf.+", line)
		if genes:
			gene=re.split(regex, line)
			begin.append(gene[1])
			end.append(gene[2])
			
	# Given the type of information it is wanted is returned
	if typ=="begin":
		return begin

	elif typ=="end":
		return end

# Function to parser from the output of prodigal.
def prodigal(files, typ="begin"):
	"""Given a prodigal output files finds the begining and the end of each gene."""
	data=open(files)
	info=data.read()
	begin=[]
	end=[]
	info=[info][0].split("\n")
	regex=re.compile("\s+")
	for line in info:
		genes=re.match("[\s]+CDS[\s]+.+", line)
		if genes:
			genes=genes.group(0)
			gene=re.split(regex, genes)
			gene=gene[2].split("..")
			try:
				new_gene=gene[0].split("complement(")[1]
				begin.append(new_gene)
			except:
				new_gene=gene[0]
				begin.append(new_gene)
			try:
				new_end=gene[1].split(")")[0]
				end.append(new_end)
			except:
				new_end=end[0]
				end.append(end[0])

	if typ=="begin":
		return begin

	elif typ=="end":
		return end

# Function to parse the EasyGene program
def easygene(files, typ):
	"""Extract the begin and the end of each gene of the easygene output data."""

	data=open(files)
	info=data.read()
	begin=[]
	end=[]
	info=[info][0].split("\n")
	for line in info:
		genes=re.match("^gi", line)
		if genes:
			gene=line.split("\t")
			new_gene=gene[3]
			begin.append(new_gene)
			new_end=gene[4]
			end.append(new_end)

	if typ=="begin":
		return begin
		
	elif typ=="end":
		return end

# Function to parse the GenemMark program
def genemark(files, typ):
	"""Extract from a genemark file the begin and end of each gene."""

	data=open(files)
	info=data.read()
	begin=[]
	end=[]
	info=[info][0].split("\n")
	regex=re.compile("\s+")
	for line in info:
		genes=re.match("^[\s]+[0-9]+[\s]+.+", line)
		if genes:
			genes=genes.group(0)
			gene=re.split(regex, genes)
			new_gene=gene[3]
			begin.append(new_gene)
			new_end=gene[4]
			end.append(new_end)

	if typ=="begin":
		return(begin)
		
	elif typ=="end":
		return end

# Function to parser the Augustus program
def augustus(files, typ):
	"""Print the begin or the end of the gene lines given the output of the augustus program."""
	data=open(files)
	info=data.read()
	begin=[]
	end=[]
	info=[info][0].split("\n")
	regex=re.compile("\s+")
	for line in info:
		genes=re.match("[\s+\w+\W+]+gene[\s+\w+\W+]+", line)
		if genes:
			genes=genes.group(0)
			gene=re.split(regex, genes)
			if gene[2]=="gene":
				new_gene=gene[3]
				begin.append(new_gene)
				new_end=gene[4]
				end.append(new_end)

	if typ=="begin":
		return begin
		
	elif typ=="end":
		return end

#This function given data of the NCBI genes in multi Fasta format get the beginning and end of the genes. 
def features(files, typ):
	"""Given a file of NCBI gene features in FASTA format it returns a list of the start codon position and the stop position, according to typ(begin, end)"""
	begin=[]
	end=[]

	data=open(files)
	info=data.read()
	regex=re.compile("\d+\.\.\d+")
	genes=re.findall(regex, info) #Select just the positions
	regex2=re.compile("\.\.")

	for positions in genes:
		gene=re.split(regex2, positions) #To split into the two numbers
		begin.append(gene[0])
		end.append(gene[1])
	
	if typ=="begin":
		return(begin)
		
	elif typ=="end":
		return end
