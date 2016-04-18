#!/usr/bin/env python
#python 3.1
#selenium 2.47.1
#Firefox 40.0.3
#09/14/2015
#Ryan O. Schenck
#XXXXXXX

#Objective: Obtain RefSeq ID's for full proteome of ICTV classified viruses of interest

#requires an input file in a tabular separated format with Family, Genus, Species, Genome_refseq ID, extra(after this it doesn't matter)

import os
import sys
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
import time
import os.path
import re
from collections import OrderedDict
import logging
import collections

#input file of genome accession numbers and descriptive name
arg1 = sys.argv[1]
arg2 = sys.argv[2]

fileOutLog = "output/" + arg2 + "_proteome_log.txt"
fileOut = 'output/' + arg2 + '_refseq_output.txt'

#set up logfile
logging.basicConfig(filename=fileOutLog, filemode='w', level=logging.INFO)
log = logging.getLogger("")

#setting up environment for outputs
outFile = open(fileOut, 'w')

#fasta parser, returns name and seq
def read_fasta(arg):
    name, seq = None, []
    for line in arg:
        line = line.rstrip()
        if line.startswith('>'):
            line = line.split('>')
            line = line[1]
            if name: yield(name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield(name,''.join(seq))

#outFile must be opened before executing
#webdriver and data extraction
def ncbi(arg):
	#counts number of samples
	with open(arg, 'r') as countme:
		sampleCount = 0
		for line in countme:
			line = line.rstrip('\n')
			split = line.split('\t')
			if split[0].startswith('Order'):
				pass
			else:
				sampleCount += 1

	#begin main section
	sampleProg = 0
	masterDict = dict()
	with open(arg, 'r') as fileIn:
		for line in fileIn:
			#progress indicator
			sampleProg = sampleProg + 1
			#percent =  repr(sampleProg) + '/' + repr(sampleCount) + ' samples processed'
			print(repr(sampleProg) + '/' + repr(sampleCount) + ' samples processed')

			#split line to obtain genome ID
			line = line.rstrip('\n')
			split = line.split('\t')
			if split[0].startswith('#'):
				pass
			else:
				genomeID = split[9]

				#try to open webdriver and get result file
				try:
					#open browser
					driver = webdriver.Chrome('/Users/ryan/Desktop/IDP_Analysis_Program/chromedriver')
					driver.get('http://www.ncbi.nlm.nih.gov')
					#select nucleotide database
					driver.find_element_by_xpath("//select[@id='database']/option[text()='Nucleotide']").click()
					#places refseq ID in search field and submits info
					driver.find_element_by_name('term').send_keys(genomeID)
					time.sleep(0.5)
					driver.find_element_by_name('term').send_keys(Keys.ENTER)
					time.sleep(1)
					#clicks the 'Send:' button in NCBI site
					driver.find_element_by_link_text('Send:').click()
					time.sleep(0.5)
					#selects the coding sequence option
					driver.find_element_by_xpath("//input[@value='codeseq']").click()
					time.sleep(0.5)
					#selects FASTA protein option from dropdown
					driver.find_element_by_xpath("//select[@id='codeseq_format']/option[text()='FASTA Protein']").click()
					time.sleep(0.5)
					#selects creat file
					driver.find_element_by_xpath("//button[@cmd='codeseq']").click()
					time.sleep(2)

					driver.quit()
				except:
					log.error('Failed in Webdriver: ' + genomeID)
					try:
						driver.quit()
					except:
						pass

				#code block to check if file was obtained from NCBI
				try:
					#need to have download location for your computer specified here
					if os.path.isfile('/Users/ryan/Downloads/sequence.txt'):
						#return a dictionary for this virus with each protein and its refseqID
						proteome = refseq_fetch('/Users/ryan/Downloads/sequence.txt')

						#adds the proteome dictionary into a new dictionary with genome refseqID as key
						masterDict.update({genomeID:proteome})
						#print(masterDict)

					else:
						log.error("Download Failed: " + genomeID)
				except Exception as e:
					log.error('Failed in refseq parser: ' + genomeID)
					log.error(e)
					try:
						os.remove('/Users/ryan/Downloads/sequence.txt')
					except:
						pass

	return masterDict

#extracts fields from fasta seq headers
def extract(s):
	start = 0
	end = s.find( ']' )
	result = s[start:end]
	return result

#extracts necessary component for NCBI_product_check
def extract2(s):
	end = s.find('"ORIGIN')
	start = end - 40
	result = s[start:end]
	return result

#check to make sure that unique values are added to dictionary in refseq_fetch otherwise you lose protein
def refseq_dictCheck(speciesDict, v1, itemCount):
	for item in speciesDict:
		if v1 in item:
			v1 = repr(itemCount) + '_' + v1
		else:
			v1 = v1
	return v1

#retrieves refseq ID for each protein in downloaded file from ncbi function
def refseq_fetch(argx):
	speciesDict = dict()
	k1 = '1212'
	v2 = '1212'
	itemCount = 0
	with open(argx, 'r') as fp:
		for name, seq in read_fasta(fp):
			title = name.replace(" ","_")
			title = title.replace("[",":")
			title = title.split(':')
			itemCount += 1
			for item in title:
				if 'protein=' in item:
					gene = extract(item)
					k1,v1 = gene.split('=')
					v1 = refseq_dictCheck(speciesDict, v1, itemCount)
				if 'protein_id=' in item:
					refseq_id = extract(item)
					k2,v2 = refseq_id.split('=')
				else:
					pass
			if '1212' not in k1 and '1212' not in v2:
				speciesDict.update({v1:v2})
			
	#deletes the downloaded file
	os.remove(argx)			
	#print(speciesDict)
	return speciesDict

#for unnamed protein in a CDS file: checks to see what the protein is based on known protein names in "note" field.
def ncbi_product_check(value, name):
	try:
		driver2 = webdriver.Chrome('/Users/ryan/Desktop/IDP_Analysis_Program/chromedriver')
		driver2.get('http://www.ncbi.nlm.nih.gov')
		#select nucleotide database
		driver2.find_element_by_xpath("//select[@id='database']/option[text()='Protein']").click()
		#places refseq ID in search field and submits info
		driver2.find_element_by_name('term').send_keys(value)
		time.sleep(0.5)
		driver2.find_element_by_name('term').send_keys(Keys.ENTER)
		time.sleep(3.5)
		
		#gets the unique identifier for the page or something
		unique = driver2.find_element_by_xpath('//meta[@name="ncbi_uidlist"]')
		unique = unique.get_attribute('content')
		#convert data from unicode to UTF-8 for handling as a string
		data = unique.encode(encoding='UTF-8', errors='strict')

		#gets the text that contains the 'note='
		myID = '"feature_' + data + '_CDS_0"]'
		stringToGo = "'//span[@id=" + myID + "'"
		body = driver2.find_element_by_xpath("//html").text
		data = body.encode(encoding='UTF-8', errors='strict')
		data = ''.join(data.split())
		myName = extract2(data)
		myName = myName.split('/note="')
		protName = myName[1]
		
		driver2.quit()
	except Exception as kf:
		log.error('Unable to perform protein check: ' + name)
		log.error(kf)
		protName = ''
		try:
			driver2.quit()
		except:
			pass

	protName = protName.lower()

	return protName

#extracts ref_seq ID's from dictionary for the virus for each virus_order function
#####ssDNA#####
def circo_order(dictin, line):
	#reformat and split the line
	line = line.rstrip('\n')
	line = line.split('\t')
	genomeID = line[9]
	lineout = '\t'.join(line)
	#extract proper key, value from ncbi function output

	outProt = OrderedDict()
	hype = list()
	for key, value in dictin[genomeID].iteritems():
		key = key.lower()
		
		#works if it has the keys (lowercase)
		if 'vp1' in key:
			outProt.update({1:value})
		elif 'vp2' in key:
			outProt.update({2:value})
		else:
			hype.append(value)
	
	#checks to make sure that items are present in the ordered dict to be added to string
	try:
		l = outProt[1]
	except:
		outProt.update({1:'.'})
	try:
		d = outProt[2]
	except:
		outProt.update({2:'.'})

	outLine = lineout + '\t' + outProt[1] + '\t' + outProt[2] + '\t' + '\t'.join(hype) + '\n'
	
	return outLine
def anello_order(dictin, line):
	#reformat and split the line
	line = line.rstrip('\n')
	line = line.split('\t')
	genomeID = line[9]
	lineout = '\t'.join(line)
	#extract proper key, value from ncbi function output

	outProt = OrderedDict()
	hype = list()
	for key, value in dictin[genomeID].iteritems():
		key = key.lower()
		#if key in ['orf1','orf2','orf3','orf4','orf1/1','orf2/2']:
		#	pass
		if 'unnamed_protein_product' in key:
			protName = ncbi_product_check(value, genomeID)
			if protName == '':
				pass
			else:
				key = protName
		elif 'hypothetical_protein' in key:
			protName = ncbi_product_check(value, genomeID)
			if protName == '':
				pass
			else:
				key = protName
		else:
			pass
		#works if it has the keys (lowercase)
		if 'orf1/1' in key:
			hype.append(value)
		elif 'orf2/2' in key:
			hype.append(value)
		elif 'orf1' in key:
			outProt.update({1:value})
		elif 'orf2' in key:
			outProt.update({2:value})
		elif 'orf3' in key:
			outProt.update({3:value})
		elif 'orf4' in key:
			outProt.update({4:value})
		else:
			hype.append(value)
	
	#checks to make sure that items are present in the ordered dict to be added to string
	try:
		l = outProt[1]
	except:
		outProt.update({1:'.'})
	try:
		d = outProt[2]
	except:
		outProt.update({2:'.'})
	try:
		l = outProt[3]
	except:
		outProt.update({3:'.'})
	try:
		d = outProt[4]
	except:
		outProt.update({4:'.'})

	outLine = lineout + '\t' + outProt[1] + '\t' + outProt[2] + '\t' + outProt[3] + '\t'  + outProt[4] + '\t' + '\t'.join(hype) + '\n'
	
	return outLine
def parvo_aveparvovirus(dictin, line):
	#reformat and split the line
	line = line.rstrip('\n')
	line = line.split('\t')
	genomeID = line[9]
	lineout = '\t'.join(line)
	#extract proper key, value from ncbi function output

	outProt = OrderedDict()
	hype = list()
	for key, value in dictin[genomeID].iteritems():
		key = key.lower()
		
		#works if it has the keys (lowercase)
		if 'vp1' in key:
			outProt.update({1:value})
		elif 'vp2' in key:
			outProt.update({2:value})
		elif 'ns1' in key:
			outProt.update({3:value})
		elif 'np1' in key:
			outProt.update({4:value})
		else:
			hype.append(value)
	
	#checks to make sure that items are present in the ordered dict to be added to string
	try:
		l = outProt[1]
	except:
		outProt.update({1:'.'})
	try:
		d = outProt[2]
	except:
		outProt.update({2:'.'})
	try:
		c = outProt[3]
	except:
		outProt.update({3:'.'})
	try:
		k = outProt[4]
	except:
		outProt.update({4:'.'})

	outLine = lineout + '\t' + outProt[1] + '\t' + outProt[2] + '\t' + outProt[3] +'\t' + outProt[4] + '\t' + '\t'.join(hype) + '\n'
	
	return outLine
def parvo_copiparvovirus(dictin, line):
	#reformat and split the line
	line = line.rstrip('\n')
	line = line.split('\t')
	genomeID = line[9]
	lineout = '\t'.join(line)
	#extract proper key, value from ncbi function output

	outProt = OrderedDict()
	hype = list()
	for key, value in dictin[genomeID].iteritems():
		key = key.lower()
		
		#works if it has the keys (lowercase)
		if 'capsid' in key:
			outProt.update({1:value})
		elif 'vp2' in key:
			outProt.update({2:value})
		elif key in ['non-structural', 'replicase']:
			outProt.update({3:value})
		else:
			hype.append(value)
	
	#checks to make sure that items are present in the ordered dict to be added to string
	try:
		l = outProt[1]
	except:
		outProt.update({1:'.'})
	try:
		d = outProt[2]
	except:
		outProt.update({2:'.'})
	try:
		c = outProt[3]
	except:
		outProt.update({3:'.'})

	outLine = lineout + '\t' + outProt[1] + '\t' + outProt[2] + '\t' + outProt[3] + '\t' + '\t'.join(hype) + '\n'
	
	return outLine
def parvo_dependoparvovirus(dictin, line):
	#reformat and split the line
	line = line.rstrip('\n')
	line = line.split('\t')
	genomeID = line[9]
	lineout = '\t'.join(line)
	#extract proper key, value from ncbi function output

	outProt = OrderedDict()
	hype = list()
	for key, value in dictin[genomeID].iteritems():
		key = key.lower()
		
		#works if it has the keys (lowercase)
		if 'vp1' in key:
			outProt.update({1:value})
		elif 'vp2' in key:
			outProt.update({2:value})
		elif 'vp3' in key:
			outProt.update({3:value})
		elif key in ['ns1', 'rep78', 'rep']:
			outProt.update({4:value})
		elif 'ns2' in key:
			outProt.update({5:value})
		else:
			hype.append(value)
	
	#checks to make sure that items are present in the ordered dict to be added to string
	try:
		l = outProt[1]
	except:
		outProt.update({1:'.'})
	try:
		d = outProt[2]
	except:
		outProt.update({2:'.'})
	try:
		c = outProt[3]
	except:
		outProt.update({3:'.'})
	try:
		k = outProt[4]
	except:
		outProt.update({4:'.'})
	try:
		v = outProt[5]
	except:
		outProt.update({5:'.'})

	outLine = lineout + '\t' + outProt[1] + '\t' + outProt[2] + '\t' + outProt[3] +'\t' + outProt[4] + '\t' + outProt[5] + '\t' + '.' + '\t' + '.' + '\t' + '.' + '\t' + '\t'.join(hype) + '\n'
	
	return outLine
def parvo_ungulateprotoparvovirus(dictin, line):
	#reformat and split the line
	line = line.rstrip('\n')
	line = line.split('\t')
	genomeID = line[9]
	lineout = '\t'.join(line)
	#extract proper key, value from ncbi function output

	outProt = OrderedDict()
	hype = list()
	for key, value in dictin[genomeID].iteritems():
		key = key.lower()
		
		#works if it has the keys (lowercase)
		if 'capsid_protein_1' in key:
			outProt.update({1:value})
		elif 'capsid_protein_2' in key:
			outProt.update({2:value})
		elif 'vp3' in key:
			outProt.update({3:value})
		elif key in ['non-structure_protein_1', 'ns1', 'nonstructural_protein_1']:
			outProt.update({4:value})
		elif 'non-structure_protein_2' in key:
			outProt.update({5:value})
		elif 'non-structure_protein_3' in key:
			outProt.update({6:value})
		else:
			hype.append(value)
	
	#checks to make sure that items are present in the ordered dict to be added to string
	try:
		l = outProt[1]
	except:
		outProt.update({1:'.'})
	try:
		d = outProt[2]
	except:
		outProt.update({2:'.'})
	try:
		c = outProt[3]
	except:
		outProt.update({3:'.'})
	try:
		k = outProt[4]
	except:
		outProt.update({4:'.'})
	try:
		v = outProt[5]
	except:
		outProt.update({5:'.'})
	try:
		v = outProt[6]
	except:
		outProt.update({6:'.'})

	outLine = lineout + '\t' + outProt[1] + '\t' + outProt[2] + '\t' + outProt[3] +'\t' + outProt[4] + '\t' + outProt[5] + '\t' + outProt[6] + '\t' + '.' + '\t' + '.' + '\t' + '\t'.join(hype) + '\n'
	
	return outLine
def parvo_tetraparvovirus(dictin, line):
	#reformat and split the line
	line = line.rstrip('\n')
	line = line.split('\t')
	genomeID = line[9]
	lineout = '\t'.join(line)
	#extract proper key, value from ncbi function output

	outProt = OrderedDict()
	hype = list()
	for key, value in dictin[genomeID].iteritems():
		key = key.lower()
		
		#works if it has the keys (lowercase)
		if 'vp1' in key:
			outProt.update({1:value})
		elif 'vp2' in key:
			outProt.update({2:value})
		elif 'vp3' in key:
			outProt.update({3:value})
		elif 'ns1' in key:
			outProt.update({4:value})
		elif 'ns2' in key:
			outProt.update({5:value})
		elif 'ns3' in key:
			outProt.update({6:value})
		else:
			hype.append(value)
	
	#checks to make sure that items are present in the ordered dict to be added to string
	try:
		l = outProt[1]
	except:
		outProt.update({1:'.'})
	try:
		d = outProt[2]
	except:
		outProt.update({2:'.'})
	try:
		c = outProt[3]
	except:
		outProt.update({3:'.'})
	try:
		k = outProt[4]
	except:
		outProt.update({4:'.'})
	try:
		v = outProt[5]
	except:
		outProt.update({5:'.'})
	try:
		v = outProt[6]
	except:
		outProt.update({6:'.'})

	outLine = lineout + '\t' + outProt[1] + '\t' + outProt[2] + '\t' + outProt[3] +'\t' + outProt[4] + '\t' + outProt[5] + '\t' + outProt[6] + '\t' + '.' + '\t' + '.' + '\t' + '\t'.join(hype) + '\n'
	
	return outLine
def parvo_tetraparvovirus2(dictin, line):
	#reformat and split the line
	line = line.rstrip('\n')
	line = line.split('\t')
	genomeID = line[9]
	lineout = '\t'.join(line)
	#extract proper key, value from ncbi function output

	outProt = OrderedDict()
	hype = list()
	for key, value in dictin[genomeID].iteritems():
		key = key.lower()
		
		#works if it has the keys (lowercase)
		if 'minor_structural_protein' in key:
			outProt.update({1:value})
		elif 'major_structural_protein' in key:
			outProt.update({2:value})
		elif 'vp3' in key:
			outProt.update({3:value})
		elif 'non-structural_protein' in key:
			outProt.update({4:value})
		elif 'ns2' in key:
			outProt.update({5:value})
		elif 'ns3' in key:
			outProt.update({6:value})
		else:
			hype.append(value)
	
	#checks to make sure that items are present in the ordered dict to be added to string
	try:
		l = outProt[1]
	except:
		outProt.update({1:'.'})
	try:
		d = outProt[2]
	except:
		outProt.update({2:'.'})
	try:
		c = outProt[3]
	except:
		outProt.update({3:'.'})
	try:
		k = outProt[4]
	except:
		outProt.update({4:'.'})
	try:
		v = outProt[5]
	except:
		outProt.update({5:'.'})
	try:
		v = outProt[6]
	except:
		outProt.update({6:'.'})

	outLine = lineout + '\t' + outProt[1] + '\t' + outProt[2] + '\t' + outProt[3] +'\t' + outProt[4] + '\t' + outProt[5] + '\t' + outProt[6] + '\t' + '.' + '\t' + '.' + '\t' + '\t'.join(hype) + '\n'
	
	return outLine
#####dsDNA####To be completed


#catches any duplicates and doesn't let the script run if duplicates in it
def dupCatcher(arg1):
	access = list()
	with open(arg1, 'r') as data:
		for line in data:
			line = line.rstrip('\n')
			line = line.split('\t')
			line = line[9]
			access.append(line)
	
	dups = [item for item, count in collections.Counter(access).items() if count > 1]
	return dups

####################################################
def main(sendToNCBI):
	#outFile = open('refseq_output.txt', 'w')

	outDict = ncbi(sendToNCBI)

	#puts output into file with proper proteome.
	log.info('Verifying Protein names and spreadsheet placement.')
	print('Verifying Protein names and spreadsheet placement.')
	with open (arg1) as gen:
		for line in gen:
			if line.startswith('#'):
				outFile.write(line)
			else:
				try:
					#pass the outDict from NCBI to the proper function for the viral family and the line from the file:
					#print(outDict)
					outLine = parvo_tetraparvovirus2(outDict, line)
					outFile.write(outLine)
				except Exception as b:
					log.error('Failed in parsing proteome')
					log.error(b)
					outFile.write(line)

	log.info("Process Complete.")
####################################################

#stops if duplicates
if len(dupCatcher(arg1)) > 0:
	print('Duplicates found in input:')
	print(dupCatcher(arg1))
else:
	log.info('No Duplicates found in input.')
	log.info('Begin proteome collection')
	main(arg1)





outFile.close()
logging.shutdown()



