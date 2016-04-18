#!/usr/bin/env python
#python 2.7
#selenium 2.47.1
#Firefox 40.0.3
#09/14/2015
#Ryan O. Schenck
#XXXXXXXX

#Objective: Obtain fasta ID's for proteome in master spreadsheet tables for data analysis

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
outName = "/Users/ryan/Desktop/IDP_Analysis_Program/output/" + arg2 + "_proteome_accessions.txt"
finalOut = "/Users/ryan/Desktop/IDP_Analysis_Program/output/" + arg2 + "_proteome_fractions.txt"

#set up logfile
logging.basicConfig(filename=fileOutLog, filemode='w', level=logging.INFO)
log = logging.getLogger("")

#setting up environment for outputs
outFile = open(outName, 'w')

#gets a list of accession numbers to get the batch entrez files
def proteome_split(arg1):
	log.info('Extracting proteome accession numbers.')

	accession_list = list()

	with open(arg1, 'r') as master:
		for line in master:
			if line.startswith('#') or line.startswith('Order'):
				pass
			else:
				line = line.rstrip('\n')
				line = line.split('\t', 10)
				#gets the proteome into a list including '.' and ''
				proteome = line[10].split('\t')

				#separates the components of the proteome
				#accession_list = list()
				for item in proteome:
					if item == '.':
						pass
					elif item == '':
						pass
					else:
						accession_list.append(item)
						outFile.write(item + '\n')
	log.info('Protein accession numbers extracted.')
	return accession_list

#to be used if you want to get a fasta file
'''
#gets a fasta file into downloads folder with all of the protein sequences
def fasta_getter(myFile):
	try:
		#open browser
		driver = webdriver.Chrome('/Users/ryan/Desktop/IDP_Analysis_Program/chromedriver')
		driver.get('http://www.ncbi.nlm.nih.gov/sites/batchentrez')
		time.sleep(1.5)
		#select protein database
		driver.find_element_by_xpath("//select[@name='db']/option[text()='Protein']").click()
		#places file into the download field
		inputFile = driver.find_element_by_name('file')
		inputFile.send_keys(myFile)
		time.sleep(2)
		
		#click retrieve button
		driver.find_element_by_xpath("//input[@value='Retrieve']").click()
		time.sleep(30)
		#time.sleep(20)
		#clicks the link to retrieve
		driver.find_element_by_xpath('//a').click()
		time.sleep(10)

		#clicks send to:
		driver.find_element_by_link_text('Send to:').click()
		time.sleep(1)
		#selects 'File'
		driver.find_element_by_xpath("//label[@for='dest_File']").click()
		time.sleep(1)
		#selects 'Fasta' from Format dropdown menu:
		driver.find_element_by_xpath("//select[@id='file_format']/option[text()='FASTA']").click()
		time.sleep(1)
		#clicks create File
		driver.find_element_by_xpath("//button[@cmd='File']").click()
		time.sleep(30)

		driver.quit()
		log.info('Fasta file downloaded.')
	except Exception as a:
		log.error('Failed at fasta_getter.')
		log.error(a)
		try:
			driver.quit()
		except:
			pass
'''
'''
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
'''

#gets the disorder scores from pondr website
def get_score(access):
	try:
		driver = webdriver.Chrome('/Users/ryan/Desktop/IDP_Analysis_Program/chromedriver')
		driver.get('http://www.pondr.com/cgi-bin/PONDR/pondr.cgi')
		
		#unclicks and clicks that appropriate fields
		driver.find_element_by_xpath("//input[@name='graphic']").click()
		driver.find_element_by_xpath("//input[@name='stats']").click()
		driver.find_element_by_xpath("//input[@name='seq']").click()
		driver.find_element_by_xpath("//input[@name='wcwraw']").click()
		time.sleep(2)

		#puts fasta file into text box
		driver.find_element_by_xpath("//input[@name='AccessionCode']").send_keys(access)
		time.sleep(5)
		
		#gets the data from website
		data = driver.find_element_by_xpath('html').text

		#convert data from unicode to UTF-8 for handling as a string
		data = data.encode(encoding='UTF-8', errors='strict')

		#puts data from website into a temporary file
		with open('/Users/ryan/Desktop/IDP_Analysis_Program/output/delete_me.txt', 'w') as temp:
			temp.write(data)

		with open('/Users/ryan/Desktop/IDP_Analysis_Program/output/delete_me.txt', 'r') as inFile:
			with open('/Users/ryan/Desktop/IDP_Analysis_Program/output/IDP_Scores.txt', 'w') as out_Final:
				for line in inFile:
					line = line.rstrip('\n')
					lineParts = line.split(' ')
					#captures the scores
					if len(lineParts) == 3:
						out_Final.write(line + '\n')
					else:
						pass
		
		os.remove('/Users/ryan/Desktop/IDP_Analysis_Program/output/delete_me.txt')
		#close the webdriver
		driver.quit()
	except Exception as c:
		log.error('Failure in PONDR website.')
		log.error(c)
		try:
			driver.quit()
		except:
			pass

#gets the proportion of disorder and returns
def get_proportion():
	with open('/Users/ryan/Desktop/IDP_Analysis_Program/output/IDP_Scores.txt', 'r') as scores:
		# get length of protein
		length = 0
		for line in scores:
			if line != '':
				length += 1

	#get fraction of disorder Pushkar et. al. 2013
	with open('/Users/ryan/Desktop/IDP_Analysis_Program/output/IDP_Scores.txt', 'r') as scores:
		sums = 0
		for line in scores:
			line = line.rstrip('\n')
			line = line.split(' ')
			score = float(line[2])
			
			if score > 0.5:
				sums += 1

	fraction = float(sums)/float(length)

	os.remove("/Users/ryan/Desktop/IDP_Analysis_Program/output/IDP_Scores.txt")

	return fraction

#===============================#
def main(arg1):
	#obtains a file and list of protein fasta files for batch entrez download
	protein_list = proteome_split(arg1)
	outFile.close()
	print('Number of Proteins: ' + repr(len(protein_list)))
	log.info('Number of Proteins: ' + repr(len(protein_list)))

	#sends the created list of protein accession numbers to NCBI batch entrez
	#fasta_getter(outName)
	
	#gets the proportion of disorder in the protein for each protein
	#open the file with all accession numbers for proteins
	with open(finalOut, 'w') as FinalOutput:
		with open(outName, 'r') as forPondr:
			for line in forPondr:
				line.rstrip('\n')
				try:
					#gets the scores file from pondr website
					get_score(line)

					#gets the fraction of disorder for that protein
					disorder = get_proportion()

					FinalOutput.write(line.rstrip('\n') + '\t' + repr(disorder) + '\n')
				except Exception as fail:
					#If a fail occurs the most likely reason is that the sequence is less than 30 AA in length
					log.error("Failed to get score: " + line)
					log.error(fail)
					FinalOutput.write(line.rstrip('\n') + '\t' + '.' + '\n')

	os.remove(outName)

	#cleanup fasta file
	#os.remove('/Users/ryan/Downloads/sequence.fasta')

#===============================#

main(arg1)

#outFile.close()
logging.shutdown()













