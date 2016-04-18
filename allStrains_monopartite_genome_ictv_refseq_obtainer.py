#!/usr/bin/env python
#python 3.1
#selenium 2.47.1
#Firefox 40.0.3
#09/14/2015
#Ryan O. Schenck
#XXXXX

#Objective: Obtain RefSeq ID's for full proteome of ICTV classified viruses of interest
#This script is designed for monopartite viruses
#Manual curation is necessary to ensure completeness.

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
import linecache

#user input (file of species), family name
arg1 = sys.argv[1]
arg2 = sys.argv[2]

fileOutLog = "output/" + arg2 + "_genomeGetter_log.txt"
fileOut = "output/" + arg2 + "_genomeGetter.txt"

#set up logfile
logging.basicConfig(filename=fileOutLog, filemode='w', level=logging.INFO)
log = logging.getLogger("")

outFile = open(fileOut, 'w')

#main body of script
def main(arg):
	log.info("NCBI Function Executed")
	
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

	#begins webdriver
	sampleProg = 0
	with open(arg, 'r') as fileIn:
		for line in fileIn:
			line = line.rstrip('\n')
			split = line.split('\t')
			if split[0].startswith('Order'):
				outFile.write(line + '\n')
			else:
				#prints progress
				sampleProg = sampleProg + 1
				percent =  repr(sampleProg) + '/' + repr(sampleCount) + ' samples processed'
				print(percent)
				
				#takes the species we are interested in
				species = split[4]
				log.info(species)
				try:
					#try to open webdriver and get result file
					driver = webdriver.Chrome('/Users/ryan/Desktop/IDP_Analysis_Program/chromedriver')
					driver.get('http://www.ncbi.nlm.nih.gov')

					#select Taxonomy database
					driver.find_element_by_xpath("//select[@id='database']/option[text()='Taxonomy']").click()

					#places refseq ID in search field and submits info
					driver.find_element_by_name('term').send_keys(species)
					time.sleep(1.5)
					
					driver.find_element_by_name('term').send_keys(Keys.ENTER)
					time.sleep(1.5)
				except:
					log.info('Failed to connect to NCBI.')
				
				'''#clicks on the taxonomy record for that species'''
				try:
					#set up string to put in driver.find_element executable
					species = species.lower()
					species = species.capitalize()

					#clicks on the taxonomy report
					driver.find_element_by_link_text(species).click()
					time.sleep(1.5)

					log.info('Taxonomy Record Found: ')
				except Exception as e:
					log.exception('NO TAXONOMY RECORD FOUND: ')
					log.error(e)

				'''#obtains taxonomic ids and builds required search information'''
				try:
					taxid, nucNum = one_taxid(driver)

					log.info("Sending to one_taxid: ")
					log.info("TaxID obtained: ")
				except:
					try:
						taxid, nucNum = multi_taxid(driver)
						log.info("Sending to multi_taxid: ")
						log.info("TaxID obtained: ")
						
					except Exception as e:
						log.info("No Taxonomic ID obtained")
						taxid = 'txid"NONE!!!!!!!!!!!!!!!!!!!!!!!!!!"[Organism:exp]AND "complete genome"[Title]'
						log.error(e)

				'''#searches for genomes using search field'''
				try:
					#try to open webdriver and get result file
					#driver2 = webdriver.Chrome('/Users/ryan/Desktop/IDP_Analysis_Program/chromedriver')
					driver.get('http://www.ncbi.nlm.nih.gov')

					#select Nucleotide database
					driver.find_element_by_xpath("//select[@id='database']/option[text()='Nucleotide']").click()

					#places refseq ID in search field and submits info
					driver.find_element_by_name('term').send_keys(taxid)
					time.sleep(1.5)
					driver.find_element_by_name('term').send_keys(Keys.ENTER)
					time.sleep(1.5)
					
					#clicks the 'Send:' button in NCBI site
					driver.find_element_by_link_text('Send to:').click()
					time.sleep(0.5)
					
					#selects the coding sequence option
					driver.find_element_by_xpath("//input[@value='File']").click()
					time.sleep(1)

					#selects create file
					driver.find_element_by_xpath("(//button[@cmd='File'])[position()=1]").click()
					if int(nucNum) < 50:
						time.sleep(4)
					elif int(nucNum) > 1000:
						time.sleep(20)
					else:
						time.sleep(6)
					
					log.info('Summary file downloaded: ')

					driver.quit()
				except Exception as c:
					log.info("Summary file NOT downloaded: ")
					log.error(c)
					driver.quit()

				'''Parses out accession number and strain name from summary file'''
				try:
					#location of downloaded summary file from above
					summary = '/Users/ryan/Downloads/nuccore_result.txt'
					#pass file location to accessGet
					
					#accessGet(summary)
					accessDict = accessGet(summary)

					#unpack dictionary strain name is item
					for item in accessDict:
						genAcc = accessDict[item]
						genAcc = genAcc.split('=rem=')
						genomeLength = genAcc[0]
						genomeAcc = genAcc[1]
						item = item.split('. ', 1)
						strainNam = item[1]
						if 'NC_' not in genomeAcc:
							if 'AC_' not in genomeAcc:
								outFile.write(line + '\t' + 'monopartite' + '\t' + strainNam + '\t' + genomeLength + '\t' + genomeAcc + '\n')
							else: 
								pass
						else:
							pass
					log.info('Strain accession numbers obtained: ')
				except Exception as d:
					outFile.write(line + '\t' + 'monopartite' + '\n')
					log.info('Unable to use accessGet(): ')
					log.error(d)



#for when there's only one strain found under species
def one_taxid(driver):
	#get to taxonomy summary page
	nuc = driver.find_element_by_xpath("(//span[@style='color:#0000A0'])[position()=2]")
	#get number of nucleotide results
	nucNum = nuc.text
	if ',' in nucNum:
		nucNum = nucNum.replace(',','')
	else:
		pass
	time.sleep(2)

	#get taxid
	mytaxid = driver.find_element_by_xpath('//td/input[@name="old_id"]').get_attribute('value')
	
	#get number of nucleotide records


	#what to add to search bar
	append = 'txid' + mytaxid + '[Organism:exp]' 'AND "complete genome"[Title]'
	
	#driver.quit()
	
	return(append, nucNum)

#gets tax id if multiple strains shown under species
def multi_taxid(driver):
	#get to taxonomy summary page
	driver.find_element_by_xpath("(//strong)[position()=2]").click()
	time.sleep(2)

	nuc = driver.find_element_by_xpath("(//span[@style='color:#0000A0'])[position()=2]")
	#get number of nucleotide results
	nucNum = nuc.text
	if ',' in nucNum:
		nucNum = nucNum.replace(',','')
	else:
		pass
	
	time.sleep(2)

	#get taxid
	mytaxid = driver.find_element_by_xpath('//td/input[@name="old_id"]').get_attribute('value')



	#what to add to search bar
	append = 'txid' + mytaxid + '[Organism:exp]' 'AND "complete genome"[Title]'

	#driver.quit()

	return(append, nucNum)

#gets the strain id + genome length + accession number
def accessGet(summary):
	#take each sample into a list with the sample name and line number separated by '=rem=' 
	with open(summary, 'r') as summaryFile:
		strainList = list()
		for k, line in enumerate(summaryFile):
			if line == '\n':
				pass
			elif 'complete genome' in line:
				line = line.rstrip('\n')
				out = line + '=rem=' + repr(k)
				strainList.append(out)
			else:
				pass
		sampleTot = len(strainList)

	accessDict = OrderedDict()

	#take the sample name and line numbers for everything in summary file
	for strain in strainList:
		#obtains strain name, genome length line number, and accession number line number
		strain = strain.split('=rem=')
		strainName = strain[0]
		sizeLine = int(strain[1]) + 2
		accessLine = int(strain[1]) + 3
		
		#puts downloaded file in cache and accesses those lines
		genSizeLine = linecache.getline(summary, sizeLine).rstrip('\n')
		accessionNumLine = linecache.getline(summary, accessLine).rstrip('\n')

		#splits and takes the proper portion of the line
		genSize = genSizeLine.split(' bp')
		accessionNum = accessionNumLine.split(' ')

		inform = genSize[0] + '=rem=' + accessionNum[0]
		
		accessDict.update({strainName:inform})

	#clear cache
	linecache.clearcache()

	#delete downloaded file that information extracted
	os.remove(summary)

	return accessDict

main(arg1)

outFile.close()
print('\nComplete.')

log.info('\nProcess Complete.\n')

logging.shutdown()



















