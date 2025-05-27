#5.25.23

#Versions
#Python 3.9.5
#Pandas 1.3.2
#Et-xml-file 1.1.0
#Requests 2.26.0

#Sources
#https://python-reference.readthedocs.io/en/latest/intro.html
#https://docs.python.org/3.11/index.html
#https://docs.python-guide.org/scenarios/scrape/
#https://pandas.pydata.org/docs/getting_started/index.html 
#https://requests.readthedocs.io/en/latest/
#https://flybase.github.io/docs/api/

#Imports
import pandas as pd

import xml.etree.ElementTree as ET
import requests
import time as time

#Initialize variables

#Read in an Excel file of deficencies to analyze
#The file should contain a list of deficencies formatted as FBab using https://flybase.org/convert/id
flies = pd.read_excel("/Users/johncsantiago/Documents/GitHub/WhartonLab/Python/DeficiencyData/renewapp/AllModifiers.xlsx", header=None)

#Number of deficencies
rows = len(flies.index)

#List of FBgns
genesFBgn = []

#List of gene names
genesName = []

#List of df for gene names
genesdf = []

#List of deficencies
dfs = []

#Initialize the list of deficencies ignoring empty cells
for i in range(rows):
    if str(flies.iat[i, 0]) != "nan" and str(flies.iat[i, 0]) != "":
        dfs.append(flies.iat[i, 0])

#Main for loop to make call to the API and extract data  
for df in dfs:
    names = []
    FBgns = []
    currentdf = []

    #Call to Flybase API
    request = requests.get("https://api.flybase.org/api/v1.0/chadoxml/"+str(df))
    
    #Import XML Data
    root = ET.fromstring(request.content)

    #Traversing the XML tree
    relationships = root.findall("./feature/feature_relationship")
    for relationship in relationships:

        #Check if relationship is a gene disrupted by deficiency (derived_computed_affected_gene)
        isgene = False

        for type_id in relationship:
            if type_id.tag=="type_id":
                 for cvterm in type_id:
                    if cvterm.tag=="cvterm":
                        for name in cvterm:
                            if name.tag == "name":
                                if name.text == "derived_computed_affected_gene":
                                    isgene = True

        #Traverse XML Tree and obtain FBgn and Name
        if isgene:
            for subject_id in relationship:
                if subject_id.tag=="subject_id":
                    for feature in subject_id:
                        if feature.tag=="feature":
                            for feature_child in feature:
                                if feature_child.tag == "uniquename":
                                    FBgns.append(feature_child.text)

                                if feature_child.tag=="feature_synonym":
                                    for synonym_id in feature_child:
                                        if synonym_id.tag =="synonym_id":
                                            for synonym in synonym_id:
                                                if synonym.tag=="synonym":
                                                    for name in synonym:
                                                        if name.tag =="name":
                                                            names.append(name.text)
                                                            
                                                        if name.tag =="name":
                                                            currentdf.append(str(df))
    #print(df, len(FBgns), len(currentdf), str(df))
    #remove duplicates
    genesFBgn+=list(FBgns)
    genesName+=list(names)
    genesdf+=list(currentdf)

    print(df, len(FBgns), len(currentdf), len(genesdf), len(genesName))

    #Delay per Flybase API request limit
    time.sleep(0.5)

#Format output as pair of FBgn and gene name removing duplicates again
final = list(map(lambda x, y, z: [x,y,z], list(genesFBgn), list(genesName), genesdf))
#final = list(map(lambda x, y: [x,y], list(set(genesFBgn)), list(set(genesName))))

#Create dateframe of output
output = pd.DataFrame(final)

#Export dataframe to Excel file
output.to_excel("Output.xlsx")
