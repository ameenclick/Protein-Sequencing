"""
Protein Sequencing Project
Name:
Roll Number:
"""

from fileinput import filename
from opcode import opname
from os import remove

from numpy import amin, math
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f = open(filename, "r")
    text = f.read().replace("\n","")
    f.close()
    return text


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    codon=""
    condons=[]
    for each in dna[startIndex:].replace("T","U"):
        codon += each
        if(len(codon) == 3):
            condons.append(codon)
            if(codon in ["UAA", "UAG", "UGA"]):
                return condons
            codon=""
    return condons


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    f=open(filename, "r")
    aminoToCondons=json.load(f)
    condonsToAmino={}
    for amino in aminoToCondons:
        for condons in aminoToCondons[amino]:
            condons=condons.replace("T","U")
            condonsToAmino[condons]=amino
    return condonsToAmino


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    proteins=[]
    for codon in codons:
            if(proteins == []):
                proteins.append("Start")
            else:
                if(codon in ["UAA", "UAG", "UGA"]):
                    proteins.append("Stop")
                    return proteins
                proteins.append(codonD[codon])
    return proteins


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    readDNA= open(dnaFilename,"r")
    dna=readDNA.read().replace("\n","")
    codonD=makeCodonDictionary(codonFilename)
    proteins=[]
    unused=0
    index=0
    while(index<len(dna)):
        if(dna[index:index+3] == "ATG"):
            rna=dnaToRna(dna,index)
            proteins.append(generateProtein(rna,codonD))
            index += (3*len(rna))
        else:
            index += 1
            unused += 1
    print("\nTotal no. of Bases:",len(dna))
    print("Unused-Base: ",unused)
    print("Total Number Of Protein Synthesized:",len(proteins))
    return proteins


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    common=[]
    for protein in proteinList1:
        if(protein in proteinList2):
            common.append(protein)
    return common


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    aminoacids=[]
    for protein in proteinList:
        for amino in protein:
            aminoacids.append(amino)
    return aminoacids


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    aaSet=set(aaList)
    aaDic={}
    for each in aaSet:
        aaDic[each]=aaList.count(each)
    return aaDic


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    aminoacidsP1=combineProteins(proteinList1)
    aminoacidsP2=combineProteins(proteinList2)
    aaDic1=aminoAcidDictionary(aminoacidsP1)
    aaDic2=aminoAcidDictionary(aminoacidsP2)
    result=[]
    totalAminoP1=len(aminoacidsP1)
    totalAminoP2=len(aminoacidsP2)
    aminoacids=set(aminoacidsP1+aminoacidsP2)
    for amino in aminoacids:
        if(amino in aaDic1):
            freq1=aaDic1[amino]/totalAminoP1
        else:
            freq1=0
        if(amino in aaDic2):
            freq2=aaDic2[amino]/totalAminoP2
        else:
            freq2=0
        if(abs(freq1-freq2)>cutoff and amino != "Start" and amino != "Stop"):
            result.append([amino,freq1,freq2])
    return result


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("\nThe following proteins occurred in both DNA Sequences:")
    allproteins=[]
    for proteins in commonalities:
        proteins.remove("Start")
        proteins.remove("Stop")
        protein= "-".join(proteins)
        allproteins.append(protein)
    allproteins=sorted(set(allproteins))
    for each in allproteins:
        print(each)
    print("\nThe following amino acids occurred at very different rates in the two DNA sequences:")
    for amino in differences:
        print(amino[0]+" : "+str(round((amino[1]*100),2))+"% in Seq1, "+str(round((amino[2]*100),2))+"% in Seq2")
    return None


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    aminoacids1=set(combineProteins(proteinList1))
    aminoacids2=set(combineProteins(proteinList2))
    aminoacids=list(aminoacids1.union(aminoacids2))
    return sorted(aminoacids)


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    aminoacids=combineProteins(proteinList)
    aminoCounts=aminoAcidDictionary(aminoacids)
    total=len(aminoacids)
    freq=[]
    for amino in labels:
        if amino not in aminoCounts:
            freq.append(0)
        else:
            freq.append(aminoCounts[amino]/total)
    return freq


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    w = 0.35  # the width of the bars

    plt.bar(xLabels, freqList1, width=-w, align='edge', label=label1, edgecolor=edgeList)
    plt.bar(xLabels, freqList2, width= w, align='edge', label=label2, edgecolor=edgeList)

    plt.xticks(rotation="horizontal")
    plt.legend()
    plt.title("Comparing Two Gene Frequency")

    plt.show()
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    result=[]
    biggestDiffs=combineProteins(biggestDiffs)
    for amino in labels:
        if amino in biggestDiffs:
            result.append("black")
        else:
            result.append("white")
    return result


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    proteins1=synthesizeProteins("data\human_p53.txt", "data/codon_table.json")
    proteins2=synthesizeProteins("data\elephant_p53.txt", "data/codon_table.json")
    commonalities=commonProteins(proteins1,proteins2)
    differences=findAminoAcidDifferences(proteins1,proteins2,0.005)
    displayTextResults(commonalities,differences)
    labels=makeAminoAcidLabels(proteins1,proteins2)
    freq1=setupChartData(labels,proteins1)
    freq2=setupChartData(labels,proteins2)
    edges=makeEdgeList(labels,differences)
    createChart(labels,freq1,"Human",freq2,"Elephant",edgeList=edges)
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()

    ## Uncomment these for Week 2 ##
    
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()


    ## Uncomment these for Week 3 ##
  
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
