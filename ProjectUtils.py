
#Problems:
	#Missing fraction calcul in function monoisotopicMass
	#in Assembling DNA Fragments, only the example gaved by Rosalind in page 45 took a fast respond, other examples must have lot in common or it will take an enormous time. 
#
###

try:
    import numpy as np
except:
    print("Please install the module by launching this in terminal : sudo pip install numpy")
    
try:
    import matplotlib.pyplot as plt
except:
    print("Please run this command on your shell command : sudo pip install matplotlib")
import random
RNA_codon_dictFullName = {'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys', 'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---', 'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg', 'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', 'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser', 'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', 'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly', 'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly', 'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'}

RNA_codon_dict = {'UUU'
: 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', 'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',
'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-', 'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',
'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', 'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S', 'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', 'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}

proteins_mass = {'A':71.03711, 'C':103.00919, 'D':115.02694, 'E':129.04259, 'F':147.06841, 'G':57.02146, 'H':137.05891, 'I':113.08406, 'K':128.09496, 'L':113.08406, 'M':131.04049, 'N':114.04293, 'P':97.05276, 'Q':128.05858, 'R':156.10111, 'S':87.03203, 'T':101.04768, 'V':99.06841, 'W':186.07931, 'Y':163.06333, '-':0.0}

data = dict()

def randomDNASequence(length):
    """random DNA sequence of desired length.
    
    This function will generate a DNA sequence known by its A_C_T_G acid
    nucleics with a length of wish.
    Enter a number in the field at button's right. 
    Please enter an integer or will return None"""
    
    try:
        length = int(length)
    except ValueError:
        return
    
    sequence = ''
    l = ['A','C','T','G']
    for i in range(length):
        sequence += l[random.randint(0,3)]

    return sequence

def countAcidNucl(sequence, acid=None, isDNA=True, fromFile=False):
    """Counting number of each type of acid nucleic in the DNA or RNA
    sequence.
    
    By passing a sequence and specify the one entered in the text field 
    is a path to a fasta file or a dna sequence.
    
    This will return only A_[TU]_G_C numbers, this will not detect if the
    sequence was invalid. 
    This will not calculate Amino Acids, Protein sequences are not
    welcome.
    ONLY RNA and DNA ARE EVALUATED.
    """
    global data #This is obligatory to change the value of the global variable 'data'
    
    if fromFile:
            sequences = CleanSequence(sequence, fromFile) #This will return a list
            counted = ['']*len(sequences)
            i=0
            for seq in sequences:
                counted[i] = {'A':seq.count("A"), 'C':seq.count("C"), 'G':seq.count("G")}
                if isDNA:
                    counted[i]['T'] = seq.count("T")
                else:
                    counted[i]['U'] = seq.count("U")
                i+=1
            
            data = counted
            return counted
                
    if acid != None:
        return sequence.count(acid)
    #Else give it all
    counted = {'A':sequence.count("A"), 'C':sequence.count("C"), 'G':sequence.count("G")}
    if isDNA:
        counted['T'] = sequence.count("T")
    else:
        counted['U'] = sequence.count("U")
     
    data = counted
    return counted
    
def CleanSequence(seq_or_file, fromFile=False):
    """Delete and sub-sequence to let only the nucliotidic bases.
    
    If you enter in the input field different sequences separated by a '>', 
    the system will consider every begining of '>' a new sequence.
    it will return a list of sequences, with a length of 0, 1 or n elements depends on the number of '>'.
    
    There is no direct call of this function in the system, used by function "sequenceIsValid"."""
    
    sequence = seq_or_file
    if fromFile == True:
        with open(seq_or_file) as _file:	
            sequence = _file.read()
            
    sequences = sequence.split('>') #In case the file contains more than one sequence, separate it
    #i=0
    for i in range(len(sequences)):		
        sequences[i] = sequences[i].split('\n',1)[-1] #split to a list of 2 elements and get the second element
        sequences[i] = sequences[i].replace('\r','')
        sequences[i] = sequences[i].replace('\n','')
     #   i+=1
    return sequences[1:] if sequences[0]=='' else sequences[0:] #To eleminate void element if the sequence is token from a file

def sequenceIsValid(seq_or_file, fromFile=False, isProtein=False, isDNA=True):
    """Checking if the sequence is valid.
    
    This function is for checking if a sequence is valid, no matter if it was a
    DNA, RNA or a protein sequence. 
    
    "seq_or_file" is an argument of the function, it can be a path to a file or
    a string containing the sequence[s]."""
    
    sequenceList = CleanSequence(seq_or_file, fromFile)
    if '' in sequenceList: #if Empty list
        print('Empty')
        return False
    listLen = len(sequenceList)
    valid = [True] * listLen #To get the non valid sequences
    if isProtein:
        values = proteins_mass.values()
        for k in range(listLen):
            for cpt in range(len(seq)):
                if seq[cpt] not in values:
                    valid[k] = False
    else:
        for i in range(listLen):
            d = countAcidNucl(sequenceList[i], isDNA=isDNA)
            count = d["A"] + d["G"] + d["C"]
            if isDNA:
                count += d["T"] 
            else:
                count += d['U']
            if count != len(sequenceList[i]):
#		     	print("ERROR! chaine non valide")
                valid[i] = False
    return valid

def Translate_DNAtoRNA(adnsequence, fromFile=False):
    """Translation of DNA to RNA.
    
    This function must have adnsequence, that it be from a file or not.
    it will not proceed if not.
    
    Argument "adnsequence" is a string which can contain a path to a file or pure
    sequence. To enter different sequences, separate them by a '>'."""
    validList = sequenceIsValid(adnsequence, fromFile)
    if validList == [True] * len(validList): #if all sequences in list are valid
        return [element.replace("T", "U") for element in CleanSequence(adnsequence, fromFile)] #returns a list
    
    print('Invalid sequence(s).', validList)

def Translate_toProtein(dnaORrna_sequence, isDNA=True, fromFile=False): 
    """transcribe of DNA to amino acids.
    
    This function will test if the sequence is valid. If so, it will 
    use a dictionnary variable declared in the Top of "ProjectUtils.py"
    script to change every 3 codons (a key in the dictionnary) with
    what it belong (the value of the key).
    
    Arguement "dnarORrna_sequence" is a string which can contain a path
    to a file or pure sequence. To enter different sequences, separate
    them by a '>'. A sequence can be a DNA or a RNA, please specify 
    by coching."""
    
    if not sequenceIsValid(dnaORrna_sequence, fromFile=fromFile, isDNA=isDNA):
        return	
    if isDNA: #==True
        rnaList = Translate_DNAtoRNA(dnaORrna_sequence, fromFile)
    else:
        rnaList = CleanSequence(dnaORrna_sequence, fromFile)
    index = i = 0
    proteins = ''
    for (i,rnaSequence) in enumerate(rnaList):
        while index+2 < len(rnaSequence):	
            value = rnaSequence[index:index+3]
            proteins += RNA_codon_dict[value]
            index += 3
	# discard last acids nucleics in sequenceif (len(rnaSequence)%3 != 0)
        rnaList[i] = proteins
    return rnaList

def ReverseDNAComplement(adnsequence, fromFile=False): 
    """Reverse complement of DNA.
    
    This function will get you the reverse of every nucliotidic base
    in the sequence based on a dictionnary initialized locally.
    
    Argument "adnsequence" is a string which can contain a path
    to a file or pure sequence. To enter different sequences, separate
    them by a '>'."""
    ref = {'A':'C', 'G':'T', 'T':'G', 'C':'A'}
    if not sequenceIsValid(adnsequence, fromFile):
        return		
    brinadn = CleanSequence(adnsequence, fromFile)
    brinadnListed = [] #Initialize
    for seq in brinadn:
        brinadnListed.append(list(seq))
	
    for (i_list, seq) in enumerate(brinadnListed):
        for i_seq in range(len(seq)):
            seq[i_seq] = ref[seq[i_seq]]
        brinadnListed[i_list] = seq
		
	#returns the reverse of every element-sequence- in the list
    return [''.join(element) for element in brinadnListed] 

def Rate_GC(seq_or_file, fromFile=False):
    """Calculate the number of GC bases in the sequence.
    
    This function will take as an argument a path to a file or string
    containing single sequence or plural."""
    
    global data
    data.clear()
    
    sequences = CleanSequence(seq_or_file, fromFile)
    tauxGC = [0] * len(sequences)
    for i in range(len(sequences)):
        count = countAcidNucl(sequences[i], 'G') + countAcidNucl(sequences[i], 'C')
        data['GC_sequence_{}'.format(i)] = count
        tauxGC[i] = str((float)(count * 100 / len(sequences[i]) )) +  "%"

    return tauxGC

def codonsFrequence(sequence, isDNA=False, fromFile=False):
    """Calculate the frequence of the codons in the sequence[s].
    
    This function will count the number of every codon in the sequence
    passed in the argument.
    
    Only DNA or RNA sequences are accepted.
    Argument "sequence" can be a path to a file or a string. To Enter
    more than one sequence, separate it with '>'."""
    
    global data
    
    if not sequenceIsValid(sequence, fromFile):
        return

    if isDNA: #==True
        sequences = Translate_DNAtoRNA(sequence, fromFile)

    index = i = 0
    codonDict = {}
    seqNumber = 1
    content=''
#    if len(sequences) > 1:
 #       seqNumber = input('Too much sequences in file, enter a seq number:')
#    if seqNumber > len(sequences):
 #       seqNumber = len(sequences) 
  #  if seqNumber <= 0:
   #     seqNumber = 1 	
    #seqNumber-=1 #for indexing the sequence in list
    for seqNumber in range(len(sequences)):
        content += 'Calculating Sequence N={}'.format(seqNumber)
        head = 'AmAcid   Codon     Number        /1000     Fraction'
        content += '\n{}'.format(head)
        
        lines = [] 
        while index+2 < len(sequences[seqNumber]):
            try:
                codonDict[sequences[seqNumber][index:index+3]] += 1 
            except KeyError:
                codonDict[sequences[seqNumber][index:index+3]] = 1
            index += 3

        references = {} #contains the repetition number of every AMINO ACID in the sequence-not codon-
#        totalCodon_foramAcid = []
        i = -1
        for (codon, number) in codonDict.items():	
            lines.append(RNA_codon_dictFullName[codon]) #get amino acid of the codon
            try:
                references[lines[len(lines)-1]] += number
#               totalCodon_foramAcid[i] = number
                
            except KeyError:
                
                references[lines[len(lines)-1]] = number
#                totalCodon_foramAcid.append(number)
                i+=1
                
            lines[len(lines)-1]+='      ' + codon #fill the column Codon
            lines[len(lines)-1]+='         ' + str(float(number)) #fill the column Number
            lines[len(lines)-1]+='         ' + str(float(number/1000)) #fill the column /1000
            #Filling the column Fraction
    #		i=0
    #		for i in range(len(lines)):
    #			lines[i] += '           ' + str(float(
        #Save	
        lines.sort()
        for line in lines:
            content += '\n{}'.format(line)
        
        content += '\n'
        
    data = references
    
    return content

def monoisotopicMass(sequence, isDNA=False, fromFile=False, isProtein=True):
    """Calculate the monoisotopic mass of the sequence.
    
    This function will sum all amino acid weights of the sequence to get
    the total weight.
    
    If the sequence passed in the argument is a DNA or RNA sequence, it
    will translate it first to a proteins sequence. Every single protein
    weight is intialized inside a dictionnary, a variable declared in top
    of the "ProjectUtils.py" script."""

    if not sequenceIsValid(sequence, fromFile):
        print('not valid/')
        return 'not valid/'		

    if isProtein: #==True
        sequences = CleanSequence(sequence, fromFile)
    else:
        sequences = Translate_toProtein(sequence, fromFile=fromFile, isDNA=isDNA)

    totalMass = [0] * len(sequences)
    for i in range(len(sequences)):
        seq = list(sequences[i]) #to split the protein sequence
        for am in seq:
            totalMass[i] += proteins_mass[am]
    return totalMass	

def RNASplicing(adnsequence, introns, fromFile=False): 
    """Splice the DNA sequence basing on list of introns.
    
    delete all sub-sequences of the dnasequence. This sub-sequences are
    gived in the enter parameters (introns) and then, after identifying
    the exons and introns, delete introns and concatenate exons.
    
    Done by using replace function of the class string. 
    
    Only DNA sequences are accepted.
    """
    if not sequenceIsValid(adnsequence, fromFile=fromFile) or not sequenceIsValid(introns, fromFile=False):
        return
    
    introns = CleanSequence(introns, fromFile=False)
    
    introns.sort(reverse=True) #From the highest length to the lowest
    
    sequences = CleanSequence(adnsequence, fromFile)
    
    for i in range(len(sequences)):
        for intron in introns:
            if intron in sequences[i]:
                sequences[i] = sequences[i].replace(intron, '') 
    
    
    return [Translate_toProtein(seq, fromFile=False, isDNA=True) for seq in sequences]

def DNAAssembler(adnsequences, fromFile=False):
    """Assemble DNA sequences.
    
    This function will determine which sequence will be assembled with
    the other first. 
    First measure, will take the sequence in the list of sequence as a base,
    then, will search for the sequence which are similar at least half length
    of the sequence and concatenate them. and so on with others.
    
    The Argument "adnsequence" can be a path to a file or a string but
    must contain more than 1 sequence. Separate the sequences with '>'."""
    
    if fromFile:
        if not sequenceIsValid(adnsequences, fromFile=fromFile, isDNA=True, isProtein=False):
#			print('Not valid sequence')
            return 'Not valid sequence'
    adnsequences = CleanSequence(adnsequences, fromFile=fromFile)
    length = len(adnsequences)	
	
    if length > 50:
        print('Too much sequences')
        return
    for seq in adnsequences:
        if len(seq) > 1000:
            print('Sequence very big')
            return
    assembledSequence = adnsequences[0]
    numbers = list(range(1,length))
    while len(numbers) != 0: #If we have sequences to compare and attach
        seq_1 = assembledSequence
        joinSeq = -1
        j = numbers[0]
        for j in numbers:			
            seq_2 = adnsequences[j]
            k = 1
            while k <= len(seq_2) and seq_2[:k] in seq_1:	
                k+=1
            if k >= len(seq_2)/2:#This the sequence to choose
                    try:
                        if joinIndex < k :
                            joinIndex = k-1
                            joinSeq = j
                    except:
                        joinIndex = k-1
                        if joinIndex < k :
                            joinSeq = j

        #Join
        if joinSeq != -1: #if there is a sequence to attach
            numbers.remove(joinSeq) #remove the attached sequence
            assembledSequence = assembledSequence.replace(adnsequences[joinSeq][:joinIndex], adnsequences[joinSeq])
#	print(assembledSequence)
    return assembledSequence

def create_plot():
    
    global data
    
    names=data.keys()
    size = data.values()
    
    y_pos = np.arange(len(names))
    
    plt.subplot(2, 1, 1)    
    plt.bar(y_pos, size, color=(0.2, 0.4, 0.6, 0.6), tick_label=tuple(names))
    
    plt.subplot(2, 1, 2)
    my_circle = plt.Circle((0,0), 0.7, color='white')
    plt.pie(size, labels=names, colors=['red','green','blue','skyblue'])
    p=plt.gcf()
    p.gca().add_artist(my_circle)
    
    plt.show()
    
    data.clear()

def Application():
    return
#    adnseq = 'ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG'
#    introns = ['ATCGGTCGAA','ATCGGTCGAGCGTGT']
#    print(RNASplicing(adnseq, introns ,fromFile=False))
    
#"bin/Sequence.fas", fromFile=True
#   print(randomDNASequence())
#	print(ReverseDNAComplement("bin/Sequence.fas", fromFile=True))
#	print(Rate_GC(randomDNASequence()))
#	print(Translate_toProtein('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA', isDNA=False))
#	codonsFrequence("bin/Sequence.fas", fromFile=True, isDNA=True)
#	print(monoisotopicMass('SKADYEK', fromFile=False, isProtein=True, isDNA=False))
	#['ATTAGACCTG','CCTGCCGGAA','AGACCTGCCG','GCCGGAATAC']	
	#[randomDNASequence(),randomDNASequence(),randomDNASequence()]
#	l = ['ATTAGACCTG','CCTGCCGGAA','AGACCTGCCG','GCCGGAATAC']	
#	DNAAssembler(l)
	#print 'ATTAGACCTG'.replace('AGACCTG','AGACCTGCCG')
#	print(ReverseDNAComplement(randomDNASequence(), fromFile=False))
    
if __name__ == '__main__':
    Application()
    
    
