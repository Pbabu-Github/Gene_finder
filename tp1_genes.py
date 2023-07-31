"""
Name: Prashanth Babu

"""




def is_stop(segment):
    l=["TAG","TAA","TGA"]
    if segment.upper()in l:
        return True
    else:
        return False
"""
List l stores the three different types of dna segment to know where to stop
it checks if the value of the dna is in the list
"""
def orf_sequence(dna):
    orf=""
    for i in range (0,len(dna),3):
        dna_seg=dna[i:(i+3)]
        if(is_stop(dna_seg)!=True):
                orf+=dna_seg
        elif(is_stop(dna_seg)!=False):
            break
    return orf
"""
orf is the string that finds the orf assuming the dna starts with the value ATG
it calls the is_stop function to decide wether it will continue or stop adding the next three segmenst of the DNA to the orf
"""
def find_orfs(dna):
    orf=[]
    i=0
    while i<=len(dna):
        dna_seg=dna[i:(i+3)]
        if dna_seg=="ATG":
            orf.append(orf_sequence(dna[i:]))
            i+=len(orf[-1])+3
        else:
            i+=3
    return orf
"""
orf= is a list that stores all the differnt orf's. The value 'i' is the index as to where the splicing has to start.
we add the( length of the orfreturned from the find_orf funtcion + 3) because we want to continue searching for orf's after that value.
The orf is added to the list if the dna segent starts with "ATG" otherwise it skips three positions and starts with three indexes after that
"""
def reverse_complement(dna):
    reverse_dna=""
    reverse_dna_complement=""
    for i in dna:
        reverse_dna=i+reverse_dna
    rev=["A","C","G","T"]
    for i in reverse_dna:
        reverse_dna_complement+=rev[(len(rev)-1)-(rev.index(i))]
    return reverse_dna_complement
"""
Accepts the string dna and reverses it and assigns that string to revers_dna .
Reverse_dna complement switches the different dna values according to their complement and the reversed dna complement is returned
"""
def gene_finder(dna):
    orf=[]
    for i in range (3):
     dna1=dna[i:]
     if (find_orfs(dna1))!=[]:
      orf+=find_orfs(dna1)
      
    reverse_dna=reverse_complement(dna)
    for i in range(3):
        reverse_dna1=reverse_dna[i:]
        if (find_orfs(reverse_dna1))!=[]:
         orf+=find_orfs(reverse_dna1)
    return orf
"""
List orf collects all the values of the orf's in the lis
dna1 takes away the first character of the string dna each time for loop runs.
The if loop checks if the list is empty and if it is not, it just adds the orf to the list
the second for loop appends all the orf values of the reversed string and takes away the first value of eachstrin
"""



def read_fasta(filename):
    """
    Read a single DNA sequence from a FASTA-formatted file
    
    For example, to read the sequence from a file named "X73525.fasta.txt"
    >>> sequence = read_fasta("X73525.fasta.txt")
    
    Args:
        filename: Filename as a string
        
    Returns: Upper case DNA sequence as a string
    """
    with open(filename, "r") as file:
        # Read (and ignore) header line
        header = file.readline()
        # Read sequence lines
        sequence = ""
        for line in file:
            sequence += line.strip()
        return sequence.upper()

def filter_orfs(orfs, min_length):
    """
    Filter ORFs to have a minimum length
    
    Args:
        orfs: List of candidate ORF sequences
        min_length: Minimum length for an ORF
    
    Returns:
        A new list containing only ORF strings longer than min_length bases
    """
    filtered_orfs = []
    for orf in orfs:
        if len(orf) > min_length:
            filtered_orfs.append(orf)
    return filtered_orfs


def write_fasta(filename, orfs):
    """
    Write list of ORFs to a FASTA formatted text file.
    
    For example, to save a list of orfs assigned to the variable my_orfs to a
    file named "genes.txt"
    >>> write_fasta("genes.txt", my_orfs)
    
    Args:
        filename: Filename as a string. Note that any existing file with this name
            will be overwritten.
        orfs: List of ORF sequences to write to the file
    """
    with open(filename, "w") as file:
        for i in range(len(orfs)):
            # A FASTA entry is a header line that begins with a '>', and then the sequence on the next line(s)
            print(">seq" + str(i), file=file)
            print(orfs[i], file=file)

           
if __name__ == "__main__":
    # If you want perform the steps for creating genes.txt as part of your
    # program, implement that code here. It will run every time you click the
    # green arrow
    pass
