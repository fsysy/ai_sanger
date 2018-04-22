from Bio import Entrez, SeqIO
Entrez.email = "fsysy@naver.com"

reference_gene_NM_number = "NM_002381.4"

def id_search(NM_number):
    handle = Entrez.esearch(db='nucleotide', term=NM_number)
    record = Entrez.read(handle)
    idlist = record["IdList"]
    if len(idlist) == 1:
        print "The id is %s"%idlist[0]
        return idlist[0]
    else:
        if len(idlist) == 0:
            print "No search NM_number in genbank. please search other names"
            return False
        else:
            print "Multple ids were found! %s"%s(",".join(idlist))
            return False

def get_sequence(Id):
    handle = Entrez.efetch(db="nucleotide", id=Id,rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    print record



gene_id = id_search("NM_002381.4")

if gene_id:
    get_sequence(gene_id)


