#This is the main program for the assignment 4: "Searching for orthologues"

#First of all we call the functions we need
require 'net/http'
require 'bio'

#Functions

#In this case we are only going to create one function

#This function will take the data from fasta format and will generate a matrix with the name and sequences of each gene
#The matrix will have the next structure: [[name_gene1, sequence_gene1], [name_gene2, sequence_gene2,...]]
#This structure will help me later in the realization of the blast
#We have to keep in mind that each gene starts with a > in fasta format

def give_me_name_and_seq(gene)
    protein_group = []
    genelist = gene.read.split(">") #Here is where we separate every gene
    genelist.each do |p|
        proteins = []
        prot = p.split("|") #Here is where we separate every protein of each gene
        protein_name = prot[0]
        if protein_name != nil
            protein_name = protein_name.delete(' ')
            protein_name = protein_name.tr('|', '')
            protein_sequence = p.split("\n")
            protein_sequence.shift #We remove the first element
            protein_sequence = protein_sequence.join #We do not want spaces among aminoacids
            proteins.push protein_name 
            proteins.push protein_sequence #Now we have the name and sequence of each protein
        end
        protein_group.push proteins  #This is the final matrix that we wanted to create
    end
    return protein_group
end

#Here, we are going to create a file to save the orthologs beetween the Arabidopsis and S.Pombe
orthologs = File.new("orthologs.txt", "w")

#Now, we define the data files we are going to use and we are going to create the databases we want by using the two last command lines
#of the next lines of code. We need to do it to make a blast, because blast is going to make a search on these databases.
DNA_data = 'TAIR10_DNA.fa'
aa_data = 'pep_aa.fa'
system "makeblastdb -in #{DNA_data} -dbtype nucl -out #{DNA_data}"
system "makeblastdb -in #{aa_data} -dbtype prot -out #{aa_data}"

#Use of the function give_me_name_and_seq. 
ARAT_hash = {}
ARAT_proteome = File.open('TAIR10_DNA.fa', 'r')
ARA_final = give_me_name_and_seq(ARAT_proteome)
ARA_final.shift
ARA_final.each do |unit|
    ARAT_hash[unit[0]] = unit[1] #We create this hash of Arabidopsis genes which contains: Gene_name => Gene_sequence
end
S_pombe_proteome = File.open('pep_aa.fa', 'r')
S_pombe_final = give_me_name_and_seq(S_pombe_proteome)
S_pombe_final.shift

#Finally, we create the BLAST: 1) a tblastn (aa --> DNA) and  2) a blastx (DNA --> aa)
###BE CAREFUL ON THE DATABASES WE USE EACH TIME### ---> Remember you had problems with this at first
factory_DNA = Bio::Blast.local('tblastn', "#{DNA_data}")
factory_aa = Bio::Blast.local('blastx', "#{aa_data}")
$stderr.puts "Searching ... Be patient, you're Oompa Loompa workers and doing everything for for you"

#Ending of the program

#We finally create a loop in which we first do the tBLASTn (aa --> DNA) and select the best hit (the first of the list) if it has a e-value lower than 1e-10.
#Then, with the best hit, we do the second BLAST: blastx (DNA --> aa) and see if the best hit, with an e-value lower than 1e-10is the same gene we used in the previous BLAST.
#Only if this occurs, both genes are added to the "orthologs.txt" file and we will see them in an output.

amount = 0
S_pombe_final.each do |sequence|
  #I have chosen as treshold a e-value of 1e-10.
  #The reason for this is because I first look at the blast information, arriving to the book called "BLAST"
  #The citation is the next: Korf, I., Yandell, M., & Bedell, J. (2003). Blast. " O'Reilly Media, Inc."
  #The book explains that is generally the case that the most similar genes between species are orthologs, and that this is used as an operational definition
  #In order to do that, it is explained later that to look for orthologs we can filter is the E-value as an initial parameter.
  #We can decide its value but it gives an example where it uses the 1e-10 value at first approximation. That is why I chosen it. Not too restrictive but enough.
  treshold1 = 1e-10
  treshold2 = 1e-10 #we need to define the treshold twice because then we are going to change its value
  puts "Using the aa sequence of S.pombe #{sequence[0]}\n" #Puts the gene that we are testing in aa sequence from S.pombe genome
  result1 = factory_DNA.query(">myseq\n#{sequence[1]}") #The result have all the hits obtein in the tblastn
  #Now we create a loop to see each of the hits obtained.
  result1.each do |res1| 
    if res1.evalue < treshold1 # We only consider that something is a hit if has an evalue greater than 1e-10
      gene1 = res1.definition #This is a long array that includes the name    
      gene_name1 = gene1.split(" | ")[0] # Selecting the first part we have the name of the gene with the best hit in Arabidopsis. This gen have a DNA sequence.
      treshold1 = res1.evalue #We define the treshold again and if we have a hit with a lower e-value we will select it as the new best hit.
      
      #Now the second part with the blastx is similar to the first one with tblastn
      puts "Using the #{gene_name1} gene of Arabidopsis\n"
      target = ARAT_hash[gene_name1] # With the hash created, we can find now the sequence we want to use 
      result2 = factory_aa.query(">myseq\n#{target}") #This result have all the hits obtein in the tblastn: DNA sequence vs aa sequence 
      result2.each do |res2|
        if res2.evalue < treshold2 #same treshold as the first treshold1
          gene2 = res2.definition
          gene_name2 = gene2.split("|")[0]
          print "Returning #{gene_name2} gene of S.pombe\n"
          treshold2 = res2.evalue
          if sequence[0] == gene_name2 #In order to check if they are orthologs 
            puts "The genes #{gene_name2} and #{gene_name1} are orthologs"
            orthologs.print "The genes #{gene_name2} and #{gene_name1} are orthologs\n" # This is to write the document
            amount += 1
            puts amount
          end
        end
      end
    end
  end
  puts "\n\n"
end


#BONUS
#In this case, to prove that the orthologs we have identified are truly orthologs, we will have to investigate FOR more than the homology
#between sequences. We can also study the function of the genes or the proccesses in which are involved. 

#It will also be interesting to do a search of phylogeny of the sequences and see IF they diverge from an event of speciation (orthologs)
#or from a duplication (not truly orthologs - paralogs indeed).The BLAST book also added that a study of the synteny will be interesting.





