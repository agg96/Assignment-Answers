#This is the program for the second part of the assignment
#The structure is similar to the first part, but here we are going to include the full chromosome coordinates

#First we introduce the requirements needed
require 'net/http' 
require 'bio'

#Second I am going to define the file I am going to create later. In this case the gff file.
gff_out = File.new("Chromosome_features.gff", "w") # "w" because to need to write the entire gff file
no_repeat = [] #Here, we will introduce the genes that do NOT have exons with the CTTCTT repeat

#Now we will define the functions we are going to use

#First function --> used to open uris and get the responses we want
def fetch(uri_str)    
  address = URI(uri_str)  
  response = Net::HTTP.get_response(address)
  case response   # with case we are testing various conditions... it is like an "if", but cleaner!
    when Net::HTTPSuccess then  # when response is of type Net::HTTPSuccess
      return response  # return that response object
    else
      raise Exception, "Something went wrong... the call to #{uri_str} failed; type #{response.class}"
      response = False
      return response  # now we are returning False
  end
end

#Second function --> used to find the target sequences
#With this function we are going to look for three main things:
#the targets located in the EXONS,
#the INITIAL POSITION of the exons to know the position of the sequence in the gene
#and the target sequence, which is the repeat motif CTTCTT.

def find_target(exon,init_position,target_motif) #seq exon, init position of the exon, cttctt or aagaag 
  target_motif = target_motif.downcase #We transform the sequence to lowercase
  target_motif_length = target_motif.length
  
  #Now we make an array with all the positions of the target sequences. Have in mind that a sequence cttcttctt counts as 2 
  target_motif_position_in_exon = (0 ... exon.length).find_all { |i| exon[i,target_motif_length] == target_motif }
  
  #IMPORTANT PART --> #Inside the exon the relative position could be 50-56, but out of the exon could be 250-256 as the exon has its own position in the genome
  #With this in mind, we generate a loop to adapt the positions of the cttctt found in an exon to the real position out of the exon
  exon_position = []
  target_motif_position_in_exon.each do |x|
    out = []
    first_position = x + init_position + 1 # The initial position is 1, not 0
    final_position = first_position + target_motif_length - 1 #The final position is one position before of the resulting sum
    out.push first_position
    out.push final_position 
    exon_position.push out 
  end
  return exon_position # Here, we will have all of the positions of the cttctt in one exon
end

#Third function --> used to transform [0,5]--> [1..6]
def transform(general_position)
    final_general_position =[]
    general_position.each do |y|
      y[0].to_s #to_s to transform it into a string
      y[1].to_s
      string = "#{y[0]}..#{y[1]}"
      final_general_position.push string
    end
    return final_general_position
end

#Once defined the functions needed, we start to write the code

#First, we open the file with the gene names and save it in an array
gene_list = File.open('ArabidopsisSubNetwork_GeneList.txt', 'r')
gene_names = gene_list.read.split() # this will read each line into an array
gene_list.close

#Until here, everything is similar to the first task, but the loop has to change in this case.

#Now we are going to examine the sequences of the genes of the previous file by retrieving them from Ensemble creating the next loop
gene_names.each do |gene_ID| 
  #Here, we find the information of the gene and save it in a BioRuby EMBL object. Then, we extract the sequence
  address = URI("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene_ID}")
  response = fetch(address)
  record = response.body
  entry = Bio::EMBL.new(record) 
  bioseq = entry.to_biosequence

  #Here, we retrieve the information of the chromosomes --> DIFFERENT PART
  chromosome= []
  chromosome=(entry.ac[0]).split(":") #Separation of the different features 
  print "#{entry.ac[0]}\n" #Here we see the information we have for the chromosomes
  chromosome_name = chromosome[1]
  chromosome_number = chromosome[2]
  chromosome_initial_position = chromosome[3]
  chromosome_final_position = chromosome[4]
  
  #Now again, the code is similar to the first task until the gff file format that we generate, which has to be different
  #We create a loop with the exons to find information in the object that tell us the position of the exon and other features of it
  general_position_plus = []
  general_position_minus = []
  target_motif ="CTTCTT"
  entry.features.each do |feature| 
    qual = feature.assoc
    position = feature.position
    qual_string = qual.to_s #BOTH ARE HASH, so we have to transform them to string
    position_string = position.to_s

    #Now we select specific parts --> the exons of the gene in both plus (+) and minus (-) strand
    if /exon_id/.match(qual_string)
        unless /[A-Z]/.match(position_string)
          #Transformation of the string [1..250] to [1,250]
          pos_aux = position.tr('complement()','')
          pos_aux = pos_aux.split("..")
          init_position = pos_aux[0].to_i #Initial position of the exon
          final_position =pos_aux[1].to_i #Final position of the exon
          #Now here, if the exon is in the minus strand, the sequence we look for have to be the complementary one to CTTCTT
          #We could have done this defining which is the complementary sequence, but we can use the next code as well.
          if /complement/.match(position_string) 
            bioseq_rev = bioseq.tr('acgt','tgca')
            bioseq_rev = bioseq_rev.reverse! #Complementary sequence
            exon = bioseq_rev [init_position..final_position]
            strand = "minus" #Tell us that is in the minus strand
          else # If not
            exon = bioseq [init_position..final_position]
            strand = "plus" #Tell us that is in the plus strand
          end
          #We have the sequence of the exon 
          #Now we are going to find the sites of the motif CTTCTT in the exon with the function find_target
          exon_position = find_target(exon,init_position,target_motif) 
          if strand == "plus"
            exon_position.each do |i|
              general_position_plus.push i
            end
          else
            exon_position.each do |j|
              general_position_minus.push j
            end
          end
          general_position_plus = general_position_plus.uniq #Positions in the plus strand  
          general_position_minus = general_position_minus.uniq# Positions in the minus strand 
        end 
    end
  end
    
   
  if general_position_plus == [] and general_position_minus == []
    puts "The gene #{gene_ID} do NOT have exons with the target_motif repeat\n\n"
    no_repeat.push gene_ID
  else
    puts "Running gene #{gene_ID}"
    puts "Target sequence positions in the exons of the plus strand are: #{general_position_plus}"
    puts "Target sequence positions in the exons of the minus strand are: #{general_position_minus}\n"

    #Use of the transform function to change [1,6] to [1..6] 
    final_general_position_plus = transform(general_position_plus)
    final_general_position_minus = transform(general_position_minus)
    
    #We add the features we want to the objets 
    final_general_position_plus.each do |z|
      f1 = Bio::Feature.new('sequence',z)
      f1.append(Bio::Feature::Qualifier.new('repeat_motif', 'CTTCTT'))
      f1.append(Bio::Feature::Qualifier.new('strand', '+'))
      bioseq.features << f1 
    end
    final_general_position_minus.each do |k|
      f2 = Bio::Feature.new('sequence',k)
      f2.append(Bio::Feature::Qualifier.new('repeat_motif', 'CTTCTT'))
      f2.append(Bio::Feature::Qualifier.new('strand', '-'))
      bioseq.features << f2 
    end
    
    #Finally, we create the gff format
    entry.features.each do |m|
      qual = m.assoc
      position = m.position
      qual_string = qual.to_s
      if /repeat_motif/.match(qual_string) 
        location = position.split("..")
        init_position = location[0].to_i #Initial position and final position
        final_position = location[1].to_i
        
        #Now, we have to define the position of the CTTCTT repeats in the full lenght chromosomes
        
        #In the minus sense
      
        #The initial position would be: The final position of the gene in the chromosome - the final position of the repeat sequence in the gene
        #The initial position of the cttctt in the chr- is the final position of the gen in the chr - initial position of the cttctt in the gen (alfa) + 
        if /"-"/.match(qual_string)
          chromosome_initial_pos_minus = chromosome_final_position.to_i - final_position + 1
          chromosome_final_pos_minus = chromosome_initial_pos_minus + 5 
          print "#{chromosome_name}:#{chromosome_number}\t#{gene_ID}\trepeat_region\t #{chromosome_initial_pos_minus}\t #{chromosome_final_pos_minus} \t-\t.\t.\n\n"
          gff_out.puts "#{chromosome_number}\t.\trepeat_region\t#{chromosome_initial_pos_minus}\t#{chromosome_final_pos_minus}\t.\t-\t.\t.\n"
          
        #In the plus sense
        
        #The initial position would be:The initial position of the gene in the chromosome + the initial position of the repeat in the gene
        #On the other hand, the final position would be: The initial position + 5, because our sequence has a lenght of 6 and if we count the edges [1,2,3,4,5,6] has a lenght of 6. 
        else
          chromosome_initial_pos_plus = chromosome_initial_position.to_i + init_position - 1 
          chromosome_final_pos_plus = chromosome_initial_pos_plus + 5
          print "#{chromosome_name}:#{chromosome_number}\t#{gene_ID}\trepeat_region\t #{chromosome_initial_pos_plus}\t #{chromosome_final_pos_plus} \t+\t.\t.\n\n"
          gff_out.puts "#{chromosome_number}\t.\trepeat_region\t#{chromosome_initial_pos_plus}\t#{chromosome_final_pos_plus}\t.\t+\t.\t.\n"
        end
      end
    end
  end 
end
