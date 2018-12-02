#This is the program for the first part of the assignment

#First we introduce the requirements needed
require 'net/http' 
require 'bio'

#Second I am going to define the file I am going to create later. In this case the gff file.
gff_out = File.new("Genes_features.gff", "w") # "w" because to need to write the entire gff file
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


#Now we are going to examine the sequences of the genes of the previous file by retrieving them from Ensemble creating the next loop
gene_names.each do |gene_ID| 
  #Here, we find the information of the gene and save it in a BioRuby EMBL object. Then, we extract the sequence
  address = URI("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene_ID}")
  response = fetch(address)
  record = response.body
  entry = Bio::EMBL.new(record) 
  bioseq = entry.to_biosequence
  
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
    if /exon_id/.match(qual_string) #Here we select the exons
        unless /[A-Z]/.match(position_string)
          #Transformation of the string [1..250] to [1,250]
          location = position.tr('complement()','')
          location = location.split("..")
          init_position = location[0].to_i #Initial position of the exon
          final_position = location[1].to_i #Final position of the exon
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
  
  # The rest of the code is for visual representation
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
        puts qual_string
        location= position.split("..")
        init_position = location[0].to_i #Initial position and final position
        final_position = location[1].to_i  
        if /"-"/.match(qual_string) # The minus strand
          print "#{gene_ID}\t.\trepeat_region\t#{init_position}\t#{final_position}\t.\t-\t.\t.\n\n"
          gff_out.puts"#{gene_ID}\t.\trepeat_region\t#{init_position}\t#{final_position}\t.\t-\t.\t.\n\n"
        else # The plus strand
          print "#{gene_ID}\t.\trepeat_region\t#{init_position}\t#{final_position}\t.\t+\t.\t.\n\n"
          gff_out.puts "#{gene_ID}\t.\trepeat_region\t#{init_position}\t#{final_position}\t.\t+\t.\t.\n\n"
        end
      end
    end
  end
  puts "\n\n"
end

#Part 4b of the task. We generate an output where we show the genes that do NOT have exons with the CTTCTT repeat
puts "#{no_repeat.length} genes do not have exons with cttctt repeat. These genes are: #{no_repeat} "