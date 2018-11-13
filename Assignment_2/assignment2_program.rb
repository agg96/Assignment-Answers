#This is the Ruby document you have to run in order to obtein the results of the assignment. 

#First, I call the objects and the functions we need in the rest of the code

require './Newgene.rb'
require './InteractionNetwork.rb'
require 'net/http'
require 'json'

#Second, and as I did in the previous assignment, I am going to create different matrix to save the data we want.

gene_name = [] # Here, we will have the names of the genes contained in the ArabidopsisSubNetwork_GeneList.txt
gene_interaction = []# Names of the genes that interact with each other
gene_position = [] #Save the position of the genes that interacts with each other
no_gene_interaction = [] #Name of the genes that do not interact with another of the list.

protein_name_matrix = []# Name of all the proteins 
no_protein_synonyms = []# Here, we have all the protein names excluding synonym names
no_protein_interaction = []# Name of proteins from the genes that do not interact with another of the list
protein_interaction = []# Names of the proteins from the genes that interact with another protein from the genes of the list

#Creation of two matrix in order to save the name of the interactors of the proteins from the genes of the list in one matrix,
#and also the name of the interactors of these previous interactors in another matrix.

interaction_level2 = [] #A group of arrays. Each array is composed by a main protein from the list and its interactors
interaction_level3 = [] #A group of arrays. Each array is composed by a main protein from the list, its interactors and the interactors of the previous interactors

networks = []
networks_genes = []
networks_object =[]
objects_interaction =[]

gene_ontology = [] # Save the data from the gene onthology (GO)
kegg_pathways = [] # Save the data from the Kegg pathway


#We open the archive that  contains the names of the genes and save it in gene_name
File.open ('ArabidopsisSubNetwork_GeneList.txt') do |f1|
    while datas =f1.gets
        datas = datas.gsub("\n", "") #gsub (pattern, replacement) --> new string #The pattern is /n in this case because each gene name is in a new line
        gene_name.push datas    
    end
end

#To obtein the information required for each gene, we create a loop:

gene_name.each do |i|
    #The code of the protein of the gen
    code_address = URI("http://togows.org/entry/ebi-uniprot/#{i}/accessions.json")
    res_code = Net::HTTP.get_response(code_address)
    protein_name = JSON.parse(res_code.body) #I wrote above protein_name_matrix to differ it
    protein_name_matrix.push protein_name
    
    #The information from de gene onthology (GO)
    go_address = URI("http://togows.org/entry/ebi-uniprot/#{i}/dr.json")
    res_go = Net::HTTP.get_response(go_address)
    data = JSON.parse(res_go.body)
    go_data = data[0]["GO"]
    gene_ontology.push go_data
    
    #The information from the Kegg pathway
    #I had some problems here because Json does not work. That is the reason of the creation of kegg_field 
    kegg_address = URI("http://togows.org/entry/kegg-enzyme/ath:#{i}") 
    res_kegg = Net::HTTP.get_response(kegg_address)
    kegg_data = res_kegg.body
    kegg_field = kegg_data.split("\n")
    kegg_path = "unknown Kegg"
    kegg_field.each do |j| #I made this loop because I could not obtein the information propperly at first.
        first_lenght = j.size
        substitution = j.sub(/PATHWAY     /,"")
        if first_lenght != substitution.size
            kegg_path = substitution
        end
    end
    kegg_pathways.push kegg_path 
end

# We can observe that genes can have more than one protein name for the same protein codificated
# We are going save the only the main one in no_protein_synonyms using the next loop and knowing that the first element taken is the main name.

(0..protein_name_matrix.size-1).each do |k|
    no_protein_synonyms.push protein_name_matrix[k][0][0]
end

# Then, we can see the protein interactions with other proteins from the list

(0..no_protein_synonyms.size-1).each do |l|
    require 'net/http'
    require 'json'
    address = URI("http://togows.org/entry/ebi-uniprot/#{no_protein_synonyms[l]}/dr.json")
    res = Net::HTTP.get_response(address)
    data = JSON.parse(res.body)
    intAct = data[0]["IntAct"]
    if intAct
        gene_interaction.push gene_name[l]
        gene_position.push l
        data[0]["IntAct"].each do |m|
        protein_interaction.push m[0]
        end
    else
        no_protein_interaction.push no_protein_synonyms[l]
        no_gene_interaction.push gene_name[l]
    end
end

#Now, we have a list of proteins that have interactions

#We are going to search for the proteins which interact with them.
#With this purpose we are going to create a matrix with the second and third levels of interactions.

#Second level of interactions
protein_interaction.each do |n|
    address = URI("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{n}")
    res_intact2 = Net::HTTP.get_response(address)
    body_intact2 = res_intact2.body
    intact2 = body_intact2.split("\n")
    interaction = []
    (0..intact2.size-1).each do |o|
        columns2 = intact2[o].split("\t")
        protein_A = columns2[0].sub(/uniprotkb:/,"")
        protein_B = columns2[1].sub(/uniprotkb:/,"")
        interaction.push protein_A    
        interaction.push protein_B   
    end
    interaction = interaction.uniq
    interaction_level2.push interaction
end

#Third level of interactions
interaction_level2.each do |array|
    interaction= []
    array.each do |p|
        address = URI("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{p}")
        res_intact3 = Net::HTTP.get_response(address)
        body_intact3 = res_intact3.body
        intact3 = body_intact3.split("\n")
        (0..intact3.size-1).each do |q|
            columns3 = intact3[q].split("\t")
            protein_A = columns3[0].sub(/uniprotkb:/,"")
            protein_B = columns3[1].sub(/uniprotkb:/,"")
            interaction.push protein_A    
            interaction.push protein_B   
        end
        interaction = interaction.uniq  
    end
    interaction_level3.push interaction  
end

# Now, we have a group of arrays composed by the main protein and its interactors of second and third level inside the last matrix

# Finally, we are going to create the objects of the class Gene.
# These objets contains the name of the gene and of the protein codified, and the Gene Ontology and Kegg information

gene_object = gene_name
(0..gene_name.size-1).each do |r|
        gene_object[r] = Gene.new(
        :gene_ID => gene_name[r],
        :protein_name => no_protein_synonyms[r], 
        :Kegg_ID_pathway => kegg_pathways[r],
        :Go_ID_Term => gene_ontology[r]
        )
end

#Here, in every position where a gene was defined to have interaction, we are going to save this gene objects
gene_position.each do |position| 
    objects_interaction.push gene_object[position]
end

# Now we obtein the protein networks by observing connections among arrays of interaction proteins
(0..interaction_level3.size-1).each do |s|
    position = []
    intact_genes = []
    network = []
    position.push s
    intact_genes.push gene_interaction[s]
    (0..interaction_level3.size-1).each do |r|
        relation = interaction_level3[s] & interaction_level3[r] 
        if relation != []
            position.push r
            intact_genes.push gene_interaction[r]
            network.push objects_interaction[r]
        end
    end
    position = position.uniq
    intact_genes = intact_genes.uniq
    network = network.uniq
    networks.push position
    networks_genes.push intact_genes
    networks_object.push network
end

#Now we are going to save the network objects
final_networks = networks_object
(0..networks_object.size-1).each do |x|
        final_networks[x] = InteractionNetwork.new(
        :genes_ID => networks_object[x],
        )
end

#Final output
output = []
(0..networks_genes.size-1).each do |y|
    network_A = networks_genes[y].shift #Use shift because we have an array and we want to show one element of it each time
    output.push network_A
end     

no =  0
yes = 0
(0..networks_genes.size-1).each do |z|
    if networks_genes[z] == []
        no +=1
    else
        puts "The gene #{output[z]} from the list have interactions with the genes #{networks_genes[z]} \n"
        yes +=1
    end
end
no_interactions = no_gene_interaction.size + no
puts "There are #{yes} interaction networks"
puts "There are #{no_interactions} genes which do not interact with other genes from the list"