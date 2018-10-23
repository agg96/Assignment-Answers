#Main program

#First, we call the tree objects
require './gene.rb'
require './seed_stock.rb'
require './cross.rb'
#Second, we are going to create matrix to save the data from the files
gene_information = [] #Saves the information of the archive gene_information.tsv
code_gene_information = [] #Saves the objects for the class gene information

seed_stock_data = [] #Saves the information of the archive seed_stock_data.tsv
code_seed_stock_data = [] #Saves the objects for the class seed stock data
seed_stock_objects =[]

cross_data = [] #Saves the information of the archive cross_data.tsv
code_cross_data = [] #Saves the objects for the class cross data
cross_data_objects1 = []
cross_data_objects2 = []

matrix_total =[] # This matrix will be used to generate an output

#Then we open the three archives

#Archive gene_information.tsv

#Open the file and save the data in gene_information
File.open ('gene_information.tsv') do |f1|
    while line =f1.gets
        @data=line.split' '
        gene_information.push @data     
    end
end

#Deletion of the first line of the data

gene_information.delete_at(0)

#Then, we save the names for each object in code_gene_information naming them with their gene ID
$length= gene_information.length
(0..$length-1).each do |i|
    code_gene_information.push gene_information[i][0]
end

#We create the objects with his attributes

(0..$length-1).each do |i|
      code_gene_information[i] = Gene.new(
      :gene_ID => gene_information [i][0],                  
      :gene_name => gene_information [i][1], 
      :mutant_phenotype => gene_information[i][2], 
      )
end

#Archive seed_stock_data
#Open the file and save the data in seed_stock_data
File.open ('seed_stock_data.tsv') do |f1|
    while line =f1.gets
        @data=line.split' '
        seed_stock_data.push @data
        matrix_total.push @data #We need this data in the matrix_total
    end
end

#Deletion of the first line of the data
seed_stock_data.delete_at(0)


#Then, we save the names for each object in code_seed_stock_data naming them as a combination between the first column (seed_stock)
#and the second colum (mutant_gene_ID)
$lenght= seed_stock_data.length
(0..$length-1).each do |a|
    code_seed_stock_data.push seed_stock_data[a][0]+" X "+seed_stock_data[a][1]
end

#We create the matrix seed_stock_objects with the objects previously created in code_gene_information.
#We keep a certain order in the objects for later adding them to the matrix "code_seed_stock_data".
#The attribute (Mutant_gene_ID) of the class StockDatabase will be an object of the class Gene
#With this the element 1 of the matrix "code_seed_stock_data" will correspond with the element 1 of the matrix code_gene_information
(0..$length-1).each do |a|
    (0..code_gene_information.length-1).each do |b|       
        if seed_stock_data [a][1] == code_gene_information[b].gene_ID
            seed_stock_objects.push code_gene_information[b]
        end
    end
end

#Now we will link the elements 

#We create the objects and then we asign the attributes to them. 
(0..$length-1).each do |c|
    code_seed_stock_data[c] = StockDatabase.new(
      :seed_Stock => seed_stock_data [c][0],
      :mutant_Gene_ID => seed_stock_objects[c], #The attribute mutant_Gene_ID is an object
      :last_Planted => seed_stock_data [c][2], 
      :storange  => seed_stock_data [c][3],
      :grams_Remaining => seed_stock_data [c][4]
      )
end

#Archive cross_data.tsv
#Open the file and save the data in cross_data
File.open ('cross_data.tsv') do |f1|
    while line =f1.gets
        @data=line.split' '
        cross_data.push @data     
    end
end

#Deletion of the first line of the data
cross_data.delete_at(0)

#Then, we save the names for each object in code_cross_data naming them as a combination between the parent1 and parent2 columns
$length= cross_data.length
(0..$length-1).each do |d| 
    code_cross_data.push cross_data[d][0]+" X "+cross_data[d][1]
end

#We compare if the field parent1 in cross data is equal to the attribute seed stock of the the objects of code_seed_stock_data
(0..$length-1).each do |x|
    (0..code_seed_stock_data.length-1).each do |y|
        if cross_data[x][0] == code_seed_stock_data[y].seed_Stock
            cross_data_objects1.push code_seed_stock_data[y]
        end
        if cross_data[x][1] == code_seed_stock_data[y].seed_Stock
            cross_data_objects2.push code_seed_stock_data[y]
        end
    end
end

#We create the objects with its atrributes.
(0..$length-1).each do |z|
      code_cross_data[z] = Cross.new(
      :parent1 => cross_data_objects1 [z],  #parent1 is an object                
      :parent2 => cross_data_objects2 [z],  #parent2 is an object
      :f2_wild => cross_data [z][2],
      :f2_P1 => cross_data [z][3],
      :f2_P2 => cross_data [z][4],
      :f2_P1P2 => cross_data [z][5],
      )
end

#-----------------------------------------
#FINAL RESULTS

#We create a loop that allows us to see if after planting 7 frams of seed still remains a quantity superior to 0 in every stock
#At the same time we put the grams of seeds remaining in our data
(0..code_seed_stock_data.length-1).each do |k|
    answer1 = code_seed_stock_data[k].planting
    if answer1 == "not true"
        code_seed_stock_data[k].grams_Remaining = code_seed_stock_data[k].grams_Remaining.to_i-7 
        matrix_total[k+1][4]= (matrix_total[k+1][4]).to_i-7
    else
        puts "WARNING: we have run out of Seed Stock #{answer1}" 
        code_seed_stock_data[k].grams_Remaining = 0
        matrix_total[k+1][4]=0
    end
end

#Finally, to check if two genes are linked we create a loop and call the chi_square function
(0..code_cross_data.length-1).each do |l|
    answer2= code_cross_data[l].chi_square
    if answer2 > 7.82
        puts "Recording: #{code_cross_data[l].parent1.mutant_Gene_ID.gene_name} is genetically linked to #{code_cross_data[l].parent2.mutant_Gene_ID.gene_name} with a chi-square score of #{answer2}"# THE USER SEE THIS
    end
end

#To see the final report I create a loop similar to the previous loop, but in this case I want to cluster all the conclusions.
puts "\nFinal Report:\n\n"

(0..code_cross_data.length-1).each do |l|
    answer2= code_cross_data[l].chi_square
    if answer2 > 7.82
        puts "#{code_cross_data[l].parent1.mutant_Gene_ID.gene_name} is linked to #{code_cross_data[l].parent2.mutant_Gene_ID.gene_name}"
        puts "#{code_cross_data[l].parent2.mutant_Gene_ID.gene_name} is linked to #{code_cross_data[l].parent1.mutant_Gene_ID.gene_name}\n\n"
    end
end

#  At the end, I create the file new_stock_file.tsv
newfile = File.new("new_stock_file.tsv", "w")

(0..matrix_total.length-1).each do |m|
    newfile.puts("#{matrix_total[m][0]} #{matrix_total[m][1]} #{matrix_total[m][2]} #{matrix_total[m][3]} #{matrix_total[m][4]}")
end
newfile.close
