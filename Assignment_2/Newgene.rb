#Here we define the new Gene class. We could have created a new Class and have included 
#this Class inside the other one we did in the previous assignment.

class Gene
    @@number_of_gene = 0
    attr_accessor :gene_ID 
    attr_accessor :protein_name
    attr_accessor :Kegg_ID_pathway
    attr_accessor :Go_ID_Term
    Wrong_ID =[]
    def initialize params={}
        @@number_of_gene+=1
        @gene_ID = params.fetch(:gene_ID, 'Unknown gene_ID')
        #Here we check the gene_Id correct format. Wrong gene_id will return the message "Wrong ID is... and then the gene ID will be rejected.
        unless @gene_ID =~ /^A[Tt]\d[Gg]\d\d\d\d\d$/
            Wrong_id.push @gene_ID
            puts "Wrong ID is #{Wrong_id}"
            @gene_ID = "null"
            @gene_ID.reject { |k| k == "null"}
        end
        @protein_name = params.fetch(:protein_name, 'Unknown protein_name')
        @Kegg_ID_pathway = params.fetch(:Kegg_ID_pathway, 'NO Kegg_ID_pathway')
        @Go_ID_Term = params.fetch(:Go_ID_Term, 'unknown gene_ontology_ID_Term')
    end
    def how_many 
      return @@number_of_gene
    end
end