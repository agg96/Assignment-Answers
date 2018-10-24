#Gene class
class Gene
    @@number_of_gene = 0
    attr_accessor :gene_ID 
    attr_accessor :gene_name
    attr_accessor :mutant_phenotype
    Wrong_id =[]
    def initialize params={}
        @@number_of_gene+=1
        @gene_ID = params.fetch(:gene_ID, 'unknow gene ID')
        #Here we check the gene_Id correct format. Wrong gene_id will return the message "Wrong ID is... and then the gene ID will be rejected.
        unless @gene_ID =~ /^A[Tt]\d[Gg]\d\d\d\d\d$/
            Wrong_id.push @gene_ID
            puts "Wrong ID is #{Wrong_id}"
            @gene_ID = "null"
            @gene_ID.reject { |k| k == "null"}
        end
        @gene_name = params.fetch(:gene_name, 'unknow gene name')
        @mutant_phenotype = params.fetch(:mutant_phenotype, 'unknow mutant phenotype')
    end
    def how_many 
      return @@number_of_gene
    end
end